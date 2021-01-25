
#' ---
#' title: Cluster Analysis of SVZ
#' author: James Choi
#' date: 2021/01/12
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)


# Load libraries
require('SingleCellExperiment')
require('Seurat')
require('ggplot2')
require('dplyr')

#' ### Extract microglia from Keren-Shaul, 2017.  
#' 

# Set file paths.
kerenshaul_counts <- '../ref/GSE98969_RAW_KerenShaul2017/'
kerenshaul_design <- '../ref/GSE98969_KerenShaul2017_experimental_design_f.txt'

# Load experimental design table, extract cells from 5XFAD/WT mice.
design <- read.table(file = kerenshaul_design,
                     header = TRUE, sep = '\t', fill = TRUE, comment.char = '')
design$Batch_desc <- gsub(pattern = ' ',
                          replacement = '_',
                          x = design$Batch_desc)

#' Original study has multiple experiments, extract only AD vs WT comparison 
#' experiment (Figure 1D). Note: authors have not provided cell-type
#' classification results in their data deposit, so we have to cluster ourselves
#' and identify microglia.  
batch_desc_keep <- c(
  'AD6m_mouse1_plate#1',
  'AD6m_mouse1_plate#2',
  'AD6m_mouse1_plate#3',
  'AD6m_mouse1_plate#4',
  'AD6m_mouse2_plate#1',
  'AD6m_mouse2_plate#2',
  'AD6m_mouse2_plate#3',
  'AD6m_mouse2_plate#4',
  'AD6m_mouse2_plate#5',
  "AD6m_mouse2_plate#6",
  "AD6m_mouse3_plate#1",
  "AD6m_mouse3_plate#2",
  "AD6m_mouse3_plate#3",
  "AD6m_mouse3_plate#4",
  "AD6m_mouse3_plate#5",
  "AD6m_mouse3_plate#6",
  "WT6m_mouse1_plate#1",
  "WT6m_mouse1_plate#2",
  "WT6m_mouse1_plate#3",
  "WT6m_mouse1_plate#4",
  "WT6m_mouse2_plate#1",
  "WT6m_mouse2_plate#2",
  "WT6m_mouse2_plate#3",
  "WT6m_mouse2_plate#4",
  "WT6m_mouse2_plate#5",
  "WT6m_mouse2_plate#6",
  "WT6m_mouse3_plate#1",
  "WT6m_mouse3_plate#2",
  "WT6m_mouse3_plate#3",
  "WT6m_mouse3_plate#4",
  "WT6m_mouse3_plate#5",
  "WT6m_mouse3_plate#6"
)
dam_keep <- with(data = design,
                 expr = Mouse_ID %in% c('5XFAD', 'C57BL/6') &
                   Batch_desc %in% batch_desc_keep &
                   Number_of_cells == 1)
dam_keep <- design[dam_keep,]
dam_keep_ampbatch <- unique(dam_keep$Amp_batch_ID)

#' Extract count matrices. Http download from GEO stored as .txt.gz. Convert to
#' sparse matrix for memory.  
plates <- grep(pattern = '.txt.gz', 
               x = list.files(kerenshaul_counts), 
               value = TRUE)
dat <- vector(mode = 'list', length = length(dam_keep_ampbatch))
names(dat) <- dam_keep_ampbatch
for (i in 1:length(plates)) {
  plate_ampbatch <- strsplit(
    x = sapply(X = strsplit(x = plates[i], split = '\\.'), FUN = `[`, 1),
    split = '_'
  )[[1]][2]
  if (plate_ampbatch %in% dam_keep_ampbatch) {
    tmp <- read.table(
      file = gzfile(
        description = paste0(kerenshaul_counts, plates[i]),
        open = 'rt'
      )
    )
    dat[[plate_ampbatch]] <- Matrix::Matrix(data = as.matrix(tmp),
                                            sparse = TRUE)
    message(paste('Loading file', plate_ampbatch, '...'))
  }
}

#' Extract cell barcodes from 5XFAD/WT mice from above. Column merge.  
dat <- lapply(
  X = dat,
  FUN = function(x) {
    x[,colnames(x) %in% dam_keep$Well_ID]
  }
)
dat <- Reduce(f = cbind, x = dat)

#' Normalize transcript abundance using ERCC spike-ins, which assumes equal 
#' concentration per cell.  
dat_sce <- SingleCellExperiment(assays = list(counts = dat))
is_spikein <- ifelse(grepl('ERCC', rownames(dat_sce)),'ERCC','endog')
dat_sce <- splitAltExps(x = dat_sce, f = is_spikein)
dat_sce <- dat_sce[,Matrix::colSums(counts(dat_sce)) > 500]
dat_sce <- scran::computeSpikeFactors(x = dat_sce, spikes = 'ERCC')
summary(sizeFactors(dat_sce))
dat_sce <- scater::logNormCounts(x = dat_sce, 
                                 size_factors = sizeFactors(dat_sce))

#' Proceed with basic clustering.  
dat_srat <- as.Seurat(x = dat_sce)
dat_srat@project.name <- 'KerenShaul'
dat_srat <- FindVariableFeatures(dat_srat, verbose = FALSE)
dat_srat <- ScaleData(dat_srat, verbose = FALSE)
dat_srat <- RunPCA(dat_srat, verbose = FALSE)

#+ ks_elbow, fig.height=2.5, fig.width=3.5, fig.cap='PCA elbow plot for Keren-Shaul, 2017'
ElbowPlot(dat_srat, ndims = 30)
dat_srat <- FindNeighbors(dat_srat, dims = 1:10, verbose = FALSE)
dat_srat <- RunTSNE(dat_srat, dims = 1:10, verbose = FALSE)
dat_srat <- FindClusters(dat_srat, resolution = 0.4, verbose = FALSE)

#+ ks_tsne, fig.height=3.5, fig.width=4, fig.cap='t-SNE result for Keren-Shaul, 2017. There will be some differences between here and their Figure 1A because different number of cells.'
DimPlot(dat_srat, label = TRUE, reduction = 'tsne') + theme_bw()

#' Next, identify microglia.  
#+ ks_microglia_tsne, fig.height=5.5, fig.width=12.5, fig.cap='Canonical microglia marker genes.'
microglia_markers <- c('Hexb','Cx3cr1','P2ry12','Tmem119','Fcrls','C1qa')
microglia_tsne <- FeaturePlot(
  object = dat_srat,
  features = microglia_markers,
  order = TRUE,
  cols = c('grey','red3'),
  combine = FALSE
)
microglia_tsne <- lapply(X = microglia_tsne, FUN = function(x) x + theme_bw())
microglia_tsne <- cowplot::plot_grid(plotlist = microglia_tsne, ncol = 3)
microglia_tsne

#+ ks_microglia_vln, fig.
VlnPlot(dat_srat, features = microglia_markers, pt.size = 0, ncol = 2)

#' Find differentially expressed genes.  
de_genes <- FindAllMarkers(
  object = dat_srat,
  only.pos = TRUE,
  logfc.threshold = 0.25
)
top_de <- de_genes %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_logFC)
knitr::kable(x = top_de)

#' Extract the microglia.  
keep_clusters <- c(0,1,2,5)
ks_micro <- dat_srat[, dat_srat$seurat_clusters %in% keep_clusters]

# Add metadata
ks_micro$source <- plyr::mapvalues(
  x = colnames(ks_micro),
  from = design$Well_ID,
  to = design$Mouse_ID,
  warn_missing = FALSE
)

#+ microglia_AD, fig.height=3.5, fig.width=4, fig.cap='t-SNE of microglia by mouse genotype (5XFAD = AD mouse model)'
DimPlot(ks_micro, group.by = 'source', pt.size = 1) + theme_bw()
saveRDS(ks_micro, file = '../ref/20210112_KerenShaul_MG.rds')

#' We successfully recapitulated the tSNE from the original Keren-Shaul paper, 
#' where a small "tail" of 5XFAD-enriched microglia are located.  
#' 
ks_micro <- readRDS(file = '')