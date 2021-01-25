
#' ---
#' title: Comparative Analysis of Microglia Injury-States
#' author: James Choi
#' date: 2021/01/18
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)


# Load libraries
# BiocManager::install('SingleCellExperiment')
require('SingleCellExperiment')
require('Seurat')
require('ggplot2')
require('dplyr')
# set.seed(123)
results_out <- '../results/20210118_MG_comparative_analysis/'
dir.create(path = results_out)

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
saveRDS(ks_micro, file = paste0(results_out, '20210112_KerenShaul_MG.rds'))
saveRDS(dat_srat, file = paste0(results_out, '20210112_KerenShaul_full.rds'))


#' We successfully recapitulated the tSNE from the original Keren-Shaul paper, 
#' where a small "tail" of 5XFAD-enriched microglia are located.  
#' 

#' ## Integrated analysis of microglia from SVZ and Keren-Shaul.  
#' 
ks_micro <- readRDS(file = '../ref/20210112_KerenShaul_MG.rds')
ks_micro$orig.ident <- 'KerenShaul'
svz <- readRDS(file = '../data/20210108_SVZ.rds')
DefaultAssay(svz) <- 'RNA'
svz_micro <- svz[,svz$celltype == 'Microglia']

#' Use genes common to both datasets (note: Keren-Shaul data aligned to mm9 
#' reference genome).  
#' 
shared_genes <- rownames(svz_micro)[which(rownames(svz_micro) %in% rownames(ks_micro))]
ks_micro <- ks_micro[shared_genes,]
svz_micro <- svz_micro[shared_genes,]
micro <- list('svz' = svz_micro,
              'ad' = ks_micro)

#' Perform integration according to Seurat tutorial.  
#' 
micro <- lapply(
  X = micro,
  FUN = FindVariableFeatures,
  nfeatures = 5000,
  verbose = FALSE
)
micro_anchors <- FindIntegrationAnchors(
  object.list = micro,
  # anchor.features = micro_feats,
  verbose = FALSE
)
micro <- IntegrateData(
  anchorset = micro_anchors,
  verbose = FALSE
)
DefaultAssay(micro) <- 'integrated'
micro <- ScaleData(micro, verbose = FALSE)
micro <- RunPCA(micro, npcs = 30, verbose = FALSE)
ElbowPlot(micro, ndims = 30)
micro <- FindNeighbors(micro, dims = 1:15, verbose = FALSE)
micro <- RunTSNE(micro, dims = 1:15, verbose = FALSE)

#+ ks_integrated_micro_tsne, fig.height=3.5, fig.width=4.25, fig.cap='t-SNE of microglia taken from SVZ (current study) and Keren-Shaul (DAM study).'
DimPlot(micro, group.by = 'orig.ident', pt.size = 1) + theme_bw()
# micro <- FindClusters(micro, resolution = 0.4, verbose = FALSE)



# Comparative analysis with Zywitza SVZ -----------------------------------

#' ## Comparative analysis with SVZ data produced by Zywitza, 2018.  
#' 

svz <- readRDS(file = '../data/20210108_SVZ.rds')
DefaultAssay(svz) <- 'RNA'

# Set file paths.
zywitza_counts <- '../ref/GSE111527_RAW_Zywitza2018/'
batches <- list.files(zywitza_counts)
batches <- batches[!grepl('KO|WT', batches)]

# Load data from WT animals.
dat <- vector(mode = 'list', length = length(batches))
names(dat) <- batches
for (i in 1:length(batches)) {
  tmp <- read.table(
    file = gzfile(
      description = paste0(zywitza_counts, batches[i]),
      open = 'rt',
    ),
    header = TRUE,
    row.names = 1
  )
  dat[[batches[i]]] <- Matrix::Matrix(data = as.matrix(tmp),
                                      sparse = TRUE)
  message(paste('Loading file', batches[i], '...'))
}
shared_genes <- lapply(
  X = dat, 
  FUN = rownames
)
shared_genes <- Reduce(intersect, shared_genes)
zyw <- lapply(
  X = dat,
  FUN = function(x) {
    x[rownames(x) %in% shared_genes,]
  }
)

#' Create count matrix and perform QC. Thresholds for QC were taken from 
#' original study. Resulting number of cells is off by ~1-10 cells from the 
#' numbers in the original study.  
#' 
zyw <- Reduce(f = cbind, x = zyw)
zyw <- CreateSeuratObject(zyw, project = 'Zywitza')
zyw <- PercentageFeatureSet(zyw, pattern = '^mt-', col.name = 'percent.mt')
keep_cells <- zyw$nCount_RNA >= 500 & 
  zyw$nFeature_RNA >= 200 &
  zyw$percent.mt <= 10
zyw <- zyw[,keep_cells]
zyw <- NormalizeData(zyw, verbose = FALSE)
zyw <- FindVariableFeatures(zyw, verbose = FALSE)
zyw <- ScaleData(zyw, verbose = FALSE)
zyw <- RunPCA(zyw, verbose = FALSE)
ElbowPlot(zyw, ndims = 30)
zyw <- FindNeighbors(zyw, dims = 1:15, verbose = FALSE)
zyw <- RunTSNE(zyw, dims = 1:15, verbose = FALSE)
zyw <- FindClusters(zyw, resolution = 0.4, verbose = FALSE)

#+ zywitza_tsne, fig.height=3.5, fig.width=4, fig.cap='t-SNE of SVZ cells from Zywitza study'
DimPlot(zyw, label = TRUE, label.size = 5) + theme_bw()
saveRDS(zyw, file = paste0(results_out, 'Zywitza_SVZ.rds'))


#' Identify common genes in both datasets (Note: Zywitza aligned to mm10 
#' reference genome).  
#' 
shared_genes <- rownames(zyw)[rownames(zyw) %in% rownames(svz)]
zyw <- zyw[shared_genes,]
svz <- svz[shared_genes,]

#' Perform integration according to Seurat tutorial.  
#' 
svz <- list('svz' = svz,
            'zy' = zyw)
svz <- lapply(
  X = svz,
  FUN = FindVariableFeatures,
  nfeatures = 4000,
  verbose = FALSE
)
svz_anchors <- FindIntegrationAnchors(
  object.list = svz,
  verbose = FALSE
)
svz <- IntegrateData(
  anchorset = svz_anchors,
  verbose = FALSE
)
DefaultAssay(svz) <- 'integrated'
svz <- ScaleData(svz, verbose = FALSE)
svz <- RunPCA(svz, npcs = 30, verbose = FALSE)
ElbowPlot(svz, ndims = 30)
svz <- FindNeighbors(svz, dims = 1:15, verbose = FALSE)
svz <- RunTSNE(svz, dims = 1:15, verbose = FALSE)
svz <- FindClusters(svz, resolution = 0.4)
svz$orig.ident <- factor(svz$orig.ident, levels = rev(c('SCtrl','S1d','S7d','Zywitza')))

#+ svz_tsne, fig.height=3.5, fig.width=4, fig.cap='t-SNE of SVZ cells from current study and SVZ cells from Zywitza. Cells are colored by cluster identity.'
DimPlot(svz, label = TRUE, label.size = 6, pt.size = 1) + theme_bw()

#+ svz_markers, fig.height=3.5, fig.width=12, fig.cap='t-SNE of marker genes for microglia and NSPCs'
FeaturePlot(svz, features = c('Cx3cr1','Ascl1','Top2a','Dcx'), order = TRUE, ncol = 4)
# p1 <- FeaturePlot(svz, features = c('Cx3cr1','Ascl1','Top2a','Dcx'), order = TRUE, ncol = 4, combine = FALSE)
# p1 <- lapply(p1, FUN = function(x) {x <- x + theme_bw() + NoLegend(); return(x)})
# p1 <- cowplot::plot_grid(plotlist = p1, ncol = 4)
# ggsave(filename = paste0(results_out, 'Zywitza_comparison_geneMarkers.svg'),
#        plot = p1, device = 'svg', height = 3.25, width = 12)


#+ svz_tsne_byStudy, fig.height=3.5, fig.width=4.25, fig.cap='t-SNE of SVZ microglia (current study) and SVZ cells from Zywitza study. Cells are colored by study of origin or time after injury (current study's groups).'
DimPlot(svz, group.by = 'orig.ident', pt.size = 1, order = TRUE) + theme_bw() +
  scale_color_manual(values = c('Zywitza' = 'grey60',
                                'SCtrl' = 'red',
                                'S1d' = 'green',
                                'S7d' = 'blue'))
p1 <- DimPlot(svz, group.by = 'orig.ident', pt.size = 1, order = TRUE) +
  theme_bw() +
  scale_color_manual(values = c('Zywitza' = 'grey60',
                                'SCtrl' = 'red',
                                'S1d' = 'green',
                                'S7d' = 'blue'))
ggsave(filename = paste0(results_out, 'Zywitza_comparision_byGroup.svg'),
       plot = p1, device = 'svg', height = 4.5, width = 6)

#+ svz_tsne_byStudy, fig.height=3.5, fig.width=4.25, fig.cap='t-SNE of SVZ microglia (current study) and SVZ cells from Zywitza study. Colored cells denote cells from current and are colored by annotated cell-type.'
DimPlot(svz, group.by = 'celltype', pt.size = 2, shuffle = TRUE) + theme_bw()

#+ count_table, fig.cap='Mapping between cluster identites from combined analysis and original cell-type annotations (from current study).'
knitr::kable(x = table(svz$celltype, svz$integrated_snn_res.0.4))
saveRDS(svz, file = '../results/2021011MG_comparative_analysis/')


#' ## Combining all three datasets  
#' 
zyw <- readRDS(file = '../ref/20210112_Zywitza_SVZ.rds')
ks_micro <- readRDS(file = '../ref/20210112_KerenShaul_MG.rds')
svz_micro <- readRDS(file = '../data/20210112_MG.rds')
zyw$orig.ident <- 'Zywitza'
ks_micro$orig.ident <- 'KerenShaul'

full_dat <- list('zyw' = zyw,
                 'ks_micro' = ks_micro,
                 'svz_micro' = svz_micro)
full_dat <- lapply(full_dat, FUN = function(x) {DefaultAssay(x) <- 'RNA'; x})
genes <- lapply(full_dat, FUN = rownames)
shared_genes <- Reduce(f = intersect, x = genes)
full_dat <- lapply(full_dat, FUN = function(x) x[rownames(x) %in% shared_genes, ])

full_dat <- lapply(
  X = full_dat,
  FUN = FindVariableFeatures,
  nfeatures = 3000,
  verbose = FALSE
)
options(future.globals.maxSize = 1000 * 1024^2)
full_dat_anchors <- FindIntegrationAnchors(full_dat, verbose = FALSE)
full_dat <- IntegrateData(anchorset = full_dat_anchors, verbose = FALSE)
DefaultAssay(full_dat) <- 'integrated'
full_dat <- ScaleData(full_dat, verbose = FALSE)
full_dat <- RunPCA(full_dat, npcs = 30, verbose = FALSE)
ElbowPlot(full_dat, ndims = 30)
full_dat <- FindNeighbors(full_dat, dims = 1:15, verbose = FALSE)
full_dat <- RunTSNE(full_dat, dims = 1:15, verbose = FALSE)
full_dat <- FindClusters(full_dat, resolution = 0.4)
DimPlot(full_dat, pt.size = 1, group.by = 'orig.ident', order = TRUE)
# DimPlot(full_dat, pt.size = 1)
DefaultAssay(full_dat) <- 'RNA'
FeaturePlot(full_dat, 'Cx3cr1', order = TRUE)

p1 <- DimPlot(full_dat, pt.size = 2, group.by = 'orig.ident', shuffle = TRUE)
p2 <- FeaturePlot(full_dat, 'Cx3cr1', order = TRUE)
p1 + p2



# DE Gene Comparison ------------------------------------------------------

mg <- readRDS(file = '../data/20210112_MG.rds')

DefaultAssay(mg) <- 'RNA'
mg[['tsne']] <- NULL
Idents(mg) <- 'SCT_snn_res.0.4'
DimPlot(mg, pt.size = 2, label = TRUE, label.size = 6) + NoLegend()

#' Data used are the following:  
#' + 

ks_dam <- read.table(file = '../ref/KerenShaul2017_TableS3_HomeostaticMG_vs_DAM.txt',
                     sep = '\t', header = TRUE)
li_p7c1 <- read.table(file = '../ref/Li2019_TableS5_MG_P7C1_DEgenes.txt',
                      sep = '\t', header = TRUE)
hammond_irm <- read.table(file = '../ref/Hammond2019_TableS1_IRM2_DEgenes.txt',
                          sep = '\t', header = TRUE)
hammond_c4 <- read.table(file = '../ref/Hammond2019_TableS1_MGcluster4_DEgenes.txt',
                         sep = '\t', header = TRUE)

mg_c1 <- FindMarkers(
  object = mg, 
  ident.1 = 1,
  ident.2 = 0,
  only.pos = TRUE,
  assay = 'RNA'
)
mg_c4 <- FindMarkers(
  object = mg, 
  ident.1 = 4,
  ident.2 = 0,
  only.pos = TRUE,
  assay = 'RNA'
)


# Keren-Shaul, 2017, DAMs vs Homeostatic Microglia
ks_dam_set <- with(
  data = ks_dam,
  expr = Fold.change..DAM.to.homeostatic.microglia. > 0 &
    DAM.FDR.p.value >= 3 &
    !is.na(DAM.FDR.p.value))
ks_dam_set <- unique(ks_dam[ks_dam_set,'Gene.name'])
length(ks_dam_set)

# Li, 2019, P7-C2 (phagocytic?) microglia vs other P7 clusters
li_p7c1_set <- with(
  data = li_p7c1,
  expr = avg_logFC > 0 &
    p_val_adj <= 0.001
)
li_p7c1_set <- unique(li_p7c1[li_p7c1_set,'gene'])
length(li_p7c1_set)

# Hammond, 2019, Injury-responsive microglia cluster 2 (IRM2) vs other microglia
# (which come from control and saline-injected)
hammond_irm_set <- with(
  data = hammond_irm,
  expr = Fold.Change > 0 & # all are positive bc theyre marker genes
    Padj..FDR. <= 0.001
)
hammond_irm_set <- unique(hammond_irm[hammond_irm_set,'Gene'])
length(hammond_irm_set)

# Hammond, 2019, P4/5 cluster 4 vs other microglia (Axon tract-associated MG)
hammond_c4_set <- with(
  data = hammond_c4,
  expr = Fold.Change > 0 & # all are positive bc theyre marker genes
    Padj..FDR. <= 0.001
)
hammond_c4_set <- unique(hammond_c4[hammond_c4_set,'Gene'])
length(hammond_c4_set)

# SVZ Microglia cluster 1
mg_c1_set <- with(
  data = mg_c1,
  expr = avg_logFC > 0 &
    p_val_adj <= 0.001
)
mg_c1_set <- unique(rownames(mg_c1[mg_c1_set,]))
length(mg_c1_set)

# SVZ Microglia cluster 4
mg_c4_set <- with(
  data = mg_c4,
  expr = avg_logFC > 0 &
    p_val_adj <= 0.001
)
mg_c4_set <- unique(rownames(mg_c4[mg_c4_set,]))
length(mg_c4_set)


require('eulerr')

mg_de_c1 <- list(
  'DAM (KerenShaul)' = ks_dam_set,
  # 'PAM (Li)' = li_p7c1_set,
  'IRM (Hammond)' = hammond_irm_set,
  # 'ATM (Hammond)' = hammond_c4_set,
  'Cluster 1 (SVZ)' = mg_c1_set
)
svg(filename = '../results/MG_C1_comparison.svg', height = 4, width = 5)
plot(euler(mg_de_c1), quantities = TRUE)
dev.off()

mg_de_c4 <- list(
  'DAM (KerenShaul)' = ks_dam_set,
  # 'PAM (Li)' = li_p7c1_set,
  'IRM (Hammond)' = hammond_irm_set,
  # 'ATM (Hammond)' = hammond_c4_set,
  'Cluster 4 (SVZ)' = mg_c4_set
)
svg(filename = '../results/MG_C4_comparison.svg', height = 4, width = 5)
plot(euler(mg_de_c4), quantities = TRUE)
dev.off()


mg_de_c4 <- list(
  'DAM (KerenShaul)' = ks_dam_set,
  'PAM (Li)' = li_p7c1_set,
  'IRM (Hammond)' = hammond_irm_set,
  'ATM (Hammond)' = hammond_c4_set,
  'Cluster 4 (SVZ)' = mg_c4_set
)
mg_de_c4 <- lapply(
  X = mg_de_c4,
  FUN = gsub,
  pattern = '-',
  replacement = '_'
)

require('UpSetR')
upset(data = fromList(input = mg_de_c4), order.by = 'freq', group.by = 'sets')
