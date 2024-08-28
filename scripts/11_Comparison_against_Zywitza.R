
# Integrated analysis with Zywitza SVZ scSeq ---------------------------------

#' ## Comparative analysis with SVZ data produced by Zywitza, 2018.  
#' 

svz <- readRDS(file = '../data/20210108_SVZ.rds')
DefaultAssay(svz) <- 'RNA'

# Set file paths.
# To download zywitza counts file:
# 1. go here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111527
# 2. at supplementary files at bottom, download GSE111527_RAW.tar
# 3. unzip .tar file (should give you folder with .gz files)
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
zyw <- readRDS(file = paste0(results_out, 'Zywitza_SVZ.rds'))
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
