
#' ---
#' title: "SVZ scRNAseq preprocessing"
#' author: "James Choi"
#' date: "`r Sys.Date()`"
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)


# Load libraries
require('Seurat')
require('dplyr')
require('ggplot2')
require('SingleCellExperiment')
set.seed(123)

# Directories for data and results
# setwd('./scripts/')
results_out <- '../results/preprocessing/'
dir.create(results_out)

# Load counts matrix
svz <- read.table(
  file = '../data/raw_counts_all_cells.csv',
  sep = ',',
  stringsAsFactors = FALSE,
  row.names = 1,
  header = TRUE
)
# Convert to sparse, round to integers (has decimals for some reason)
svz <- Matrix::Matrix(data = as.matrix(svz), sparse = TRUE)
svz <- round(svz)
dat_sce <- SingleCellExperiment(assay = list('counts' = svz))
svz <- CreateSeuratObject(counts = svz, project = 'SVZ', assay = 'RNA')

#' First, we plot some sequencing metrics such as UMIs per cell and number of
#' unique genes detected per cell.
#+ detection_scatter, fig.height=4, fig.width=4.5, fig.cap='Scatter plot of sequencing metrics by cell.'
qc <- svz@meta.data
qc <- qc[sample(x = nrow(qc), size = nrow(qc)),]
plot(x = qc$nCount_RNA,
     y = qc$nFeature_RNA,
     xlab = 'total UMI',
     ylab = '# genes detected',
     main = 'Detection rates by sample',
     sub = 'Black = SCtrl; Red = S1d; Green = S7d',
     col = sapply(X = as.character(qc$orig.ident), FUN = function(x) switch(x, 'SCtrl' = 'black', 'S1d' = 'red3', 'S7d' = 'green')))


#' From the above, we see that there are a few cells with outlier numbers of 
#' total UMI and genes detected. These are likely to be doublets. Beyond that, 
#' there are no obvious biases based on sample of origin.  

#' Next we examine cell-level QC metrics.  
#' 
#' * nCount_RNA = number of UMIs per cell.
#' * nFeature_RNA = number of unique genes detected per cell.
#' * percent.mt = percent of UMIs mapping to a mitochondrial protein-encoding 
#' gene.
#' * log10_total_counts = log10(nCount_RNA)
#+ qc_vln, fig.height=3.5, fig.width=6, fig.cap='Violin plot of various per-cell QC metrics.'
svz <- PercentageFeatureSet(svz, pattern = '^mt-', col.name = 'percent.mt')
svz$log10_total_counts <- log10(svz$nCount_RNA)
qc_plot <- svz@meta.data %>%
  select(-c('log10_total_counts')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  ggplot(mapping = aes(x = orig.ident, y = value)) +
  geom_violin(scale = 'width') +
  ggbeeswarm::geom_quasirandom(size = 1) +
  facet_wrap(. ~ variable, scales = 'free_y') +
  labs(title = 'Quality control metrics') +
  theme(axis.title = element_blank(),
        plot.title = element_text(size = 18, face = 'bold'),
        axis.line.y = element_line(size = 1),
        axis.ticks.y = element_line(size = 1))
qc_plot


#' In general, the samples share very similar distributions across quality 
#' control metric - there are a few outlier cells (see nCount_RNA). However, we 
#' note that the nCount_RNA and nFeature_RNA appear to have bimodal
#' distributions. There is also a long tail of cells with unusually high 
#' percent.mt values (mean of `r round(mean(svz$percent.mt), 1)`%). These cells
#' collectively should be considered suspect and removed from downstream.  

umi_lo_threshold <- 1000
# umi_hi_threshold <- median(svz$log10_total_counts) + 2.6*mad(svz$log10_total_counts, constant = 1)
feat_hi_threshold <- 6000
# log_umi_hi_threshold <- median(svz$log10_total_counts) + 3*mad(x = svz$log10_total_counts, constant = 1)
percent_mt_threshold <- 25

bad_cells <- svz$nCount_RNA < umi_lo_threshold | 
  # svz$log10_total_counts > log_umi_hi_threshold | 
  svz$nFeature_RNA > feat_hi_threshold |
  svz$percent.mt > percent_mt_threshold
table(bad_cells)
svz$outlier <- bad_cells
svz <- svz[,!svz$outlier]


#+ variable_features, fig.height=2.5, figh.width=4.5
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}
DefaultAssay(svz) <- 'RNA'
summary(svz$nFeature_RNA)
svz <- NormalizeData(svz, verbose = FALSE)
svz <- CellCycleScoring(
  object = svz,
  s.features = intersect(firstup(cc.genes$s.genes), rownames(svz)),
  g2m.features = intersect(firstup(cc.genes$g2m.genes), rownames(svz))
)
svz$CC.difference <- svz$S.Score - svz$G2M.Score
svz <- SCTransform(svz, 
                   assay = 'RNA', 
                   variable.features.n = 2000, 
                   vars.to.regress = 'CC.difference', 
                   verbose = FALSE)
# VariableFeaturePlot(svz) + theme_bw()

#+ pca_elbow, fig.height=2.5, fig.width=4
svz <- RunPCA(svz, npcs = 40, verbose = FALSE)
# ElbowPlot(svz, ndims = 40) + theme_bw()

#+ pca_iteration, fig.height=6, fig.width=15, fig.cap="UMAPs using varying number of PCs"
# test_pcs <- 10:19
# umaps <- vector(mode = "list", length = length(test_pcs))
# for (i in 1:length(test_pcs)) {
#   svz <- FindNeighbors(svz, dims = 1:test_pcs[i], verbose = FALSE)
#   svz <- FindClusters(svz, resolution = 1.5, verbose = FALSE)
#   svz <- RunUMAP(svz, dims = 1:test_pcs[i], verbose = FALSE)
#   umaps[[i]] <- DimPlot(svz, pt.size = 1, label = TRUE, label.size = 6) +
#     theme_bw() + labs(title = paste("# PCs:", test_pcs[i]))
# }
# umaps <- cowplot::plot_grid(plotlist = umaps, ncol = length(test_pcs)/2)
# umaps

#' UMAPs with varying number of PCs and high resolution shows that the small
#' group of cells successfully clusters out using lower number of PCs (10-12). 
#' To better choose PCs, we look at gene loadings for each component:  
#' 

#+ dimloadings, fig.height=4, fig.width=8
# VizDimLoadings(svz, dims = 9:12, nfeatures = 30, ncol = 4)

#' We also inspect the distrbution of cells along thesse components:  
#' 
# DimPlot(svz, pt.size = 2, reduction = 'pca', dims = c(10,11))
# DimPlot(svz, pt.size = 2, reduction = 'pca', dims = c(10,12))
# DimPlot(svz, pt.size = 2, reduction = 'pca', dims = c(11,12))

#' Seems that up to PC 11 has good distribution of cells, whereas PC 12 has a
#' weird skewed distribution.  
#' 


#+ svz_umap, fig.height=3, fig.width=4
svz <- FindNeighbors(svz, dims = 1:11, verbose = FALSE)
svz <- RunUMAP(svz, dims = 1:11, verbose = FALSE)
svz <- FindClusters(svz, resolution = 2, verbose = FALSE)
DimPlot(svz, pt.size = 2, label = TRUE, label.size = 6) + theme_bw()
