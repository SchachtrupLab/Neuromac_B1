
#' ---
#' title: Examining MG variation in gene expression
#' author: James Choi
#' date: 2020-12-22
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center')

require('Seurat')
require('dplyr')
require('ggplot2')

# setwd('./scripts/')
mg <- readRDS(file = '../data/20200625_microglia.rds')


DefaultAssay(mg) <- 'RNA'

mg <- RunPCA(mg, npcs = 30, verbose = FALSE)
#+ pca_plots, fig.height=5, fig.width=12
pca_plots <- vector(mode = 'list', length = 10)
for (i in 1:10) {
  pca_plots[[i]] <- PCAPlot(
    object = mg,
    pt.size = 1,
    dims = c(i, 1)
  ) + NoLegend()
}
pca_plots <- cowplot::plot_grid(plotlist = pca_plots, ncol = 5)
pca_plots

#' PC_3, PC_5, PC_6, PC_8, PC_9, and PC_10 are dominated by few cells, but the
#' most blatant offender is PC_3. Which genes have the greatest weights for 
#' these components?
#'

#+ fig.height=10, fig.width=14
VizDimLoadings(
  object = mg,
  dims = 1:10,
  ncol = 5,
  balanced = TRUE
)

#' There are some components with genes that are typically not associated with
#' microglia. For example, Otx2 and Snap91 are neuron-associated genes and have
#' large weight values in component 2. Also, we notice that PC_1, which contains
#' the greatest amount of variance, is dominated by "Gm-" genes. Not
#' particularly informative.

mg <- FindNeighbors(mg, dims = 1:10, verbose = FALSE)
mg <- FindClusters(mg, verbose = FALSE)
mg <- RunUMAP(mg, dims = 1:10, verbose = FALSE)
DimPlot(mg, pt.size = 2)


#' Next, we inspect PCA results and dim loadings for when data are normalized 
#' using SCTransform.

DefaultAssay(mg) <- 'RNA'
mg <- SCTransform(mg, variable.features.n = 2000, verbose = FALSE)
mg <- RunPCA(mg, npcs = 30,verbose = FALSE)
mg <- FindNeighbors(mg, dims = 1:10, verbose = FALSE)
mg <- FindClusters(mg, verbose = FALSE)
mg <- RunUMAP(mg, dims = 1:10, verbose = FALSE)
DimPlot(mg, pt.size = 2)

tmp <- FindAllMarkers(mg, only.pos = TRUE)

mg <- readRDS(file = '../data/20200625_microglia.rds')
log_dist <- dist(mg[['pca']]@cell.embeddings[,1:10])


mg <- SCTransform(
  object = mg,
  variable.features.n = 2000,
  verbose = FALSE
)
DefaultAssay(mg) <- 'SCT'
mg <- RunPCA(mg, npcs = 30, verbose = FALSE)

sct_dist <- dist(mg[['pca']]@cell.embeddings[,1:10])

log_tree <- hclust(d = log_dist, method = 'ward.D2')
log_tree$labels <- seq_along(log_tree$labels)
log_dend <- as.dendrogram(log_tree)
sct_tree <- hclust(d = sct_dist, method = 'ward.D2')
sct_tree$labels <- seq_along(sct_tree$labels)
sct_dend <- as.dendrogram(sct_tree)

require('dendextend')
log_dend
leaf_groups <- mg$orig.ident
leaf_cols <- c('SCtrl' = 'black', 
               'S1d' = 'indianred',
               'S7d' = 'dodgerblue')
labels_colors(log_dend) <- leaf_cols[leaf_groups][order.dendrogram(log_dend)]
labels_colors(sct_dend) <- leaf_cols[leaf_groups][order.dendrogram(sct_dend)]

#+ figure.height=4, fig.width=12
par(mfrow = c(1,2))
plot(log_dend, 
     main = 'Dendrogram of original MG data',
     sub = 'SCtrl = black; S1d = indianred; S7d = dodgerblue')
plot(sct_dend, 
     main = 'Dendrogram of SCTransform MG data',
     sub = 'SCtrl = black; S1d = indianred; S7d = dodgerblue')