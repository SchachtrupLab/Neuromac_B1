
#' ---
#' title: Cluster Analysis of SVZ
#' author: James Choi
#' date: 2021/01/08
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)

# install.packages('dendextend')
setwd(dir = './scripts/')
require('Seurat')
require('ggplot2')
require('dplyr')
require('dendextend')
set.seed(123)

#' ## Background
#' It is apparent, based on initial cluster analyses, that we have more than two
#' cell-types. In addition to microglia and NSCs, we have an undefined 
#' population of Flt1(+) cells and Foxj1(+) cells. In order to have confidence
#' that downstream analysis of microglia and NSCs are not driven by 
#' contaminating cell-types, we repeat cluster analysis of all SVZ cells to 
#' ensure proper classification of cell-types.  
#' 
#' NOTE(!!!): This script is a work in progress - figures on gene expression, as
#' done in 20200511_PreliminaryClustering_allCells, needs to be added.  
#' 


# This Seurat object contains all cells passed quality control. Cluster analysis
# can be repeated.
svz <- readRDS(file = '../data/20200511_SVZ.rds')

#+ variable_features, fig.height=2.5, figh.width=4.5
DefaultAssay(svz) <- 'RNA'
summary(svz$nFeature_RNA)
svz <- NormalizeData(svz, verbose = FALSE)
svz <- SCTransform(svz, assay = 'RNA', variable.features.n = 2000, 
                   vars.to.regress = 'CC.difference', verbose = FALSE)
VariableFeaturePlot(svz) + theme_bw()

#+ pca_elbow, fig.height=2.5, fig.width=4
svz <- RunPCA(svz, npcs = 40, verbose = FALSE)
ElbowPlot(svz, ndims = 40) + theme_bw()

#+ pca_iteration, fig.height=6, fig.width=15, fig.cap="UMAPs using varying number of PCs"
test_pcs <- 10:19
umaps <- vector(mode = "list", length = length(test_pcs))
for (i in 1:length(test_pcs)) {
  svz <- FindNeighbors(svz, dims = 1:test_pcs[i], verbose = FALSE)
  svz <- FindClusters(svz, resolution = 1.5, verbose = FALSE)
  svz <- RunUMAP(svz, dims = 1:test_pcs[i], verbose = FALSE)
  umaps[[i]] <- DimPlot(svz, pt.size = 1, label = TRUE, label.size = 6) +
  theme_bw() + labs(title = paste("# PCs:", test_pcs[i]))
}
umaps <- cowplot::plot_grid(plotlist = umaps, ncol = length(test_pcs)/2)
umaps

#' UMAPs with varying number of PCs and high resolution shows that the small
#' group of cells successfully clusters out using lower number of PCs (10-12). 
#' To better choose PCs, we look at gene loadings for each component:  
#' 

#+ dimloadings, fig.height=4, fig.width=8
VizDimLoadings(svz, dims = 9:12, nfeatures = 30, ncol = 4)

#' We also inspect the distrbution of cells along thesse components:  
#' 
DimPlot(svz, pt.size = 2, reduction = 'pca', dims = c(10,11))
DimPlot(svz, pt.size = 2, reduction = 'pca', dims = c(10,12))
DimPlot(svz, pt.size = 2, reduction = 'pca', dims = c(11,12))

#' Seems that up to PC 11 has good distribution of cells, whereas PC 12 has a
#' weird skewed distribution.  
#' 


#+ svz_umap, fig.height=3, fig.width=4
svz <- FindNeighbors(svz, dims = 1:11, verbose = FALSE)
svz <- RunUMAP(svz, dims = 1:11, verbose = FALSE)
svz <- FindClusters(svz, resolution = 2, verbose = FALSE)
DimPlot(svz, pt.size = 2, label = TRUE, label.size = 6) + theme_bw()

#' Identify marker genes for each cluster to help classify as MG or NSC.  
#' 

Idents(svz) <- 'SCT_snn_res.2'
gene_markers <- FindAllMarkers(svz, assay = 'RNA', only.pos = TRUE)
top_markers <- gene_markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = -p_val_adj) %>%
  ungroup()
knitr::kable(x = top_markers, caption = 'Top 3 DE genes per cluster.')


#+ marker_genes, fig.width=15, fig.height=6
known_markers <- c('Ptprc','Itgam','Cx3cr1','Tmem119','P2ry12',
                   'Gfap','Dcx','Ascl1','Dlx1','Mki67')
known_markers_umaps <- vector(mode = 'list', length = length(known_markers))
for(i in 1:length(known_markers)) {
  known_markers_umaps[[i]] <- FeaturePlot(
    object = svz, 
    features = known_markers[i],
    cols = c('grey','red3'),
    pt.size = 1,
    order = TRUE
  ) +
    theme_bw() +
    labs(title = known_markers[i])
}
known_markers_umaps <- cowplot::plot_grid(
  plotlist = known_markers_umaps,
  ncol = 5
)
known_markers_umaps

#+ svz_celltype_umap, fig.height=3, fig.width=4.25
svz$celltype <- plyr::mapvalues(
  x = svz$SCT_snn_res.2,
  from = 0:12,
  to = c('NSC','Microglia','Microglia','Microglia','Microglia','NSC',
         'Microglia','NSC','NSC','Microglia','Microglia','Microglia', 
         'Unknown/Ependymal')
)
Idents(svz) <- 'celltype'
DimPlot(svz, pt.size = 2, label = TRUE, label.size = 6) + theme_bw()
# p1 <- DimPlot(svz, pt.size = 1.5, label = TRUE, label.size = 4) + theme_bw() + 
#   NoLegend()
# ggsave(filename = '../results/svz_umap.svg', plot = p1, device = 'svg', 
#        height = 3, width = 3.25)

saveRDS(svz, file = '../data/20210108_SVZ.rds')


#' We can also compare the data with a previously published study on NSCs vs 
#' ependymal cells in V-SVZ by [Shah, Stratton, et al.](https://doi.org/10.1016/j.cell.2018.03.063).
#' Specifically, we refer to Figure 3. Acta2 and Foxj1 are unique to ependymal
#' cells. Interestingly, they also note a small Flt1(+) population that we also
#' identified in the initial clustering of all SVZ cells.

#+ ependymal_genes, fig.height=3, fig.width=17.5
ependymal_genes <- c('Foxj1','Cfap126','Ccdc153', 'Rarres2','Tmem212')
ependymal_umaps <- FeaturePlot(svz, features = ependymal_genes, combine = FALSE)
ependymal_umaps <- lapply(
  X = ependymal_umaps,
  FUN = function(x) {
    x <- x + theme_bw()
    return(x)
  }
)
ependymal_umaps <- cowplot::plot_grid(plotlist = ependymal_umaps, ncol = 5)
ependymal_umaps
