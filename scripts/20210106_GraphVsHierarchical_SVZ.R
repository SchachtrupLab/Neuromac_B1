
#' ---
#' title: Graph-based vs Hierarchical clustering of SVZ
#' author: James Choi
#' date: 2021/1/6
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)

## If you don't have these packages installed already, run the following four 
## lines to install them.
# install.packages('Seurat')
# install.packages('tidyverse')
# install.packages('ggplot2')
# install.packages('dendextend')
require('Seurat')
require('ggplot2')
require('dplyr')
require('dendextend')

#' ### Background
#' Goal of this analysis is to determine how SCTransform-normalized cluster 
#' results compare to log-normalized cluster results in classifying MG vs NSC. 
#' We want to do this because our cluster results for MG show unknown clusters,
#' which may be possible NSC contamination. NSCs also show a small group of 
#' cells that are distinct from NSC (A, B, C cells), and this may potentially 
#' be MG. 
#' 

#' First, build dendrogram of cells in log-normalized PCA-space and compare to
#' graph-based clustering results.  
#' 

svz <- readRDS(file = '../data/20200511_SVZ.rds')

#+ log_umap, fig.height=3.5, fig.width=4
DefaultAssay(svz) <- 'RNA'
svz$graph_log_clust <- svz$RNA_snn_res.0.8
svz$graph_log_clust <- ifelse(
  test = svz$RNA_snn_res.0.8 %in% c(0,1,2,6),
  yes = 1,
  no = 2
)
DimPlot(svz, group.by = 'RNA_snn_res.0.8', pt.size = 2, label = TRUE, 
        label.size = 6) + 
  theme_bw() +
  labs(title = 'graph-based, log-normalized')
p1 <- DimPlot(svz,
              group.by = 'graph_log_clust',
              pt.size = 1.5,
              label = TRUE, 
              label.size = 5) + 
  theme_bw() +
  labs(title = 'graph-based, log-norm')

#+ log_dendrogram, fig.height=3.5, fig.width=5
svz_dist <- dist(x = svz[['pca']]@cell.embeddings[,1:11])
svz_tree <- hclust(d = svz_dist, method = 'ward.D2')
svz_dend <- as.dendrogram(svz_tree, hang = 1)
svz_dend <- dendsort::dendsort(d = svz_dend)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
tmp_cols <- gg_color_hue(length(levels(svz$RNA_snn_res.0.8)))
names(tmp_cols) <- levels(svz$RNA_snn_res.0.8)
labels_colors(svz_dend) <- 
  tmp_cols[svz$RNA_snn_res.0.8][order.dendrogram(svz_dend)]

par(mar = c(7.5, 2.5, 2.5, 0))
plot(svz_dend,
     main = 'Hierarchical clustering of SVZ (log-normalized)')
par(mar = c(5,4,4,2))

svz_log_dclus <- cutree(tree = svz_dend, h = 350)
svz$hier_log_clust <- 
  svz_log_dclus[match(rownames(svz@meta.data), names(svz_log_dclus))]
p2 <- DimPlot(svz,
              group.by = 'hier_log_clust',
              pt.size = 1.5,
              label = TRUE, 
              label.size = 5) +
  theme_bw() +
  labs(title = 'hierarchical, log-norm')

#' Next we test how SCT-normalized data affects downstream clusteirng. 
#' We choose 2000 variable genes even with SCTransform for consistency.  
#' 
svz <- NormalizeData(svz, verbose = FALSE)
svz <- SCTransform(svz, 
                   vars.to.regress = 'CC.difference', 
                   variable.features.n = 3000,
                   verbose = FALSE)
head(VariableFeatures(svz), 50)
tail(VariableFeatures(svz), 50)
svz <- RunPCA(svz, npcs = 40, verbose = FALSE)

#+ elbow, fig.height=3, fig.width=4
ElbowPlot(svz, ndims = 40) + theme_bw()

#+ dimloadings, fig.height=15, fig.width=16
VizDimLoadings(svz, dims = 1:15, ncol = 5)

#+ sct_umap, fig.height=3.5, fig.width=4
svz <- FindNeighbors(svz, dims = 1:15, verbose = FALSE)
svz <- FindClusters(svz, verbose = FALSE)
svz <- RunUMAP(svz, dims = 1:15, verbose = FALSE)
svz$graph_sct_clust <- svz$SCT_snn_res.0.8
svz$graph_sct_clust <- ifelse(
  test = svz$SCT_snn_res.0.8 %in% c(0,1,4,6,7),
  yes = 1,
  no = 2
)
DimPlot(svz, pt.size = 1, group.by = 'SCT_snn_res.0.8', label = TRUE, 
        label.size = 6) + 
  theme_bw() + 
  labs(title = 'graph-based, sct-normalized')
p3 <- DimPlot(svz, 
              group.by = 'graph_sct_clust',
              pt.size = 1.5,
              label = TRUE, 
              label.size = 5) + 
  theme_bw() + 
  labs(title = 'graph-based, SCT norm')

#' Graph-based clustering of SCT-normalized data shows roughly the same number
#' of clusters that we previously identified, with the exception of an
#' additional MG cluster.  
#' 

#+ sct_dendrogram, fig.height=3.5, fig.width=5
svz_dist <- dist(x = svz[['pca']]@cell.embeddings[,1:15])
svz_tree <- hclust(d = svz_dist, method = 'ward.D2')

svz_dend <- as.dendrogram(svz_tree, hang = 1)
svz_dend <- dendsort::dendsort(d = svz_dend)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
tmp_cols <- gg_color_hue(length(levels(svz$SCT_snn_res.0.8)))
names(tmp_cols) <- levels(svz$SCT_snn_res.0.8)
labels_colors(svz_dend) <- 
  tmp_cols[svz$SCT_snn_res.0.8][order.dendrogram(svz_dend)]

par(mar = c(7.5, 2.5, 2.5, 0))
plot(svz_dend,
     main = 'Hierarchical clustering of SVZ (SCTransform)')
par(mar = c(5,4,4,2))

svz_sct_dclus <- cutree(tree = svz_dend, h = 400)
svz$hier_sct_clust <- 
  svz_sct_dclus[match(rownames(svz@meta.data), names(svz_sct_dclus))]
p4 <- DimPlot(svz,
              group.by = 'hier_sct_clust',
              pt.size = 1.5,
              label = TRUE, 
              label.size = 5) +
  theme_bw() +
  labs(title = 'hierarchical, SCT-norm')

#' Hierarchical clustering of the cells in SCT-normalized PCA-space shows two
#' large branches, presumably MG vs NSC. We extract these results to compare 
#' with previous graph-based results of log-normazlied and sct-normalized data. 
#'   


#+ clust_comparison, fig.height=7.5, fig.width=9
cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)


#' Cluster results show that, even with hierarchical clustering, 100% accurate 
#' classification of MG vs NSC cannot be achieved.  
#' 


#' ### Conclusions
#' There is mixing of cell-types regardless of normalization/clustering method.
#' 