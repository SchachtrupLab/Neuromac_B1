
#' ---
#' title: Revisiting SVZ clustering
#' author: James Choi
#' date: 2020-12-21
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center')

######## 20201221 exploratory analysis - revisiting clustering #########


require('Seurat')
require('dplyr')
require('ggplot2')

# setwd('./scripts/')
results_out <- '../results/20201221_revisitingClustering/'
dir.create(path = results_out)

svz <- readRDS(file = '../data/20200511_SVZ.rds')


#' ## Reclustering all SVZ cells

#' Inspect the number of uniquely expressed genes per cell to determine whether
#' the number of variable features we selected (2000) is reasonable. I think we
#' want to aim for ~median(nFeature_RNA).

#+ fig.height=3, fig.width=6
Idents(svz) <- ''
VlnPlot(svz, features = c('nFeature_RNA','nCount_RNA','percent.mt'), ncol = 3)
summary(svz$nFeature_RNA)

#' We test SCTransform method of normalization + variable feature selection. 
#' Purpose of this is to evaluate how cluster results compare to initial result,
#' given that there are some microglial clusters that cannot be validated (i.e.
#' Developmental MG [Krt genes], and Unknown MG ["GM-"]).
DefaultAssay(svz) <- 'RNA'
log_feats <- svz[['RNA']]@var.features

#+ warning=FALSE
svz <- SCTransform(
  object = svz, 
  vars.to.regress = 'CC.difference',
  variable.features.n = 2000,
  verbose = FALSE
)
sct_feats <- svz[['SCT']]@var.features

#' There are substantial differences in variable gene selection depending on 
#' normalization method.
shared_feats <- intersect(log_feats, sct_feats)
sct_feats[!sct_feats %in% shared_feats]
log_feats[!log_feats %in% shared_feats]

#' Perform usual Seurat analysis.
#+ fig.width=4.5, fig.height=3
svz <- RunPCA(svz, npcs = 30, verbose = FALSE)
ElbowPlot(svz, ndims = 30)
svz <- RunUMAP(svz, dims = 1:10, verbose = FALSE)
svz <- FindNeighbors(svz, dims = 1:10, verbose = FALSE)
svz <- FindClusters(svz, resolution = 0.8, verbose = FALSE) #default

#+ fig.width=10, fig.height=9
p1 <- DimPlot(svz, label = TRUE, pt.size = 2)
p2 <- DimPlot(svz, label = TRUE, group.by = 'orig.ident', pt.size = 2)
p3 <- DimPlot(svz, label = TRUE, group.by = 'RNA_snn_res.0.8', pt.size = 2)
cowplot::plot_grid(p1, p2, p3, ncol = 2)


#' By eye, NSC clusters more or less agree well with previous cluster results. 
#' Major differences arise from cells which I presume to be homeostatic-MG. 
#' Identify differentially expressed genes per cluster to see if "Gm" cluster
#' is reproduced.
sct_de <- FindAllMarkers(
  object = svz, 
  only.pos = TRUE,
  assay = 'RNA',
  verbose = FALSE
)
sct_de_top <- sct_de %>%
  group_by(cluster) %>%
  filter(p_val_adj < 1e-03) %>%
  top_n(wt = avg_logFC, n = 5)
print(data.frame(sct_de_top))

table(svz$SCT_snn_res.0.8, svz$orig.ident)



#' ## Reclustering Microglia cells with SCTransform

#+ fig.height = 4.5, fig.width=5.5
mg <- readRDS('../data/20200625_microglia.rds')

#+ fig.height=4, fig.width=4.5, fig.align='center'
DimPlot(mg, label = TRUE, pt.size = 2, label.size = 8)


#+ warnings=FALSE
DefaultAssay(mg) <- 'RNA'
mg <- SCTransform(
  object = mg,
  variable.features.n = 2000,
  verbose = FALSE
)
log_feats <- mg[['RNA']]@var.features
sct_feats <- mg[['SCT']]@var.features
shared_feats <- intersect(log_feats, sct_feats)
# Genes that are selected only using SCTransform
sct_feats[!sct_feats %in% shared_feats]
# Genes that are selected only using log-normalization
log_feats[!log_feats %in% shared_feats]


#+ fig.width=4.5, fig.height=3
mg <- RunPCA(mg, npcs = 30, verbose = FALSE)
ElbowPlot(mg, ndims = 30)
mg <- RunUMAP(mg, dims = 1:10, verbose = FALSE)
mg <- FindNeighbors(mg, dims = 1:10, verbose = FALSE)
mg <- FindClusters(mg, resolution = 0.8, verbose = FALSE)

#+ fig.width=10, fig.height=9
p1 <- DimPlot(mg, label = TRUE, pt.size = 2, label.size = 8)
p2 <- DimPlot(mg, label = TRUE, group.by = 'orig.ident', pt.size = 2,
              label.size = 8)
p3 <- DimPlot(mg, label = TRUE, group.by = 'RNA_snn_res.0.8', pt.size = 2,
              label.size = 8)
cowplot::plot_grid(p1, p2, p3, ncol = 2)

#' Clusters using SCTransform are different from log-normalization. Log-norm 
#' clusters 0 and 1 are now more split into 3. Log-norm cluster 2 is split into
#' 2 clusters. Now, identify DE genes to determine if new clusters are distinct.
mg_sct_de <- FindAllMarkers(
  object = mg,
  only.pos = TRUE,
  verbose = FALSE
)
mg_sct_top <- mg_sct_de %>%
  group_by(cluster) %>%
  filter(p_val_adj < 1e-3) %>%
  top_n(5, wt = avg_logFC)
data.frame(mg_sct_top)


#' Some of the DE results suggest that the Activated-MG can be further split 
#' into two distinct clusters. Meanwhile, homeostatic-MG clusters 0, 2, and 3 
#' are still characterized by "Gm-" genes, some "mt-" genes, and other genes for
#' which we cannot assign any functional relevance. For now, we focus on the two
#' activated-MG clusters to determine if they are visually distinct.
#+ fig.width=7, fig.height=15, warnings=FALSE, fig.align='center'
act_mg_genes <- mg_sct_top %>%
  filter(cluster %in% c(1, 5)) %>%
  ungroup() %>%
  select(gene) %>%
  .[['gene']]
DefaultAssay(mg) <- 'RNA'
act_mg_genes_umap <- FeaturePlot(
  object = mg, 
  features = act_mg_genes,
  pt.size = 3,
  combine = FALSE,
  reduction = 'umap',
  order = TRUE
)
act_mg_genes_umap <- lapply(
  X = act_mg_genes_umap,
  FUN = function(x) {
    x + 
      theme_bw() +
      scale_color_gradientn(
        colours = rev(RColorBrewer::brewer.pal(n = 8, name = 'Spectral'))
      )
  }
)
act_mg_genes_umap <- cowplot::plot_grid(plotlist = act_mg_genes_umap, ncol = 2,
                                        byrow = FALSE)
act_mg_genes_umap


#' Based on feature plots of the top genes, there is some visible distinction 
#' between the two activated microglia clusters. How are the clusters spread 
#' over time?
#+ fig.height=3, fig.width=4, fig.align='center'
prop <- round(prop.table(
  table(mg$SCT_snn_res.0.8, mg$orig.ident), 
  margin = 2)*100,1) %>%
  as.data.frame() %>%
  ggplot(mapping = aes(x = Var2, y = Freq)) +
  geom_bar(mapping = aes(fill = Var1), position = 'stack', stat = 'identity') +
  theme_bw()
prop




#' ## Reclustering Microglia cells with SCTransform (BY INJURY TIME-POINT)

#' We repeat above analysis but this time perform SCTransform on each time-point
#' individually. The idea is that different preps from different injury times 
#' will have different ambient RNA profiles, giving rise to batch effects that
#' can be modeled out by using SCTransform's relative normalization procedure.

#+ warnings=FALSE
# The argument "batch_var" allows building variance~mean parameter models with
# batch_var included as interaction term (whatever that means).
mg <- readRDS('../data/20200625_microglia.rds')
DefaultAssay(mg) <- 'RNA'
mg <- SCTransform(
  object = mg,
  variable.features.n = 2000,
  batch_var = 'orig.ident',
  verbose = FALSE
)
log_feats <- mg[['RNA']]@var.features
sct_feats <- mg[['SCT']]@var.features
shared_feats <- intersect(log_feats, sct_feats)
# Genes that are selected only using SCTransform
sct_feats[!sct_feats %in% shared_feats]
# Genes that are selected only using log-normalization
log_feats[!log_feats %in% shared_feats]


#+ fig.width=4.5, fig.height=3
mg <- RunPCA(mg, npcs = 30, verbose = FALSE)
ElbowPlot(mg, ndims = 30)
mg <- RunUMAP(mg, dims = 1:10, verbose = FALSE)
mg <- FindNeighbors(mg, dims = 1:10, verbose = FALSE)
mg <- FindClusters(mg, resolution = 0.8, verbose = FALSE)

#+ fig.width=10, fig.height=9
p1 <- DimPlot(mg, label = TRUE, pt.size = 2,
              label.size = 8)
p2 <- DimPlot(mg, label = TRUE, group.by = 'orig.ident', pt.size = 2,
              label.size = 8)
p3 <- DimPlot(mg, label = TRUE, group.by = 'RNA_snn_res.0.8', pt.size = 2,
              label.size = 8)
cowplot::plot_grid(p1, p2, p3, ncol = 2)


#' Now, clusters are really mixed up. Need to characterize each cluster from 
#' scratch because there is no obvious mapping between log-norm clusters and 
#' these clusters.
mg_sct_de <- FindAllMarkers(
  object = mg,
  only.pos = TRUE,
  verbose = FALSE
)
mg_sct_top <- mg_sct_de %>%
  group_by(cluster) %>%
  filter(p_val_adj < 1e-3) %>%
  top_n(5, wt = avg_logFC)
print(data.frame(mg_sct_top))


#' Some of the DE results suggest that the Activated-MG can be further split 
#' into two distinct clusters. Meanwhile, homeostatic-MG clusters 0, 2, and 3 
#' are still characterized by "Gm-" genes, some "mt-" genes, and other genes for
#' which we cannot assign any functional relevance. For now, we focus on the two
#' activated-MG clusters to determine if they are visually distinct.
#+ fig.width=13, fig.height=9, warnings=FALSE
act_mg_genes <- mg_sct_top %>%
  top_n(n = 2, wt = avg_logFC) %>%
  ungroup() %>%
  select(gene) %>%
  .[['gene']]
DefaultAssay(mg) <- 'RNA'
act_mg_genes_umap <- FeaturePlot(
  object = mg, 
  features = act_mg_genes,
  pt.size = 2,
  combine = FALSE,
  reduction = 'umap',
  order = TRUE
)
act_mg_genes_umap <- lapply(
  X = act_mg_genes_umap,
  FUN = function(x) {
    x + 
      theme_bw() +
      scale_color_gradientn(
        colours = rev(RColorBrewer::brewer.pal(n = 8, name = 'Spectral'))
      )
  }
)
act_mg_genes_umap <- cowplot::plot_grid(plotlist = act_mg_genes_umap, ncol = 4)
act_mg_genes_umap


#' These figures display some gradation of expression of some genes, but the
#' defining marker genes for each cluster are certainly less interpretable than
#' all previous resultly. For example, one cluster only as a single DE gene. 
#' This approach is sub-optimal, discard and work with other approaches. 
#' Conclude that including batch_var does not improve results, but cannot rule
#' out that batch effects are present. At this point, unclear if this can even
#' be diagnosed, because S1d MG are certainly distinct from MG at other times.



#' ## Hierarchical clustering of MG

#' Part of the reason graph-based clustering is commonly implemented in scSeq 
#' studies is computational efficiency for large datasets. Since we have 
#' relatively few cells here, and since we seek to identify distinct clusters, 
#' we can opt to test hierarchical clustering. Previous results have revealed
#' some distinct clusters in MG (e.g. activated-MG), but we also have some that
#' are not obviously distinct (e.g. Unknown-MG ["Gm-" genes]).


#' First step is variable feature selection. We've already noted substantial 
#' differences in selection outputs using log-norm vs SCTransform. 

mg <- readRDS('../data/20200625_microglia.rds')
DefaultAssay(mg) <- 'RNA'

log_dist <- dist(mg[['pca']]@cell.embeddings[,1:10])

mg <- SCTransform(
  object = mg,
  variable.features.n = 2000,
  verbose = FALSE
)
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

