

#' ---
#' title: Hierarchical clustering of SVZ
#' author: James Choi
#' date: 2021/1/4
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)

#' ## Background  
#' 
#' Purpose of this analysis is to revisit the cluster analysis of the NSCs. As
#' Suvra and I discussed early December 2020, a trajectory analysis of the NSCs
#' would be informative to study NSC fate after stroke. In particular, we would
#' like to know: does stroke affect NSC differentiation dynamics? To answer this
#' question with a trajectory analysis, we need to a) normalize counts, b) find
#' variable genes, c) embed cells in dimensionally-reduced space (ie PCA). In 
#' some previous studies (https://doi.org/10.1016/j.cels.2016.08.011), and even
#' in the slingshot R package tutorial, trajectories are inferred from this 
#' reduced space without having to use specialized tools like Monocle3.  


# Libraries
# install.packages('dendsort')
require('Seurat')
require('ggplot2')
require('dplyr')
require('dendextend')



#' ### Exploratory PCA 
#' 
nsc <- readRDS(file = '../data/20200625_nsc.rds')

#+ original_umap, fig.height=4, fig.width=5.5
DimPlot(nsc, pt.size = 3, label.size = 6, label = TRUE) + theme_bw()


#' First, we examine the previous cluster results. Of note, we regressed 
#' "CC.difference" in the ScaleData() step.  
DefaultAssay(nsc) <- 'RNA'

head(VariableFeatures(nsc), 50)
tail(VariableFeatures(nsc), 50)

summary(nsc$nFeature_RNA)


#+ original_pca, fig.height=5, fig.width=12
pca_plots <- vector(mode = 'list', length = 10)
for (i in 1:10) {
  pca_plots[[i]] <- PCAPlot(
    object = nsc,
    pt.size = 1,
    dims = c(i, 1)
  ) +
    NoLegend()
}
pca_plots_log <- cowplot::plot_grid(plotlist = pca_plots, ncol = 5)
pca_plots_log


#' The PCA plots show us that PC1, which contains most variance, is dominated by
#' ~10 cells. PC3 may represent spectrum of NSC differentiation. PC4 and 5 have
#' low spread for cells low in PC1. PC6 and 7 look very similar, and PC9 looks 
#' like a mirror image, which indicates they are measuring similar properties.  

#+ log_dim_loadings, fig.height=12, fig.width=14
VizDimLoadings(nsc, dims = 1:10, balanced = TRUE, ncol = 5, nfeatures = 40)


#' Based on gene weightings, PC1 may determine ependymal cells (contamination?).
#' PC3 has cell-cycle at + end and Dcx at - end. PC4 and 5 have familiar genes
#' at + end, but unclear what it means - may be contamination (e.g. Cebpb is a
#' myeloid gene marker).  
#' 


#' ### Exploratory PCA & clustering with alternative count normalization  
#' 
#' As an alternative, we test how results differ with SCTransform, which 
#' performs relative normalization of the counts. We select an equal number of
#' features.  
#' 

DefaultAssay(nsc) <- 'RNA'
nsc <- SCTransform(nsc, 
                   variable.features.n = 2000,
                   vars.to.regress = 'CC.difference',
                   verbose = FALSE)

head(VariableFeatures(nsc), 50)
tail(VariableFeatures(nsc), 50)


#+ pca_plots, fig.height=5, fig.width=12
nsc <- RunPCA(nsc, npcs = 30, verbose = FALSE)
pca_plots <- vector(mode = 'list', length = 10)
for (i in 1:10) {
  pca_plots[[i]] <- PCAPlot(
    object = nsc,
    pt.size = 1,
    dims = c(i, 1)
  ) +
    NoLegend()
}
pca_plots_sct <- cowplot::plot_grid(plotlist = pca_plots, ncol = 5)
pca_plots_sct

#' PC1 and 2 together seem to capture variation along the NSC differentiation 
#' trajectory (based on our previous discussions on the putative identities of
#' clusters 0-4). PCs 4-10 pretty much look similar... with what appears to be 
#' minor differences in the spread of cells in clusters 0, 2, 3 along the axis. 
#'  


#' The conclusion here is that normalization method affects downstream PCA. As
#' we have already seen, the first few components capture the most variance and
#' PC's beyond 5 are almost identical. We can visualize contribution of variance
#' with an elbow plot:  

#+ elbow_sct, fig.height=2.5, fig.width=4
ElbowPlot(nsc, ndims = 30) + theme_bw()

#+ sct_dimloadings, fig.height=12, fig.width=14
VizDimLoadings(nsc, dims = 1:10, ncol = 5, nfeatures = 40)


#' How do cluster results compare when using SCTransform?
#+ sct_umap, fig.height=4, fig.width=5.5
nsc <- FindNeighbors(nsc, dims = 1:10, verbose = FALSE)
nsc <- FindClusters(nsc, verbose = FALSE)
nsc <- RunUMAP(nsc, dims = 1:10, verbose = FALSE)
DimPlot(nsc, label = TRUE, label.size = 6, pt.size = 2) + theme_bw()

#' Based on UMAP, there are two larger groups of cells: 1) a larger group of 
#' cells that appears to form a continuum, and 2) a smaller group of cells. 
#' Using a graph-based clustering approach, at default Seurat resolution, 
#' the small group of clusters togethers with some cells of the larger group and
#' is not part of the continuum. Assuming the UMAP is accurate (since even some 
#' trajectory analyses e.g. Monocle3 and PAGA depend on UMAP embeddings), can 
#' we achieve separation of cluster 1 into two smaller clusters?  

#+ higher_res, fig.height=4, fig.width=5
nsc <- FindClusters(nsc, resolution = 2, verbose = FALSE)
DimPlot(nsc, label = TRUE, label.size = 6, pt.size = 2) + theme_bw()

#' Even at high resolutions, the smaller group of cells continues to cluster 
#' with a portion of the larger group of cells (namely, the "bridge" segment)
#' i.e. the UMAP separation and cluster results do not align very well.  
#' 


#' ### Hierarchical clustering  
#' 
#' How do cluster results differ if we use a different clustering approach e.g.
#' hierarchical clustering? First, we test how graph-based clustering results
#' compare to hierarchical clustering results using log-normalized data (as 
#' done in initial analysis). We use the cell coordinates in PCA-space (same
#' number of dimensions used for graph-based) and perform hierarchical 
#' clustering using ward.D2 method, which aims to minimize within-cluster
#' variance during agglomerization.  
#' 

nsc <- readRDS(file = '../data/20200625_nsc.rds')

# Note here that PC's used in previous analysis were determined using JackStraw
# method. Jackstraw method is currently not compatible with SCTransform, which 
# forces uses to use elbow-heuristic.  
# 
log_pcDims <- nsc@commands$FindNeighbors.RNA.pca$dims
log_pcDims

nsc_dist <- dist(x = nsc[['pca']]@cell.embeddings[,log_pcDims])
nsc_tree <- hclust(d = nsc_dist, method = 'ward.D2')

nsc_tree$labels <- seq_along(nsc_tree$labels)
nsc_dend <- dendsort::dendsort(d = as.dendrogram(nsc_tree, hang = 1))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
tmp_cols <- gg_color_hue(4)
names(tmp_cols) <- factor(0:3)
labels_colors(nsc_dend) <- 
  tmp_cols[nsc$RNA_snn_res.0.8][order.dendrogram(nsc_dend)]

#+ hierarchical_dend, fig.height=5, fig.width=7.5
plot(nsc_dend,
     main = 'Hierarchical clustering of NSCs (log-normalized)',
     sub = 'Labels colored by original cluster results.')


#' The dendrogram shows some differences in classification compared to graph-
#' based clustering. Most notably, cluster 1 in the graph-based approach (green)
#' is shared across two larger branches of the dendrogram, which spli very high
#' in the three. We also notice that there is some overlap between cluster 3
#' (purple) and cluster cluster 1 (green), and some overlap between cluster 3 (
#' purple) and cluster 2 (blue). These overlaps are not surprising.  
#' 


#' We have shown that hierarchical clustering (in PCA-space) can provide more
#' resolved cluster results. But as mentioned before, normalization method (and
#' downstream feature selection) can influence cluster results as well. One 
#' particular observation from the previous PCA plots (see: Exploratory PCA) is
#' that there is overdispersion towards one end of the PC pair plots.  
#' 

#' Examine PC5 and focus on Nr4a1, Rhoj:
#+ pc5, fig.height=6, fig.width=7.5
FeaturePlot(nsc,
            features = c('log10_total_counts','nFeature_RNA','Nr4a1', 'Rhoj'),
            dims = c(5,1),
            cols = c('grey','red3'),
            reduction = 'pca',
            pt.size = 2,
            order = TRUE)
nsc$blank <- 'blank'
#+ pc5_gene_vln, fig.height=3.5, fig.width=5.5
VlnPlot(nsc, 
        features = c('Nr4a1','Rhoj'),
        assay = 'RNA',
        slot = 'counts',
        group.by = 'blank',
        pt.size = 1)


#' Two of the top (+)-weighted genes for PC5, Nr4a1 and Rhoj, are dominated by
#' unusually high counts in a single cell or very low overall abundance.  
#' 


#' ### Normalization and Feature Selection
#' 
#' Generally, having genes with large weights due to a few number of cells is
#' not undesirable, since these can be rare cell populations. However, given the
#' low number of sequenced cells to begin with, we want relatively highly 
#' expressed genes to ensure that our cluster results are robust to 1-3 cells 
#' driving variation. From the results so far, we see that there is room for 
#' improvement in normalization and/or feature selection.   
#' 

nsc <- readRDS(file = '../data/20200625_nsc.rds')

#' Using the default Seurat pipeline method, we select the top 2000 variable 
#' features. See ?FindVariableFeatures for more details on exact method.  
#' 
#+ default_selection, fig.height=3.5, fig.width=5
DefaultAssay(nsc) <- 'RNA'
nsc <- NormalizeData(nsc, verbose = FALSE)
nsc <- FindVariableFeatures(nsc, verbose = FALSE)
meta_features <- nsc[['RNA']]@meta.features
min(meta_features$vst.mean[meta_features$vst.variable])
min(meta_features$vst.variance.standardized[meta_features$vst.variable])
VariableFeaturePlot(nsc) + theme_bw()

#' To remove genes with low average expression (i.e. genes with low counts
#' across some cells or with high counts in very few cells), we filter based on
#' the values from the fit model.  
#' 

#+ manual_selection_1, fig.height=3.5, fig.width=5
# 1.24493 is the minimum gene variance under default method + params.
is_var <- with(data = nsc[['RNA']]@meta.features,
               expr = vst.mean > 0.1 & vst.variance.standardized > 1.24493)
var_genes <- rownames(nsc[['RNA']]@meta.features)[is_var]
nsc[['RNA']]@meta.features$vst.variable <- is_var
VariableFeatures(nsc) <- var_genes 
VariableFeaturePlot(nsc)

#' This results in 1135 variable features selected. Proceed to downstream.  
#' 

#+ pca_dims, fig.height=3.5, fig.width=4.5
nsc <- ScaleData(nsc, vars.to.regress = 'CC.difference', verbose = FALSE)
nsc <- RunPCA(nsc, verbose = FALSE)
ElbowPlot(nsc, ndims = 30) + theme_bw()

#+ manual_variable_clusters, fig.height=4, fig.width=5.5
# For now, we do not concern with determining # of PCs because we want to
# isolate effect of changing feature selection method on downstream.
# nsc <- JackStraw(nsc, dims = 20, verbose = FALSE)
# nsc <- ScoreJackStraw(nsc, dims = 1:20)
# sig_pcs <- which(nsc[['pca']]@jackstraw$overall.p.values[,2] < 1e-03)
# JackStrawPlot(nsc, dims = 1:20)
nsc <- FindNeighbors(nsc, dims = 1:5, verbose = FALSE)
nsc <- FindClusters(nsc, verbose = FALSE)
nsc <- RunUMAP(nsc, dims = 1:5, verbose = FALSE)
DimPlot(nsc, pt.size = 2, label = TRUE, label.size = 6) + theme_bw()

#' We again reproduce results similar to our original result under default 
#' feature selection method. What happens if we lower the variance threshold in
#' order to capture more features with the same minimum average expression?  
#' 

#+ manual_selection_2, fig.height=3.5, fig.width=5
nsc <- FindVariableFeatures(nsc, verbose = FALSE)
is_var <- with(data = nsc[['RNA']]@meta.features,
               expr = vst.mean > 0.1 & vst.variance.standardized > 1) # set to 1
var_genes <- rownames(nsc[['RNA']]@meta.features)[is_var]
nsc[['RNA']]@meta.features$vst.variable <- is_var
VariableFeatures(nsc) <- var_genes 
VariableFeaturePlot(nsc)

#' With a minimum standardized variance of 1, we now have 2944 variable genes.
#' Proceed with downstream.  
#' 

#+ lower_gene_var_elbow, fig.height=3.5, fig.width=4.5
nsc <- ScaleData(nsc, vars.to.regress = 'CC.difference', verbose = FALSE)
dim(nsc[['RNA']]@scale.data)
nsc <- RunPCA(nsc, npcs = 30, verbose = FALSE)
ElbowPlot(nsc, ndims = 30) + theme_bw()

#+ lower_gene_var_clusters, fig.height=4, fig.width=5.5
nsc <- FindNeighbors(nsc, dims = 1:10, verbose = FALSE)
nsc <- FindClusters(nsc, verbose = FALSE)
nsc <- RunUMAP(nsc, dims = 1:10, verbose = FALSE)
DimPlot(nsc, pt.size = 2, label = TRUE, label.size = 6) + theme_bw()


#' For completion, we try removing cc.genes from the variable features to test
#' how it affects downstream results compared to setting "vars.to.regress = 
#' 'CC.difference'".  
#' 

#+ manual_selection_cc, fig.height=3.5, fig.width=5
nsc <- FindVariableFeatures(nsc, verbose = FALSE)
is_var <- with(data = nsc[['RNA']]@meta.features,
               expr = vst.mean > 0.1 & vst.variance.standardized > 1) # set to 1
var_genes <- rownames(nsc[['RNA']]@meta.features)[is_var]
firstup <- function(x) {
  x <- tolower(x)
  substr(x,1,1) <- toupper(substr(x,1,1))
  return(x)
}
var_genes <- var_genes[!var_genes %in% firstup(unlist(cc.genes.updated.2019))]
is_var <- rownames(nsc[['RNA']]@meta.features) %in% var_genes
nsc[['RNA']]@meta.features$vst.variable <- is_var
VariableFeatures(nsc) <- var_genes 
VariableFeaturePlot(nsc)

nsc <- ScaleData(nsc, verbose = FALSE)
dim(nsc[['RNA']]@scale.data)
nsc <- RunPCA(nsc, npcs = 30, verbose = FALSE)

#+ lower_gene_var_cc_clusters, fig.height=4, fig.width=5.5
nsc <- FindNeighbors(nsc, dims = 1:10, verbose = FALSE)
nsc <- FindClusters(nsc, verbose = FALSE)
nsc <- RunUMAP(nsc, dims = 1:10, verbose = FALSE)
DimPlot(nsc, pt.size = 2, label = TRUE, label.size = 6) + theme_bw()

#' Removing cc.genes from variable features manually does not drastically change
#' results compared to setting vars.to.regress argument. If anything, there is
#' even greater separation between clusters of cells.  
#' 


#' As we demonstrated earlier in this analysis during the exploratory PCA, 
#' SCTransform produces a continuum of UMAP cell clusters that recapitulates the
#' differentiation trajectory we are interested in. First, let's examine the 
#' mean-variance relationship (following script from here: 
#' https://rawgit.com/ChristophH/sctransform/supp_html/supplement/variance_stabilizing_transformation.html).  
#' 
raw_count <- nsc[['RNA']]@counts
gene_attr <- data.frame(
  mean = Matrix::rowMeans(raw_count),
  detection_rate = Matrix::rowMeans(raw_count > 0),
  var = apply(X = raw_count, MARGIN = 1, FUN = var)
)
gene_attr$log_mean <- log10(gene_attr$mean)
gene_attr$log_var <- log10(gene_attr$var)
rownames(gene_attr) <- rownames(raw_count)
cell_attr <- data.frame(
  n_umi = Matrix::colSums(raw_count),
  n_gene = Matrix::colSums(raw_count > 0)
)
rownames(cell_attr) <- colnames(raw_count)

#' The mean-variance relationship has properties similar to those described in
#' the script linked above.  
#' 

#+ mean_var, fig.height=3.5, fig.width=4
ggplot(gene_attr, aes(log_mean, log_var)) +
  geom_point(alpha = 0.3, shape = 16) + 
  geom_density_2d(size = 0.3) + 
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw()

#' Mean-detection rate relationship also holds.  
#' 

#+ mean_detectionRate, fig.height=3.5, fig.width=4
x = seq(from = -3, to = 2, length.out = 1000)
poisson_model <- data.frame(log_mean = x, detection_rate = 1 - dpois(0, lambda = 10^x))
ggplot(gene_attr, aes(log_mean, detection_rate)) + 
  geom_point(alpha = 0.3, shape = 16) + 
  geom_line(data = poisson_model, color = "red") + 
  theme_bw()

#' Based on this, SCTransform normalization can be used. Whether it is the best
#' method, I am still unsure. Below, we repeat SCTransform and cluster analysis
#' as above, but also explore how selecting different # of PCs affects clusters.
#'   
#' 

#+ sct_varfeature, fig.height=3.5, fig.width=5
DefaultAssay(nsc) <- 'RNA'
nsc <- NormalizeData(nsc, verbose = FALSE)
nsc <- SCTransform(nsc, 
                   variable.features.n = 2000, 
                   vars.to.regress = 'CC.difference',
                   verbose = FALSE)
VariableFeaturePlot(nsc) + theme_bw()

#+ sct_elbow, fig.height=3, fig.width=4
nsc <- RunPCA(nsc, npcs = 50, verbose = FALSE)
ElbowPlot(nsc, ndims = 50)

#+ sct_loadings, fig.height=14, fig.width=14
VizDimLoadings(nsc, dims = 1:15, ncol = 5)

#+ sct_umap_2, fig.height=4, fig.width=5.5
nsc <- FindNeighbors(nsc, dims = 1:10, verbose = FALSE)
nsc <- FindClusters(nsc, verbose = FALSE)
nsc <- RunUMAP(nsc, dims = 1:10, verbose = FALSE)
DimPlot(nsc, pt.size = 2, label = TRUE, label.size = 6) + theme_bw()


#' How do results compare with different numbers of PCs?  
#' 

#+ pc_umap_iterate, fig.height=6, fig.width=18.5
test_pcs <- 10:19
umaps <- vector(mode = 'list', length = length(test_pcs))
for(i in 1:length(test_pcs)) {
  nsc <- FindNeighbors(nsc, dims = 1:test_pcs[i], verbose = FALSE)
  nsc <- FindClusters(nsc, verbose = FALSE)
  nsc <- RunUMAP(nsc, dims = 1:test_pcs[i], verbose = FALSE)
  umaps[[i]] <- 
    DimPlot(nsc, pt.size = 1, label = TRUE, label.size = 6) + 
    theme_bw() +
    labs(title = paste('# PCs:', test_pcs[i]))
}
umaps <- cowplot::plot_grid(plotlist = umaps, ncol = length(test_pcs)/2)
umaps


#' ### Checkpoint conclusions
#' 
#' Based on the results so far, SCTransform appears to better recapitulate the
#' continuous trajectory of differentiation along B cell --> neuroblast path. 
#' This suggests that normalization/feature selection should be revisited. There
#' is also a group of cells that, despite appearing distinct in UMAPs, is not 
#' separable through graph-based clustering.  
#' 
#' 
#' It is also possible that microglia vs NSC cell-type classification is also 
#' not entirely accurate since it appears graph-based clustering cannot resolve 
#' certain groups of cells.  


