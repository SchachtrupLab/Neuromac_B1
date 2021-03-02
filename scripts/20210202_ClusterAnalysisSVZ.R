
#' ---
#' title: "Cluster Analysis of SVZ"
#' author: "James Choi"
#' date: "`r Sys.Date()`"
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)

if (!grepl('scripts', getwd())) {
  setwd('./scripts/')
}
results_out <- '../results/ClusterAnalysisSVZ/'
dir.create(path = results_out)

require('Seurat')
require('ggplot2')
require('dplyr')
require('dendextend')


# This Seurat object contains all cells passed quality control.  
svz <- readRDS(file = '../data/20210202_SVZ.rds')


#+ variable_features, fig.height=2.5, fig.width=4.5
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}
DefaultAssay(svz) <- 'RNA'
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
VariableFeaturePlot(svz) + theme_bw()


#+ pca_elbow, fig.height=2.5, fig.width=4
svz <- RunPCA(svz, npcs = 40, verbose = FALSE)
ElbowPlot(svz, ndims = 40) + theme_bw()


#+ pca_iteration, fig.height=6, fig.width=17, fig.cap="UMAPs using varying number of PCs"
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


#+ dimloadings, fig.height=4.5, fig.width=8
VizDimLoadings(svz, dims = 9:12, nfeatures = 30, ncol = 4)


#' We also inspect the distrbution of cells along thesse components:  
#+ pca_plots, fig.height=3.5, fig.width=12
p1 <- DimPlot(svz, pt.size = 2, reduction = 'pca', dims = c(10,11))
p2 <- DimPlot(svz, pt.size = 2, reduction = 'pca', dims = c(10,12))
p3 <- DimPlot(svz, pt.size = 2, reduction = 'pca', dims = c(11,12))
p1+p2+p3

#' Seems that up to PC 11 has good distribution of cells, whereas PC 12 has a
#' weird skewed distribution.  


#+ svz_umap, fig.height=3.5, fig.width=4.5
svz <- FindNeighbors(svz, dims = 1:11, verbose = FALSE)
svz <- RunUMAP(svz, dims = 1:11, verbose = FALSE)
svz <- FindClusters(svz, resolution = 2, verbose = FALSE)
Idents(svz) <- 'SCT_snn_res.2'
p1 <- DimPlot(svz, pt.size = 2, label = TRUE, label.size = 6) + theme_bw()
p1 <- FetchData(svz, vars = c('SCT_snn_res.2', 'UMAP_1', 'UMAP_2')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(mapping = aes(fill = SCT_snn_res.2),
             color = 'black', pch = 21, size = 3) +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave(filename = paste0(results_out, 'SVZ_res-2-clusters_umap.svg'),
       plot = p1, device = 'svg', height = 3.5, width = 4.5)
p1


#' Identify marker genes for each cluster to help classify as MG or NSC.  
Idents(svz) <- 'SCT_snn_res.2'
gene_markers <- FindAllMarkers(svz, assay = 'RNA', only.pos = TRUE)
top_markers <- gene_markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = -p_val_adj) %>%
  ungroup()
knitr::kable(x = top_markers, caption = 'Top 3 DE genes per cluster.')


#' Generate a heatmap of the top markers to better identify overlapping and 
#' distinct clusters.  
#+ marker_heatmap, fig.height=7.75, fig.width=5.5, fig.cap='Heatmap of top 3 markers per cluster'
top_markers <- gene_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = -p_val_adj)
top_markers <- unique(top_markers$gene)

Idents(svz) <- 'SCT_snn_res.2'
DefaultAssay(svz) <- 'SCT'
my_cols <- rev(RColorBrewer::brewer.pal(n = 11, name = 'RdBu'))
avg_exp <- ScaleData(object = svz[['SCT']]@data, features = top_markers)
marker_heatmap <- cbind(data.frame(t(avg_exp)), 'ident' = svz@active.ident) %>%
  reshape2::melt(id.vars = 'ident') %>%
  group_by(ident, variable) %>%
  summarise(avg_exp = mean(value)) %>%
  mutate(variable = factor(variable, levels = rev(levels(variable)))) %>%
  ggplot(mapping = aes(x = ident, y = variable, fill = avg_exp)) +
  geom_tile(color = 'black') +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient2(low = "#053061",
                       high = "#B2182B",
                       limits = c(NA, 3),
                       na.value ="#B2182B") +
  theme(axis.text.x = element_text(angle = 0),
        axis.title = element_blank()) +
  guides(fill = guide_colorbar(
    title = 'Expression\nz-score',
    barwidth = 1.25, 
    frame.colour = 'black', 
    frame.linewidth = 0.5, 
    ticks.colour = 'black',
    ticks.linewidth = 0.5))
ggsave(filename = paste0(results_out, 'SVZ_res-2-clusters_marker_heatmap.svg'),
       plot = marker_heatmap, device = 'svg', height = 7.75, width = 5.5)
marker_heatmap



#' Generate violin plots of common microglia and myeloid-associated genes. Gene
#' products of some of these genes are commonly targeted epitopes for antibodies
#' (e.g. Adgre == F4/80, Ptprc = CD11b, Mrc1 = CD206).  
#+ myeloid_marker_violin, fig.height=4.5, fig.width=11, fig.cap='Violin plots of various myeloid marker genes pulled from various previously published studies (contact JSC for exact references).'
myeloid_genes <- c('P2ry12','Tmem119','Cx3cr1','Csf1r','Adgre1','Ptprc','Ly6c2',
                   'Ccr2','Mrc1','Lyve1','Mki67','Top2a')
svz$seurat_clusters <- svz$SCT_snn_res.2
Idents(svz) <- 'seurat_clusters'
DefaultAssay(svz) <- 'RNA'
expr_dat <- FetchData(
  object = svz,
  vars = c('seurat_clusters', myeloid_genes),
  slot = 'data'
)
myeloid_vln <- expr_dat %>%
  reshape2::melt(id.vars = 'seurat_clusters') %>%
  ggplot(mapping = aes(x = seurat_clusters, y = value)) +
  geom_violin(mapping = aes(fill = seurat_clusters), scale = 'width') +
  facet_wrap(variable ~ ., scales = 'free_y') +
  xlab(label = 'Cluster') + 
  ylab(label = 'log-normalized expression') + 
  theme_bw() + 
  theme(legend.position = 'none',
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(size = 12, hjust = 0))
ggsave(
  filename = paste0(results_out, 'SVZ_res-2-clusters_myeloid-genes_vln.svg'),
  plot = myeloid_vln,
  device = 'svg',
  height = 4.5, width = 11
)
myeloid_vln

#' Plot UMAPs of myeloid gene expression to better visualize spread of detection
#' over all cells.  
#+ myeloid_marker_umap, fig.height=8, fig.width=13.5, fig.cap='Myeloid gene expression overlaid on SVZ UMAPs.'
myeloid_genes <- c('P2ry12','Tmem119','Cx3cr1','Csf1r','Adgre1','Ptprc','Ly6c2',
                   'Ccr2','Mrc1','Lyve1','Mki67','Top2a')
svz$seurat_clusters <- svz$SCT_snn_res.2
Idents(svz) <- 'seurat_clusters'
DefaultAssay(svz) <- 'RNA'
myeloid_umap <- FeaturePlot(
  object = svz,
  features = myeloid_genes,
  cols = c('grey','red3'),
  slot = 'data',
  combine = FALSE
)
myeloid_umap <- lapply(
  X = myeloid_umap,
  FUN = function(x) {
    x + theme_bw()
  }
)
myeloid_umap <- cowplot::plot_grid(plotlist = myeloid_umap)
ggsave(filename = paste0(results_out, 'SVZ_myeloid-genes_umap.svg'),
       plot = myeloid_umap, device = 'svg', height = 8, width = 13.5)
myeloid_umap



#' Generate violin plots of previously described NSPC genes from this 
#' [paper](https://www.ncbi.nlm.nih.gov/pubmed/30777863).  
#+ nsc_marker_violin, fig.height=5, fig.width=11, fig.cap='Violin plots of various NSPC genes associated with quiescence, activation, priming.'
nsc_genes <- c('Gfap','Sox2','Dcx','Ascl1','Sox9','Id3','Id4','Foxm1','Egfr',
               'Mki67','Top2a','Dlx1','Dlx2','Fos')
svz$seurat_clusters <- svz$SCT_snn_res.2
Idents(svz) <- 'seurat_clusters'
DefaultAssay(svz) <- 'RNA'
expr_dat <- FetchData(
  object = svz,
  vars = c('seurat_clusters', nsc_genes),
  slot = 'data'
)
nsc_vln <- expr_dat %>%
  reshape2::melt(id.vars = 'seurat_clusters') %>%
  ggplot(mapping = aes(x = seurat_clusters, y = value)) +
  geom_violin(mapping = aes(fill = seurat_clusters), scale = 'width') +
  facet_wrap(variable ~ ., scales = 'free_y') +
  xlab(label = 'Cluster') + 
  ylab(label = 'log-normalized expression') + 
  theme_bw() + 
  theme(legend.position = 'none',
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(size = 12, hjust = 0))
ggsave(
  filename = paste0(results_out, 'SVZ_res-2-clusters_nsc-genes_vln.svg'),
  plot = nsc_vln,
  device = 'svg',
  height = 5, width = 11
)
nsc_vln

#' Plot UMAPs of NSC gene expression to better visualize spread of detection
#' over all cells.  
#+ myeloid_marker_umap, fig.height=11, fig.width=13.5, fig.cap='NSC gene expression overlaid on SVZ UMAPs.'
nsc_genes <- c('Gfap','Sox2','Dcx','Ascl1','Sox9','Id3','Id4','Foxm1','Egfr',
               'Mki67','Top2a','Dlx1','Dlx2','Fos')
svz$seurat_clusters <- svz$SCT_snn_res.2
Idents(svz) <- 'seurat_clusters'
DefaultAssay(svz) <- 'RNA'
nsc_umap <- FeaturePlot(
  object = svz,
  features = nsc_genes,
  cols = c('grey','red3'),
  slot = 'data',
  combine = FALSE
)
nsc_umap <- lapply(
  X = nsc_umap,
  FUN = function(x) {
    x + theme_bw()
  }
)
nsc_umap <- cowplot::plot_grid(plotlist = nsc_umap)
ggsave(filename = paste0(results_out, 'SVZ_nsc-genes_umap.svg'),
       plot = nsc_umap, device = 'svg', height = 11, width = 13.5)
nsc_umap



#' We can also compare the data with a previously published study on NSCs vs 
#' ependymal cells in V-SVZ by [Shah, Stratton, et al.](https://doi.org/10.1016/j.cell.2018.03.063).
#' Specifically, we refer to Figure 3. Acta2 and Foxj1 are unique to ependymal
#' cells. Interestingly, they also note a small Flt1(+) population that we also
#' identified in the initial clustering of all SVZ cells.

#+ ependymal_genes, fig.height=4.75, fig.width=9
ependymal_genes <- c('Foxj1','Cfap126','Ccdc153', 'Rarres2','Tmem212')
ependymal_umaps <- FeaturePlot(
  object = svz, 
  features = ependymal_genes, 
  cols = c('grey','red3'),
  combine = FALSE
)
ependymal_umaps <- lapply(
  X = ependymal_umaps,
  FUN = function(x) {
    x <- x + theme_bw()
    return(x)
  }
)
ependymal_umaps <- cowplot::plot_grid(plotlist = ependymal_umaps, ncol = 3)
ggsave(filename = paste0(results_out, 'SVZ_ependymal-gene_umap.svg'),
       plot = ependymal_umaps, device = 'svg', height = 4.75, width = 9)
ependymal_umaps


#+ svz_celltype_umap, fig.height=3, fig.width=5.25
svz$celltype <- plyr::mapvalues(
  x = svz$SCT_snn_res.2,
  from = 0:12,
  to = c('NSC','Microglia','Microglia','Microglia','Microglia','NSC',
         'Microglia','NSC','NSC','Microglia','Microglia','Microglia', 
         'Unknown/Ependymal')
)
Idents(svz) <- 'celltype'
p1 <- DimPlot(svz, pt.size = 2, label = TRUE, label.size = 4) + theme_bw()
ggsave(filename = paste0(results_out, 'SVZ_celltype-umap.svg'),
       plot = p1, height = 3, width = 5, device = 'svg')
p1


saveRDS(svz, file = '../data/20210202_SVZ.rds')


sessionInfo()