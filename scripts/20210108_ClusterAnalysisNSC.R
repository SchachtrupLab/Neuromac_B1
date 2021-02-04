
#' ---
#' title: Cluster Analysis of NSC
#' author: James Choi
#' date: 2021/1/8
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)

require('Seurat')
require('ggplot2')
require('dplyr')
require('dendextend')
# DO NOT set.seed() !!

results_out <- '../results/ClusterAnalysisNSC/'
dir.create(path = results_out)
svz <- readRDS(file = '../data/20210202_SVZ.rds')

#' ## Cluster Analysis of NSCs
# subset nsc
nsc <- svz[,svz$celltype == 'NSC']
DefaultAssay(nsc) <- 'RNA'

#+ variable_features, fig.height=3, fig.width=5.25
nsc <- NormalizeData(nsc, verbose = FALSE)
summary(nsc$nFeature_RNA)
nsc <- SCTransform(nsc,
                   assay = 'RNA',
                   variable.features.n = 2000,
                   vars.to.regress = 'CC.difference',
                   verbose = FALSE)
VariableFeaturePlot(nsc) + theme_bw()

#+ elbow, fig.height=2.5, fig.width=3
ElbowPlot(nsc, ndims = 50) + theme_bw()

#+ iterate_pca, fig.height=6, fig.width=17.5
nsc <- RunPCA(nsc, npcs = 50, verbose = FALSE)
test_pcs <- 6:15
umaps <- vector(mode = "list", length = length(test_pcs))
for (i in 1:length(test_pcs)) {
  nsc <- FindNeighbors(nsc, dims = 1:test_pcs[i], verbose = FALSE)
  nsc <- FindClusters(nsc, verbose = FALSE)
  nsc <- RunUMAP(nsc, dims = 1:test_pcs[i], n.neighbors = 30, verbose = FALSE)
  umaps[[i]] <- DimPlot(nsc, pt.size = 1, label = TRUE, label.size = 6) +
    theme_bw() + labs(title = paste("# PCs:", test_pcs[i]))
}
umaps <- cowplot::plot_grid(plotlist = umaps, ncol = length(test_pcs)/2)
umaps

#+ dimloadings_sct, fig.height=8, fig.width=14
VizDimLoadings(nsc, dims = 6:15, nfeatures = 20, ncol = 5)

#+ pca_plots, fig.height=2.75, fig.width=9
p1 <- DimPlot(nsc, pt.size = 2, reduction = 'pca', dims = c(8,11)) + theme_bw()
p2 <- DimPlot(nsc, pt.size = 2, reduction = 'pca', dims = c(8,13)) + theme_bw()
p3 <- DimPlot(nsc, pt.size = 2, reduction = 'pca', dims = c(11,13)) + theme_bw()
pca_plots <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
pca_plots


#+ nsc_umap, fig.height=3.5, fig.width=4.25
nsc <- FindNeighbors(nsc, dims = 1:8, verbose = FALSE)
nsc <- RunUMAP(nsc, dims = 1:8, verbose = FALSE)
nsc <- FindClusters(nsc, verbose = FALSE)
DimPlot(nsc, pt.size = 2, label = TRUE, label.size = 6) + theme_bw()


#+ hierarchical_dend, fig.height=5, fig.width=7.5
nsc_dist <- dist(x = nsc[['pca']]@cell.embeddings[,1:8])
nsc_tree <- hclust(d = nsc_dist, method = 'ward.D2')
nsc_dend <- dendsort::dendsort(d = as.dendrogram(nsc_tree, hang = 1))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
tmp_cols <- gg_color_hue(length(levels(nsc$SCT_snn_res.0.8)))
names(tmp_cols) <- levels(nsc$SCT_snn_res.0.8)
labels_colors(nsc_dend) <- 
  tmp_cols[nsc$SCT_snn_res.0.8][order.dendrogram(nsc_dend)]
plot(nsc_dend,
     main = 'Hierarchical clustering of NSCs',
     sub = 'Labels colored by original cluster results.')
abline(h = 80)

#+ cluster_comparison, fig.height=3.5, fig.width=8
nsc_dclus <- cutree(tree = nsc_dend, h = 80)
nsc$hier_clust <- paste('H', nsc_dclus[match(rownames(nsc@meta.data),
                                             names(nsc_dclus))],
                        sep = '_')
nsc$hier_clust <- factor(nsc$hier_clust)

p1 <- DimPlot(nsc, 
              group.by = 'hier_clust', 
              pt.size = 1.5, 
              label = TRUE, 
              label.size = 5) +
  theme_bw() +
  labs(title = 'Hierarchical cluster results\noverlaid on UMAP')
p2 <- DimPlot(nsc, 
              group.by = 'SCT_snn_res.0.8', 
              pt.size = 1.5, 
              label = TRUE, 
              label.size = 5) +
  theme_bw() +
  labs(title = 'SNN-graph cluster results\noverlaid on UMAP')
cowplot::plot_grid(p1, p2, ncol = 2)

#' Inspect distribution of cells/clusters over injury time:  
#' 
#+ umap_over_time, fig.height=2.75, fig.width=8.5
Idents(nsc) <- 'SCT_snn_res.0.8'
DimPlot(nsc, pt.size = 2, split.by = 'orig.ident') + 
  theme_bw() +
  theme(strip.text = element_text(size = 12))


#' Identify marker genes:  
#' 

Idents(nsc) <- 'SCT_snn_res.0.8'
markers <- FindAllMarkers(
  object = nsc,
  assay = 'RNA',
  only.pos = TRUE
)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = -p_val_adj)
knitr::kable(x = top_markers)  


# Convert gene-names to Ensembl IDs
ensembl <- useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', 
                      dataset = 'mmusculus_gene_ensembl')
conversion_input <- rownames(nsc[['RNA']]@counts)
conversion_input <- plyr::mapvalues(conversion_input, 'Sepp1', 'Selenop')
markers$gene <- plyr::mapvalues(markers$gene, 'Sepp1', 'Selenop')
ensembl_conversion <- getBM(
  attributes = c('mgi_symbol','ensembl_gene_id','external_gene_name'),
  filters = 'external_gene_name',
  values = conversion_input,
  mart = ensembl
)
markers$ensembl_id <- plyr::mapvalues(
  x = markers$gene,
  from = ensembl_conversion$mgi_symbol,
  to = ensembl_conversion$ensembl_gene_id,
  warn_missing = FALSE
)
table(grepl('ENS', markers$ensembl_id))
write.table(markers, file = paste0(results_out, 'NSC_subcluster_markers.csv'),
            sep = ',', row.names = FALSE)

# p1 <- DimPlot(nsc, pt.size = 2, label = TRUE, label.size = 6) + theme_bw() +
#   NoLegend()
# ggsave(filename = '../results/NSC_umap_byCluster.svg', plot = p1, device = 'svg',
#        height = 3, width = 3.25)
# p1 <- DimPlot(nsc, pt.size = 2, group.by = 'orig.ident') + theme_bw() +
#   NoLegend()
# ggsave(filename = '../results/NSC_umap_byTime.svg', plot = p1, device = 'svg',
#        height = 3, width = 3.25)


saveRDS(nsc, file = '../data/20210108_NSC.rds')


# nsc <- FindNeighbors(nsc, dims = 1:11, verbose = FALSE)
# nsc <- RunUMAP(nsc, dims = 1:11, verbose = FALSE)
# nsc <- FindClusters(nsc, verbose = FALSE)
# DimPlot(nsc, pt.size = 2, label = TRUE, label.size = 6) + theme_bw()
# DimPlot(nsc, pt.size = 2, split.by = 'orig.ident') + theme_bw()
