
#' ---
#' title: Cluster Analysis of Microglia
#' author: James Choi
#' date: 2021/1/12
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)

require('Seurat')
require('ggplot2')
require('dplyr')
require('dendextend')

setwd('./scripts')
results_out <- '../results/ClusterAnalysisMG/'
dir.create(path = results_out)
svz <- readRDS(file = '../data/20210202_SVZ.rds')

#' ## Cluster Analysis of Microglia
# subset microglia
mg <- svz[,svz$celltype == 'Microglia']
DefaultAssay(mg) <- 'RNA'

#+ variable_features, fig.height=3, fig.width=5.25
mg <- NormalizeData(mg, verbose = FALSE)
summary(mg$nFeature_RNA)
mg <- SCTransform(mg,
                   assay = 'RNA',
                   variable.features.n = 2000,
                   vars.to.regress = 'CC.difference',
                   verbose = FALSE)
VariableFeaturePlot(mg) + theme_bw()

#+ elbow, fig.height=2.5, fig.width=3
ElbowPlot(mg, ndims = 50) + theme_bw()

#+ iterate_pca, fig.height=6, fig.width=17.5
mg <- RunPCA(mg, npcs = 50, verbose = FALSE)
test_pcs <- 4:13
umaps <- vector(mode = "list", length = length(test_pcs))
for (i in 1:length(test_pcs)) {
  mg <- FindNeighbors(mg, dims = 1:test_pcs[i], verbose = FALSE)
  mg <- FindClusters(mg, verbose = FALSE)
  mg <- RunUMAP(mg, dims = 1:test_pcs[i], n.neighbors = 30, verbose = FALSE)
  umaps[[i]] <- DimPlot(mg, pt.size = 1, label = TRUE, label.size = 6) +
    theme_bw() + labs(title = paste("# PCs:", test_pcs[i]))
}
umaps <- cowplot::plot_grid(plotlist = umaps, ncol = length(test_pcs)/2)
umaps

#+ dimloadings_sct, fig.height=8, fig.width=14
VizDimLoadings(mg, dims = test_pcs, nfeatures = 20, ncol = 5)

#+ pca_plots, fig.height=2.75, fig.width=9
p1 <- DimPlot(mg, pt.size = 2, reduction = 'pca', dims = c(2,3)) + theme_bw()
p2 <- DimPlot(mg, pt.size = 2, reduction = 'pca', dims = c(4,5)) + theme_bw()
p3 <- DimPlot(mg, pt.size = 2, reduction = 'pca', dims = c(6,7)) + theme_bw()
pca_plots <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
pca_plots


#+ mg_umap, fig.height=3.5, fig.width=4.25
pcs <- 1:7
mg <- FindNeighbors(mg, dims = pcs, verbose = FALSE)
mg <- RunUMAP(mg, dims = pcs, verbose = FALSE)
mg <- FindClusters(mg, resolution = 0.4, verbose = FALSE)
DimPlot(mg, pt.size = 2, label = TRUE, label.size = 6) + theme_bw()


#+ hierarchical_dend, fig.height=5, fig.width=7.5
mg_dist <- dist(x = mg[['pca']]@cell.embeddings[,pcs])
mg_tree <- hclust(d = mg_dist, method = 'ward.D2')
mg_dend <- dendsort::dendsort(d = as.dendrogram(mg_tree, hang = 1))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
tmp_cols <- gg_color_hue(length(levels(mg$SCT_snn_res.0.4)))
names(tmp_cols) <- levels(mg$SCT_snn_res.0.4)
labels_colors(mg_dend) <- 
  tmp_cols[mg$SCT_snn_res.0.4][order.dendrogram(mg_dend)]
plot(mg_dend,
     main = 'Hierarchical clustering of Microglia',
     sub = 'Labels colored by original cluster results.')
abline(h = 115)

#+ cluster_comparison, fig.height=3.5, fig.width=8
mg_dclus <- cutree(tree = mg_dend, h = 115)
mg$hier_clust <- paste('H', mg_dclus[match(rownames(mg@meta.data),
                                             names(mg_dclus))],
                        sep = '_')
mg$hier_clust <- factor(mg$hier_clust)

p1 <- DimPlot(mg, 
              group.by = 'hier_clust', 
              pt.size = 1.5, 
              label = TRUE, 
              label.size = 5) +
  theme_bw() +
  labs(title = 'Hierarchical cluster results\noverlaid on UMAP')
p2 <- DimPlot(mg, 
              group.by = 'SCT_snn_res.0.4', 
              pt.size = 1.5, 
              label = TRUE, 
              label.size = 5) +
  theme_bw() +
  labs(title = 'SNN-graph cluster results\noverlaid on UMAP')
cowplot::plot_grid(p1, p2, ncol = 2)

#' Inspect distribution of cells/clusters over injury time:  
#' 

#+ umap_over_time, fig.height=2.75, fig.width=8.5
Idents(mg) <- 'SCT_snn_res.0.4'
DimPlot(mg, pt.size = 2, split.by = 'orig.ident') + 
  theme_bw() +
  theme(strip.text = element_text(size = 12))


#' Identify marker genes:  
#' 
ensembl_conversion <- read.table(file = '../')
Idents(mg) <- 'SCT_snn_res.0.4'
markers <- FindAllMarkers(
  object = mg,
  assay = 'RNA',
  # test.use = 'MAST',
  only.pos = TRUE
)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = -p_val_adj)
top_markers$p_val <- as.character(signif(top_markers$p_val, digits = 3))
top_markers$p_val_adj <- as.character(signif(top_markers$p_val_adj, digits = 3))
knitr::kable(x = top_markers, digits = 3)
table("# of DE genes by cluster" = markers$cluster[markers$p_val_adj < 1e-03])

# Convert gene-names to Ensembl IDs
ensembl <- useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', 
                      dataset = 'mmusculus_gene_ensembl')
conversion_input <- rownames(mg[['RNA']]@counts)
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
write.table(markers, file = paste0(results_out, 'MG_subcluster_markers.csv'),
            sep = ',', row.names = FALSE)


# p1 <- DimPlot(mg, pt.size = 2, label = TRUE, label.size = 6) + theme_bw() +
#   NoLegend()
# ggsave(filename = '../results/MG_umap_byCluster.svg', plot = p1, device = 'svg',
#        height = 3, width = 3.25)
# p1 <- DimPlot(mg, pt.size = 2, group.by = 'orig.ident') + theme_bw() +
#   NoLegend()
# ggsave(filename = '../results/MG_umap_byTime.svg', plot = p1, device = 'svg',
#        height = 3, width = 3.25)


saveRDS(mg, file = '../data/20210112_MG.rds')




# Marker sum exp umaps ----------------------------------------------------

mg <- readRDS('../data/20210112_MG.rds')

DimPlot(mg, pt.size = 2, label = TRUE, label.size = 6)

Idents(mg) <- 'SCT_snn_res.0.4'
markers <- FindAllMarkers(
  object = mg,
  assay = 'RNA',
  # test.use = 'MAST',
  only.pos = TRUE
)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = -p_val_adj)

markers_umap <- vector(mode = 'list', 
                       length = length(levels(top_markers$cluster)))
names(markers_umap) <- levels(top_markers$cluster)
for (c in levels(top_markers$cluster)) {
  top_genes <- top_markers$gene[top_markers$cluster == c]
  exp_sum <- Matrix::colSums(mg[['RNA']]@data[top_genes,])
  plot_dat <- data.frame(cbind(exp_sum, mg[['umap']]@cell.embeddings))
  markers_umap[[c]] <- plot_dat %>%
    arrange(exp_sum) %>%
    ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(mapping = aes(fill = exp_sum),
               color = 'black',
               pch = 21,
               size = 3) +
    labs(title = paste0('Cluster ', c),
         subtitle = paste0('(', paste(top_genes, collapse = ', '), ')')) +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = 'RdYlBu'))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12),
          panel.background = element_rect(fill = NA, colour = 'black'),
          panel.border = element_rect(fill = NA, color = 'black'),
          legend.title = element_text(angle = 90),
          plot.title = element_text(face = 'bold')) +
    guides(fill = guide_colorbar(
      title = 'Sum of expression',
      ticks.colour = 'black',
      frame.colour = 'black',
      title.position = 'left'))
}
markers_umap <- cowplot::plot_grid(plotlist = markers_umap, ncol = 3)
markers_umap
ggsave(filename = '../results/MG_markers_sumexp_umap.svg',
       plot = markers_umap, device = 'svg', height = 6, width = 11)



# Plot specific genes on UMAP -------------------------------------------------

# CHANGE THIS HERE: 
plot_genes <- data.frame(
  'subtype' = c('HMG','HMG','HMG','AMG','AMG','AMG'),
  'gene' = c('P2ry12','Tmem119','Cx3cr1','C1qa','Lamp1','Cd68')
)

markers_umap <- vector(mode = 'list', length = length(unique(plot_genes$subtype)))
names(markers_umap) <- unique(plot_genes$subtype)
for (c in unique(plot_genes$subtype)) {
  top_genes <- plot_genes$gene[plot_genes$subtype == c]
  exp_sum <- Matrix::colSums(mg[['RNA']]@data[top_genes,])
  plot_dat <- data.frame(cbind(exp_sum, mg[['umap']]@cell.embeddings))
  markers_umap[[c]] <- plot_dat %>%
    arrange(exp_sum) %>%
    ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(mapping = aes(fill = exp_sum),
               color = 'black',
               pch = 21,
               size = 3) +
    labs(title = paste0('Cluster ', c),
         subtitle = paste0('(', paste(top_genes, collapse = ', '), ')')) +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = 'RdYlBu'))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 12),
          panel.background = element_rect(fill = NA, colour = 'black'),
          panel.border = element_rect(fill = NA, color = 'black'),
          legend.title = element_text(angle = 90),
          plot.title = element_text(face = 'bold')) +
    guides(fill = guide_colorbar(
      title = 'Sum of expression',
      ticks.colour = 'black',
      frame.colour = 'black',
      title.position = 'left'))
}
markers_umap <- cowplot::plot_grid(plotlist = markers_umap, ncol = 2)
markers_umap


# Marker violins ----------------------------------------------------------

mg <- readRDS('../data/20210112_MG.rds')

DimPlot(mg, pt.size = 2, label = TRUE, label.size = 6)

mg$subtype <- mg$SCT_snn_res.0.4
Idents(mg) <- 'subtype'
markers <- FindAllMarkers(
  object = mg,
  assay = 'RNA',
  # test.use = 'MAST',
  only.pos = TRUE
)
top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = -p_val_adj)

markers_vln <- vector(mode = 'list',
                      length = length(levels(top_markers$cluster)))
names(markers_vln) <- levels(top_markers$cluster)
for (c in levels(top_markers$cluster)) {
  top_genes <- top_markers$gene[top_markers$cluster == c]
  markers_vln[[c]] <- FetchData(
    object = mg,
    vars = c(top_genes, 'subtype'),
    slot = 'data'
  ) %>%
    reshape2::melt(id.vars = 'subtype') %>%
    ggplot(mapping = aes(x = subtype, y = value)) +
    geom_violin(mapping = aes(fill = subtype), 
                color = 'black', 
                scale = 'width') +
    # ggbeeswarm::geom_quasirandom(size = 0.5, alpha = 0.5) +
    facet_wrap(. ~ variable, 
               scales = 'free_y', 
               ncol = 1, 
               strip.position = 'left') +
    scale_y_continuous(position = 'right') +
    labs(title = paste('Cluster', c, 'DE genes:')) +
    ylab(label = 'Normalized expression') +
    theme(panel.background = element_rect(fill = NA, color = 'black'),
          panel.border = element_rect(fill = NA, color = 'black'),
          legend.position = 'none',
          strip.background = element_rect(fill = NA, color = NA),
          strip.text.y.left = element_text(angle = 0, size = 12),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y.right = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12))
}
markers_vln <- cowplot::plot_grid(plotlist = markers_vln, ncol = 5)
markers_vln
ggsave(filename = '../results/MG_markers_vln.svg', device = 'svg',
       plot = markers_vln, height = 3, width = 15)



# Plot specific genes on Violin -------------------------------------------

# CHANGE THIS HERE: 
plot_genes <- data.frame(
  'subtype' = c('HMG','HMG','HMG','AMG','AMG','AMG'),
  'gene' = c('P2ry12','Tmem119','Cx3cr1','C1qa','Lamp1','Cd68')
)
mg$subtype <- mg$SCT_snn_res.0.4
Idents(mg) <- 'subtype'
DefaultAssay(mg) <- 'RNA'

markers_vln <- vector(mode = 'list', length = length(unique(plot_genes$subtype)))
names(markers_vln) <- unique(plot_genes$subtype)
for (c in unique(plot_genes$subtype)) {
  top_genes <- plot_genes$gene[plot_genes$subtype == c]
  markers_vln[[c]] <- FetchData(
    object = mg,
    vars = c(top_genes, 'subtype'),
    slot = 'data'
  ) %>%
    reshape2::melt(id.vars = 'subtype') %>%
    ggplot(mapping = aes(x = subtype, y = value)) +
    geom_violin(mapping = aes(fill = subtype), 
                color = 'black', 
                scale = 'width') +
    # ggbeeswarm::geom_quasirandom(size = 0.5, alpha = 0.5) +
    facet_wrap(. ~ variable, 
               scales = 'free_y', 
               ncol = 1, 
               strip.position = 'left') +
    scale_y_continuous(position = 'right') +
    labs(title = paste('Cluster', c, 'DE genes:')) +
    ylab(label = 'Normalized expression') +
    theme(panel.background = element_rect(fill = NA, color = 'black'),
          panel.border = element_rect(fill = NA, color = 'black'),
          legend.position = 'none',
          strip.background = element_rect(fill = NA, color = NA),
          strip.text.y.left = element_text(angle = 0, size = 12),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y.right = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 12))
}
markers_vln <- cowplot::plot_grid(plotlist = markers_vln, ncol = 2)
markers_vln



# Dot plot ----------------------------------------------------------------

my_genes <- c('Tlr4','C1qa','Lamp1','P2ry12','Cd68')
DotPlot(mg, features = my_genes, group.by = 'orig.ident') # by time-point
DotPlot(mg, features = my_genes) # by cluster 
