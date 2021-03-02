

require('Seurat')
require('ggplot2')
require('dplyr')


if (!grepl('scripts', getwd())) {
  setwd('./scripts/')
}
results_out <- '../results/figs/'
dir.create(path = results_out)

svz <- readRDS(file = '../data/20210202_SVZ.rds')
mg <- readRDS(file = '../data/20210112_MG.rds')
nsc <- readRDS(file = '../data/20210108_NSC.rds')


# Preliminary cluster (not labeled) ---------------------------------------

p1 <- FetchData(svz, vars = c('SCT_snn_res.2', 'UMAP_1', 'UMAP_2')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(mapping = aes(fill = SCT_snn_res.2),
             color = 'black', pch = 21, size = 3) +
  theme_bw() +
  theme(legend.title = element_blank())
p1
ggsave(filename = paste0(results_out, 'PreliminaryCluster_notLabeled.svg'),
       plot = p1, device = 'svg', height = 3.5, width = 4.5)



# Preliminary Cluster depecting categorization genes ----------------------

# Myeloid genes
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
  order = TRUE,
  combine = FALSE
)
myeloid_umap <- lapply(
  X = myeloid_umap,
  FUN = function(x) {
    x + theme_bw()
  }
)
myeloid_umap <- cowplot::plot_grid(plotlist = myeloid_umap)
ggsave(filename = paste0(results_out, 'PreliminaryCluster_myeloid-genes_umap.svg'),
       plot = myeloid_umap, device = 'svg', height = 8, width = 13.5)


# NSPC genes
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
nsc_umap <- cowplot::plot_grid(plotlist = nsc_umap, ncol = 5)
ggsave(filename = paste0(results_out, 'PreliminaryCluster_nsc-genes_umap.svg'),
       plot = nsc_umap, device = 'svg', height = 7.5, width = 15)


ependymal_genes <- c('Foxj1','Cfap126','Ccdc153', 'Rarres2','Tmem212')
svz$seurat_clusters <- svz$SCT_snn_res.2
Idents(svz) <- 'seurat_clusters'
DefaultAssay(svz) <- 'RNA'
ependymal_umap <- FeaturePlot(
  object = svz,
  features = ependymal_genes,
  cols = c('grey','red3'),
  slot = 'data',
  order = TRUE,
  combine = FALSE
)
ependymal_umap <- lapply(
  X = ependymal_umap,
  FUN = function(x) {
    x + theme_bw()
  }
)
ependymal_umap <- cowplot::plot_grid(plotlist = ependymal_umap, ncol = 3)
ggsave(filename = paste0(results_out, 'PreliminaryCluster_ependymal-genes_umap.svg'),
       plot = ependymal_umap, device = 'svg', height = 4.75, width = 9)



# Preliminary cluster celltype categorization -----------------------------

Idents(svz) <- 'celltype'
p1 <-  FetchData(svz, vars = c('celltype', 'UMAP_1', 'UMAP_2')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(mapping = aes(fill = celltype),
             color = 'black', pch = 21, size = 3) +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave(filename = paste0(results_out, 'PreliminaryCluster_celltype-umap.svg'),
       plot = p1, device = 'svg', height = 3.5, width = 5.5)


# Preliminary celltype depecting categorization genes ----------------------

# Myeloid genes
myeloid_genes <- c('P2ry12','Tmem119','Cx3cr1','Csf1r','Adgre1','Ptprc','Ly6c2',
                   'Ccr2','Mrc1','Lyve1','Mki67','Top2a')
svz$seurat_clusters <- svz$celltype
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
  xlab(label = 'Cell-type') + 
  ylab(label = 'log-normalized expression') + 
  theme_bw() + 
  theme(legend.position = 'none',
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(size = 12, hjust = 0),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
ggsave(
  filename = paste0(results_out, 'PreliminaryCluster_celltype_myeloid-genes_vln.svg'),
  plot = myeloid_vln,
  device = 'svg',
  height = 5, width = 8
)

# NSPC genes
nsc_genes <- c('Gfap','Sox2','Dcx','Ascl1','Sox9','Id3','Id4','Foxm1','Egfr',
               'Mki67','Top2a','Dlx1','Dlx2','Fos')
svz$seurat_clusters <- svz$celltype
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
  xlab(label = 'Cell-type') + 
  ylab(label = 'log-normalized expression') + 
  theme_bw() + 
  theme(legend.position = 'none',
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(size = 12, hjust = 0),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
ggsave(
  filename = paste0(results_out, 'PreliminaryCluster_celltype_nsc-genes_vln.svg'),
  plot = nsc_vln,
  device = 'svg',
  height = 6, width = 8
)

# ependymal genes
ependymal_genes <- c('Foxj1','Cfap126','Ccdc153', 'Rarres2','Tmem212')
svz$seurat_clusters <- svz$celltype
Idents(svz) <- 'seurat_clusters'
DefaultAssay(svz) <- 'RNA'
expr_dat <- FetchData(
  object = svz,
  vars = c('seurat_clusters', ependymal_genes),
  slot = 'data'
)
ependymal_vln <- expr_dat %>%
  reshape2::melt(id.vars = 'seurat_clusters') %>%
  ggplot(mapping = aes(x = seurat_clusters, y = value)) +
  geom_violin(mapping = aes(fill = seurat_clusters), scale = 'width') +
  facet_wrap(variable ~ ., scales = 'free_y') +
  xlab(label = 'Cell-type') + 
  ylab(label = 'log-normalized expression') + 
  theme_bw() + 
  theme(legend.position = 'none',
        strip.background = element_rect(color = NA, fill = NA),
        strip.text = element_text(size = 12, hjust = 0),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1))
ggsave(
  filename = paste0(results_out, 'PreliminaryCluster_celltype_ependymal-genes_vln.svg'),
  plot = ependymal_vln,
  device = 'svg',
  height = 4, width = 6
)



# Microglia subcluster (by time + cluster) --------------------------------

Idents(mg) <- 'SCT_snn_res.0.4'
p1 <- FetchData(mg, vars = c('SCT_snn_res.0.4', 'orig.ident','UMAP_1', 'UMAP_2')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(mapping = aes(fill = SCT_snn_res.0.4),
             color = 'black', pch = 21, size = 3) +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave(filename = paste0(results_out, 'MG_subcluster_umap.svg'), 
       plot = p1, device = 'svg', height = 3, width = 4)

p1 <- FetchData(mg, vars = c('SCT_snn_res.0.4', 'orig.ident','UMAP_1', 'UMAP_2')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(mapping = aes(fill = SCT_snn_res.0.4),
             color = 'black', pch = 21, size = 3) +
  facet_grid(. ~ orig.ident) + 
  theme_bw() +
  theme(legend.title = element_blank(),
        strip.text = element_text(size = 12))
ggsave(filename = paste0(results_out, 'MG_subcluster-byTime_umap.svg'), 
       plot = p1, device = 'svg', height = 3, width = 9)




# MG DE genes heatmap -----------------------------------------------------

Idents(mg) <- 'SCT_snn_res.0.4'
DefaultAssay(mg) <- 'SCT'

mg_markers <- read.csv(file = '../results/ClusterAnalysisMG/MG_subcluster_markers.csv')
top_genes <- mg_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = -p_val_adj) %>%
  top_n(n = 10, wt = avg_logFC)
top_genes <- unique(top_genes$gene)

my_cols <- rev(RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))
avg_exp <- ScaleData(object = mg[['SCT']]@data, features = top_genes)
marker_heatmap <- cbind(data.frame(t(avg_exp)), 'ident' = mg@active.ident) %>%
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
                       limits = c(NA, 2.5),
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
marker_heatmap
ggsave(filename = paste0(results_out, 'MG_10-DE-genes_heatmap.svg'),
       plot = marker_heatmap, device = 'svg', height = 7, width = 3.5)




# MG Cluster specific Summed expression (marcos) --------------------------

Idents(mg) <- 'SCT_snn_res.0.4'
DefaultAssay(mg) <- 'SCT'

mg_markers <- read.csv(file = '../results/ClusterAnalysisMG/MG_subcluster_markers.csv')
top_genes <- mg_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = -p_val_adj) %>%
  top_n(n = 5, wt = avg_logFC)
# top_genes <- unique(top_genes$gene)

markers_umap <- vector(mode = 'list', 
                       length = length(unique(top_genes$cluster)))
names(markers_umap) <- unique(mg_markers$cluster)
for (c in unique(mg_markers$cluster)) {
  tmp_genes <- top_genes$gene[top_genes$cluster == c]
  exp_sum <- Matrix::colSums(mg[['RNA']]@data[tmp_genes,])
  plot_dat <- data.frame(cbind(exp_sum, mg[['umap']]@cell.embeddings))
  markers_umap[[c+1]] <- plot_dat %>%
    arrange(exp_sum) %>%
    ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(mapping = aes(fill = exp_sum),
               color = 'black',
               pch = 21,
               size = 3) +
    labs(title = paste0('Cluster ', c),
         subtitle = paste0('(', paste(tmp_genes, collapse = ', '), ')')) +
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
ggsave(filename = paste0(results_out, 'MG_10-DE-genes_sumExp-umap.svg'),
       plot = markers_umap, device = 'svg', height = 6, width = 11)




# MG population numbers ---------------------------------------------------


Idents(mg) <- 'SCT_snn_res.0.4'
counts <- data.frame(table(mg$SCT_snn_res.0.4, mg$orig.ident))
names(counts) <- c('subtype','orig.ident','count')
numbers <- counts %>%
  ggplot(mapping = aes(x = orig.ident, y = count, group = subtype)) + 
  geom_bar(aes(fill = subtype), stat = 'identity', color = 'black', size = 0.75) +
  scale_y_continuous() +
  scale_x_discrete() + 
  ylab('# of cells') + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.line = element_line(size = 1),
        panel.background = element_rect(fill = NA),
        axis.text.y = element_text(size = 12))
props <- counts %>% 
  ggplot(mapping = aes(x = orig.ident, y = count, group = subtype)) + 
  geom_bar(aes(fill = subtype), stat = 'identity', color = 'black', size = 0.75, position = 'fill') +
  scale_y_continuous() + 
  scale_x_discrete() +
  ylab('% of cells') + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.line = element_line(size = 1), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(size = 12), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12))
legend <- cowplot::get_legend(plot = props)
pop.dynamics <- cowplot::plot_grid(numbers + theme(legend.position = 'none'), props + theme(legend.position = 'none'), legend, ncol = 3, rel_widths = c(1, 1, 0.25))
ggsave(filename = paste0(results_out, 'MG_population_dynamics.svg'),
       plot = pop.dynamics, device = 'svg', height = 3, width = 6)




# NSPC subcluster (by time + cluster) --------------------------------

Idents(nsc) <- 'SCT_snn_res.0.8'
p1 <- FetchData(nsc, vars = c('SCT_snn_res.0.8', 'orig.ident','UMAP_1', 'UMAP_2')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(mapping = aes(fill = SCT_snn_res.0.8),
             color = 'black', pch = 21, size = 3) +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave(filename = paste0(results_out, 'NSC_subcluster_umap.svg'), 
       plot = p1, device = 'svg', height = 3, width = 4)

p1 <- FetchData(nsc, vars = c('SCT_snn_res.0.8', 'orig.ident','UMAP_1', 'UMAP_2')) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(mapping = aes(fill = SCT_snn_res.0.8),
             color = 'black', pch = 21, size = 3) +
  facet_grid(. ~ orig.ident) + 
  theme_bw() +
  theme(legend.title = element_blank(),
        strip.text = element_text(size = 12))
ggsave(filename = paste0(results_out, 'NSC_subcluster-byTime_umap.svg'), 
       plot = p1, device = 'svg', height = 3, width = 9)




# NSPC DE genes heatmap -----------------------------------------------------

Idents(nsc) <- 'SCT_snn_res.0.8'
DefaultAssay(nsc) <- 'SCT'

nsc_markers <- read.csv(file = '../results/ClusterAnalysisNSC/NSC_subcluster_markers.csv')
top_genes <- nsc_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = -p_val_adj) %>%
  top_n(n = 10, wt = avg_logFC)
top_genes <- unique(top_genes$gene)

my_cols <- rev(RColorBrewer::brewer.pal(n = 9, name = 'RdBu'))
avg_exp <- ScaleData(object = nsc[['SCT']]@data, features = top_genes)
marker_heatmap <- cbind(data.frame(t(avg_exp)), 'ident' = nsc@active.ident) %>%
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
                       limits = c(NA, 2.5),
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
marker_heatmap
ggsave(filename = paste0(results_out, 'NSC_10-DE-genes_heatmap.svg'),
       plot = marker_heatmap, device = 'svg', height = 6, width = 3.5)




# NSC population dynamics -------------------------------------------------

Idents(nsc) <- 'SCT_snn_res.0.8'
counts <- data.frame(table(nsc$SCT_snn_res.0.8, nsc$orig.ident))
names(counts) <- c('subtype','orig.ident','count')
numbers <- counts %>%
  ggplot(mapping = aes(x = orig.ident, y = count, group = subtype)) + 
  geom_bar(aes(fill = subtype), stat = 'identity', color = 'black', size = 0.75) +
  scale_y_continuous() +
  scale_x_discrete() + 
  ylab('# of cells') + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.line = element_line(size = 1),
        panel.background = element_rect(fill = NA),
        axis.text.y = element_text(size = 12))
props <- counts %>% 
  ggplot(mapping = aes(x = orig.ident, y = count, group = subtype)) + 
  geom_bar(aes(fill = subtype), stat = 'identity', color = 'black', size = 0.75, position = 'fill') +
  scale_y_continuous() + 
  scale_x_discrete() +
  ylab('% of cells') + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.line = element_line(size = 1), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(size = 12), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12))
legend <- cowplot::get_legend(plot = props)
pop.dynamics <- cowplot::plot_grid(numbers + theme(legend.position = 'none'), props + theme(legend.position = 'none'), legend, ncol = 3, rel_widths = c(1, 1, 0.25))
ggsave(filename = paste0(results_out, 'NSC_population_dynamics.svg'),
       plot = pop.dynamics, device = 'svg', height = 3, width = 6)


# MG Marker violins ----------------------------------------------------------

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
ggsave(filename = paste0(results_out, 'MG_markers_vln.svg'), 
       device = 'svg', plot = markers_vln, height = 3, width = 15)


# NSPC Marker violins ----------------------------------------------------------

nsc$subtype <- nsc$SCT_snn_res.0.8
Idents(nsc) <- 'subtype'
markers <- FindAllMarkers(
  object = nsc,
  assay = 'RNA',
  # test.use = 'MAST',
  only.pos = TRUE
)
top_markers <- markers %>%
  filter(gene %in% c('Stmn2','Tubb3','Dcx','Basp1','Aldoc','Gja1','Gfap','Ttyh1','Ascl1','Egfr','Thrsp','Hes6','Top2a','Rrm2','Fbxo5','Mki67')) %>%
  group_by(cluster) %>%
  top_n(n = 4, wt = -p_val_adj)

markers_vln <- vector(mode = 'list',
                      length = length(levels(top_markers$cluster)))
names(markers_vln) <- levels(top_markers$cluster)
for (c in levels(top_markers$cluster)) {
  top_genes <- top_markers$gene[top_markers$cluster == c]
  markers_vln[[c]] <- FetchData(
    object = nsc,
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
markers_vln <- markers_vln[[1]] + markers_vln[[2]] + markers_vln[[3]] + markers_vln[[4]] + patchwork::plot_layout(ncol = 4)
# markers_vln <- cowplot::plot_grid(plotlist = markers_vln, ncol = 4)
markers_vln
ggsave(filename = paste0(results_out, 'nsc_markers_vln.svg'), 
       device = 'svg', plot = markers_vln, height = 3, width = 12)

