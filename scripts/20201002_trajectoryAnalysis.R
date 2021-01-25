
######### 20201002 Trajectory analysis of microglia and NSCs ##########

# Objective: Determine whether trajectories can be built from data and identify
# gene expression dynamics + important transcriptional regulators.


# Data Prep -----------------------------------------------------------------


# Libraries
library('Seurat')
library('ggplot2')
library('dplyr')

# Directories and data
data_path <- 'data/'
results_path <- 'results/'
new_results_path <- paste0(results_path, '20201002_trajectoryAnalysis/')
dir.create(path = new_results_path)
svz <- readRDS(file = paste0(data_path, '20200918_svz.rds'))
mg <- readRDS(file = paste0(data_path, '20200625_microglia.rds'))
nsc <- readRDS(file = paste0(data_path, '20200625_nsc.rds'))


# Set celltype/subcluster across datasets
mg$subcluster <- plyr::mapvalues(x = mg$RNA_snn_res.0.8,
                                 from = c(0,1,2,3),
                                 to = c('Homeo_MG','Unknown_MG','Activated_MG', 'Dev_MG'))
nsc$subcluster <- plyr::mapvalues(x = nsc$RNA_snn_res.0.8,
                                  from = c(0,1,2,3),
                                  to = c('Neuroblast','B_cell', 'Prolif_C_cell', 'Activated_C_cell'))
subcluster_labels <- c(as.character(mg$subcluster), as.character(nsc$subcluster))
names(subcluster_labels) <- c(names(mg$subcluster), names(nsc$subcluster))
svz$subcluster <- colnames(x = svz)
svz$subcluster <- plyr::mapvalues(x = svz$subcluster,
                                  from = names(subcluster_labels),
                                  to = subcluster_labels)
svz$subcluster <- ifelse(svz$subcluster %in% colnames(svz), yes = 'Unknown', no = svz$subcluster)
svz$subcluster <- factor(svz$subcluster,
                         levels = c('Homeo_MG','Unknown_MG','Activated_MG', 'Dev_MG','Neuroblast','B_cell', 'Prolif_C_cell', 'Activated_C_cell', 'Unknown'))
DimPlot(svz, group.by = 'subcluster', label = TRUE, label.size = 6, pt.size = 2, repel = TRUE)


# ggplot2 aesthetics
my_theme <- theme(panel.background = element_rect(fill = NA, color = 'black'),
                  panel.border = element_rect(fill = NA, color = 'black'),
                  axis.title = element_text(size = 14, color = 'black'),
                  axis.ticks = element_blank(),
                  axis.text = element_blank())



# Microglia Trajectory analysis using Monocle3 ------------------------------

# installation
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
library('monocle3')


# Convert object type
mg_sce <- as.SingleCellExperiment(x = mg)
rowData(mg_sce)[,'gene_short_name'] <- rownames(mg_sce)
mg_sce <- new_cell_data_set(expression_data = counts(mg_sce),
                          cell_metadata = colData(mg_sce),
                          gene_metadata = rowData(mg_sce))


# normalize, variable genes, pca, 
mg_sce <- preprocess_cds(cds = mg_sce, num_dim = 10)
plot_pc_variance_explained(mg_sce)


# nonlinear dimreduce
mg_sce <- reduce_dimension(mg_sce, preprocess_method = 'PCA', reduction_method = 'UMAP', umap.fast_sgd = FALSE)


# initial UMAP (via monocle 3 algorithm parameters)
umap_mg_sce <- plot_cells(cds = mg_sce, 
                          color_cells_by = 'subcluster',
                          cell_size = 2,
                          group_label_size = 8) + my_theme; umap_mg_sce
ggsave(filename = paste0(new_results_path, 'mg_umap_monocle.tiff'),
       plot = umap_mg_sce, device = 'tiff', height = 5, width = 5.5)

# import previous seurat result to learn path
tmp_umap <- mg@reductions$umap@cell.embeddings
reducedDims(mg_sce)$UMAP <- tmp_umap


# Cluster the cells. Iterate through resolutions for greatest modularity.
set.seed(123) # !!! Cluster result is stochastically initialized!!!
mods <- c()
for(i in seq(0.01, 0.3, 0.02)) {
  tmp <- cluster_cells(mg_sce, reduction_method = 'UMAP', resolution = i)
  mods <- c(mods, tmp@clusters$UMAP$cluster_result$optim_res$modularity)
}
mod_plot <- data.frame('res' = seq(0.01, 0.3, 0.02),
                       'mod' = mods) %>%
  ggplot(mapping = aes(x = res, y = mod)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.025)) +
  scale_x_continuous(breaks = seq(0.01, 0.3, 0.02)) +
  xlab(label = 'cluster_cells() resolution') +
  ylab(label = 'Modularity') +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black')); mod_plot
ggsave(filename = paste0(new_results_path, 'mg_modularity_optimization.tiff'),
       plot = mod_plot, device = 'tiff', height = 3, width = 7.5)
  

mg_sce <- cluster_cells(mg_sce, reduction_method = 'UMAP', resolution = 0.07)
mg_sce@clusters$UMAP <- colData(mg_sce)['subcluster']
plot_cells(mg_sce, color_cells_by = 'partition', cell_size = 3)
plot_cells(mg_sce, color_cells_by = 'cluster', cell_size = 3)


# generate trajectory + pseudotime (programmatically, see Monocle 3 tutorial)
mg_sce <- learn_graph(mg_sce)
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin='SCtrl'){
  cell_ids <- which(colData(cds)[,'orig.ident'] == time_bin)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  return(root_pr_nodes)
}
mg_sce <- order_cells(mg_sce, root_pr_nodes = get_earliest_principal_node(mg_sce, time_bin = 'SCtrl'))

# Visualize results
tmp_plot <- function(cds, grp, show_graph = FALSE) {
  tmp <- plot_cells(cds, 
                    color_cells_by = grp,
                    cell_size = 2, 
                    graph_label_size = 0, 
                    show_trajectory_graph = show_graph,
                    label_groups_by_cluster = FALSE,
                    group_label_size = 6,
                    trajectory_graph_segment_size = 2,
                    trajectory_graph_color = 'black') +
    my_theme
  return(tmp)
}
p1 <- tmp_plot(mg_sce, grp = 'subcluster', show_graph = FALSE); p1
p2 <- tmp_plot(mg_sce, grp = 'orig.ident', show_graph = FALSE); p2
p3 <- tmp_plot(mg_sce, grp = 'pseudotime', show_graph = TRUE) +
  guides(color = guide_colorbar(title = 'Pseudotime',
                                frame.colour = 'black', 
                                frame.linewidth = 1, 
                                ticks.colour = 'black', 
                                ticks.linewidth = 1,
                                barwidth = 1)); p3
tmp <- cowplot::plot_grid(NULL, p3, NULL, ncol = 1, rel_heights = c(0.5, 1, 0.5))
monocle_umaps <- cowplot::plot_grid(cowplot::plot_grid(p1, p2, ncol = 1), tmp, ncol = 2, rel_widths = c(0.8, 1)); monocle_umaps
ggsave(filename = paste0(new_results_path, 'mg_monocle_pseudotime.tiff'),
       plot = monocle_umaps, device = 'tiff', height = 6.5, width = 8)


# Finding modules of co-regulated genes using graph-autocorrelation
mg_graph_cor <- graph_test(mg_sce, neighbor_graph = 'principal_graph', cores = 1) # graph test
mg_graph_cor <- mg_graph_cor[order(mg_graph_cor$q_value, decreasing = FALSE),] # sort by significance
mg_graph_de <- mg_graph_cor$gene_short_name[mg_graph_cor$q_value < 0.05] # subset sig genes
mg_graph_de <- mg_graph_de[!is.na(mg_graph_de)] 
gene_module_df <- find_gene_modules(mg_sce[mg_graph_de,], resolution = c(10^seq(-5, -1, 1)))


# Plot module expression by cluster, cluster the modules
cell_group_df <- tibble::tibble('cell' = row.names(colData(mg_sce)), 
                                'cell_group' = colData(mg_sce)$subcluster)
agg_mat <- aggregate_gene_expression(mg_sce, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
module_heatmap <- pheatmap::pheatmap(agg_mat, 
                                     scale = "column", 
                                     clustering_method = "ward.D2",
                                     fontsize = 14)
ggsave(filename = paste0(new_results_path, 'mg_gene_module_heatmap.tiff'),
       plot = module_heatmap, device = 'tiff', height = 6, width = 4.5)


# Visualize gene modules in UMAP
mod_plot <- function(cds, genes) {
  tmp <- plot_cells(cds, 
                    genes = genes, 
                    cell_size = 3, 
                    trajectory_graph_segment_size = 2, 
                    trajectory_graph_color = 'red',
                    graph_label_size = 0) + 
    my_theme + 
    scale_color_viridis_c(option = 'A') +
    theme(strip.text = element_text(size = 16, color = 'black'),
          legend.text = element_text(size = 12, color = 'black')) +
    guides(color = guide_colorbar(title = 'Module\nscore',
                                  frame.colour = 'black', 
                                  frame.linewidth = 1, 
                                  ticks.colour = 'black', 
                                  ticks.linewidth = 1,
                                  barwidth = 1))
  return(tmp)
}
my_mods <- c(1,3,5,8)
tmp_plots <- list()
for(i in 1:length(my_mods)) {
  tmp_plots[[i]] <- mod_plot(mg_sce, genes = gene_module_df %>% filter(module %in% c(my_mods[i])))
}
sample_modules <- cowplot::plot_grid(plotlist = tmp_plots, ncol = 2); sample_modules
ggsave(filename = paste0(new_results_path, 'mg_sample_gene_modules_umap.tiff'),
       plot = sample_modules,
       height = 7, width = 8.5, device = 'tiff')
# plot_genes_in_pseudotime(mg_sce[rownames(mg_sce) %in% mg_graph_de[1:10]], min_expr = 0.5, color_cells_by = 'subcluster')

# Save data
saveRDS(mg_sce, file = 'results/20200927_trajectoryAnalysis/mg_sce_monocle.rds')



# Microglia Trajectory analysis using Slingshot -------------------------------------


# Libraries and data
BiocManager::install('slingshot')
library('slingshot')
library('SingleCellExperiment')
mg_sce <- as.SingleCellExperiment(mg)
mg_sce <- mg_sce[rownames(mg_sce) %in% VariableFeatures(mg),]


# Load PCA coordinates + subclusters (original subcluster analysis used 10 PCs)
mg_pca <- slot(slot(mg, 'reductions')[['pca']], 'cell.embeddings')
mg_umap <- slot(slot(mg, 'reductions')[['umap']], 'cell.embeddings')
npcs <- 1:10
mg_pca <- mg_pca[,npcs]
mg_umap <- mg_umap


# Calculate lineages using subcluster (MST built from the 4 subclusters)
mg_cols <- c('Activated_MG' = 'red', 'Unknown_MG' = 'orange', 'Dev_MG' = 'green', 'Homeo_MG' = 'blue')
plot(mg_umap, col = mg_cols[mg_sce$subcluster], pch=16, asp = 1)
reducedDim(mg_sce, type = 'UMAP') <- mg_umap
mg_sce <- slingshot(data = mg_sce, 
                    clusterLabels = 'subcluster', 
                    reducedDim = 'UMAP',
                    start.clus = 'Unknown_MG')
mg_sds <- SlingshotDataSet(mg_sce)


# Visualize the lineages
{
  tiff(filename = paste0(new_results_path, 'mg_lineages_slingshot.tiff'), height = 600, width = 650)
  par(mar = c(5,5,2,2))
  plot(reducedDims(mg_sce)$UMAP, 
       col = mg_cols[mg_sce$subcluster],
       pch=16, asp = 1, cex = 2, cex.lab = 2, cex.axis = 0.001)
  lines(SlingshotDataSet(mg_sce), lwd = 3, col='black')
  dev.off()
}



# Load tradeSeq for examining temporal changes in gene expression
BiocManager::install('tradeSeq', version = 'devel')
library('tradeSeq')
set.seed(123)


# Fit tradeSeq models and do Association Test b/w expression ~ pseudotime
mg_gam <- fitGAM(counts = counts(mg_sce), sds = mg_sds)
mg_at <- associationTest(mg_gam) # test within lineage, exp ~ pseudotime
mg_diffEnd <- diffEndTest(mg_gam) # test b/w lineages, exp ~ pseudotime(late) + lineage
mg_earlyDE <- earlyDETest(mg_gam) # test b/w lineages, exp ~ pseudotime(early) + lineage
mg_patternDE <- patternTest(mg_gam) 


# Visualize some results
tmp_theme <- my_theme + theme(legend.title = element_text(size = 14, color = 'black'),
                              legend.text = element_text(size = 12, color = 'black'),
                              plot.title = element_text(size = 16, color = 'black', face = 'bold'))
my_genes <- c('Actg1','Aqp5','Cd63', 'Lypd2', 'Ifitm1','Mbp','Plp1','Sdc1')
smooth_plots <- vector(mode = 'list', length = length(my_genes))
for(i in 1:length(my_genes)) {
  smooth_plots[[i]] <- plotSmoothers(gene = my_genes[i],
                                     mg_gam, counts = counts(mg_sce),
                                     lwd = 2, size = 1.5, border = TRUE) + 
    tmp_theme +
    labs(title = my_genes[i])
}
smooth_plots <- cowplot::plot_grid(plotlist = smooth_plots, ncol = 4)
smooth_plots
ggsave(filename = paste0(new_results_path, 'geneCurves_mg_slingshotLineages.tiff'),
       plot = smooth_plots, device = 'tiff', height = 5.5, width = 14)


# visualize top genes by p-val
topgenes <- rownames(mg_patternDE[order(mg_patternDE$pvalue), ])[1:250]
pst_order <- order(mg_sce$slingPseudotime_1, na.last = NA)
heat_data <- logcounts(mg_sce)[topgenes, pst_order]
heat_clus <- mg_sce$subcluster[pst_order]
heatmap(as.matrix(heat_data), Colv = NA,
        ColSideColors = RColorBrewer::brewer.pal(9,"Set1")[heat_clus])




# NSC Trajectory analysis using Monocle3 ------------------------------


library('monocle3')
nsc <- readRDS(file = paste0(data_path, '20200625_nsc.rds'))

# Convert object type
nsc_sce <- as.SingleCellExperiment(x = nsc)
rowData(nsc_sce)[,'gene_short_name'] <- rownames(nsc_sce)
nsc_sce <- new_cell_data_set(expression_data = counts(nsc_sce),
                            cell_metadata = colData(nsc_sce),
                            gene_metadata = rowData(nsc_sce))


# normalize, variable genes, pca, 
nsc_sce <- preprocess_cds(cds = nsc_sce, num_dim = 5)
plot_pc_variance_explained(nsc_sce)


# nonlinear dimreduce
nsc_sce <- reduce_dimension(nsc_sce, preprocess_method = 'PCA', reduction_method = 'UMAP', umap.fast_sgd = FALSE)


# initial UMAP (via monocle 3 algorithm parameters)
umap_nsc_sce <- plot_cells(cds = nsc_sce, 
                          color_cells_by = 'subcluster',
                          cell_size = 2,
                          group_label_size = 8) + my_theme; umap_nsc_sce
ggsave(filename = paste0(new_results_path, 'nsc_umap_monocle.tiff'),
       plot = umap_nsc_sce, device = 'tiff', height = 5, width = 5.5)

# import previous seurat result to learn path
tmp_umap <- nsc@reductions$umap@cell.embeddings
reducedDims(nsc_sce)$UMAP <- tmp_umap


# Cluster the cells. Iterate through resolutions for greatest modularity.
set.seed(123) # !!! Cluster result is stochastically initialized!!!
mods <- c()
for(i in seq(0.01, 0.3, 0.02)) {
  tmp <- cluster_cells(nsc_sce, reduction_method = 'UMAP', resolution = i)
  mods <- c(mods, tmp@clusters$UMAP$cluster_result$optim_res$modularity)
}
mod_plot <- data.frame('res' = seq(0.01, 0.3, 0.02),
                       'mod' = mods) %>%
  ggplot(mapping = aes(x = res, y = mod)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  scale_y_continuous(breaks = seq(0, 1, 0.025)) +
  scale_x_continuous(breaks = seq(0.01, 0.3, 0.02)) +
  xlab(label = 'cluster_cells() resolution') +
  ylab(label = 'Modularity') +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black')); mod_plot
ggsave(filename = paste0(new_results_path, 'nsc_modularity_optimization.tiff'),
       plot = mod_plot, device = 'tiff', height = 3, width = 7.5)


nsc_sce <- cluster_cells(nsc_sce, reduction_method = 'UMAP', resolution = 0.1)
plot_cells(nsc_sce, color_cells_by = 'partition', cell_size = 3)
plot_cells(nsc_sce, color_cells_by = 'cluster', cell_size = 3)


# generate trajectory + pseudotime (programmatically, see Monocle 3 tutorial)
nsc_sce <- learn_graph(nsc_sce)

# We do not use the earliest time-point helper function because we have two partitions
# PARTITION 1
nsc_sce <- order_cells(nsc_sce)

# Visualize results
tmp_plot <- function(cds, grp, show_graph = FALSE) {
  tmp <- plot_cells(cds, 
                    color_cells_by = grp,
                    cell_size = 2, 
                    graph_label_size = 0, 
                    show_trajectory_graph = show_graph,
                    label_groups_by_cluster = FALSE,
                    group_label_size = 6,
                    trajectory_graph_segment_size = 2,
                    trajectory_graph_color = 'black') +
    my_theme
  return(tmp)
}
p1 <- tmp_plot(nsc_sce, grp = 'subcluster', show_graph = FALSE); p1
p2 <- tmp_plot(nsc_sce, grp = 'orig.ident', show_graph = FALSE); p2
nsc_sce <- order_cells(nsc_sce)
p3 <- tmp_plot(nsc_sce, grp = 'pseudotime', show_graph = TRUE) +
  guides(color = guide_colorbar(title = 'Pseudotime',
                                frame.colour = 'black', 
                                frame.linewidth = 1, 
                                ticks.colour = 'black', 
                                ticks.linewidth = 1,
                                barwidth = 1)); p3
nsc_sce <- order_cells(nsc_sce)
p4 <- tmp_plot(nsc_sce, grp = 'pseudotime', show_graph = TRUE) +
  guides(color = guide_colorbar(title = 'Pseudotime',
                                frame.colour = 'black', 
                                frame.linewidth = 1, 
                                ticks.colour = 'black', 
                                ticks.linewidth = 1,
                                barwidth = 1)); p4
tmp <- cowplot::plot_grid(p3, p4, ncol = 1)
monocle_umaps <- cowplot::plot_grid(cowplot::plot_grid(p1, p2, ncol = 1), tmp, ncol = 2, rel_widths = c(0.8, 1)); monocle_umaps
ggsave(filename = paste0(new_results_path, 'nsc_monocle_pseudotime.tiff'),
       plot = monocle_umaps, device = 'tiff', height = 6.5, width = 8)


# Finding modules of co-regulated genes using graph-autocorrelation
nsc_graph_cor <- graph_test(nsc_sce, neighbor_graph = 'principal_graph', cores = 1) # graph test
nsc_graph_cor <- nsc_graph_cor[order(nsc_graph_cor$q_value, decreasing = FALSE),] # sort by significance
nsc_graph_de <- nsc_graph_cor$gene_short_name[nsc_graph_cor$q_value < 0.05] # subset sig genes
nsc_graph_de <- nsc_graph_de[!is.na(nsc_graph_de)] 
gene_module_df <- find_gene_modules(nsc_sce[nsc_graph_de,], resolution = c(10^seq(-5, -1, 1)))


# Plot module expression by cluster, cluster the modules
cell_group_df <- tibble::tibble('cell' = row.names(colData(nsc_sce)), 
                                'cell_group' = colData(nsc_sce)$subcluster)
agg_mat <- aggregate_gene_expression(nsc_sce, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
module_heatmap <- pheatmap::pheatmap(agg_mat, 
                                     scale = "column", 
                                     clustering_method = "ward.D2",
                                     fontsize = 14)
ggsave(filename = paste0(new_results_path, 'nsc_gene_module_heatmap.tiff'),
       plot = module_heatmap, device = 'tiff', height = 7.5, width = 4)


# Visualize gene modules in UMAP
mod_plot <- function(cds, genes) {
  tmp <- plot_cells(cds, 
                    genes = genes, 
                    cell_size = 3, 
                    trajectory_graph_segment_size = 2, 
                    trajectory_graph_color = 'red',
                    graph_label_size = 0) + 
    my_theme + 
    scale_color_viridis_c(option = 'A') +
    theme(strip.text = element_text(size = 16, color = 'black'),
          legend.text = element_text(size = 12, color = 'black')) +
    guides(color = guide_colorbar(title = 'Module\nscore',
                                  frame.colour = 'black', 
                                  frame.linewidth = 1, 
                                  ticks.colour = 'black', 
                                  ticks.linewidth = 1,
                                  barwidth = 1))
  return(tmp)
}
mod13 <- mod_plot(nsc_sce, genes = gene_module_df %>% filter(module %in% c(13)))
mod20 <- mod_plot(nsc_sce, genes = gene_module_df %>% filter(module %in% c(20)))
mod11 <- mod_plot(nsc_sce, genes = gene_module_df %>% filter(module %in% c(11)))
mod8 <- mod_plot(nsc_sce, genes = gene_module_df %>% filter(module %in% c(8)))
sample_modules <- cowplot::plot_grid(mod13, mod20, mod11, mod8, ncol = 2); sample_modules
ggsave(filename = paste0(new_results_path, 'nsc_sample_gene_modules_umap.tiff'),
       plot = sample_modules,
       height = 7, width = 8.5, device = 'tiff')
# plot_genes_in_pseudotime(nsc_sce[rownames(nsc_sce) %in% nsc_graph_de[1:10]], min_expr = 0.5, color_cells_by = 'subcluster')

# Save data
saveRDS(nsc_sce, file = 'results/20201002_trajectoryAnalysis/nsc_sce_monocle.rds')



# Microglia Trajectory analysis using Slingshot -------------------------------------


# Libraries and data
library('slingshot')
library('SingleCellExperiment')
nsc_sce <- as.SingleCellExperiment(nsc)
nsc_sce <- nsc_sce[rownames(nsc_sce) %in% VariableFeatures(nsc),]


# Load PCA coordinates + subclusters (original subcluster analysis used 10 PCs)
nsc_pca <- slot(slot(nsc, 'reductions')[['pca']], 'cell.embeddings')
nsc_umap <- slot(slot(nsc, 'reductions')[['umap']], 'cell.embeddings')
npcs <- 1:10
nsc_pca <- nsc_pca[,npcs]
nsc_umap <- nsc_umap


# Calculate lineages using subcluster (MST built from the 4 subclusters)
nsc_cols <- c('B_cell' = 'red', 'Activated_C_cell' = 'orange', 'Prolif_C_cell' = 'green', 'Neuroblast' = 'blue')
plot(nsc_umap, col = nsc_cols[nsc_sce$subcluster], pch=16, asp = 1)
reducedDim(nsc_sce, type = 'UMAP') <- nsc_umap
nsc_sce <- slingshot(data = nsc_sce, 
                     clusterLabels = 'subcluster', 
                     reducedDim = 'UMAP',
                     start.clus = 'Activated_C_cell')
nsc_sds <- SlingshotDataSet(nsc_sce)


# Visualize the lineages
{ # together
  tiff(filename = paste0(new_results_path, 'nsc_lineages_slingshot.tiff'), height = 600, width = 650)
  par(mar = c(5,5,2,2))
  plot(reducedDims(nsc_sce)$UMAP, 
       col = nsc_cols[nsc_sce$subcluster],
       pch=16, asp = 1, cex = 2, cex.lab = 2, cex.axis = 0.001)
  lines(SlingshotDataSet(nsc_sce), lwd = 3, col='black')
  dev.off()
}

{ #separately
  tiff(filename = paste0(new_results_path, 'nsc_slingshot_lineage1_umap.tiff'), height = 600, width = 650)
  plot(reducedDims(nsc_sce)$UMAP, col = nsc_cols[nsc_sce$subcluster], pch = 16, asp = 1)
  lines(nsc_sds@curves$curve1, lwd= 2, col = 'black')
  dev.off()
  
  tiff(filename = paste0(new_results_path, 'nsc_slingshot_lineage2_umap.tiff'), height = 600, width = 650)
  plot(reducedDims(nsc_sce)$UMAP, col = nsc_cols[nsc_sce$subcluster], pch = 16, asp = 1)
  lines(nsc_sds@curves$curve2, lwd= 2, col = 'black')
  dev.off()
}


# Load tradeSeq for examining temporal changes in gene expression
library('tradeSeq')
set.seed(123)


# Fit tradeSeq models and do Association Test b/w expression ~ pseudotime
nsc_gam <- fitGAM(counts = counts(nsc_sce), sds = nsc_sds)
nsc_at <- associationTest(nsc_gam) # test within lineage, exp ~ pseudotime
nsc_diffEnd <- diffEndTest(nsc_gam) # test b/w lineages, exp ~ pseudotime(late) + lineage
nsc_earlyDE <- earlyDETest(nsc_gam) # test b/w lineages, exp ~ pseudotime(early) + lineage
nsc_patternDE <- patternTest(nsc_gam)


# Visualize some results
tmp_theme <- my_theme + theme(legend.title = element_text(size = 14, color = 'black'),
                              legend.text = element_text(size = 12, color = 'black'),
                              plot.title = element_text(size = 16, color = 'black', face = 'bold'))
my_genes <- c('Dcx','Fxyd1','Fgfr1','Sox11','Ncam1', 'Ntrk2','Mt1','Mt2','Plpp3','Tubb3')
smooth_plots <- vector(mode = 'list', length = length(my_genes))
for(i in 1:length(my_genes)) {
  smooth_plots[[i]] <- plotSmoothers(gene = my_genes[i],
                                     nsc_gam, counts = counts(nsc_sce),
                                     lwd = 2, size = 1.5, border = TRUE) + 
    tmp_theme +
    labs(title = my_genes[i])
}
smooth_plots <- cowplot::plot_grid(plotlist = smooth_plots, ncol = 5)
smooth_plots
ggsave(filename = paste0(new_results_path, 'geneCurves_nsc_slingshotLineages.tiff'),
       plot = smooth_plots, device = 'tiff', height = 5.5, width = 16)



# visualize top genes by p-val
topgenes <- rownames(nsc_patternDE[order(nsc_patternDE$pvalue), ])[1:250]
pst_order <- order(nsc_sce$slingPseudotime_1, na.last = NA)
heat_data <- t(scale(t(logcounts(nsc_sce)[topgenes, pst_order])))
heat_clus <- nsc_sce$subcluster[pst_order]


# heat_data <- nsc[['RNA']]@scale.data[topgenes, pst_order]
data.frame(t(heat_data), 'subcluster' = heat_clus) %>%
  tibble::rownames_to_column(var = 'barcode') %>%
  reshape2::melt(id.vars = c('barcode', 'subcluster')) %>%
  ggplot() + 
  geom_raster(mapping = aes(x = barcode, y = variable, fill = value))
heatmap(t(scale(t(as.matrix(heat_data)))), Colv = NA,
        ColSideColors = RColorBrewer::brewer.pal(9,"Set1")[heat_clus],
        col = RColorBrewer::brewer.pal(n = 9, name = 'Spectral'))


# Cluster the genes by similar temporal pattern
library('clusterExperiment')
npoints <- 20
patterns <- clusterExpressionPatterns(nsc_gam, nPoints = npoints, genes = VariableFeatures(nsc), nReducedDims = 5)
clusterLabels <- primaryCluster()


save(list = ls(), file = paste0(new_results_path, 'session_data.RData'))
