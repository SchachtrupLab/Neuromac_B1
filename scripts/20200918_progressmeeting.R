

######## Data for progress meeting Sept 18th 2020 ########

# Progress report for meeting dated 2020-09-18 between Suvra and James PIs
# Christian Schachtrup (Freiburg) and Jae Lee (Miami).


# libraries
library('Seurat')
library('ggplot2')
library('dplyr')


# Directories and data
data_path <- 'data/'
results_path <- 'results/'
new_results_path <- paste0(results_path, '20200918_progressmeeting/')
dir.create(path = new_results_path)

# Import
svz <- readRDS(file = paste0(data_path, '20200511_SVZ.rds'))
mg <- readRDS(file = paste0(data_path, '20200625_microglia.rds'))
nsc <- readRDS(file = paste0(data_path, '20200625_nsc.rds'))

# Set celltype/subcluster across datasets
svz$full_clusters <- svz$RNA_snn_res.0.8
svz$celltype <- plyr::mapvalues(x = svz$RNA_snn_res.0.8,
                                from = 0:7,
                                to = c('Microglia',
                                       'Microglia',
                                       'Microglia',
                                       'NSC',
                                       'NSC',
                                       'NSC',
                                       'Microglia',
                                       'Unknown'))
mg$subcluster <- paste('mg', mg$RNA_snn_res.0.8, sep = '')
nsc$subcluster <- paste('nsc', nsc$RNA_snn_res.0.8, sep = '')
subcluster_labels <- c(mg$subcluster, nsc$subcluster)
svz$subcluster <- colnames(x = svz)
svz$subcluster <- plyr::mapvalues(x = svz$subcluster,
                                  from = names(subcluster_labels),
                                  to = subcluster_labels)
svz$subcluster <- ifelse(svz$subcluster %in% colnames(svz), yes = 'Unknown', no = svz$subcluster)
svz$subcluster <- factor(svz$subcluster, levels = sort(unique(svz$subcluster)))
DimPlot(svz, group.by = 'subcluster', label = TRUE, label.size = 6, pt.size = 2)
saveRDS(svz, file = 'data/20200918_svz.rds')



# Cell type annotation ----------------------------------------------------

# According to Dulken et al. Cell Reports 2017.

# quiescent-like          :Egfr(-)
# early activated         :Egfr(+) Cdk1(-)
# mid activated           :Egfr(+) Cdk1(+) Dlx2(lo)
# late activated          :Egfr(+) Cdk1(+) Dlx2(hi)
# NPC-like                :Dlx2(+) Dck(+)
# Reports unable to distinguish B cells/astrocytes but have FACS info to gate
# for pseudotime w/o astrocytes.

# neuroblasts (A cells)             : Dcx
# transit amplifying (C cells)      : Mki67, 
# activated NSC, neurogenic (A )    : Ascl1
# quiescent, early (B cells)        : Notch2
# B cell, astrocyte

FeaturePlot(svz, features = c('Egfr','Cdk1','Dlx2','Dck'))
FeaturePlot(svz, features = c('Ascl1','Id3','Notch2','Gfap'))

# Gjb6 == Cx30 == Connexin 30
more_genes <- c('P2ry12','Tmem119','Cx3cr1','Csf1r','Nes','Gjb6','Dcx','Thbs4','Ascl1','Egfr')
extended_plots <- FeaturePlot(object = svz, features = more_genes, slot = 'data', combine = FALSE, pt.size = 1.5, order = TRUE)
names(extended_plots) <- more_genes
extended_plots <- lapply(X = seq_along(extended_plots), 
                         tmp_plot = extended_plots,
                         FUN = function(x, tmp_plot) {
                           max.exp <- ceiling(max(tmp_plot[[x]]$data[[4]])*10)/10
                           tmp_plot[[x]] + 
                             geom_text(mapping = aes(x = -7.5, y = 17.5, label = names(tmp_plot)[x]),
                                       fontface = 'italic',
                                       size = 4.5,
                                       hjust = 1) +
                             theme(plot.title = element_blank(), 
                                   panel.border = element_rect(fill = NA, color = 'black'), 
                                   panel.background = element_rect(fill = NA), 
                                   axis.text = element_blank(), 
                                   axis.title = element_blank(), 
                                   axis.ticks = element_blank(), 
                                   legend.key = element_rect(fill = NA),
                                   legend.text = element_text(size = 10),
                                   legend.position = c(0.8, 0.175), 
                                   legend.key.size = unit(0.25, units = 'cm'), 
                                   legend.spacing.x = unit(0.1, units = 'cm')) +
                             # scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11,'Spectral')),
                             #                       breaks = c(0, max.exp),
                             #                       limits = c(0, max.exp)) + 
                             scale_color_gradient2(mid = 'grey80', 
                                                   high = 'red3', 
                                                   breaks = c(0, max.exp),
                                                   limits = c(0, max.exp)) +
                             guides(color = guide_colorbar(barwidth = 0.5, 
                                                           frame.colour = 'black', 
                                                           frame.linewidth = 1,
                                                           ticks.colour = 'black',
                                                           ticks.linewidth = 1,
                                                           barheight = 1.75
                             )
                             )
                         }
)
extended_plots <- cowplot::plot_grid(plotlist = extended_plots, ncol = 5)

# display and save
extended_plots
ggsave(filename = paste0(new_results_path, 'microglia_nsc_canonicalmarkers.tiff'), plot = extended_plots, device = 'tiff', height = 5.25, width = 14)






# Trajectory analysis using Monocle3 --------------------------------------

# Import data (start with NSCs)
nsc <- readRDS(file = paste0(data_path, '20200625_nsc.rds'))

# installation
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

# libraries
library('monocle3')
library('Seurat')
library('dplyr')
library('ggplot2')


# Convert object type
nsc <- as.SingleCellExperiment(x = nsc)
rowData(nsc)[,'gene_short_name'] <- rownames(nsc)
nsc2 <- new_cell_data_set(expression_data = counts(nsc),
                          cell_metadata = colData(nsc),
                          gene_metadata = rowData(nsc))

# normalize, variable genes, pca, 
nsc2 <- preprocess_cds(cds = nsc2, num_dim = 5)
plot_pc_variance_explained(nsc2)

# nonlinear dimreduce
nsc2 <- reduce_dimension(nsc2, preprocess_method = 'PCA', )

# plotting
umap_nsc <- plot_cells(cds = nsc2, 
                       color_cells_by = 'RNA_snn_res.0.8',
                       cell_size = 2,
                       group_label_size = 8)
umap_nsc <- umap_nsc + theme(axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             axis.title = element_text(size = 16, color = 'black'))
umap_nsc
ggsave(filename = paste0(new_results_path, 'nsc_umap_alternative.tiff'),
       plot = umap_nsc, device = 'tiff', height = 5, width = 5.5)

# plot by gene
plot_cells(nsc2, genes = c('Egfr','Cdk1','Dlx2','Dck'), cell_size = 2) + 
  theme(strip.text = element_text(size = 18, color = 'black'),
        axis.text = element_blank(),
        panel.border = element_rect(size = 1, fill = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) + 
  guides(color = guide_colorbar(barwidth = 1, 
                                barheight = 6,
                                frame.colour = 'black', 
                                frame.linewidth = 1,
                                ticks.colour = 'black',
                                ticks.linewidth = 1))


# cluster cells
nsc2 <- cluster_cells(nsc2, reduction_method = 'UMAP', k = 15, resolution = 0.1, verbose = TRUE) # force finer clustering
plot_cells(nsc2, color_cells_by = 'cluster', cell_size = 2)
nsc2 <- cluster_cells(nsc2, reduction_method = 'UMAP')
plot_cells(nsc2, color_cells_by = 'cluster', cell_size = 2)


# import previous seurat result to learn path
nsc <- readRDS(file = paste0(data_path, '20200625_nsc.rds'))
tmp_umap <- nsc@reductions$umap@cell.embeddings
reducedDims(nsc2)$UMAP <- tmp_umap
nsc2 <- cluster_cells(nsc2, reduction_method = 'UMAP', k = 15, resolution = 0.)
nsc2 <- learn_graph(nsc2, use_partition = TRUE)
prelim_result <- plot_cells(nsc2,
                            color_cells_by = 'cluster',
                            label_branch_points = FALSE,
                            label_leaves = FALSE,
                            label_groups_by_cluster = FALSE,
                            cell_size = 2) +
  theme(strip.text = element_text(size = 18, color = 'black'),
        axis.text = element_blank(),
        panel.border = element_rect(size = 1, fill = NA),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))
prelim_result



# Trajectory analysis using Slingshot -------------------------------------




nsc2_nsc3_deg <- FindMarkers(svz, ident.1 = 'nsc2', ident.2 = 'nsc3', logfc.threshold = 0.5)
tmp <- FindMarkers(svz, ident.1 = 'nsc2', ident.2 = 'nsc3', logfc.threshold = 0.5)


mg_markers <- readRDS(file = paste0(results_path, '20200625_subclustering_microglia/microglia_de_markers.rds'))
nsc_markers <- readRDS(file = paste0(results_path, '20200625_subclustering_nsc/nsc_de_markers.rds'))

DimPlot(svz, group.by = 'orig.ident', pt.size = 2)

DimPlot(svz, label = TRUE, pt.size = 2)
DimPlot(mg, label = TRUE, pt.size = 2)
DimPlot(nsc, label = TRUE, pt.size = 2)

mg$subcluster <- mg$RNA_snn_res.0.8
nsc$subcluster <- nsc$RNA_snn_res.0.8



svz$
Idents(svz) <- svz$seurat_clusters

FeaturePlot(svz, 'Id3', cols = c('grey85','red'), pt.size = 2, order = TRUE)





# SCENIC ------------------------------------------------------------------

# Following tutorial here: https://rawcdn.githack.com/aertslab/SCENIC/701cc7cc4ac762b91479b3bd2eaf5ad5661dd8c2/inst/doc/SCENIC_Setup.html

## Required
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) 

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "pheatmap", "R2HTML", "Rtsne"))


if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC")


# Testing SCENIC using scRNAseq data derived from subventricular zone. Cells 
# were enriched for Nestin(+) expression using cytometry from mice that receive 
# a photothrombotic ischemic injury. Cells were collected before, 1 day, and 7 
# days post-injury (SCtrl, S1d, S7d).
