
#' ---
#' title: Trajectory Analysis of NSC
#' author: James Choi
#' date: 2021/1/20
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)

require('Seurat')
require('ggplot2')
require('dplyr')
require('dendextend')
# BiocManager::install('slingshot')
require('slingshot')
require('SingleCellExperiment')

nsc <- readRDS(file = '../data/20210108_NSC.rds')

#+ subtype_umap, fig.height=7, fig.width=7
DefaultAssay(nsc) <- 'SCT'
nsc$subtype <- plyr::mapvalues(
  x = nsc$SCT_snn_res.0.8,
  from = 0:3,
  to = c('Neuroblast','B_cell','Activated_C_cell', 'Prolif_C_cell')
)
Idents(nsc) <- 'subtype'
p1 <- DimPlot(nsc, pt.size = 2, label = TRUE, label.size = 6, repel = TRUE) + 
  theme_bw() + NoLegend()
p2 <- DimPlot(nsc, pt.size = 1.5, split.by = 'orig.ident') + theme_bw() +
  NoLegend()
cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(1,0.6))


# Load PCA coordinates + subclusters (original subcluster analysis used 8 PCs)
nsc_sce <- as.SingleCellExperiment(nsc, assay = 'RNA')
nsc_pca <- slot(slot(nsc, 'reductions')[['pca']], 'cell.embeddings')
nsc_umap <- slot(slot(nsc, 'reductions')[['umap']], 'cell.embeddings')
npcs <- 1:8
nsc_pca <- nsc_pca[,npcs]
nsc_umap <- nsc_umap

# Calculate lineages using subcluster (MST built from the 4 subclusters)
nsc_cols <- c('B_cell' = 'red', 'Activated_C_cell' = 'orange', 'Prolif_C_cell' = 'green', 'Neuroblast' = 'blue')
# plot(nsc_umap, col = nsc_cols[nsc_sce$subtype], pch=16, asp = 1)
reducedDim(nsc_sce, type = 'UMAP') <- nsc_umap
nsc_sce <- slingshot(data = nsc_sce, 
                     clusterLabels = 'subtype', 
                     reducedDim = 'UMAP',
                     start.clus = 'B_cell')
nsc_sds <- SlingshotDataSet(nsc_sce)

# par(mar = c(5,5,2,2))
# plot(reducedDims(nsc_sce)$UMAP, 
#      col = nsc_cols[nsc_sce$subtype],
#      pch=16, asp = 1, cex = 2, cex.lab = 2, cex.axis = 0.001)
# lines(SlingshotDataSet(nsc_sce), lwd = 3, col='black')
# par(mar = c(5,5,4,2))


#+ trajectory_umap, fig.height=3.5, fig.width=9
nsc$pseudotime <- nsc_sce$slingPseudotime_1
p1 <- DimPlot(nsc, pt.size = 2) +
  geom_line(data = data.frame(nsc_sds@curves$curve1$s), 
            mapping = aes(x = UMAP_1, y = UMAP_2),
            lwd = 1.5) + 
  theme_bw() +
  labs(title = 'NSC subtype')
p2 <- FeaturePlot(nsc, pt.size = 2, features = 'pseudotime') + 
  scale_color_viridis_c() +
  geom_line(data = data.frame(nsc_sds@curves$curve1$s), 
            mapping = aes(x = UMAP_1, y = UMAP_2),
            lwd = 1.5) + 
  theme_bw() +
  labs(title = 'Pseudotime')
trajectory_umap <- p1 + p2
trajectory_umap

#+ pt_density, fig.height=3.5, fig.width=5.5
plot_dat <- cbind(nsc@meta.data, nsc[['umap']]@cell.embeddings)
pt_density <- plot_dat %>%
  ggplot(mapping = aes(x = pseudotime)) +
  geom_density(mapping = aes(color = orig.ident), adjust = 1, lwd = 1) +
  geom_point(mapping = aes(x = pseudotime, y = 0, fill = subtype), 
             color = 'black', pch = 21, size = 5) +
  theme_bw() +
  labs(title = 'Distribution of cells along pseudotime') +
  ylab(label = 'Density') + 
  xlab(label = 'Pseudotime') +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  guides(fill = guide_legend(title = 'NSC subtype'),
         color = guide_legend(title = 'Injury time-point'))
pt_density

#+ some_genes, fig.height=6, fig.width=6.5
genes <- c('Gfap','Sox2','Ascl1','Dcx')
pt_genes <- FetchData(nsc, vars = c(genes, 'pseudotime','orig.ident')) %>%
  reshape2::melt(id.vars = c('pseudotime','orig.ident')) %>%
  ggplot(mapping = aes(x = pseudotime, y = value)) +
  geom_point(size = 1) +
  geom_smooth(lwd = 1.5) +
  facet_wrap(. ~ variable, scales = 'free_y') + 
  theme_bw() +
  ylab(label = 'Normalized Expression') +
  theme(strip.text = element_text(size = 12))
pt_genes



# Figures (2021-01-20) ----------------------------------------------------

# Run this first:
nsc$subtype <- plyr::mapvalues(
  x = nsc$SCT_snn_res.0.8,
  from = 0:3,
  to = c('Neuroblast','B_cell','Activated_C_cell', 'Prolif_C_cell')
)
Idents(nsc) <- 'subtype'

# UMAPs of NSCs by time
umap_byTime <- vector(mode = 'list', length = length(levels(nsc$orig.ident)))
for (t in 1:length(levels(nsc$orig.ident))) {
  this_time <- levels(nsc$orig.ident)[t]
  tmp_time <- ifelse(test = nsc$orig.ident == this_time, yes = 1, no = 0)
  umap_dat <- as.data.frame(cbind(nsc[['umap']]@cell.embeddings, tmp_time))
  umap_dat[['tmp_time']] <- as.character(umap_dat[['tmp_time']])
  umap_byTime[[t]] <- umap_dat %>%
    arrange(tmp_time) %>%
    ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(mapping = aes(color = tmp_time), size = 3) +
    scale_color_manual(values = c('0' = 'grey', '1' = 'red2')) +
    # geom_text(mapping = aes(x = 7, y = 6.5, label = this_time), size = 6) +
    theme(panel.background = element_rect(fill = NA, color = 'black'),
          panel.border = element_rect(color = 'black', fill = NA),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none')
}
umap_byTime <- cowplot::plot_grid(plotlist = umap_byTime, ncol = 1)
umap_byTime
ggsave(filename = '../results/NSC_trajectory_byTime.svg', device = 'svg',
       plot = umap_byTime, height = 7, width = 6)

# UMAPs of NSCs by pseudotime
umap_byPseudotime <- FeaturePlot(nsc, features = 'pseudotime', pt.size = 3) + 
  scale_color_viridis_c() +
  geom_line(data = data.frame(nsc_sds@curves$curve1$s),
            mapping = aes(x = UMAP_1, y = UMAP_2),
            lwd = 1.5) +
  labs(title = 'Pseudotime', hjust = 1) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0),
        legend.position = 'bottom',
        legend.justification = 'center',
        legend.direction = 'horizontal',
        legend.text = element_blank()) +
  guides(color = guide_colorbar(ticks = FALSE))
ggsave(filename = '../results/NSC_trajectory_byPseudotime.svg', device = 'svg',
       plot = umap_byPseudotime, height = 3.25, width = 6)

# UMAPs by NSC subtype
umap_byType <- DimPlot(nsc, pt.size = 3) +
  geom_line(data = data.frame(nsc_sds@curves$curve1$s), 
            mapping = aes(x = UMAP_1, y = UMAP_2),
            lwd = 1.5) + 
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0),
        legend.position = 'bottom',
        legend.justification = 'center',
        legend.direction = 'horizontal',
        legend.text = element_text(size = 14)) +
  labs(title = 'NSPC Subtype') +
  guides(color = guide_legend(override.aes = list('size' = 5)))
ggsave(filename = '../results/NSC_trajectory_bySubtype.svg', device = 'svg',
       plot = umap_byType, height = 3.25, width = 6)

# Plotting individual cells in trajectory
umap_byGene <- FeaturePlot(nsc, pt.size = 3, cols = c('grey','red2'),
                           features = c('Dcx','Ascl1','Gfap','Sox2'),
                           order = TRUE, combine = FALSE)
umap_byGene <- lapply(
  X = umap_byGene,
  FUN = function(x) {
    x + 
      scale_color_gradientn(colors = 
                              RColorBrewer::brewer.pal(n = 9, name = 'OrRd')) +
      theme(panel.background = element_rect(fill = NA, color = 'black'),
            panel.border = element_rect(fill = NA, color = 'black'),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.position = 'none',
            legend.text = element_blank()) +
      guides(color = guide_colorbar(ticks = FALSE))
  }
)
umap_byGene <- cowplot::plot_grid(plotlist = umap_byGene, ncol = 2)
ggsave(filename = '../results/NSC_trajectory_byGene.svg', device = 'svg',
       plot = umap_byGene, height = 4, width = 7.5)



# Misc --------------------------------------------------------------------


# BiocManager::install('tradeSeq')
# library('tradeSeq')
# set.seed()
# 
# 
# # Fit tradeSeq models and do Association Test b/w expression ~ pseudotime
# nsc_gam <- fitGAM(counts = counts(nsc_sce), sds = nsc_sds)
# nsc_at <- associationTest(nsc_gam) # test within lineage, exp ~ pseudotime
# nsc_diffEnd <- diffEndTest(nsc_gam) # test b/w lineages, exp ~ pseudotime(late) + lineage
# nsc_earlyDE <- earlyDETest(nsc_gam) # test b/w lineages, exp ~ pseudotime(early) + lineage
# nsc_patternDE <- patternTest(nsc_gam)
# 
# 
# # Visualize some results
# tmp_theme <- my_theme + theme(legend.title = element_text(size = 14, color = 'black'),
#                               legend.text = element_text(size = 12, color = 'black'),
#                               plot.title = element_text(size = 16, color = 'black', face = 'bold'))
# my_genes <- c('Dcx','Fxyd1','Fgfr1','Sox11','Ncam1', 'Ntrk2','Mt1','Mt2','Plpp3','Tubb3')
# smooth_plots <- vector(mode = 'list', length = length(my_genes))
# for(i in 1:length(my_genes)) {
#   smooth_plots[[i]] <- plotSmoothers(gene = my_genes[i],
#                                      nsc_gam, counts = counts(nsc_sce),
#                                      lwd = 2, size = 1.5, border = TRUE) + 
#     tmp_theme +
#     labs(title = my_genes[i])
# }
# smooth_plots <- cowplot::plot_grid(plotlist = smooth_plots, ncol = 5)
# smooth_plots
# ggsave(filename = paste0(new_results_path, 'geneCurves_nsc_slingshotLineages.tiff'),
#        plot = smooth_plots, device = 'tiff', height = 5.5, width = 16)
# 
# 
# 
# # visualize top genes by p-val
# topgenes <- rownames(nsc_patternDE[order(nsc_patternDE$pvalue), ])[1:250]
# pst_order <- order(nsc_sce$slingPseudotime_1, na.last = NA)
# heat_data <- t(scale(t(logcounts(nsc_sce)[topgenes, pst_order])))
# heat_clus <- nsc_sce$subcluster[pst_order]
# 
# 
# # heat_data <- nsc[['RNA']]@scale.data[topgenes, pst_order]
# data.frame(t(heat_data), 'subcluster' = heat_clus) %>%
#   tibble::rownames_to_column(var = 'barcode') %>%
#   reshape2::melt(id.vars = c('barcode', 'subcluster')) %>%
#   ggplot() + 
#   geom_raster(mapping = aes(x = barcode, y = variable, fill = value))
# heatmap(t(scale(t(as.matrix(heat_data)))), Colv = NA,
#         ColSideColors = RColorBrewer::brewer.pal(9,"Set1")[heat_clus],
#         col = RColorBrewer::brewer.pal(n = 9, name = 'Spectral'))
# 
# 
# # Cluster the genes by similar temporal pattern
# library('clusterExperiment')
# npoints <- 20
# patterns <- clusterExpressionPatterns(nsc_gam, nPoints = npoints, genes = VariableFeatures(nsc), nReducedDims = 5)
# clusterLabels <- primaryCluster()