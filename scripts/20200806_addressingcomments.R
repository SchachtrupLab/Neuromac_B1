
####### Addressing comments from Suvra's progress meeting #######

library('Seurat')
library('dplyr')
library('ggplot2')

svz <- readRDS(file = 'data/20200511_SVZ.rds')
microglia <- readRDS(file = 'data/20200625_microglia.rds')
nsc <- readRDS(file = 'data/20200625_nsc.rds')

results_path <- 'results/20200806_address_comments/'
dir.create(path = results_path)


# Need for figure clearly defining cell-types ----------------------------------

#testing
DimPlot(microglia, pt.size = 2, label = TRUE)
FeaturePlot(nsc, features = c('Mki67'), pt.size = 2) + 
  scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 9, name = 'Spectral')))

# Classify cell-types
svz$orig_clusters <- svz$RNA_snn_res.0.8
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
Idents(svz) <- 'celltype'
DimPlot(svz, label = TRUE, label.size = 6)

# Rename microglia subcluster to better characterize.
microglia$subcluster <- microglia$seurat_clusters
microglia$subcluster <- plyr::mapvalues(x = microglia$subcluster,
                                        from = c(0:3),
                                        to = c('Homeostatic\nmicroglia',
                                               'Unknown\nmicroglia',
                                               'Inflammatory\nmicroglia',
                                               'Development\nmicroglia'))
Idents(microglia) <- 'subcluster'

# Rename NSC subclusters to better characterize.
nsc$subcluster <- nsc$seurat_clusters
nsc$subcluster <- plyr::mapvalues(x = nsc$subcluster,
                                  from = c(0:3),
                                  to = c('Neuroblast',
                                         'C-cell',
                                         'Activated\nB-cell',
                                         'B-cell'))
Idents(nsc) <- 'subcluster'


# Map cell-types onto UMAP with all cells
joined_subcluster <- c(colnames(microglia), colnames(nsc), colnames(svz)[svz$celltype == 'Unknown'])
names(joined_subcluster) <- c(as.character(microglia$subcluster), as.character(nsc$subcluster), rep('Unknown', sum(svz$celltype == 'Unknown')))
svz$subcluster <- colnames(svz)
svz$subcluster <- plyr::mapvalues(x = svz$subcluster,
                                  from = joined_subcluster,
                                  to = names(joined_subcluster))
Idents(svz) <- 'subcluster'
labelled_umap <- DimPlot(svz, label = TRUE, label.size = 6, repel = TRUE, pt.size = 2) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        axis.ticks = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.text = element_blank())
ggsave(filename = paste0(results_path, 'celltype_labelled_umap.tiff'), plot = labelled_umap,
       device = 'tiff', height = 6, width = 7.5)


# Here insert plotting marker genes





# Concerns with displaying gene expression -------------------------------------


# Checking for double positive gene detection
bcl2_dcx_double <- FetchData(object = svz, vars = c('subcluster', 'Bcl2', 'Dcx', 'orig.ident'), slot = 'data') %>%
  group_by(subcluster, orig.ident) %>%
  summarise('Bcl2(+) / Dcx(-)' = round(mean(Bcl2 > 0 & Dcx == 0) * 100, 1),
            'Bcl2(-) / Dcx(+)' = round(mean(Dcx > 0 & Bcl2 == 0) * 100, 1),
            'Bcl2(+) / Dcx(+)' = round(mean(c(Bcl2 > 0 & Dcx > 0)) * 100, 1)) %>%
  reshape2::melt(is.vars = c('subcluster', 'orig.ident')) %>%
  ggplot(mapping = aes(x = orig.ident, y = value)) +
  geom_bar(mapping = aes(fill = variable), stat = 'identity', color = 'black', position = 'stack') +
  facet_grid(. ~ subcluster) +
  xlab(label = 'Days after injury') +
  ylab(label = 'Percent detection') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        strip.background = element_rect(fill = NA, color = 'black'),
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        strip.text = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 14, color = 'black'),
        panel.grid = element_line(color = 'grey80', linetype = 'dashed')) +
  guides(fill = guide_legend(title = ''))
bcl2_dcx_double

bcl2_dcx_counts <- FetchData(object = svz, vars = c('subcluster', 'Bcl2', 'Dcx', 'orig.ident'), slot = 'data') %>%
  group_by(subcluster, orig.ident) %>%
  summarise('Bcl2(+) / Dcx(-)' = sum(Bcl2 > 0 & Dcx == 0),
            'Bcl2(-) / Dcx(+)' = sum(Dcx > 0 & Bcl2 == 0),
            'Bcl2(+) / Dcx(+)' = sum(Bcl2 > 0 & Dcx > 0)) %>%
  reshape2::melt(is.vars = c('subcluster', 'orig.ident')) %>%
  ggplot(mapping = aes(x = orig.ident)) +
  geom_bar(mapping = aes(y = value, fill = variable), stat = 'identity', color = 'black', position = 'stack') +
  # geom_text(mapping = aes(label = value, y = value), position = 'stack') +
  facet_grid(. ~ subcluster) +
  xlab(label = 'Days after injury') +
  ylab(label = 'Number of cells') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        strip.background = element_rect(fill = NA, color = 'black'),
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        strip.text = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 14, color = 'black')) +
  guides(fill = guide_legend(title = ''))
bcl2_dcx_counts 

ggsave(filename = paste0(results_path, 'bcl2_dcx_double_positive_percent.tiff'), plot = bcl2_dcx_double,
       device = 'tiff', height = 3, width = 15)




