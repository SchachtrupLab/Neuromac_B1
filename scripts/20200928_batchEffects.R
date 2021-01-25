
######## 20200928 Testing for plate-based batch effects ########

library('Seurat')
library('ggplot2')
library('dplyr')
svz <- readRDS(file = 'data/20200918_svz.rds')

# Extract unique barcodes
id <- colnames(svz)
# Split to extract Tx group and plate #
id <- strsplit(x = id, split = '_')
grp <- sapply(X = id, FUN = `[`, 1)
plate <- sapply(X = id, FUN = `[`, 2)
# Add metadata and plot
svz$plate <- plate
DimPlot(svz, group.by = 'plate')
# Check for subcluster biases across plates
tmp <- as.data.frame(table(svz$subcluster, svz$plate))
colnames(tmp) <- c('subcluster','Plate','n_cells')  
tmp_plot <- tmp %>%
  ggplot(mapping = aes(x = subcluster, y = n_cells)) + 
  geom_bar(mapping = aes(fill = Plate), stat = 'identity', position = 'stack', color = 'black') +
  # geom_text(mapping = aes(label = n_cells)) +
  labs(title = 'Detecting plate-based batch effects') +
  ylab(label = '# of cells') +
  xlab(label = 'Cell-type') +
  scale_y_continuous(breaks = seq(0, 300, 25)) +
  theme(panel.background = element_rect(color = 'black', fill = NA),
        panel.border = element_rect(color = 'black', fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = 'black'),
        plot.title = element_text(size = 16, color = 'black'),
        axis.text.y = element_text(size = 14, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 12, color = 'black'))
  
# Save results
dir.create(path = 'results/20200928_batchEffects')
ggsave(filename = 'results/20200928_batchEffects/subcluster_plate_composition.tiff', device = 'tiff', height = 3.5, width = 5)


svz <- svz %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = 'plate') %>%
  RunPCA()
svz <- svz %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:10)


mg$plate <- plyr::mapvalues(x = colnames(mg),
                            from = colnames(svz),
                            to = svz$plate)
mg <- mg %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = 'plate') %>%
  RunPCA()
mg <- mg %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:10)
plate_mg_umap <- DimPlot(mg, group.by = 'plate', split.by = 'orig.ident', pt.size = 3)
ggsave(filename = 'results/20200928_batchEffects/plate_mg_umap.tiff', device = 'tiff', height = 3.5, width = 8.5)
