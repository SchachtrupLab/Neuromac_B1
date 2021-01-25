
library('Seurat')
library('ggplot2')

svz <- readRDS(file = 'data/20200918_svz.rds')

# Get tsne coordinates and cell-type labels
tsne_coord <- svz[['tsne']]@cell.embeddings
svz_celltype <- svz@meta.data[['celltype']]
svz_data <- data.frame(
  tsne_coord,
  svz_celltype
)

# Plot
svz_tsne <- ggplot(data = svz_data,
                   mapping = aes(x = tSNE_1, y = tSNE_2)) +
  geom_point(mapping = aes(fill = svz_celltype),
             size = 4,
             shape = 21,
             color = 'black') +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 16, color = 'black'),
        legend.text = element_text(size = 14, color = 'black'),
        legend.key = element_rect(fill = NA)) +
  guides(fill = guide_legend(title = 'Cell-type',
                             override.aes = list(size = 4)))

# SVG files can be written using base R (ie it's built-in)
svg(filename = 'testSVG_SVZtSNE.svg', width = 6.75, height = 5)
svz_tsne
dev.off()
