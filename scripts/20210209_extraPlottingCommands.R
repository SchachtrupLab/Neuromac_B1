DimPlot(mg, pt.size = 2, label = TRUE, label.size = 6, split.by = 'orig.ident') + theme_bw()


# Saving images in SVG ----------------------------------------------------

DefaultAssay(mg) <- 'RNA'
p1 <- FeaturePlot(
  object = mg, 
  features = c('Lamp1','Apoe','Cd63'),
  cols = c('grey','red3'),
  pt.size = 2
)
ggsave(filename = '../results/test.svg', plot = p1, device = 'svg', 
       height = 5, width = 6.5)



# Sample violin plots -----------------------------------------------------


VlnPlot(
  object = mg,
  features = c('Lamp1','Apoe','Cd63'),
  pt.size = 0,
  group.by = 'orig.ident'
)

VlnPlot(
  object = mg,
  features = c('Lamp1','Apoe','Cd63'),
  pt.size = 0.5,
  split.by = 'orig.ident'
)



# Differential expression testing -----------------------------------------

# Question: Is APOE differentially expressed by cluster 4 cells at 7D versus all
# other microglia from other time-points combined?
mg$tmp_ident <- paste(mg$SCT_snn_res.0.4, mg$orig.ident, sep = '_')
Idents(mg) <- 'tmp_ident'
de_test <- FindMarkers(
  object = mg,
  features = 'Apoe',
  ident.1 = '4_S7d'
)

# Checking double positive gene expression rate ------------------------------

DefaultAssay(mg) <- 'RNA'
tmp <- FetchData(
  object = mg,
  vars = c('Cx3cr1','Tmem119','orig.ident'), # change this
  slot = 'counts'
) %>%
  mutate('double_pos' = Cx3cr1 > 0 & Tmem119 > 0) # change this too

prop.table(table(tmp$orig.ident, tmp$double_pos), margin = 1) * 100




# Combine multiple marker genes -------------------------------------------

Idents(mg) <- 'SCT_snn_res.0.4'
DimPlot(mg, pt.size = 2, label = TRUE, label.size = 6)

# Or choose your own genes of interest:
my_genes <- c('P2ry12','Tmem119','Hexb','Cx3cr1')
my_dat <- Matrix::colSums(mg[['RNA']]@data[my_genes,])
my_dat <- data.frame(cbind(my_dat, mg[['umap']]@cell.embeddings))
my_umap <- my_dat %>%
  arrange(my_dat) %>%
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(mapping = aes(fill = my_dat),
             color = 'black',
             pch = 21,
             size = 4) +
  labs(title = paste(my_genes, collapse = ', ')) + 
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
my_umap

ggsave(filename = paste0('../results/test.svg'),
       plot = my_umap, device = 'svg', height = 3, width = 3.5)