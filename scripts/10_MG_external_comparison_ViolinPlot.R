
#### External MG comparison (2021-10-11) -------------------

require('SingleCellExperiment')
require('Seurat')
require('ggplot2')
require('dplyr')

results_out <- '../results/5_MG_External_Comparison/'


# Import full Zywitza dataset
zyw <- readRDS(file = paste0(results_out, 'Zywitza_SVZ.rds'))

# Extract Zywitza MG
zyw_mg <- zyw[,zyw$RNA_snn_res.0.4 == 4]
zyw_mg$celltype <- 'Zyw_MG'
Idents(zyw_mg) <- 'celltype'

# Import SVZ MG
mg <- readRDS(file = '../data/mg.rds')
DefaultAssay(mg) <- 'RNA'


### 2021-10-12 response to sequencing depth email ----------
p1 <- VlnPlot(zyw_mg, features = 'nFeature_RNA', pt.size = 0) +
  scale_y_continuous(limits = c(0, 4000),
                     breaks = seq(0, 4000, 1000)) +
  ggtitle(label = 'Zywitza MG') +
  ylab(label = 'Unique genes detected per cell') +
  xlab(label = '')
p2 <- VlnPlot(mg, features = 'nFeature_RNA', group.by = 'orig.ident', 
              pt.size = 0) +
  ggtitle(label = 'SVZ MG (current study)') +
  scale_y_continuous(limits = c(0, 4000),
                     breaks = seq(0, 4000, 1000)) +
  ylab(label = 'Unique genes detected per cell') +
  xlab(label = '')
p3 <- (p1 | p2) & patchwork::plot_layout(widths = c(0.4, 1))
ggsave(filename = paste0(results_out, 'Zywitza_sequencingDepthComparison.tiff'),
       plot = p3, height = 4, width = 6, device = 'tiff')
# ---------------------------------------------------------


# Find shared genes between datasets
shared_genes <- rownames(zyw_mg)[rownames(zyw_mg) %in% rownames(mg)]
zyw_mg <- zyw_mg[shared_genes,]
mg <- mg[shared_genes,]

# Remove development + unknown
mg <- mg[, !mg$subtype %in% c('Development-associated', 'Unknown')]

# Chunk for generating combined dataset + UMAP
zmg <- list('mg' = mg,
            'zyw_mg' = zyw_mg)
zmg <- lapply(
  X = zmg,
  FUN = FindVariableFeatures,
  nfeatures = 500,
  verbose = FALSE
)
zmg_anchors <- FindIntegrationAnchors(
  object.list = zmg,
  verbose = FALSE
)
zmg <- IntegrateData(
  anchorset = zmg_anchors,
  verbose = FALSE
)
DefaultAssay(zmg) <- 'integrated'
zmg <- ScaleData(zmg, verbose = FALSE)
zmg <- RunPCA(zmg, npcs = 30, verbose = FALSE)
ElbowPlot(zmg, ndims = 30)
zmg <- FindNeighbors(zmg, dims = 1:8, verbose = FALSE)
zmg <- RunUMAP(zmg, dims = 1:8, verbose = FALSE)
zmg <- FindClusters(zmg, resolution = 0.4)

# UMAP of combined data
DimPlot(zmg, group.by = 'subtype', shuffle = TRUE, pt.size = 2)

# Setting up microglia subtype order
zmg$subtype <- plyr::mapvalues(
  zmg$subtype,
  from = NA,
  to = 'Zywitza'
)
zmg$subtype <- factor(x = zmg$subtype,
                      levels = c('Homeostatic',
                                 'Activated',
                                 'Disease-associated',
                                 'Zywitza'))

# Violin plot for comparing specific genes
DefaultAssay(zmg) <- 'RNA'
genes <- c('P2ry12','Tmem119','Hexb','Cx3cr1', # homeostatic
           'Apoe','Axl','Cd52','Lgals3bp', # DAM
           'Srgn','Tlr7','Socs3','Id2') # activated
VlnPlot(zmg, features = genes, ncol = 4, pt.size = 0,
        group.by = 'subtype') &
  theme(axis.title = element_blank())


# Compact violin plot
DefaultAssay(zmg) <- 'RNA'
genes <- c('P2ry12','Tmem119','Hexb','Cx3cr1', # homeostatic
           'Apoe','Axl','Cd52','Lgals3bp',
           'Srgn','Tlr7','Socs3','Id2') 
compact_vln <- FetchData(object = zmg, 
                 vars = c(genes, 'subtype'),
                 slot = 'data') %>% 
  reshape2::melt(id.vars = c('subtype')) %>% 
  ggplot(mapping = aes(x = subtype, y = value)) +
  geom_violin(mapping = aes(fill = subtype), scale = 'width') +
  facet_wrap(. ~ variable, scales = 'free_y', ncol = 1,
             strip.position = 'left') +
  theme(axis.text.y.left = element_blank(),
        panel.background = element_rect(color = 'black', fill = NA),
        strip.text = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA),
        panel.spacing.x=unit(0, "lines"),panel.spacing.y=unit(0, "lines"))
compact_vln
ggsave(filename = paste0(results_out, 'Zywitza_MG-gene-comparison_vln.tiff'),
       plot = compact_vln, height = 4.5, width = 3, device = 'tiff')
ggsave(filename = paste0(results_out, 'Zywitza_MG-gene-comparison_vln.svg'),
       plot = compact_vln, height = 4.5, width = 3, device = 'svg')
