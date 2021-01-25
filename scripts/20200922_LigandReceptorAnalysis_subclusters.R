

######## 20200922 Ligand-receptor co-expression analysis (subclusters) ########

getwd()

# libraries
library('Seurat')
library('ggplot2')
library('dplyr')
library('ggrepel')
source('scripts/LRfunctions.R')
source('scripts/VolcanoPlot.R')


# Directories and data
data_path <- 'data/'
results_path <- 'results/'
new_results_path <- paste0(results_path, '20200922_LigandReceptorAnalysis_subclusters/')
lr_ref <- 'ref/fantom_PairsLigRec_mouse.csv'
lr_ref <- read.csv(file = lr_ref)
dir.create(path = new_results_path)


# Import
svz <- readRDS(file = paste0(data_path, '20200918_svz.rds'))
mg <- readRDS(file = paste0(data_path, '20200625_microglia.rds'))
nsc <- readRDS(file = paste0(data_path, '20200625_nsc.rds'))


# Set celltype/subcluster across datasets
mg$subcluster <- plyr::mapvalues(x = mg$RNA_snn_res.0.8,
                                 from = c(0,1,2,3),
                                 to = c('H0-MG','H1-MG','Activated MG', 'Unknown MG'))
nsc$subcluster <- plyr::mapvalues(x = nsc$RNA_snn_res.0.8,
                                  from = c(0,1,2,3),
                                  to = c('Neuroblast','B cell', 'Prolif. C cell', 'Activated C cell'))
subcluster_labels <- c(as.character(mg$subcluster), as.character(nsc$subcluster))
names(subcluster_labels) <- c(names(mg$subcluster), names(nsc$subcluster))
svz$subcluster <- colnames(x = svz)
svz$subcluster <- plyr::mapvalues(x = svz$subcluster,
                                  from = names(subcluster_labels),
                                  to = subcluster_labels)
svz$subcluster <- ifelse(svz$subcluster %in% colnames(svz), yes = 'Unknown', no = svz$subcluster)
svz$subcluster <- factor(svz$subcluster,
                         levels = c('H0-MG','H1-MG','Activated MG', 'Unknown MG','Neuroblast','B cell', 'Prolif. C cell', 'Activated C cell', 'Unknown'))
DimPlot(svz, group.by = 'subcluster', label = TRUE, label.size = 6, pt.size = 2, repel = TRUE)


# Calculate LR scores
Idents(svz) <- 'subcluster'
DefaultAssay(svz) <- 'RNA'
my_setup <- setupLR(seurat_object = svz, lr_ref = lr_ref, split_by = 'orig.ident')
my_results <- calculateLR(setup = my_setup)
write.csv(x = my_results, file = paste0(results_path, 'LR_results_subcluster.csv'), quote = FALSE)


# Identify LR genes that are also DE between mg2 (S1d mg) vs mg0 (SCtrl mg).
# This line will give you the 62 genes that are 1) DE between mg0 and mg2, and
# 2) included in the ligand-receptor database.
de_lr <- FindMarkers(svz, ident.1 = 'mg2', ident.2 = 'mg0', features = my_setup$lr_genes) %>%
  filter(p_val < 0.05)


# View specific results
subset_lr <- my_results %>%
  filter(Ligand == 'P2ry12') %>%
  # filter(Ligand %in% my_setup$lr_genes) %>%
  # filter()
  filter(Ligand_cell == 'mg2') %>%
  filter(pval < 0.05) 
View(subset_lr)

# plot specific ligand-receptor pairs
# e.g. all interactions with receptor P2ry12
plotLR(results = my_results, receptors = 'P2ry12', pval_threshold = 0.05)


# Generating the LR plots for each pair.
p1 <- plotLR(results = my_results, l_cells = c('mg0'), r_cells = c('nsc0'), pval_threshold = 0.05)
p2 <- plotLR(results = my_results, l_cells = c('mg0'), r_cells = c('nsc1'), pval_threshold = 0.05)
p3 <- plotLR(results = my_results, l_cells = c('mg0'), r_cells = c('nsc2'), pval_threshold = 0.05)
p4 <- plotLR(results = my_results, l_cells = c('mg0'), r_cells = c('nsc3'), pval_threshold = 0.05)
mg0_plot <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)

p1 <- plotLR(results = my_results, l_cells = c('mg2'), r_cells = c('nsc0'), pval_threshold = 0.05)
p2 <- plotLR(results = my_results, l_cells = c('mg2'), r_cells = c('nsc1'), pval_threshold = 0.05)
p3 <- plotLR(results = my_results, l_cells = c('mg2'), r_cells = c('nsc2'), pval_threshold = 0.05)
p4 <- plotLR(results = my_results, l_cells = c('mg2'), r_cells = c('nsc3'), pval_threshold = 0.05)
mg2_plot <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)



# LR results and observations ---------------------------------------------

# 1) IGF2 is down after injury and does not return to pre-injury levels 
# Refs:
#   https://pubmed.ncbi.nlm.nih.gov/31484925/
#   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3338187/
genes <- c('Igf2', 'Igf1', 'Insr','Igf1r','Igf2r')
Idents(svz) <- 'subcluster'
p1 <- VlnPlot(svz, features = genes, pt.size = 0.5, ncol = 2)
p2 <- VlnPlot(svz, features = genes, pt.size = 0.5, ncol = 2, group.by = 'celltype', split.by = 'orig.ident')
tmp_plot <- cowplot::plot_grid(p1, p2, ncol = 2)
ggsave(filename = paste0(new_results_path, 'Igf_family_vln.tiff'), plot = p, height = 6, width = 9, device = 'tiff')


# 2) ???



