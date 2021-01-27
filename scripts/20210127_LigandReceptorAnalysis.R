
#' ---
#' title: Ligand-Receptor Analysis of SVZ after PT
#' author: James Choi
#' date: `r Sys.Date()`
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)


require('Seurat')
require('dplyr')
require('ggplot2')
source('20210125_LRfunctions.R')

setwd(dir = './scripts/')
results_out <- '../results/20210127_LigandReceptorAnalysis/'
dir.create(results_out)
svz <- readRDS(file = '../data/20210108_SVZ.rds')
nsc <- readRDS(file = '../data/20210108_NSC.rds')
mg <- readRDS(file = '../data/20210112_MG.rds')
lr_ref <- read.csv(file = '../ref/fantom_PairsLigRec_mouse.csv')

DimPlot(svz, label = TRUE)

nsc$subtype <- paste('NSC', nsc$SCT_snn_res.0.8, sep = '.')
mg$subtype <- paste('MG', mg$SCT_snn_res.0.4, sep = '.')
unknown <- as.character(svz$celltype[svz$celltype == 'Unknown/Ependymal'])
names(unknown) <- colnames(svz)[svz$celltype == 'Unknown/Ependymal']

subtype <- c(nsc$subtype,
             mg$subtype,
             unknown)
svz$subtype <- plyr::mapvalues(
  x = colnames(svz),
  from = names(subtype),
  to = subtype
)

#+ celltype_subtype_umap, fig.height=3.5, fig.width=8.5, fig.cap='Left: UMAP of cells by cell-type. Right: UMAP of cells by sub-type'
p1 <- DimPlot(svz, label = TRUE, label.size = 4, group.by = 'celltype') + 
  theme_bw()
p2 <- DimPlot(svz, label = TRUE, label.size = 4, group.by = 'subtype') + 
  theme_bw()
p1 + p2


#' Do ligand-receptor analysis of interactions between cell-types at different 
#' time-points.  
Idents(svz) <- 'celltype'
svz_setup <- setupLR(
  seurat_object = svz,
  lr_ref = lr_ref,
  split_by = 'orig.ident'
)
svz_results <- calculateLR(
  setup = svz_setup,
  resample = 10000,
  adjust_pval = TRUE
)
write.table(x = svz_results,
            file = paste0(results_out, 'SVZ_celltype_LRanalysis.csv'),
            sep = ',',
            quote = FALSE,
            row.names = FALSE)

#' Do ligand-receptor analysis of interactions between sub-types at different
#' time-points.  
Idents(svz) <- 'subtype'
svz_setup <- setupLR(
  seurat_object = svz,
  lr_ref = lr_ref,
  split_by = 'orig.ident'
)
svz_results <- calculateLR(
  setup = svz_setup,
  resample = 10000,
  adjust_pval = TRUE
)
write.table(x = svz_results,
            file = paste0(results_out, 'SVZ_subtype_LRanalysis.csv'),
            sep = ',',
            quote = FALSE,
            row.names = FALSE)



# LR visualization --------------------------------------------------------

celltype_results <- read.table(
  file = paste0(results_out, 'SVZ_celltype_LRanalysis.csv'),
  sep = ',',
  header = TRUE
)
celltype_results$split_by <- factor(celltype_results$split_by, 
                                    levels = c('SCtrl','S1d','S7d'))
subtype_results <- read.table(
  file = paste0(results_out, 'SVZ_subtype_LRanalysis.csv'),
  sep = ',',
  header = TRUE
)
subtype_results$split_by <- factor(subtype_results$split_by, 
                                    levels = c('SCtrl','S1d','S7d'))

# # Visual some of the pathways highlighted by Suvra:
# -B2m-TFrc
# -C1qa-Cr1l
# -Apoe-Lrp8
# -Cd14-Itgb1
# -Lgals3bp-Itgb1
# -Pros1-Tyro3
# -Igf2-Igf1r

celltype_plot <- plotLR(
  results = celltype_results,
  ligands = c('B2m','C1qa','Apoe','Cd14','Lgals3bp','Pros1','Igf2'),
  receptors = c('Tfrc','Cr1l','Lrp8','Itgb1','Tyro3','Igf1r'),
  l_cells = c('Microglia','NSC'),
  r_cells = c('Microglia','NSC'),
  resample = 10000
)
celltype_plot
ggsave(filename = paste0(results_out, 'celltype_selectPairs_LRplot.tiff'),
       plot = celltype_plot, device = 'tiff', height = 4.25, width = 7.5)
ggsave(filename = paste0(results_out, 'celltype_selectPairs_LRplot.svg'),
       plot = celltype_plot, device = 'svg', height = 4.25, width = 7.5)


subtype_plot <- plotLR(
  results = svz_results,
  ligands = c('B2m','C1qa','Apoe','Cd14','Lgals3bp','Pros1','Igf2'),
  receptors = c('Tfrc','Cr1l','Lrp8','Itgb1','Tyro3','Igf1r'),
  resample = 10000
)
subtype_plot
ggsave(filename = paste0(results_out, 'subtype_selectPairs_LRplot.tiff'),
       plot = subtype_plot, device = 'tiff', height = 21, width = 8)
ggsave(filename = paste0(results_out, 'subtype_selectPairs_LRplot.svg'),
       plot = subtype_plot, device = 'svg', height = 21, width = 8)
