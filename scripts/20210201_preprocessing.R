
#' ---
#' title: "SVZ scRNAseq preprocessing"
#' author: "James Choi"
#' date: "Last run: `r Sys.Date()`"
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)


# Load libraries
require('Seurat')
require('dplyr')
require('ggplot2')
require('SingleCellExperiment')
set.seed(123)

# Directories for data and results
# setwd('./scripts/')
results_out <- '../results/preprocessing/'
dir.create(results_out)

# Load counts matrix
svz <- read.table(
  file = '../data/raw_counts_all_cells.csv',
  sep = ',',
  stringsAsFactors = FALSE,
  row.names = 1,
  header = TRUE
)
# Convert to sparse, round to integers (has decimals for some reason)
svz <- Matrix::Matrix(data = as.matrix(svz), sparse = TRUE)
svz <- round(svz)
dat_sce <- SingleCellExperiment(assay = list('counts' = svz))
svz <- CreateSeuratObject(counts = svz, project = 'SVZ', assay = 'RNA')

#' First, we plot some sequencing metrics such as UMIs per cell and number of
#' unique genes detected per cell.
#+ detection_scatter, fig.height=4, fig.width=4.5, fig.cap='Scatter plot of sequencing metrics by cell.'
qc <- svz@meta.data
qc <- qc[sample(x = nrow(qc), size = nrow(qc)),]
plot(
  x = qc$nCount_RNA,
  y = qc$nFeature_RNA,
  xlab = 'total UMI',
  ylab = '# genes detected',
  main = 'Detection rates by sample',
  sub = 'Black = SCtrl; Red = S1d; Green = S7d',
  col = sapply(
    X = as.character(qc$orig.ident),
    FUN = function(x) switch(x, 'SCtrl'='black','S1d'='red3','S7d'='green')
  )
)


#' From the above, we see that there are a few cells with outlier numbers of 
#' total UMI and genes detected. These are likely to be doublets. Beyond that, 
#' there are no obvious biases based on sample of origin.  

#' Next we examine cell-level QC metrics.  
#' 
#' * nCount_RNA = number of UMIs per cell.
#' * nFeature_RNA = number of unique genes detected per cell.
#' * percent.mt = percent of UMIs mapping to a mitochondrial protein-encoding 
#' gene.
#' * log10_total_counts = log10(nCount_RNA)
#+ qc_vln, fig.height=3.5, fig.width=6, fig.cap='Violin plot of various per-cell QC metrics.'
svz <- PercentageFeatureSet(svz, pattern = '^mt-', col.name = 'percent.mt')
svz$log10_total_counts <- log10(svz$nCount_RNA)
qc_plot <- svz@meta.data %>%
  select(-c('log10_total_counts')) %>%
  reshape2::melt(id.vars = 'orig.ident') %>%
  ggplot(mapping = aes(x = orig.ident, y = value)) +
  geom_violin(scale = 'width') +
  ggbeeswarm::geom_quasirandom(size = 1) +
  facet_wrap(. ~ variable, scales = 'free_y') +
  labs(title = 'Quality control metrics') +
  theme(axis.title = element_blank(),
        plot.title = element_text(size = 18, face = 'bold'),
        axis.line.y = element_line(size = 1),
        axis.ticks.y = element_line(size = 1))
qc_plot


#' In general, the samples share very similar distributions across quality 
#' control metric - there are a few outlier cells (see nCount_RNA). However, we 
#' note that the nCount_RNA and nFeature_RNA appear to have bimodal
#' distributions. There is also a long tail of cells with unusually high 
#' percent.mt values (mean of `r round(mean(svz$percent.mt), 1)`%). These cells
#' collectively should be considered suspect and removed from downstream.  

umi_lo_threshold <- 1000
feat_hi_threshold <- 6000 # removes 2 cells
percent_mt_threshold <- 25

bad_cells <- svz$nCount_RNA < umi_lo_threshold | 
  svz$nFeature_RNA > feat_hi_threshold |
  svz$percent.mt > percent_mt_threshold
table(bad_cells)
svz$outlier <- bad_cells

#' Re-examine the disitrbution of cells colored by whether considered outlier:  
meta_vars <- c('orig.ident','outlier','nFeature_RNA','nCount_RNA', 'percent.mt')
qc_plot <- svz@meta.data[meta_vars] %>%
  reshape2::melt(id.vars = c('orig.ident','outlier')) %>%
  ggplot(mapping = aes(x = orig.ident, y = value)) +
  geom_violin(scale = 'width') +
  ggbeeswarm::geom_quasirandom(mapping = aes(color = ifelse(outlier, 'True','False')), size = 0.35) +
  scale_color_manual(values = c('True' = 'red', 'False' = 'black')) +
  facet_wrap(. ~ variable, scales = 'free_y') +
  labs(title = 'Quality control metrics') +
  theme(axis.title = element_blank(),
        plot.title = element_text(size = 18, face = 'bold'),
        axis.line.y = element_line(size = 1),
        axis.ticks.y = element_line(size = 1),
        legend.key = element_rect(fill = NA),
        legend.position = 'bottom') +
  guides(color = guide_legend(title = 'Is outlier?'))
qc_plot

#' Filter the cells:  
table('High quality cells:' = !svz$outlier,
      'Injury time-point:' = svz$orig.ident)

#' Number of high-quality cells: `r sum(!svz$outlier)`
svz <- svz[,!svz$outlier]
saveRDS(object = svz, file = '../data/20210202_SVZ.rds')

sessionInfo()
