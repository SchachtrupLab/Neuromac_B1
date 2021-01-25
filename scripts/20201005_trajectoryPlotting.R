

BiocManager::install('slingshot')
BiocManager::install('tradeSeq', version = 'devel')

library('Seurat')
library('dplyr')
library('ggplot2')
library('slingshot')
library('SingleCellExperiment')
library('tradeSeq')

load(file = 'session_data.RData')


# Plot individual genes along pseudotime
my_genes <- c('Dcx','Fxyd1','Fgfr1','Sox11','Ncam1', 'Ntrk2','Mt1','Mt2','Plpp3','Tubb3') # Change genes here


# Run the following chunk all together
tmp_theme <- my_theme + theme(legend.title = element_text(size = 14, color = 'black'),
                              legend.text = element_text(size = 12, color = 'black'),
                              plot.title = element_text(size = 16, color = 'black', face = 'bold'))
smooth_plots <- vector(mode = 'list', length = length(my_genes))
for(i in 1:length(my_genes)) {
  # Change "nsc_gam" to "mg_gam", and "counts = counts(nsc_sce)" to "counts = counts(mg_sce)"
  # if you want to plot along the MG trajectories.
  smooth_plots[[i]] <- plotSmoothers(gene = my_genes[i],
                                     nsc_gam,
                                     counts = counts(nsc_sce),
                                     lwd = 2, size = 1.5, border = TRUE) + 
    tmp_theme +
    labs(title = my_genes[i])
}
smooth_plots <- cowplot::plot_grid(plotlist = smooth_plots, ncol = 5)
smooth_plots