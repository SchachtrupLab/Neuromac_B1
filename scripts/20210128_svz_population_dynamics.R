
#' ---
#' title: "SVZ population dynamics"
#' author: "James Choi"
#' date: "`r Sys.Date()`"
#' output: pdf_document
#' ---

#+ setup, include=TRUE
knitr::opts_chunk$set(warning=FALSE, message=FALSE, fig.align='center', 
                      tidy.opts=list(width.cutoff=60), tidy=TRUE)

require('Seurat')
require('ggplot2')
require('dplyr')
require('dendextend')

# setwd('./scripts')
results_out <- '../results/population_dynamics/'
dir.create(results_out)

mg <- readRDS(file = '../data/20210112_MG.rds')
nsc <- readRDS(file = '../data/20210108_NSC.rds')



# Microglia population dynamics -------------------------------------------

Idents(mg) <- 'SCT_snn_res.0.4'
counts <- data.frame(table(mg$SCT_snn_res.0.4, mg$orig.ident))
names(counts) <- c('subtype','orig.ident','count')
numbers <- counts %>%
  ggplot(mapping = aes(x = orig.ident, y = count, group = subtype)) + 
  geom_bar(aes(fill = subtype), stat = 'identity', color = 'black', size = 0.75) +
  scale_y_continuous() +
  scale_x_discrete() + 
  ylab('# of cells') + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.line = element_line(size = 1),
        panel.background = element_rect(fill = NA),
        axis.text.y = element_text(size = 12))
props <- counts %>% 
  ggplot(mapping = aes(x = orig.ident, y = count, group = subtype)) + 
  geom_bar(aes(fill = subtype), stat = 'identity', color = 'black', size = 0.75, position = 'fill') +
  scale_y_continuous() + 
  scale_x_discrete() +
  ylab('% of cells') + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.line = element_line(size = 1), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(size = 12), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12))
legend <- cowplot::get_legend(plot = props)
pop.dynamics <- cowplot::plot_grid(numbers + theme(legend.position = 'none'), props + theme(legend.position = 'none'), legend, ncol = 3, rel_widths = c(1, 1, 0.25))

#+ mg_dynamics, fig.height=3, fig.width=6, fig.cap='Left: bar graph of microglia subtype counts across injury time-points. Right: bar graph of NSC subtype proportions across injury time-points.'
pop.dynamics
# ggsave(filename = paste0(results_out, 'microglia_population_dynamics.tiff'),
#        plot = pop.dynamics, device = 'tiff', height = 3, width = 6)




# NSC population dynamics -------------------------------------------------

Idents(nsc) <- 'SCT_snn_res.0.8'
counts <- data.frame(table(nsc$SCT_snn_res.0.8, nsc$orig.ident))
names(counts) <- c('subtype','orig.ident','count')
numbers <- counts %>%
  ggplot(mapping = aes(x = orig.ident, y = count, group = subtype)) + 
  geom_bar(aes(fill = subtype), stat = 'identity', color = 'black', size = 0.75) +
  scale_y_continuous() +
  scale_x_discrete() + 
  ylab('# of cells') + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12), 
        axis.line = element_line(size = 1),
        panel.background = element_rect(fill = NA),
        axis.text.y = element_text(size = 12))
props <- counts %>% 
  ggplot(mapping = aes(x = orig.ident, y = count, group = subtype)) + 
  geom_bar(aes(fill = subtype), stat = 'identity', color = 'black', size = 0.75, position = 'fill') +
  scale_y_continuous() + 
  scale_x_discrete() +
  ylab('% of cells') + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12),
        axis.line = element_line(size = 1), 
        panel.background = element_rect(fill = NA), 
        axis.text.y = element_text(size = 12), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12))
legend <- cowplot::get_legend(plot = props)
pop.dynamics <- cowplot::plot_grid(numbers + theme(legend.position = 'none'), props + theme(legend.position = 'none'), legend, ncol = 3, rel_widths = c(1, 1, 0.25))

#+ nsc_dynamics, fig.height=3, fig.width=6, fig.cap='Left: bar graph of NSC subtype counts across injury time-points. Right: bar graph of NSC subtype proportions across injury time-points.'
pop.dynamics
# ggsave(filename = paste0(results_out, 'NSC_population_dynamics.tiff'),
#        plot = pop.dynamics, device = 'tiff', height = 3, width = 6)


sessionInfo()