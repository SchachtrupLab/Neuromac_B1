##### Preliminary Data for Schachtrup Collaboration ###########################

# Hypothesis: TGF-B1 expression by microglia in SCI activates fibroblasts

library('Seurat')
library('ggplot2')
library('dplyr')
library('reshape2')
library('cowplot')
source('../sci_cellAtlas_scRNAseq/R/Functions/functions.R')

# Set up seurat object
sci.full <- readRDS('../sci_cellAtlas_scRNAseq/R/Data/Fig1_sci.rds')
# vascular <- readRDS('../sci_cellAtlas_scRNAseq/R/Data/Fig3_vascular.rds')
DefaultAssay(sci.full) <- 'SCT'
sci.full$Figure_Full <- factor(x = sci.full$Figure_Full)
Idents(sci.full) <- 'Figure_Full'
sci <- sci.full[,!sci.full$Figure_Full %in% c('Low_Quality','U1-Myeloid','U1-Vascular','U2-Myeloid',
                                              'U2-Vascular','U3-Myeloid','U4-Myeloid','U5-Myeloid',
                                              'B Cell','T Cell','Neuron')]
Idents(sci) <- 'broad_cellType'
DimPlot(sci, label = TRUE)

sci$broad_cellType <- plyr::mapvalues(
  x = sci$Figure_Full,
  from = c('Apop-Microglia','DA-Microglia','Div-Microglia','Foamy-Microglia','H-Microglia',
           'Monocyte','Trans-Monocyte',
           'Apoe-Macrophage','BA-Macrophage','Dendritic Cell','Foamy-Macrophage',
           'A-Endothelial','C-Endothelial','IFN-Vascular','SMC','V-Endothelial','Tip Cell',
           'OPC','Div-OPC','Pre-Oligo',
           'H1-Ependymal','H2-Ependymal','Astroependymal','Trans-Ependymal',
           'Oligodendrocyte','Fibroblast','Div-Myeloid','Astrocyte','Pericyte','Neutrophil','Ph-Neutrophil'),
  to = c('Microglia','Microglia','Microglia','Microglia','Microglia',
         'Monocyte','Monocyte',
         'Macrophage','Macrophage','Macrophage','Macrophage',
         'Endothelial','Endothelial','Endothelial','Endothelial','Endothelial', 'Endothelial',
         'OPC','OPC','OPC',
         'Ependymal','Ependymal','Ependymal','Ependymal',
         'Oligodendrocyte','Fibroblast','Myeloid_Derived_Cell','Astrocyte','Pericyte','Neutrophil','Neutrophil')
)
sci$broad_cellType <- factor(sci$broad_cellType, levels = c('Microglia','Monocyte','Macrophage','Myeloid_Derived_Cell',
                                                            'Endothelial','Pericyte','Fibroblast','Oligodendrocyte',
                                                            'OPC','Astrocyte','Ependymal','Neutrophil','B Cell','T Cell','Neuron'))
Idents(sci) <- 'broad_cellType'
DimPlot(sci, label = TRUE)

DimPlot(sci, label = TRUE)
saveRDS(sci, '../sci_cellAtlas_scRNAseq/R/Data/Fig1_sci_adjusted.rds')

# Specific genes of interest
my.genes <-  c('Tgfb1','Tgfbr1','Tgfbr2','Smad2','Smad3')


# Overall UMAP ------------------------------------------------------------

umap <- DimPlot(object = sci, label = TRUE, label.size = 5, repel = TRUE) +
  theme(axis.ticks = element_blank(),
        axis.line = element_line(size = 1),
        axis.text = element_blank(),
        legend.position = 'none')
ggsave(filename = 'TGFB-signaling_SCI_overall-umap.png',
       plot = umap,
       height = 6, width = 7)

# Representative UMAPs -----------------------------------------------------------------
my.umaps <- {
  tmp.umap <- vector(mode = 'list',
                     length = length(my.genes) * length(levels(sci$Time)))
  counter <- 1
  for(t in levels(sci$Time)) {
    tmp.data <- sci[,sci$Time == t]
    for(g in my.genes) {
      tmp.umap[[counter]] <- FeaturePlot(object = tmp.data,
                                         features = g,
                                         cols = c('grey','red3'),
                                         min.cutoff = 0,
                                         pt.size = 1) +
        theme(axis.ticks = element_blank(),
              axis.line = element_line(size = 1),
              axis.text = element_blank(),
              legend.position = c(0.9, 0.2),
              legend.key.size = unit(0.4, units = 'cm'),
              legend.text = element_text(size = 10),
              legend.spacing.x = unit(0.1, 'cm'),
              plot.title = element_text(hjust = 0)) +
        labs(title = t)
      counter <- counter + 1
    }
  }; counter <- 1; rm(tmp.data, t, g)
  
  my.umaps <- vector(mode = 'list', 
                     length = length(my.genes))
  for(g in 1:length(my.genes)) {
    tmp.list <- seq(from = g, 
                    by = length(my.genes),
                    length.out = length(levels(sci$Time)))
    tmp.list <- tmp.umap[tmp.list]
    my.umaps[[g]] <- plot_grid(plotlist = tmp.list, 
                               ncol = length(levels(sci$Time))/2)
    tmp.title <- ggdraw() + draw_label(label = my.genes[g], fontface = 'bold', size = 20)
    my.umaps[[g]] <- plot_grid(tmp.title, my.umaps[[g]], ncol = 1, rel_heights = c(0.1,1))
  }; rm(tmp.list, g)
  my.umaps
}


# VlnPlot expression --------------------------------------------------------------------------
my.vlnplots <- {
  # Generate vlns by gene, time
  tmp.vln <- vector(mode = 'list', 
                    length = length(my.genes) * length(levels(sci$Time)))
  counter <- 1
  for(t in levels(sci$Time)) {
    tmp.data <- sci[,sci$Time == t]
    for(g in my.genes) {
      tmp.vln[[counter]] <- VlnPlot(object = tmp.data,
                                    features = g,
                                    pt.size = 0) +
        theme(legend.position = 'none', 
              axis.line = element_line(size = 1),
              axis.title = element_blank(),
              axis.ticks = element_line(size = 1),
              plot.title = element_text(hjust = 0)) +
        labs(title = t)
      counter <- counter + 1
    }
  }; counter <- 1; rm(tmp.data, t, g)
  # Aggregate data for each gene across times
  my.vlnplots <- vector(mode = 'list', 
                        length = length(my.genes))
  for(g in 1:length(my.genes)) {
    tmp.list <- seq(from = g, 
                    by = length(my.genes),
                    length.out = length(levels(sci$Time)))
    tmp.list <- tmp.vln[tmp.list]
    tmp.title <- list(ggdraw() + draw_label(label = my.genes[g],
                                            fontface = 'bold',
                                            size = 20))
    tmp.list <- append(tmp.title, tmp.list)
    my.vlnplots[[g]] <- plot_grid(plotlist = tmp.list,
                                  ncol = 1, 
                                  rel_heights = c(0.15, rep(1, length(levels(sci$Time)))))
  }; rm(tmp.list, tmp.title, g)
  my.vlnplots
}


# Dotplot expression -----------------------------------------------------------------------
my.dotplots <- {
  percent <- function(x) {
    return(sum(x > 0)/length(x))
  }
  MyDotPlot <- function(gene, data) {
    tmp.mean <- data %>%
      group_by(ident, Time) %>%
      summarise_at(c(gene), mean)
    tmp.pct <- data %>%
      group_by(ident, Time) %>%
      summarise_at(c(gene), percent)
    tmp.mean <- melt(data = tmp.mean, id.vars = c('ident','Time'))
    tmp.pct <- melt(data = tmp.pct, id.vars = c('ident','Time'))
    tmp.gg <- merge(x = tmp.mean, y = tmp.pct, by = c('ident','Time','variable'))
    colnames(tmp.gg) <- c('ident','Time','gene','mean','pct')
    tmp.gg[['gene']] <- factor(tmp.gg[['gene']], levels = rev(levels(tmp.gg[['gene']])))
    tmp.gg[['ident']] <- factor(tmp.gg[['ident']], levels = rev(levels(tmp.gg[['ident']])))
    
    tmp.gg <- tmp.gg %>%
      ggplot(mapping = aes(x = Time, y = ident, size = pct, color = mean)) +
      geom_point(mapping = aes(color = mean, size = pct)) +
      scale_color_gradient2(mid = 'blue3', high = 'red3') +
      scale_size_continuous(range = c(0, 10)) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 14),
            plot.title = element_text(size = 22),
            axis.title = element_blank(),
            legend.text = element_text(size = 11),
            legend.title = element_text(size = 14)) +
      guides(color = guide_colorbar(title = 'Mean\nExpression'),
             size = guide_legend(title = 'Percent\nExpression')) +
      # xlab(label = 'Time') +
      # ylab(label = 'Cell Type') +
      labs(title = gene)
    
    return(tmp.gg)
  }
  my.data <- FetchData(object = sci, vars = c(my.genes, 'ident','Time'), slot = 'data')
  my.dotplots <- lapply(
    X = my.genes,
    data = my.data,
    FUN = function(x, data) {MyDotPlot(gene = x, data = data)}
  ); rm(my.data)
  my.dotplots
}



# Save Expression results ------------------------------------------------------------

for(i in 1:length(my.genes)) {
  ggsave(filename = paste('TGFB-signaling_SCI', my.genes[i], 'Umap.png', sep = '_'),
         plot = my.umaps[[i]],
         device = 'png',
         height = 11.5, width = 10)
  ggsave(filename = paste('TGFB-signaling_SCI', my.genes[i], 'Dotplot.png', sep = '_'),
         plot = my.dotplots[[i]],
         device = 'png',
         height = 6, width = 5.5)
  ggsave(filename = paste('TGFB-signaling_SCI', my.genes[i], 'Vlnplot.png', sep = '_'),
         plot = my.vlnplots[[i]],
         device = 'png',
         height = 11.5, width = 5.5)
}



# TGFB Ligand-Receptor Interactions ---------------------------------------

source('../sci_cellAtlas_scRNAseq/R/Functions/20190921LR_plotting.R')

lr <- readRDS('../sci_cellAtlas_scRNAseq/R/Data/LR_results_broadcelltype.rds')
lr.ref <- read.table(file = 'fantom_PairsLigRec.txt',
                     header = TRUE,
                     sep = '\t')
lr.ref.known <- lr.ref[lr.ref$Pair.Source == 'known' & lr.ref$Pair.Evidence == 'literature supported',]
lr.ref.known <- lr.ref.known$Pair.Name

lr.subset <- lr[lr$Pair_name %in% lr.ref.known,]

# TGFB1 in MICROGLIA %%%%%%%%%%%%%%%%%%%
microglia.tgfb1.nonsig <- {
  myPalette <- colorRampPalette(colors = rev(brewer.pal(11, 'RdBu')))
  tmp <- selectLR(score.results.df = lr.subset,
                  ligand = 'Tgfb1',
                  ligand.cell = 'Microglia',
                  significant.only = FALSE)
  sc <- scale_fill_gradientn(colors = myPalette(20), limits = c(min(tmp$data$Score), max(tmp$data$Score)))
  tmp <- TilePlotLR_simple(select.out = tmp,
                           title.text = 'TGFB1 Signaling by Microglia',
                           subtitle.text = 'All results',
                           title.main.font.scale = 1,
                           title.font.scale = 0.75,
                           uniform.scale = TRUE,
                           mid.point = -1) +
    facet_grid(Ligand_Cell~Time, drop = TRUE, switch = 'y') +
    xlab('Ligand-Receptor Pairs') + 
    ylab('Receptor Cell') +
    sc
  tmp
}
microglia.tgfb1.sig <- {
  tmp <- selectLR(score.results.df = lr.subset,
                  ligand = 'Tgfb1',
                  ligand.cell = 'Microglia',
                  significant.only = TRUE)
  tmp <- TilePlotLR_simple(select.out = tmp,
                           title.text = 'TGFB1 Signaling by Microglia',
                           subtitle.text = 'P-value < 0.05 results only',
                           title.main.font.scale = 1,
                           title.font.scale = 0.75) +
    facet_grid(Ligand_Cell~Time, drop = TRUE, switch = 'y') +
    xlab('Ligand-Receptor Pairs') + 
    ylab('Receptor Cell') +
    sc
  tmp
}
ggsave(filename = 'TGFB1_microglia_nonsignificant_LR.png',
       plot = microglia.tgfb1.nonsig,
       device = 'png',
       height = 10, width = 13)
ggsave(filename = 'TGFB1_microglia_significant_LR.png',
       plot = microglia.tgfb1.sig,
       device = 'png',
       height = 10, width = 13)


# TGFBR2 in FIBROBLASTS %%%%%%%%%%%%%%%%%%%
fibroblast.tgfbr2.nonsig <- {
  myPalette <- colorRampPalette(colors = rev(brewer.pal(11, 'RdBu')))
  tmp <- selectLR(score.results.df = lr.subset,
                  receptor = 'Tgfbr2',
                  receptor.cell = 'Fibroblast',
                  significant.only = FALSE)
  sc <- scale_fill_gradientn(colors = myPalette(20), limits = c(min(tmp$data$Score), max(tmp$data$Score)))
  tmp <- TilePlotLR_simple(select.out = tmp,
                           title.text = 'TGFBR2 Signaling by Fibroblasts',
                           subtitle.text = 'All results',
                           title.main.font.scale = 0.8,
                           title.font.scale = 0.6,
                           uniform.scale = TRUE,
                           mid.point = -1) +
    facet_grid(Receptor_Cell~Time, drop = TRUE) +
    xlab('Ligand-Receptor Pairs') + 
    ylab('Ligand Cell') +
    sc
  tmp
}
fibroblast.tgfbr2.sig <- {
  myPalette <- colorRampPalette(colors = rev(brewer.pal(11, 'RdBu')))
  tmp <- selectLR(score.results.df = lr.subset,
                  receptor = 'Tgfbr2',
                  receptor.cell = 'Fibroblast',
                  significant.only = TRUE)
  tmp <- TilePlotLR_simple(select.out = tmp,
                           title.text = 'TGFBR2 Signaling by Fibroblasts',
                           subtitle.text = 'P-value < 0.05 results only',
                           title.main.font.scale = 0.8,
                           title.font.scale = 0.6,
                           uniform.scale = TRUE,
                           mid.point = -1) +
    facet_grid(Receptor_Cell~Time, drop = TRUE) +
    xlab('Ligand-Receptor Pairs') + 
    ylab('Ligand Cell') +
    sc
  tmp
}
ggsave(filename = 'TGFBR2_fibroblast_nonsignificant_LR.png',
       plot = fibroblast.tgfbr2.nonsig,
       device = 'png',
       height = 10, width = 11.5)
ggsave(filename = 'TGFBR2_fibroblast_significant_LR.png',
       plot = fibroblast.tgfbr2.sig,
       device = 'png',
       height = 10, width = 11.5)


# TGFBR2 in Pericytes vs Fibroblasts %%%%%%%%%%%%%%%%%%%%%
pericyte.tgfbr2.nonsig <- {
  myPalette <- colorRampPalette(colors = rev(brewer.pal(11, 'RdBu')))
  tmp <- selectLR(score.results.df = lr.subset,
                  receptor = 'Tgfbr2',
                  receptor.cell = c('Fibroblast','Pericyte'),
                  organize.y = 'Receptor_Cell',
                  significant.only = FALSE)
  sc <- scale_fill_gradientn(colors = myPalette(20), limits = c(min(tmp$data$Score), max(tmp$data$Score)))
  tmp <- TilePlotLR_simple(select.out = tmp,
                           title.text = 'TGFBR2 Signaling by Fibroblasts and Pericytes',
                           subtitle.text = 'All results',
                           title.main.font.scale = 0.8,
                           title.font.scale = 0.6,
                           uniform.scale = TRUE,
                           mid.point = -1) +
    facet_grid(Receptor_Cell~Time, drop = TRUE) +
    xlab('Ligand-Receptor Pairs') +
    ylab('Ligand Cell') +
    sc
  tmp
}
pericyte.tgfbr2.sig <- {
  myPalette <- colorRampPalette(colors = rev(brewer.pal(11, 'RdBu')))
  tmp <- selectLR(score.results.df = lr.subset,
                  receptor = 'Tgfbr2',
                  receptor.cell = c('Fibroblast','Pericyte'),
                  organize.y = 'Receptor_Cell',
                  significant.only = TRUE)
  tmp <- TilePlotLR_simple(select.out = tmp,
                           title.text = 'TGFBR2 Signaling by Fibroblasts and Pericytes',
                           subtitle.text = 'P-value < 0.05 results only',
                           title.main.font.scale = 0.8,
                           title.font.scale = 0.6,
                           uniform.scale = TRUE,
                           mid.point = -1) +
    facet_grid(Receptor_Cell~Time, drop = TRUE) +
    xlab('Ligand-Receptor Pairs') +
    ylab('Ligand Cell') +
    sc
  tmp
}
ggsave(filename = 'TGFBR2_pericyte-vs-fibroblast_nonsignificant_LR.png',
       plot = pericyte.tgfbr2.nonsig,
       device = 'png',
       height = 15.5, width = 11.5)
ggsave(filename = 'TGFBR2_pericyte-vs-fibroblast_significant_LR.png',
       plot = pericyte.tgfbr2.sig,
       device = 'png',
       height = 15.5, width = 11.5)


# Additional Gene expressions of interest ---------------------------------

Idents(sci) <- 'broad_cellType'
my.genes <- c('Inha','Inhba','Inhbb')
my.plot <- FeaturePlot(object = sci, features = my.genes, split.by = 'Time',
                       cols = c('grey','red3'), pt.size = 0.75)
ggsave(filename = 'Inhibin_umap.png',
       plot = my.plot,
       device = 'png',
       height = 10.5, width = 15.75)

my.plots <- vector(mode = 'list', length = length(my.genes)*length(levels(sci$Time)))
counter <- 1
for(t in levels(sci$Time)) {
  tmp.data <- sci[,sci$Time == t]
  for(g in my.genes) {
    my.plots[[counter]] <- FeaturePlot(object = tmp.data,
                                       features = g,
                                       cols = c('grey','red3'),
                                       min.cutoff = 0,
                                       pt.size = 1) +   
      theme(axis.ticks = element_blank(),
            axis.line = element_line(size = 1),
            axis.text = element_blank(),
            legend.position = c(0.9, 0.2),
            legend.key.size = unit(0.4, units = 'cm'),
            legend.text = element_text(size = 10),
            legend.spacing.x = unit(0.1, 'cm'),
            plot.title = element_text(hjust = 0)) +
      labs(title = t)
    counter <- counter + 1
  }
}
for(g in my.genes) {
  my.plots[[g]] <- 
}

for(g in my.genes) {
  tmp.umap[[counter]] <- FeaturePlot(object = tmp.data,
                                     features = g,
                                     cols = c('grey','red3'),
                                     min.cutoff = 0,
                                     pt.size = 1) +
    theme(axis.ticks = element_blank(),
          axis.line = element_line(size = 1),
          axis.text = element_blank(),
          legend.position = c(0.9, 0.2),
          legend.key.size = unit(0.4, units = 'cm'),
          legend.text = element_text(size = 10),
          legend.spacing.x = unit(0.1, 'cm'),
          plot.title = element_text(hjust = 0)) +
    labs(title = t)
  counter <- counter + 1
}