# Take subset of LR data.frame for plotting
selectLR <- function(
  score.results.df,
  timepoints = c('Uninjured','1dpi','3dpi','7dpi'),
  receptor.cell = NULL,
  receptor = NULL,
  ligand.cell = NULL,
  ligand = NULL,
  pair.name = NULL,
  significant.only = FALSE,
  organize.y = 'Ligand_Cell',
  sort.order.x = 'Ligand'
) {
  require('dplyr')
  
  params <- list()
  params$data <- score.results.df
  # Functionality to add: DotPlot to add p-val info (low priority)
  params$style <- list('dims' = 1,
                       'p' = 'TilePlotLR')
  
  if(length(ligand.cell) == 0) {
    ligand.cell <- unique(score.results.df[['Ligand_Cell']])
  } 
  # else if(any(ligand.cell %in% names(LR_Interactions_Functions_objects$cell))) {
  #   ligand.cell <- union(ligand.cell[-which(ligand.cell %in% names(LR_Interactions_Functions_objects$cell))],
  #                        unlist(LR_Interactions_Functions_objects$cell[ligand.cell[ligand.cell %in% names(LR_Interactions_Functions_objects$cell)]]))
  # }
  if(length(receptor.cell) == 0) {
    receptor.cell <- unique(score.results.df[['Receptor_Cell']])
  } 
  # else if(any(receptor.cell %in% names(LR_Interactions_Functions_objects$cell))) {
  #   receptor.cell <- union(receptor.cell[-which(receptor.cell %in% names(LR_Interactions_Functions_objects$cell))],
  #                          unlist(LR_Interactions_Functions_objects$cell[receptor.cell[receptor.cell %in% names(LR_Interactions_Functions_objects$cell)]]))
  # }
  if(length(timepoints) == 0) {
    stop('No timepoint selected. Please choose a timepoint.')
  }
  
  if(length(pair.name) != 0) {
    pair.name.tmp <- pair.name
  } else if(all(c(length(pair.name), length(ligand), length(receptor)) == 0)) {
    pair.name.tmp <- levels(params$data$Pair_name)
  } else {
    if(length(ligand) == 0) {ligand <- levels(score.results.df$Ligand)}
    if(length(receptor) == 0) {receptor <- levels(score.results.df$Receptor)}
    pair.name.tmp <- union(pair.name, 
                           intersect(params$data$Pair_name[params$data$Ligand %in% ligand],
                                     params$data$Pair_name[params$data$Receptor %in% receptor]))
  }
  
  params$data <- params$data %>% 
    {if(significant.only) filter(., Pval < 0.05) else .} %>%
    filter(Ligand_Cell %in% ligand.cell) %>%
    filter(Receptor_Cell %in% receptor.cell) %>%
    filter(Pair_name %in% pair.name.tmp) %>%
    filter(Time %in% timepoints)
  
  if(!nrow(params$data) > 0) {stop('No interaction score available.')}
  
  params$data$Pair_name <- factor(params$data$Pair_name, 
                                  levels = unique(params[['data']][['Pair_name']][order(params[['data']][sort.order.x])]))
  
  
  if(all(c(length(ligand.cell), length(receptor.cell)) > 1)) {
    params$style$dims <- 2
    params$x.axis <- 'Pair_name'
    params$y.axis <- c('Ligand_Cell','Receptor_Cell')[which(c('Ligand_Cell','Receptor_Cell') != organize.y)]
    params$y.facet <- organize.y
    params$switch <- list('y.position' = ifelse(organize.y == 'Ligand_Cell', 'right', 'left'),
                          'switch' = switch((organize.y == 'Ligand_Cell')+1, NULL, 'y'))
  } else if(length(ligand.cell) == 1) {
    params$x.axis <- 'Pair_name'
    params$y.axis <- 'Receptor_Cell'
    params$y.facet <- 'Ligand_Cell'
    params$switch <- list('y.position' = 'right',
                          'switch' = 'y')
  } else if(length(receptor.cell) == 1) {
    params$x.axis <- 'Pair_name'
    params$y.axis <- 'Ligand_Cell'
    params$y.facet <- 'Receptor_Cell'
    params$switch <- list('y.position' = 'left',
                          'switch' = NULL)
  }
  
  return(params)
}



TilePlotLR_simple <- function(
  select.out,
  title.text = '',
  subtitle.text = '',
  title.font.scale = 1,
  title.main.font.scale = 1,
  axis.text.font.scale = 1,
  strip.text.override = NULL,
  uniform.scale = FALSE,
  mid.point = 0,
  uniform.scale.lo = NULL,
  uniform.scale.hi = NULL
) {
  require('ggplot2')
  require('RColorBrewer')

  select.out$data[[select.out$y.axis]] <- factor(
    x = select.out$data[[select.out$y.axis]],
    levels = rev(levels(select.out$data[[select.out$y.axis]]))
  )
  myPalette <- colorRampPalette(colors = rev(brewer.pal(11, 'RdBu')))
  if(!uniform.scale) {
    hi = quantile(select.out$data$Score, 1, na.rm = TRUE)
    lo = min(select.out$data$Score, na.rm = TRUE)
    sc <- scale_fill_gradientn(colors = myPalette(20), limits = c(lo, hi))
  } else {
    vals <- seq(from = min(select.out$data$Score), to = max(select.out$data$Score), length.out = 20)
    sc <- scale_fill_gradientn(colors = myPalette(20), values = rev(vals))
  }
  tmp <- select.out$data %>%
    ggplot(aes_string(y = select.out$y.axis, x = select.out$x.axis, fill = 'Score')) +
    geom_tile(color = 'black', size = 0.5) +
    facet_grid(reformulate('Time', select.out$y.facet), 
               switch = select.out$switch[[2]],
               drop = FALSE) +
    scale_y_discrete(position = select.out$switch[[1]],
                     drop = FALSE) +
    xlab(label = 'Ligand_Receptor Pairs') +
    ylab(label = select.out$y.axis) +
    labs(title = title.text, subtitle = subtitle.text) +
    theme(panel.grid.major.y = element_line(color = 'grey50'),
          panel.grid.major.x = element_line(color = 'grey85'),
          panel.background = element_rect(fill = 'grey85'),
          strip.text.y = element_text(size = title.font.scale*34),
          strip.text.x = element_text(size = title.font.scale*34),
          axis.title.y = element_text(size = (title.main.font.scale^2)*28),
          axis.title.x = element_text(size = (title.main.font.scale^2)*28),
          axis.text.y = element_text(size = axis.text.font.scale*20),
          axis.text.x = element_text(size = axis.text.font.scale*14, angle = 90, hjust = 1, vjust = 0.35),
          plot.title = element_text(size = (title.main.font.scale^2)*40),
          plot.subtitle = element_text(size = (title.main.font.scale^2)*28),
          legend.background = element_rect(fill = 'grey85', size = 1, linetype = 'solid'),
          legend.key.size = unit(3, 'line'),
          legend.text = element_text(size = title.font.scale*28), 
          legend.title = element_text(size = (title.main.font.scale)^2 * 34),
          plot.margin = unit(c(2,2,2,2), 'cm')) + sc
  
  if(!is.null(strip.text.override)) {tmp <- tmp + theme(strip.text.y = element_text(size = strip.text.override))}
  return(tmp)
}
