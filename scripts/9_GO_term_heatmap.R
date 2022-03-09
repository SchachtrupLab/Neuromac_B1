require('dplyr')
require('ggplot2')
require('Seurat')

setwd('scripts')
results_out <- '../results/8_GOanalysis/'
dir.create(path = results_out)

mg <- readRDS(file = '../data/mg.rds')

dam <- readxl::read_xlsx(path = '../data/DAM_100.xlsx')
act <- readxl::read_xlsx(path = '../data/Activated_100.xlsx')
hom <- readxl::read_xlsx(path = '../data/Homeostatic_100.xlsx')

dam$subtype <- 'DAM'
act$subtype <- 'Activated MG'
hom$subtype <- 'Homeostatic MG'

results <- Reduce(f = rbind, x = list(dam, act, hom))
results$subtype <- factor(
  x = results$subtype,
  levels = c('Homeostatic MG',
             'Activated MG',
             'DAM')
)


# Make plot like in Figure 3b of: https://www.nature.com/articles/s41467-021-26059-4
head(results)

# Using ggplot2
cols <- viridis::viridis_pal(option = 'A')
cols <- rev(cols(100))
GOheatmap <- results %>% 
  filter(source == 'GO:BP') %>%
  group_by(subtype) %>% 
  top_n(n = 15, wt = negative_log10_of_adjusted_p_value) %>% 
  mutate('tmp_term' = paste(subtype, term_name, sep = ' ')) %>%
  group_by(subtype) %>% 
  mutate('tmp_term' = reorder(tmp_term, negative_log10_of_adjusted_p_value)) %>% 
  mutate(tmp_term = factor(tmp_term, levels = rev(levels(tmp_term)))) %>%
  mutate('order' = as.numeric(tmp_term)) %>% 
  ggplot(mapping = aes(x = subtype, y = reorder(term_name, order))) +
  geom_tile(mapping = aes(fill = negative_log10_of_adjusted_p_value),
            color = 'black') +
  scale_fill_gradientn(colors = cols, 
                       na.value = 'black',
                       limits = c(0, NA),
                       breaks = seq(0, 30, 5)) +
  ylab(label = 'GO:Biological Process') +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid.minor = element_line(color = 'black'),
        panel.border = element_rect(fill = NA, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.text = element_blank(),
        legend.title = element_text(angle = 90)) +
  guides(fill = guide_colorbar(title = '-log10(adjusted p-value)',
                               ticks.colour = 'black',
                               frame.colour = 'black',
                               title.position = 'left'))
ggsave(filename = paste0(results_out, 'MGsubtype_GOheatmap.svg'),
       plot = GOheatmap, device = 'svg', height = 8, width = 8)
ggsave(filename = paste0(results_out, 'MGsubtype_GOheatmap.tiff'),
       plot = GOheatmap, device = 'tiff', height = 8, width = 8)
GOheatmap




# 2021-10-11: You can ignore this portion. We will use ggplot2 for
# plotting the GO terms heatmap.
# Using Complexheatmap
top_terms <- results %>% 
  filter(source == 'GO:BP') %>%
  group_by(subtype) %>% 
  top_n(n = 15, wt = negative_log10_of_adjusted_p_value) %>% 
  select(term_name, subtype, negative_log10_of_adjusted_p_value) %>% 
  arrange(subtype)
tmp <- reshape2::acast(
  data = top_terms, 
  formula = term_name ~ subtype,
  value.var = 'negative_log10_of_adjusted_p_value'
)
ComplexHeatmap::Heatmap(
  matrix = tmp,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = scales::viridis_pal(option = 'A')(100),
  na_col = 'black'
)
