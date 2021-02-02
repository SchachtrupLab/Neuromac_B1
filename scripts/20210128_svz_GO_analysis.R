
#' ---
#' title: "SVZ subtype GO analysis"
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
require('biomaRt')
require('topGO')
require('org.Mm.eg.db')


# setwd('./scripts')
results_out <- '../results/gene_ontology/'
dir.create(results_out)

mg <- readRDS(file = '../data/20210112_MG.rds')
nsc <- readRDS(file = '../data/20210108_NSC.rds')



# Ensembl + TopGO setup ---------------------------------------------------

ensembl <- useEnsembl(biomart = 'ENSEMBL_MART_ENSEMBL', 
                      dataset = 'mmusculus_gene_ensembl')

# Use biomaRt to convert MGI symbols to Ensembl IDs
conversion_input <- rownames(mg[['RNA']]@counts)
conversion_input <- plyr::mapvalues(conversion_input, 'Sepp1', 'Selenop')
ensembl_conversion <- getBM(
  attributes = c('mgi_symbol','ensembl_gene_id','external_gene_name'),
  filters = 'external_gene_name',
  values = conversion_input,
  mart = ensembl
)

# Get Ensembl IDs for all genes
global_set <- rownames(mg[['RNA']]@counts)
global_set <- plyr::mapvalues(x = global_set, from = 'Sepp1', to = 'Selenop')
global_set <- plyr::mapvalues(
  x = global_set,
  from = ensembl_conversion$external_gene_name,
  to = ensembl_conversion$ensembl_gene_id,
  warn_missing = FALSE
)
table('Has Ensembl ID:' = grepl('ENSMUS', global_set))


# Microglia GO analysis ------------------------------------------------------

# Compute DE genes
Idents(mg) <- 'SCT_snn_res.0.4'
mg_markers <- FindAllMarkers(
  object = mg,
  only.pos = TRUE,
  assay = 'RNA',
  slot = 'data'
)

# Number of significant DE genes per cluster
ndeg <- mg_markers %>%
  filter(p_val_adj <= 0.001) %>%
  count('Cluster' = cluster)
knitr::kable(x = ndeg, caption = 'Number of DE genes per cluster.')

# MGI symbol to Ensembl ID conversion
mg_markers$gene <- plyr::mapvalues(
  x = mg_markers$gene, 
  from = c('Sepp1'), 
  to = c('Selenop')
)
mg_markers$ensembl <- plyr::mapvalues(
  x = mg_markers$gene, 
  from = ensembl_conversion$external_gene_name,
  to = ensembl_conversion$ensembl_gene_id,
  warn_missing = FALSE
)


#' Following genes do not have an Ensembl ID and were excluded from GO analysis:
knitr::kable(x = table('Has Ensembl ID:' = grepl('ENSMUS', mg_markers$ensembl)))
no_match <- unique(mg_markers$gene[!grepl('ENSMUS', mg_markers$ensembl)])
mg_markers <- mg_markers[!mg_markers$gene %in% no_match,]
no_match


# Perform GO for each microglia subtype
mg_markers$cluster <- as.numeric(as.character(mg_markers$cluster))
subtype <- levels(mg$SCT_snn_res.0.4)
microglia_GO <- vector(mode = 'list', length = length(subtype))

for(i in 1:length(subtype)) {
  # upregulated
  upreg <- with(
    data = mg_markers,
    expr = avg_logFC > 0 & cluster == subtype[i]
  )
  upreg <- mg_markers[upreg, c('gene','p_val_adj','ensembl')]
  gene_pval <- c(upreg$p_val_adj, rep(1, sum(!global_set %in% upreg$ensembl)))
  names(gene_pval) <- 
    c(upreg$ensembl, global_set[!global_set %in% upreg$ensembl])
  
  GOdata <- new("topGOdata", 
                ontology = "BP", 
                allGenes = gene_pval, 
                geneSel = function(p) p <= 1e-03, 
                description = "",
                annot = annFUN.org, 
                mapping="org.Mm.eg.db", 
                ID="Ensembl")
  GOres <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  upreg_result <- GenTable(GOdata, pvalue = GOres, topNodes = 50)
  upreg_result$cluster <- subtype[i]
  microglia_GO[[i]] <- upreg_result
}


# Write results as csv
for(i in 1:length(subtype)) {
  tmp_name <- 
    paste0(results_out, 'microglia_','cluster',subtype[i],'_GOterms.csv')
  write.csv(x = microglia_GO[[i]], 
            file = tmp_name,
            quote = FALSE,
            row.names = FALSE)
}


# Visualization
GO_plots <- vector(mode = 'list', length = length(subtype))
for(i in 1:length(subtype)) {
  tmp <- microglia_GO[[i]]
  tmp <- tmp[,c("GO.ID","Term","pvalue")]
  tmp$Term <- gsub(" [a-z]*\\.\\.\\.$", "", tmp$Term)
  tmp$Term <- gsub("\\.\\.\\.$", "", tmp$Term)
  tmp$Term <- paste(tmp$GO.ID, tmp$Term, sep=", ")
  tmp$Term <- factor(tmp$Term, levels=rev(tmp$Term))
  tmp$pvalue <- as.numeric(tmp$pvalue)
  GO_plots[[i]] <- tmp[1:20,] %>%
    ggplot(mapping = aes(x = Term, y = -log10(pvalue))) +
    geom_bar(position = 'dodge', stat = 'identity') +
    xlab("Biological process") +
    ylab("Enrichment") +
    labs(title = paste('Cluster', i-1, 'GO Enrichment')) +
    scale_y_continuous(
      breaks = round(seq(0, max(-log10(tmp$pvalue)), by = 2), 1)) +
    theme_bw(base_size = 18) +
    theme(legend.position = 'none',
          legend.background = element_rect(),
          plot.title = element_text(angle = 0, size = 12, vjust = 1),
          axis.text.x = element_text(angle = 0, size = 10, hjust = 1),
          axis.text.y = element_text(angle = 0, size = 10, vjust = 0.5),
          axis.title = element_text(size = 12),
          legend.key = element_blank(), #removes the border
          legend.key.size = unit(1, "cm"), # Sets overall ize of the legend
          legend.text = element_text(size = 12), # Text size
          title = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size = 2.5))) +
    coord_flip()
}
GO_plot_final <- cowplot::plot_grid(plotlist = GO_plots, ncol = 3)

#+ mg_GO, fig.height=8, fig.width=18, fig.cap='GO biological process enrichment for each microglia subtype. Enrichment = -log10(p-value).'
GO_plot_final
ggsave(filename = paste0(results_out, 'microglia_GOterms_bargraph.tiff'),
       plot = GO_plot_final,
       device = 'tiff',
       height = 8,
       width = 18)



# NSC GO analysis ---------------------------------------------------------

# Compute DE genes
Idents(nsc) <- 'SCT_snn_res.0.8'
nsc_markers <- FindAllMarkers(
  object = nsc,
  only.pos = TRUE,
  assay = 'RNA',
  slot = 'data'
)

# Number of significant DE genes per cluster
ndeg <- nsc_markers %>%
  filter(p_val_adj <= 0.001) %>%
  count('Cluster' = cluster)
knitr::kable(x = ndeg, caption = 'Number of DE genes per cluster.')

# MGI symbol to Ensembl ID conversion
nsc_markers$gene <- plyr::mapvalues(
  x = nsc_markers$gene, 
  from = c('Sepp1'), 
  to = c('Selenop')
)
nsc_markers$ensembl <- plyr::mapvalues(
  x = nsc_markers$gene, 
  from = ensembl_conversion$external_gene_name,
  to = ensembl_conversion$ensembl_gene_id,
  warn_missing = FALSE
)


#' Following genes do not have an Ensembl ID and were excluded from GO analysis:
knitr::kable(x = table('Has Ensembl ID:' = grepl('ENSMUS', nsc_markers$ensembl)))
no_match <- unique(nsc_markers$gene[!grepl('ENSMUS', nsc_markers$ensembl)])
nsc_markers <- nsc_markers[!nsc_markers$gene %in% no_match,]
no_match


# Perform GO for each nsc subtype
nsc_markers$cluster <- as.numeric(as.character(nsc_markers$cluster))
subtype <- levels(nsc$SCT_snn_res.0.8)
nsc_GO <- vector(mode = 'list', length = length(subtype))

for(i in 1:length(subtype)) {
  # upregulated
  upreg <- with(
    data = nsc_markers,
    expr = avg_logFC > 0 & cluster == subtype[i]
  )
  upreg <- nsc_markers[upreg, c('gene','p_val_adj','ensembl')]
  gene_pval <- c(upreg$p_val_adj, rep(1, sum(!global_set %in% upreg$ensembl)))
  names(gene_pval) <- 
    c(upreg$ensembl, global_set[!global_set %in% upreg$ensembl])
  
  GOdata <- new("topGOdata", 
                ontology = "BP", 
                allGenes = gene_pval, 
                geneSel = function(p) p <= 1e-03, 
                description = "",
                annot = annFUN.org, 
                mapping="org.Mm.eg.db", 
                ID="Ensembl")
  GOres <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  upreg_result <- GenTable(GOdata, pvalue = GOres, topNodes = 50)
  upreg_result$cluster <- subtype[i]
  nsc_GO[[i]] <- upreg_result
}


# Write results as csv
for(i in 1:length(subtype)) {
  tmp_name <- paste0(results_out, 'nsc_', 'cluster', subtype[i], '_GOterms.csv')
  write.csv(x = nsc_GO[[i]], 
            file = tmp_name,
            quote = FALSE,
            row.names = FALSE)
}


# Visualization
GO_plots <- vector(mode = 'list', length = length(subtype))
for(i in 1:length(subtype)) {
  tmp <- nsc_GO[[i]]
  tmp <- tmp[,c("GO.ID","Term","pvalue")]
  tmp$Term <- gsub(" [a-z]*\\.\\.\\.$", "", tmp$Term)
  tmp$Term <- gsub("\\.\\.\\.$", "", tmp$Term)
  tmp$Term <- paste(tmp$GO.ID, tmp$Term, sep=", ")
  tmp$Term <- factor(tmp$Term, levels=rev(tmp$Term))
  tmp$pvalue <- as.numeric(tmp$pvalue)
  GO_plots[[i]] <- tmp[1:20,] %>%
    ggplot(mapping = aes(x = Term, y = -log10(pvalue))) +
    geom_bar(position = 'dodge', stat = 'identity') +
    xlab("Biological process") +
    ylab("Enrichment") +
    labs(title = paste('Cluster', i-1, 'GO Enrichment')) +
    scale_y_continuous(breaks = round(seq(0, max(-log10(tmp$pvalue)), by = 2), 1)) +
    theme_bw(base_size = 18) +
    theme(legend.position = 'none',
          legend.background = element_rect(),
          plot.title = element_text(angle = 0, size = 12, vjust = 1),
          axis.text.x = element_text(angle = 0, size = 10, hjust = 1),
          axis.text.y = element_text(angle = 0, size = 10, vjust = 0.5),
          axis.title = element_text(size = 12),
          legend.key = element_blank(),     #removes the border
          legend.key.size = unit(1, "cm"),      #Sets overall area/size of the legend
          legend.text = element_text(size = 12),  #Text size
          title = element_text(size = 12)) +
    guides(colour = guide_legend(override.aes = list(size = 2.5))) +
    coord_flip()
}
GO_plot_final <- cowplot::plot_grid(plotlist = GO_plots, ncol = 2)

#+ nsc_GO, fig.height=8, fig.width=12, fig.cap='GO biological process enrichment for each NSC subtype. Enrichment = -log10(p-value).'
GO_plot_final
ggsave(filename = paste0(results_out, 'nsc_GOterms_bargraph.tiff'),
       plot = GO_plot_final,
       device = 'tiff',
       height = 8,
       width = 12)
