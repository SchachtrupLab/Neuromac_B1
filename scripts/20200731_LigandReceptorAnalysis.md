---
title: "Ligand-Receptor Interaction Analysis"
author: "James Choi"
date: "2020-07-31"
output:  
  html_document:
    keep_md: true
    code_folding: show
    toc: true
---




```r
library('Seurat')
library('ggplot2')
library('dplyr')
library('ggrepel')
source('LRfunctions.R')
source('VolcanoPlot.R')
```


```r
svz <- readRDS(file = '../data/20200511_SVZ.rds')
lr_ref <- '../ref/fantom_PairsLigRec_mouse.csv'
lr_ref <- read.csv(file = lr_ref)
results_path <- '../results/20200731_LigandReceptorAnalysis/'
dir.create(path = results_path)
```


```r
svz$full_clusters <- svz$RNA_snn_res.0.8
svz$celltype <- plyr::mapvalues(x = svz$RNA_snn_res.0.8,
                                from = 0:7,
                                to = c('Microglia',
                                       'Microglia',
                                       'Microglia',
                                       'NSC',
                                       'NSC',
                                       'NSC',
                                       'Microglia',
                                       'Unknown'))
Idents(svz) <- 'celltype'
DefaultAssay(svz)
my_setup <- setupLR(seurat_object = svz, lr_ref = lr_ref, split_by = 'orig.ident')
my_results <- calculateLR(setup = my_setup)
write.csv(x = my_results, file = paste0(results_path, 'LR_results_celltype.csv'), quote = FALSE)
```

```
## [1] "RNA"
##   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  15%  |                                                                              |=============                                                         |  19%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  26%  |                                                                              |=====================                                                 |  30%  |                                                                              |=======================                                               |  33%  |                                                                              |==========================                                            |  37%  |                                                                              |=============================                                         |  41%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=================================================                     |  70%  |                                                                              |====================================================                  |  74%  |                                                                              |======================================================                |  78%  |                                                                              |=========================================================             |  81%  |                                                                              |============================================================          |  85%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
```


```r
# Changes in microglia gene expression after injury
microglia_1d_ctrl <- FindMarkers(object = svz,
                                 subset.ident = 'Microglia',
                                 group.by = 'orig.ident',
                                 ident.1 = 'S1d',
                                 ident.2 = 'SCtrl')
microglia_7d_1d <- FindMarkers(object = svz,
                               subset.ident = 'Microglia',
                               group.by = 'orig.ident',
                               ident.1 = 'S7d',
                               ident.2 = 'S1d')
microglia_1d_ctrl_volcano <- VolcanoPlot(de_results = microglia_1d_ctrl,
                                         label_by_fc = 0.75,
                                         label_by_pval = 0.001,
                                         label_size = 5) +
  labs(title = 'S1d vs SCtrl',
       subtitle = 'Log2(fold-change) > 0 indicates higher expression in S1d') +
  theme(plot.title = element_text(size = 16, color = 'black'))
microglia_7d_1d_volcano <- VolcanoPlot(de_results = microglia_7d_1d,
                                       label_by_fc = 0.75,
                                       label_by_pval = 0.001,
                                       label_size = 5) +
  labs(title = 'S7d vs S1d',
       subtitle = 'Log2(fold-change) > 0 indicates higher expression in S7d') +
  theme(plot.title = element_text(size = 16, color = 'black'))
microglia_volcano <- cowplot::plot_grid(microglia_1d_ctrl_volcano, microglia_7d_1d_volcano, ncol = 2)

# Changes in NSC gene expression after injury
nsc_1d_ctrl <- FindMarkers(object = svz,
                           subset.ident = 'NSC',
                           group.by = 'orig.ident',
                           ident.1 = 'S1d',
                           ident.2 = 'SCtrl')
nsc_7d_1d <- FindMarkers(object = svz,
                         subset.ident = 'NSC',
                         group.by = 'orig.ident',
                         ident.1 = 'S7d',
                         ident.2 = 'S1d')
nsc_1d_ctrl_volcano <- VolcanoPlot(de_results = nsc_1d_ctrl,
                                   label_by_fc = 0.5,
                                   label_by_pval = 0.001,
                                   label_size = 5,
                                   use_adj_pval = FALSE) +
  labs(title = 'S1d vs SCtrl',
       subtitle = 'Log2(fold-change) > 0 indicates higher expression in S1d') +
  theme(plot.title = element_text(size = 16, color = 'black'))
nsc_7d_1d_volcano <- VolcanoPlot(de_results = nsc_7d_1d,
                                 label_by_fc = 0.5,
                                 label_by_pval = 0.001,
                                 label_size = 5,
                                 use_adj_pval = FALSE) +
  labs(title = 'S7d vs S1d',
       subtitle = 'Log2(fold-change) > 0 indicates higher expression in S7d') +
  theme(plot.title = element_text(size = 16, color = 'black'))
nsc_volcano <- cowplot::plot_grid(nsc_1d_ctrl_volcano, nsc_7d_1d_volcano, ncol = 2)
de_genes_over_time <- cowplot::plot_grid(microglia_volcano, nsc_volcano, ncol = 1)
ggsave(filename = paste0(results_path, 'DEgenes_byTime_volcanoplot.tiff'), 
       plot = de_genes_over_time,
       device = 'tiff',
       height = 10, width = 12)
de_genes_over_time
```

<img src="20200731_LigandReceptorAnalysis_files/figure-html/DE_tests-1.png" style="display: block; margin: auto;" />



```r
top_genes_1 <- rownames(microglia_1d_ctrl[order(abs(microglia_1d_ctrl$avg_logFC), decreasing = TRUE),])[1:100]
top_genes_2 <- rownames(microglia_7d_1d[order(abs(microglia_7d_1d$avg_logFC), decreasing = TRUE),])[1:100]
top_genes_microglia <- union(top_genes_1, top_genes_2)

top_genes_1 <- rownames(nsc_1d_ctrl[order(abs(nsc_1d_ctrl$avg_logFC), decreasing = TRUE),])[1:100]
top_genes_2 <- rownames(nsc_7d_1d[order(abs(nsc_7d_1d$avg_logFC), decreasing = TRUE),])[1:100]
top_genes_nsc <- union(top_genes_1, top_genes_2)

top_genes <- union(top_genes_microglia, top_genes_nsc)
tmp_lr_results <- my_results[my_results$Ligand %in% top_genes | my_results$Receptor %in% top_genes,]

LR_plots1 <- plotLR(results = tmp_lr_results,
                    l_cells = 'Microglia',
                    r_cells = 'NSC',
                    split_along_y = TRUE)
ggsave(filename = paste0(results_path, 'LRplots_byDEgenes.tiff'),
       plot = LR_plots1,
       device = 'tiff',
       height = 4.5, width = 9)
```


```r
my_results_known <- my_results[my_results$Pair_name %in% lr_ref$Pair.Name[lr_ref$Pair.Source == 'known'],]
my_results_known <- my_results_known[my_results_known$pval < 0.05,]
my_results_known <- my_results_known[order(my_results_known$Score, decreasing = TRUE),]

LR_plots2 <- plotLR(results = my_results_known,
                    l_cells = 'Microglia',
                    r_cells = 'NSC',
                    split_along_y = TRUE)
ggsave(filename = paste0(results_path, 'LRplots_allFantomKnown.tiff'),
       plot = LR_plots2,
       device = 'tiff',
       height = 4.5, width = 9)
```




```r
sessionInfo()
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 18363)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] tibble_3.0.1  ggrepel_0.8.2 dplyr_1.0.0   ggplot2_3.3.2 Seurat_3.2.0 
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-148          RcppAnnoy_0.0.16      RColorBrewer_1.1-2   
##  [4] httr_1.4.1            sctransform_0.2.1     tools_4.0.2          
##  [7] R6_2.4.1              irlba_2.3.3           rpart_4.1-15         
## [10] KernSmooth_2.23-17    uwot_0.1.8            mgcv_1.8-31          
## [13] lazyeval_0.2.2        colorspace_1.4-1      withr_2.2.0          
## [16] tidyselect_1.1.0      gridExtra_2.3         compiler_4.0.2       
## [19] plotly_4.9.2.1        labeling_0.3          scales_1.1.1         
## [22] lmtest_0.9-37         spatstat.data_1.4-3   ggridges_0.5.2       
## [25] pbapply_1.4-2         rappdirs_0.3.1        spatstat_1.64-1      
## [28] goftest_1.2-2         stringr_1.4.0         digest_0.6.25        
## [31] spatstat.utils_1.17-0 rmarkdown_2.3         pkgconfig_2.0.3      
## [34] htmltools_0.5.0       fastmap_1.0.1         htmlwidgets_1.5.1    
## [37] rlang_0.4.6           shiny_1.5.0           farver_2.0.3         
## [40] generics_0.0.2        zoo_1.8-8             jsonlite_1.7.0       
## [43] ica_1.0-2             magrittr_1.5          patchwork_1.0.1      
## [46] Matrix_1.2-18         Rcpp_1.0.4.6          munsell_0.5.0        
## [49] abind_1.4-5           ape_5.4               reticulate_1.16      
## [52] lifecycle_0.2.0       stringi_1.4.6         yaml_2.2.1           
## [55] MASS_7.3-51.6         Rtsne_0.15            plyr_1.8.6           
## [58] grid_4.0.2            parallel_4.0.2        listenv_0.8.0        
## [61] promises_1.1.1        crayon_1.3.4          deldir_0.1-28        
## [64] miniUI_0.1.1.1        lattice_0.20-41       cowplot_1.0.0        
## [67] splines_4.0.2         tensor_1.5            knitr_1.29           
## [70] pillar_1.4.4          igraph_1.2.5          future.apply_1.6.0   
## [73] reshape2_1.4.4        codetools_0.2-16      leiden_0.3.3         
## [76] glue_1.4.1            evaluate_0.14         data.table_1.12.8    
## [79] vctrs_0.3.1           png_0.1-7             httpuv_1.5.4         
## [82] gtable_0.3.0          RANN_2.6.1            purrr_0.3.4          
## [85] polyclip_1.10-0       tidyr_1.1.0           future_1.17.0        
## [88] xfun_0.15             rsvd_1.0.3            mime_0.9             
## [91] xtable_1.8-4          later_1.1.0.1         survival_3.1-12      
## [94] viridisLite_0.3.0     cluster_2.1.0         globals_0.12.5       
## [97] fitdistrplus_1.1-1    ellipsis_0.3.1        ROCR_1.0-11
```
