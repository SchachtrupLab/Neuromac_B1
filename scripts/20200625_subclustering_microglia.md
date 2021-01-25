---
title: "Subclustering Microglia"
author: "James Choi"
date: "2020-06-25"
output:  
  html_document:
    keep_md: true
    code_folding: hide
    toc: true
---



### Background

We build upon the results of the previous cluster analysis of SVZ microglia and NSCs in a photothrombotic injury model of stroke. In this analysis, we take the extract the microglia from the full SVZ dataset and perform an additional round of clustering with the overall goal of uncovering potential subsets of microglia with distinct transcriptional signatures that were not identified in the initial SVZ analysis.  


### Loading libraries


```r
dir.create(path = results_outpath)
library('SingleCellExperiment')
library('Seurat')
library('dplyr')
library('ggplot2')
```


### Loading data

Based on marker gene expression (see previous analysis report), clusters 0, 1, 2, and 6 correspond to microglia while clusters 3, 4, and 5 correspond to NSCs.  



```r
svz <- readRDS(file = data_inpath)
DimPlot(svz)
```

<img src="20200625_subclustering_microglia_files/figure-html/load_data-1.png" style="display: block; margin: auto;" />


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
DimPlot(svz, label = TRUE, label.size = 6)
```

<img src="20200625_subclustering_microglia_files/figure-html/set_identities-1.png" style="display: block; margin: auto;" />

We extract the microglia and perform the standard single-cell clustering procedure as implemented in Seurat.  


```r
microglia <- svz[,svz$celltype == 'Microglia']
DefaultAssay(microglia) <- 'RNA'
```


### Normalization and PCA


```r
microglia <- microglia %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c('CC.difference')) %>%
  RunPCA(npcs = 50)
```


```r
elbow <- ElbowPlot(microglia, ndims = 50)
elbow <- cowplot::plot_grid(NULL, elbow, NULL, rel_widths = c(0.3, 1, 0.3), ncol = 3)
microglia <- JackStraw(microglia, dims = 30)
microglia <- ScoreJackStraw(microglia, dims = 1:30)
jackstraw <- JackStrawPlot(microglia, dims = 1:30)
jackstraw <- cowplot::plot_grid(NULL, jackstraw, NULL, rel_widths = c(0.3, 1, 0.3), ncol = 3)

npcs <- which(microglia@reductions$pca@jackstraw$overall.p.values[,2] < 1e-5)
npcs <- 1:10
total_var <- sum(matrixStats::rowVars(x = microglia[['RNA']]@scale.data))
pc_ev <- microglia[['pca']]@stdev^2
explained_var <- round(pc_ev/total_var * 100, digits = 1)

color_by <- 'orig.ident'
tmp <- cbind(microglia[['pca']]@cell.embeddings, microglia@meta.data)
biplot <- vector(mode = 'list', length = length(npcs))
for(i in 1:(length(npcs)-1)) {
  biplot[[i]] <- vector(mode = 'list', length = length(npcs))
  for(j in (i+1):length(npcs)) {
    dim_i <- paste0('PC_', npcs[i])
    dim_j <- paste0('PC_', npcs[j])
    biplot[[i]][[j]] <- tmp %>%
      ggplot(mapping = aes_string(y = dim_i, x = dim_j)) +
      geom_point(mapping = aes_string(color = color_by), size = 1) +
      geom_hline(yintercept = 0, linetype = 'dashed') +
      geom_vline(xintercept = 0, linetype = 'dashed') +
      theme(panel.background = element_rect(fill = NA, color = 'black'),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = 'none')
  }
  biplot[[i]][[i]] <- cowplot::ggdraw() + cowplot::draw_label(x = 0.5, y = 0.5, paste0(dim_i, '\n', explained_var[npcs[i]]))
  biplot[[i]] <- cowplot::plot_grid(plotlist = biplot[[i]], ncol = length(biplot[[i]]))
}
biplot[[length(npcs)]] <- vector(mode = 'list', length = length(npcs))
biplot[[length(npcs)]][[length(npcs)]] <- cowplot::ggdraw() + cowplot::draw_label(x = 0.5, y = 0.5, paste0('PC_', npcs[length(npcs)], '\n', explained_var[npcs[length(npcs)]]))
biplot[[length(npcs)]] <- cowplot::plot_grid(plotlist = biplot[[length(npcs)]], ncol = length(biplot[[length(npcs)]]))
biplot <- cowplot::plot_grid(plotlist = biplot, ncol = 1)
tmp <- tmp %>%
      ggplot(mapping = aes_string(x = dim_i, y = dim_j)) +
      geom_point(mapping = aes_string(color = color_by), size = 1) +
  guides(color = guide_legend(title = 'Cluster', override.aes = list(size = 4)))
tmp_legend <- cowplot::get_legend(tmp)
biplot <- cowplot::plot_grid(biplot, tmp_legend, rel_widths = c(1, 0.05))
pca_plots <- cowplot::plot_grid(elbow, jackstraw, biplot, ncol = 1, rel_heights = c(0.5, 0.7, 1.5))
pca_plots
```

<img src="20200625_subclustering_microglia_files/figure-html/pca_dim_selection-1.png" style="display: block; margin: auto;" />

```r
# ggsave(filename = paste0(results_outpath, 'PCA_biplot.tiff'), plot = tmp_plots, height = 14, width = 14.5)
```


```r
VizDimLoadings(microglia, dims = npcs, ncol = round(length(npcs)/2))
```

<img src="20200625_subclustering_microglia_files/figure-html/vizdimloadings-1.png" style="display: block; margin: auto;" />


### Dimensional Reduction and Clustering  



```r
# distance matrix in PCA-space
microglia_dist <- dist(x = microglia[['pca']]@cell.embeddings[,npcs])

# Assign values of k.param to iterate
k_range <- seq(from = 5, to = 35, by = 5) # sqrt(ncol(microglia)) == 25.6907
names(k_range) <- paste('k', k_range, sep ='_')
k_results <- vector(mode = 'list', length = length(k_range))
names(k_results) <- names(k_range)

# Cluster data through k.params using a set louvain resolution
for(k in 1:length(k_range)) {
  tmp <- FindNeighbors(object = microglia, k.param = k_range[k], dims = npcs)
  tmp <- FindClusters(tmp, resolution = 0.8, algorithm = 3)
  k_results[[names(k_range)[k]]] <- tmp@active.ident
}
k_results <- do.call(cbind, k_results)

# Compute silhouette coefficients for clustering results
silhouette_val <- vector(mode = 'list', length = ncol(k_results))
names(silhouette_val) <- colnames(k_results)
for(k in 1:ncol(k_results)) {
  silhouette_val[[colnames(k_results)[k]]] <- cluster::silhouette(x = k_results[,k], dist = microglia_dist)
}

# Plot silhouette coefficients
tiff(filename = paste0(results_outpath, 'SilhouettePlots_kparamTesting.tiff'), height = 2400, width = 2440, res = 220)
plot_dim <- ceiling(sqrt(length(silhouette_val)))
{
  par(mfrow = c(plot_dim, plot_dim))
  for(k in 1:length(silhouette_val)) {
    cat(plot(silhouette_val[[k]], border = NA, main = paste(names(silhouette_val)[k], '; res = 0.8'), cex = 1.5))
  }
}
dev.off()
```


```r
microglia <- FindNeighbors(object = microglia, dims = npcs, k.param = 20)
```


```r
# distance matrix in PCA-space
microglia_dist <- dist(x = microglia[['pca']]@cell.embeddings[,npcs])

res_range <- seq(from = 0.2, to = 1.0, by = 0.1) # Exclude 0.1 because results in single cluster
names(res_range) <- paste('res', res_range, sep ='_')
res_results <- vector(mode = 'list', length = length(res_range))
names(res_results) <- names(res_range)

# Cluster data through k.params using a set louvain resolution
for(r in 1:length(res_range)) {
  tmp <- FindClusters(microglia, resolution = res_range[r])
  res_results[[names(res_range)[r]]] <- tmp@active.ident
}
res_results <- do.call(cbind, res_results)

# Compute silhouette coefficients for clustering results
silhouette_val <- vector(mode = 'list', length = ncol(res_results))
names(silhouette_val) <- colnames(res_results)
for(r in 1:ncol(res_results)) {
  silhouette_val[[colnames(res_results)[r]]] <- cluster::silhouette(x = res_results[,r], dist = microglia_dist)
}

# Plot silhouette coefficients
tiff(filename = paste0(results_outpath, 'SilhouettePlots_resTesting.tiff'), height = 2400, width = 2440, res = 220)
plot_dim <- ceiling(sqrt(length(silhouette_val)))
{
  par(mfrow = c(ceiling(length(res_range)/plot_dim), plot_dim))
  for(r in 1:length(silhouette_val)) {
    cat(plot(silhouette_val[[r]], border = NA, main = paste('k.param = 20;', names(silhouette_val)[r]), cex = 1.5))
  }
}
dev.off()
```


```r
microglia <- FindNeighbors(object = microglia, dims = npcs, k.param = 20)
microglia <- FindClusters(object = microglia, resolution = 0.8)
```


```r
microglia <- RunUMAP(microglia, dims = npcs)
plot1 <- DimPlot(microglia, group.by = 'orig.ident', reduction = 'umap', pt.size = 2) + 
  theme(axis.line = element_line(size = 0.5),
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_text(size = 10), 
        panel.border = element_rect(color = 'black', fill = NA)) +
  labs(title = 'By time')
plot2 <- DimPlot(microglia, group.by = 'seurat_clusters', label = TRUE, label.size = 6, reduction = 'umap', pt.size = 2) + 
  theme(axis.line = element_line(size = 0.5),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_text(size = 10), 
        panel.border = element_rect(color = 'black', fill = NA)) +
  labs(title = 'By microglial subcluster')
plot3 <- DimPlot(microglia, group.by = 'full_clusters', label = TRUE, label.size = 6, reduction = 'umap', pt.size = 2) +
  theme(axis.line = element_line(size = 0.5),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_text(size = 10), 
        panel.border = element_rect(color = 'black', fill = NA)) +
  labs(title = 'By full dataset cluster')
plot <- cowplot::plot_grid(plot1, plot2, plot3, ncol = 3, rel_widths = c(1,0.925, 0.925))
plot
```

<img src="20200625_subclustering_microglia_files/figure-html/dimensional_reduction-1.png" style="display: block; margin: auto;" />

```r
ggsave(filename = paste0(results_outpath, 'microglia_umap.tiff'), plot = plot, device = 'tiff', height = 4.5, width = 16)
```

--------------------------------------------------------------------------------

### Microglia sub-cluster DE genes


```r
microglia_markers <- FindAllMarkers(object = microglia, assay = 'RNA', slot = 'data', test.use = 'wilcox', only.pos = TRUE)
saveRDS(object = microglia_markers, file = paste0(results_outpath, 'microglia_de_markers.rds'))
write.csv(file = paste0(results_outpath, 'microglia_de_markers.csv'), x = microglia_markers[microglia_markers$p_val_adj < 0.001,], quote = FALSE)
```


```r
markers <- readRDS(file = paste0(results_outpath, 'microglia_de_markers.rds'))
markers <- markers[markers$p_val_adj < 0.001, c(6,7,2,5,1,3,4)]
markers[c('p_val_adj','avg_logFC','p_val','pct.1','pct.2')] <- signif(markers[c('p_val_adj','avg_logFC','p_val','pct.1','pct.2')]
, digits = 3)
DT::datatable(data = markers, rownames = FALSE, extensions = c('Scroller','Buttons'), options = list(scroller = TRUE, deferRender = TRUE, scrollY = 350, fixedColumns = TRUE, autoWidth = TRUE, pageLength = 10, dom = 'Bfrtip', buttons = list(extend = 'collection', buttons = c('csv', 'excel', 'pdf'))))
```

<!--html_preserve--><div id="htmlwidget-f843a680a1fcd62a0dd1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f843a680a1fcd62a0dd1">{"x":{"filter":"none","extensions":["Scroller","Buttons"],"data":[["0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3"],["Gpr34","Marcks","Cx3cr1","P2ry12","Ttr","Tmsb4x","Gm42418","mt-Rnr2","Hexb","Sepp1","Actb","Enpp2","mt-Rnr1","Tmem119","Dpysl2","Olfml3","Selplg","Itgam","mt-Nd3","Qk","mt-Nd5","Basp1","Lrrc58","Lpcat2","Cst3","Rasgrp3","mt-Nd4","Col27a1","Sall1","Plxdc2","Gnai2","Siglech","Sparc","Bin1","Kctd12","Glul","Cyfip1","Aif1","Ccni","Gm10800","Gm21738","Gm37158","Gm26870","Gm21833","Gm37746","Gm10715","A430073D23Rik","Gm10719","Gm17555","Gm37954","Gm42519","Gm10717","Gm29055","Gm26905","Gm43238","RP23-341H2.4","Gm21750","Gm37900","Gm37903","Cd300lf","Srgn","Id2","Cd63","Ctsb","Lilr4b","Rbm3","Fth1","Mt1","Mt2","Ms4a6c","B2m","Klf13","Msr1","Ctsl","Ftl1","Ier3","C5ar1","Ms4a6d","C1qb","Ifi30","Cndp2","Ccdc86","Il4ra","Nfkbia","Lamp1","Fcgr3","Il1r2","Klf9","Gpr65","Tlr7","Pfkfb3","Tlr2","Socs3","Cd83","C3ar1","Plaur","Tnf","Smim3","Sgpl1","Ccl3","AI607873","Arhgdia","Eif4ebp1","Mapkapk2","Gm3839","H3f3b","Rps19","Cd14","Galc","Fcgr4","Sema6b","Nccrp1","Aqp3","Elf3","Lypd2","Egln3","Ly6g6c","Fam129b","Lad1","Mal","Emp1","Gsto1","Upk3bl","Aqp5","S100a10","Gsta4","Calml3","S100a6","Dsp","Mgarp","Anxa2","Krt6a","Evpl","Irf6","S100a11","Krt13","Mal2","Slurp1","Krt12","Cbr2","Krt5","Aldh3a1","Psca","Sfn","Hspb1","Gm9573","Muc4","Hspa1a","Asprv1","Adh7","Dsc2","Ifitm1","Pdzk1ip1","Aldh3b2","Tmprss11a","Pkp1","2310007B03Rik","Dapl1","Upk1b","Rab25","Krt16","S100a14","Tacstd2","Ppl","Rnf152","Trim29","Adh1","AI661453","Col17a1","Hopx","Ppp1r13l","Fam46b","Lgals3","Ybx3","Prss22","Ifi202b","Cdo1","Anxa1","Scel","Perp","Tmprss11g","Anxa8","Klf5","Cwh43","Capns1","Gsn","2200002D01Rik","Ifit1bl1","Pim1","Cyp4a12a","Actg1","Prss27","Sdc1","Fgfbp1","St6galnac1","Serpinb11","Dsg1a","Hbegf","Sptbn2","H60c","Zfp750","Plekhn1","Pax6","Plet1","Dsc3","Pard3","Ahnak","Cd44","Osbpl5","Capn5","Prkcdbp","Rapgefl1","Cdh1","AA986860","Fdft1","Sqle","Gprc5a","Cltb","Dgka","Cd24a","Net1","Lypd3","2610528A11Rik","Sox15","Fosl2","Fos","Krt4","Pdlim1","Slc5a1","St8sia6","Fosb","S100a16","Chit1","Sytl1","Eps8l2","Myo6","Muc15","Lmo7","Ablim1","Serpinb5","Actn4","Rarg","Ezr","Krt15","Ptgr1","Foxc1","Nr1d1","Ldlr","Cdkn1a","Wfs1","Akr1b7","Csta1","Epb41l1","Ephb3","Pglyrp1","Ehf","Aim1","Mall","Ptprf","Rhod","Ptrf","Prdx5","Notch3","Lama3","Cav2","Grhl1","Vsig10","Plec","Yap1","Ucp2","Pof1b","Slc5a8","Cldn4","Clic3","Mreg","Gltp","Phlda1","Enpp3","Dusp1","Spint2","Maff","Pkp3","Sdr39u1","Ptpn13","Sprr1a","Lynx1","Cxcl2","Dusp5","Cldn7","Sptan1","Galnt18","Scnn1a","Ckmt1","Aldh1a1","Kit","Defb1","Ier3","Otx1","Sik1","Alcam","Oasl1","Cldn23","Dnajb1","Hlf","Idi1","Amotl2","Ppp1r3c","Ovol1","Sc5d","Sorbs2","Arpc2","Kif21a","Rora","mt-Cytb","Tmem154","Arhgap23","Foxq1","Anxa11","Plekha6","Bace2","Atf3","Mboat1","Myl6","Dock1","Scin","Eif5a","Mast4","Muc20","Myh14","Urah","Ephx3","Ctnna1","Ifit3","Sh3rf1","Hras","Prss23","Scd2","Cysrt1","Serpinb1a","Cyr61"],[0.739,0.597,0.664,0.695,1.66,0.469,0.763,0.518,0.416,0.5,0.387,1.2,0.536,0.486,0.526,0.562,0.529,0.542,0.414,0.493,0.335,0.458,0.449,0.497,0.393,0.497,0.284,0.586,0.458,0.484,0.463,0.484,0.42,0.535,0.411,0.431,0.452,0.496,0.479,0.43,0.49,0.601,0.441,0.591,0.534,0.484,0.548,0.487,0.548,0.483,0.65,0.433,0.448,0.645,0.566,0.661,0.527,0.476,0.557,0.801,1.14,1.07,0.845,0.724,0.782,0.585,0.59,0.526,0.519,0.798,0.496,0.685,0.643,0.702,0.589,0.822,0.754,0.735,0.348,0.559,0.526,0.54,0.645,0.666,0.39,0.524,0.721,0.483,0.64,0.613,0.773,0.589,0.647,0.567,0.6,0.381,0.731,0.421,0.503,1.16,0.36,0.439,0.466,0.369,0.531,0.332,0.399,0.636,0.301,0.587,0.332,1.71,1.49,1.4,2.05,1.61,2.04,1.57,1.35,1.88,2.31,2.68,1.42,1.5,1.47,2.25,2.19,2.26,1.83,1.96,1.6,2.59,1.35,0.898,2.09,2.4,1.51,2.01,2.52,2.11,2.45,2.46,1.5,1.68,1.95,1.4,1.34,1.32,2.37,1.16,1.27,1.99,1.25,0.758,0.824,0.944,0.977,1.63,1.18,1.11,0.641,0.779,1.13,1.01,0.832,0.891,1.73,0.692,0.68,1.47,0.878,0.65,1.18,1.3,0.611,0.859,0.784,0.685,0.724,1.26,0.706,1.05,0.589,0.797,1.42,1.39,0.576,0.756,1.12,1.71,1.31,0.806,0.875,0.658,0.572,0.444,0.846,1.04,0.776,0.633,0.502,0.413,0.777,1.24,0.73,0.542,1.01,1.09,0.503,0.469,0.861,0.511,0.749,0.609,0.524,0.657,0.581,0.988,0.712,1.01,0.795,0.594,1.25,0.445,0.713,1.16,1.24,0.919,0.461,0.783,0.816,0.861,0.442,0.263,0.656,0.75,0.457,0.448,0.661,0.551,0.953,0.468,0.734,1.84,0.745,0.397,0.583,0.634,0.968,0.598,0.317,0.266,0.614,0.735,0.874,0.531,0.422,0.464,0.586,0.484,0.449,0.907,0.274,0.356,0.277,0.468,0.342,0.697,0.46,0.697,0.334,0.3,0.351,0.433,0.256,0.666,0.572,0.415,0.749,0.678,0.332,0.382,0.282,0.267,0.608,0.55,0.725,0.384,0.599,0.587,0.325,0.418,0.343,0.503,0.283,0.404,0.644,0.278,0.268,0.44,0.311,0.345,0.889,0.452,0.512,0.437,0.367,0.401,0.543,0.362,0.691,0.586,0.385,0.439,0.357,0.445,0.292,0.447,0.294,0.395,0.822,0.311,0.785,0.505,0.426,0.683,0.685,0.314,0.411,0.407,0.374,0.789,0.52,0.272,0.649,0.388,0.874,0.612,0.272,0.463],[1.19e-013,9.67e-013,2.08e-011,3.75e-011,6.76e-010,8.69e-009,1.49e-008,5.42e-008,1.79e-007,2.58e-007,3.81e-007,4.06e-007,4.14e-007,4.83e-007,1.26e-006,4.02e-006,4.65e-006,4.94e-006,1.29e-005,1.36e-005,1.65e-005,2.03e-005,2.45e-005,2.56e-005,7.83e-005,9.47e-005,0.000108,0.000131,0.000145,0.000145,0.000146,0.000322,0.000434,0.00048,0.000589,0.000617,0.000694,0.000745,0.000781,3.03e-014,1.77e-010,7.89e-010,2e-009,4.08e-009,1.24e-008,1.22e-007,2.04e-007,3.03e-007,4.57e-007,1.91e-006,2.23e-006,2.9e-006,4.27e-006,6.68e-006,8.45e-006,9.21e-006,0.0001,0.000625,0.000844,5.55e-018,1.82e-017,2.65e-016,2.05e-015,1.8e-013,2.86e-013,6.6e-012,3.96e-011,4.58e-011,1.34e-010,2.54e-010,7.93e-010,1.51e-009,3.21e-009,7.41e-009,8.68e-009,1.44e-008,8.73e-008,8.87e-008,1.04e-007,1.29e-007,1.52e-007,5.46e-007,1.01e-006,1.41e-006,1.91e-006,2.47e-006,3.15e-006,3.64e-006,3.96e-006,5.93e-006,6.14e-006,7.95e-006,1.1e-005,2.62e-005,4.57e-005,8.28e-005,0.0001,0.000127,0.000207,0.00021,0.00023,0.000263,0.000294,0.000315,0.000345,0.000684,0.000813,0.000836,0.000884,0.000909,0.000932,9.08e-031,6.24e-029,1.3e-028,5.92e-028,1.11e-026,1.2e-026,1.37e-026,2.43e-026,6.56e-026,7.15e-026,1.38e-025,4.7e-025,5.22e-025,5.81e-025,7.1e-025,7.47e-025,7.87e-025,7.96e-025,1.35e-024,6.68e-024,8.02e-024,1.18e-023,2.15e-023,2.97e-023,3.77e-023,9.55e-023,2.07e-022,2.18e-022,2.19e-022,2.63e-022,1.48e-021,1.6e-021,1.68e-021,1.73e-021,2.94e-021,8.42e-021,2.37e-020,3.93e-020,3.99e-020,6.45e-020,9.07e-020,9.1e-020,1.13e-019,3.85e-019,4.43e-019,5.4e-019,6.33e-019,2.38e-018,5.13e-018,8.04e-018,9.04e-018,1.22e-017,1.27e-017,1.5e-017,2.09e-017,1.49e-016,1.77e-016,2.55e-016,3.19e-016,3.64e-016,5.04e-016,1.71e-015,1.81e-015,2.44e-015,3.71e-015,4.82e-015,4.86e-015,8.21e-015,1e-014,1.03e-014,1.18e-014,2.26e-014,3.94e-014,4.18e-014,5.27e-014,5.41e-014,5.48e-014,9.49e-014,1.19e-013,1.32e-013,1.52e-013,2.4e-013,6.11e-013,6.93e-013,7.25e-013,1.42e-012,2.26e-012,2.8e-012,2.91e-012,3.25e-012,3.45e-012,3.64e-012,3.78e-012,3.87e-012,4.59e-012,6.92e-012,8.35e-012,1.94e-011,2.02e-011,2.58e-011,3.61e-011,4.09e-011,7.04e-011,7.13e-011,7.14e-011,7.22e-011,8.23e-011,9.48e-011,9.51e-011,1.02e-010,1.27e-010,1.3e-010,1.45e-010,1.98e-010,2.24e-010,2.33e-010,4.51e-010,5.8e-010,6.68e-010,7.25e-010,1e-009,1.2e-009,1.87e-009,1.96e-009,2.4e-009,3.65e-009,3.93e-009,4.05e-009,4.27e-009,5.76e-009,5.86e-009,1.07e-008,1.44e-008,1.46e-008,1.56e-008,2.92e-008,3.44e-008,4.04e-008,4.24e-008,5.26e-008,6e-008,8.46e-008,8.55e-008,9.97e-008,1.06e-007,1.51e-007,1.55e-007,1.6e-007,1.87e-007,2.45e-007,2.92e-007,3.26e-007,3.44e-007,3.73e-007,4.09e-007,8.73e-007,8.88e-007,1.17e-006,1.23e-006,1.78e-006,1.78e-006,1.95e-006,1.99e-006,2.59e-006,2.62e-006,2.71e-006,3.13e-006,3.49e-006,4.43e-006,4.88e-006,6.52e-006,7.8e-006,1.13e-005,1.38e-005,1.43e-005,1.59e-005,1.66e-005,1.66e-005,1.75e-005,1.87e-005,1.97e-005,2.01e-005,2.77e-005,2.78e-005,4.86e-005,5.66e-005,5.78e-005,6.14e-005,6.15e-005,6.34e-005,6.82e-005,6.9e-005,7.09e-005,8.07e-005,8.29e-005,9.06e-005,9.67e-005,9.9e-005,0.000103,0.000111,0.000113,0.000126,0.000143,0.000157,0.000184,0.00021,0.000211,0.000215,0.000325,0.00036,0.000364,0.000364,0.000379,0.000382,0.000393,0.000398,0.000415,0.00042,0.000451,0.000489,0.000538,0.000561,0.000563,0.00059,0.000615,0.000623,0.000713,0.000772,0.000882],[3.48e-018,2.83e-017,6.08e-016,1.1e-015,1.98e-014,2.54e-013,4.37e-013,1.58e-012,5.22e-012,7.55e-012,1.11e-011,1.19e-011,1.21e-011,1.41e-011,3.68e-011,1.18e-010,1.36e-010,1.44e-010,3.78e-010,3.98e-010,4.84e-010,5.95e-010,7.15e-010,7.49e-010,2.29e-009,2.77e-009,3.17e-009,3.83e-009,4.24e-009,4.24e-009,4.28e-009,9.41e-009,1.27e-008,1.4e-008,1.72e-008,1.8e-008,2.03e-008,2.18e-008,2.28e-008,8.87e-019,5.17e-015,2.31e-014,5.86e-014,1.19e-013,3.63e-013,3.56e-012,5.97e-012,8.87e-012,1.33e-011,5.57e-011,6.53e-011,8.49e-011,1.25e-010,1.95e-010,2.47e-010,2.69e-010,2.94e-009,1.83e-008,2.47e-008,1.62e-022,5.31e-022,7.74e-021,5.99e-020,5.27e-018,8.37e-018,1.93e-016,1.16e-015,1.34e-015,3.92e-015,7.44e-015,2.32e-014,4.42e-014,9.38e-014,2.17e-013,2.54e-013,4.21e-013,2.55e-012,2.59e-012,3.04e-012,3.77e-012,4.45e-012,1.6e-011,2.95e-011,4.12e-011,5.58e-011,7.21e-011,9.21e-011,1.06e-010,1.16e-010,1.73e-010,1.8e-010,2.32e-010,3.21e-010,7.65e-010,1.34e-009,2.42e-009,2.92e-009,3.7e-009,6.06e-009,6.15e-009,6.73e-009,7.7e-009,8.6e-009,9.22e-009,1.01e-008,2e-008,2.38e-008,2.44e-008,2.58e-008,2.66e-008,2.72e-008,2.66e-035,1.83e-033,3.79e-033,1.73e-032,3.23e-031,3.5e-031,4.01e-031,7.11e-031,1.92e-030,2.09e-030,4.04e-030,1.37e-029,1.53e-029,1.7e-029,2.07e-029,2.18e-029,2.3e-029,2.33e-029,3.94e-029,1.95e-028,2.35e-028,3.45e-028,6.28e-028,8.69e-028,1.1e-027,2.79e-027,6.06e-027,6.36e-027,6.4e-027,7.7e-027,4.32e-026,4.69e-026,4.91e-026,5.06e-026,8.61e-026,2.46e-025,6.94e-025,1.15e-024,1.17e-024,1.89e-024,2.65e-024,2.66e-024,3.31e-024,1.13e-023,1.29e-023,1.58e-023,1.85e-023,6.96e-023,1.5e-022,2.35e-022,2.64e-022,3.57e-022,3.72e-022,4.4e-022,6.1e-022,4.36e-021,5.16e-021,7.47e-021,9.32e-021,1.06e-020,1.47e-020,5e-020,5.29e-020,7.14e-020,1.09e-019,1.41e-019,1.42e-019,2.4e-019,2.93e-019,3e-019,3.45e-019,6.6e-019,1.15e-018,1.22e-018,1.54e-018,1.58e-018,1.6e-018,2.78e-018,3.48e-018,3.86e-018,4.43e-018,7.02e-018,1.79e-017,2.02e-017,2.12e-017,4.14e-017,6.61e-017,8.17e-017,8.51e-017,9.52e-017,1.01e-016,1.07e-016,1.11e-016,1.13e-016,1.34e-016,2.02e-016,2.44e-016,5.67e-016,5.9e-016,7.56e-016,1.06e-015,1.2e-015,2.06e-015,2.08e-015,2.09e-015,2.11e-015,2.41e-015,2.77e-015,2.78e-015,2.97e-015,3.72e-015,3.79e-015,4.24e-015,5.8e-015,6.55e-015,6.82e-015,1.32e-014,1.7e-014,1.95e-014,2.12e-014,2.93e-014,3.5e-014,5.47e-014,5.72e-014,7.01e-014,1.07e-013,1.15e-013,1.18e-013,1.25e-013,1.68e-013,1.71e-013,3.14e-013,4.21e-013,4.28e-013,4.56e-013,8.53e-013,1.01e-012,1.18e-012,1.24e-012,1.54e-012,1.76e-012,2.47e-012,2.5e-012,2.91e-012,3.1e-012,4.41e-012,4.53e-012,4.67e-012,5.48e-012,7.17e-012,8.53e-012,9.52e-012,1.01e-011,1.09e-011,1.2e-011,2.55e-011,2.6e-011,3.42e-011,3.6e-011,5.21e-011,5.21e-011,5.71e-011,5.8e-011,7.57e-011,7.66e-011,7.92e-011,9.14e-011,1.02e-010,1.3e-010,1.43e-010,1.91e-010,2.28e-010,3.29e-010,4.03e-010,4.18e-010,4.65e-010,4.84e-010,4.86e-010,5.11e-010,5.47e-010,5.75e-010,5.88e-010,8.1e-010,8.12e-010,1.42e-009,1.65e-009,1.69e-009,1.8e-009,1.8e-009,1.85e-009,1.99e-009,2.02e-009,2.07e-009,2.36e-009,2.42e-009,2.65e-009,2.83e-009,2.9e-009,3.02e-009,3.23e-009,3.3e-009,3.7e-009,4.18e-009,4.58e-009,5.37e-009,6.13e-009,6.16e-009,6.29e-009,9.5e-009,1.05e-008,1.06e-008,1.06e-008,1.11e-008,1.12e-008,1.15e-008,1.16e-008,1.21e-008,1.23e-008,1.32e-008,1.43e-008,1.57e-008,1.64e-008,1.65e-008,1.72e-008,1.8e-008,1.82e-008,2.09e-008,2.26e-008,2.58e-008],[0.881,0.963,0.948,0.91,0.731,1,0.993,1,0.978,0.918,1,0.53,1,0.881,0.612,0.896,0.955,0.828,0.993,0.881,1,0.903,0.985,0.843,1,0.694,1,0.56,0.701,0.776,0.925,0.858,0.963,0.627,0.873,0.933,0.806,0.634,0.657,1,1,0.992,1,0.992,0.992,0.992,0.975,0.983,0.983,1,0.959,1,1,0.95,0.835,0.926,0.942,0.934,0.826,0.605,0.912,0.807,0.904,0.982,0.43,0.763,1,0.895,0.719,0.518,1,0.772,0.325,0.974,0.921,0.623,0.491,0.544,1,0.649,0.553,0.412,0.544,0.886,0.956,0.974,0.272,0.43,0.421,0.746,0.596,0.482,0.614,0.351,0.763,0.254,0.36,0.316,0.658,0.474,0.211,0.789,0.351,0.474,0.868,0.939,0.807,0.807,0.254,0.377,0.202,0.788,0.75,0.731,0.923,0.769,0.923,0.923,0.673,0.904,0.981,1,0.788,0.769,0.846,1,0.962,0.981,0.808,0.923,0.885,1,0.692,0.481,0.962,0.942,0.692,0.923,1,0.923,1,1,0.635,0.75,0.885,0.692,0.654,0.769,0.962,0.558,0.577,0.827,0.596,0.462,0.404,0.5,0.519,0.712,0.635,0.519,0.308,0.404,0.558,0.558,0.462,0.462,0.769,0.365,0.404,0.712,0.462,0.346,0.654,0.788,0.346,0.442,0.385,0.462,0.385,0.596,0.346,0.5,0.308,0.365,0.942,0.731,0.288,0.327,0.731,0.75,1,0.269,0.404,0.327,0.327,0.231,0.308,0.346,0.404,0.346,0.269,0.269,0.404,0.712,0.404,0.365,0.519,0.5,0.327,0.308,0.404,0.327,0.385,0.25,0.385,0.327,0.288,0.731,0.481,0.712,0.442,0.327,0.654,0.308,0.538,0.75,0.385,0.385,0.212,0.423,0.346,0.481,0.192,0.173,0.327,0.481,0.269,0.231,0.442,0.327,0.731,0.346,0.615,0.558,0.385,0.269,0.346,0.288,0.769,0.308,0.173,0.154,0.308,0.385,0.327,0.327,0.269,0.231,0.288,0.231,0.308,0.827,0.192,0.192,0.192,0.269,0.173,0.346,0.25,0.981,0.154,0.154,0.231,0.192,0.135,0.712,0.327,0.192,0.673,0.385,0.212,0.231,0.173,0.173,0.269,0.25,0.385,0.212,0.231,0.558,0.192,0.212,0.25,0.212,0.192,0.173,0.731,0.135,0.135,0.25,0.173,0.173,0.5,0.288,0.269,0.212,0.192,0.192,0.327,0.192,0.942,0.308,0.327,1,0.327,0.269,0.154,0.231,0.154,0.173,0.404,0.173,0.923,0.308,0.212,0.827,0.404,0.173,0.212,0.25,0.192,0.519,0.346,0.192,0.442,0.231,0.731,0.212,0.135,0.25],[0.53,0.871,0.815,0.822,0.387,0.99,0.997,1,0.958,0.836,0.997,0.216,1,0.69,0.303,0.728,0.892,0.686,0.927,0.756,1,0.746,0.937,0.645,0.99,0.404,1,0.286,0.415,0.523,0.84,0.679,0.976,0.408,0.735,0.871,0.655,0.404,0.397,1,0.893,0.927,1,0.967,0.897,0.927,0.883,0.947,0.89,0.937,0.823,0.987,0.987,0.833,0.747,0.767,0.827,0.813,0.647,0.13,0.599,0.345,0.55,0.919,0.078,0.368,0.889,0.456,0.277,0.166,0.967,0.423,0.055,0.814,0.691,0.264,0.169,0.231,0.997,0.293,0.202,0.124,0.228,0.645,0.818,0.879,0.052,0.137,0.143,0.453,0.293,0.186,0.29,0.101,0.489,0.052,0.114,0.085,0.381,0.208,0.036,0.485,0.114,0.186,0.691,0.739,0.531,0.609,0.062,0.143,0.036,0.098,0.098,0.092,0.241,0.13,0.236,0.238,0.084,0.238,0.412,0.737,0.133,0.13,0.19,0.455,0.339,0.477,0.173,0.274,0.247,0.699,0.108,0.03,0.461,0.317,0.114,0.374,0.713,0.344,0.566,0.669,0.087,0.157,0.312,0.127,0.103,0.163,0.466,0.065,0.076,0.252,0.087,0.038,0.024,0.051,0.057,0.16,0.108,0.062,0.008,0.027,0.079,0.079,0.043,0.046,0.233,0.022,0.033,0.173,0.051,0.019,0.144,0.255,0.022,0.049,0.033,0.054,0.033,0.125,0.024,0.073,0.016,0.03,0.61,0.222,0.014,0.022,0.23,0.257,0.873,0.011,0.046,0.024,0.024,0.005,0.022,0.033,0.049,0.033,0.014,0.014,0.049,0.217,0.049,0.038,0.103,0.095,0.03,0.024,0.057,0.03,0.049,0.014,0.049,0.033,0.022,0.271,0.089,0.238,0.076,0.033,0.209,0.027,0.117,0.276,0.054,0.057,0.008,0.07,0.043,0.098,0.005,0.003,0.038,0.103,0.022,0.014,0.081,0.041,0.287,0.046,0.176,0.173,0.065,0.024,0.049,0.033,0.336,0.038,0.005,0.003,0.041,0.07,0.049,0.046,0.03,0.019,0.035,0.019,0.041,0.442,0.011,0.011,0.011,0.03,0.008,0.06,0.027,0.889,0.005,0.005,0.022,0.014,0.003,0.304,0.054,0.014,0.266,0.081,0.019,0.024,0.011,0.011,0.038,0.033,0.084,0.022,0.027,0.184,0.016,0.022,0.033,0.022,0.016,0.014,0.309,0.005,0.005,0.035,0.014,0.014,0.173,0.049,0.043,0.024,0.019,0.019,0.065,0.019,0.875,0.06,0.065,1,0.065,0.046,0.011,0.033,0.011,0.016,0.114,0.016,0.71,0.062,0.027,0.547,0.117,0.016,0.027,0.041,0.022,0.192,0.081,0.022,0.144,0.035,0.423,0.03,0.008,0.043]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>cluster<\/th>\n      <th>gene<\/th>\n      <th>avg_logFC<\/th>\n      <th>p_val_adj<\/th>\n      <th>p_val<\/th>\n      <th>pct.1<\/th>\n      <th>pct.2<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scroller":true,"deferRender":true,"scrollY":350,"fixedColumns":true,"autoWidth":true,"pageLength":10,"dom":"Bfrtip","buttons":{"extend":"collection","buttons":["csv","excel","pdf"]},"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6]}],"order":[],"orderClasses":false}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->


### Cluster Dendrograms


```r
Idents(microglia) <- 'seurat_clusters' 
avg <- AverageExpression(object = microglia, assays = 'RNA', slot = 'scale.data')[[1]]
microglia.dist <- dist(x = t(avg))
microglia.tree <- hclust(d = microglia.dist, method = 'ward.D2')
microglia.dend <- dendsort::dendsort(as.dendrogram(microglia.tree, hang = 0.1))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
microglia.cols <- gg_color_hue(n = length(levels(microglia$seurat_clusters)))
dendextend::labels_colors(microglia.dend) <- microglia.cols[microglia.tree$labels][order.dendrogram(microglia.dend)]
dendextend::labels_colors(microglia.dend) <- microglia.cols
dendextend::labels_cex(microglia.dend) <- 0.75
tiff(filename = paste0(results_outpath, 'Cluster_dendrogram.tiff'), height = 630, width = 630)
{par(mar = c(2,2,2,2), cex = 3, lwd = 2); dendextend::plot_horiz.dendrogram(microglia.dend, side = FALSE, main = 'Cluster dendrogram')}
dev.off()
{par(mar = c(2,2,2,2), cex = 2, lwd = 2, cex.main = 1); dendextend::plot_horiz.dendrogram(microglia.dend, side = FALSE, main = 'Cluster dendrogram')}
```

<img src="20200625_subclustering_microglia_files/figure-html/dendrogram_subcluster-1.png" style="display: block; margin: auto;" />

```r
rm(avg, microglia.dist, microglia.tree, microglia.dend)
```

```
## png 
##   2
```


```r
counts <- data.frame(table(microglia$seurat_clusters, microglia$orig.ident))
names(counts) <- c('seurat_clusters','orig.ident','count')
numbers <- counts %>%
  ggplot(mapping = aes(x = orig.ident, y = count, group = seurat_clusters)) + 
  geom_bar(aes(fill = seurat_clusters), stat = 'identity', color = 'black', size = 0.75) +
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
  ggplot(mapping = aes(x = orig.ident, y = count, group = seurat_clusters)) + 
  geom_bar(aes(fill = seurat_clusters), stat = 'identity', color = 'black', size = 0.75, position = 'fill') +
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

pop.dynamics
```

<img src="20200625_subclustering_microglia_files/figure-html/microglia_population_dynamics-1.png" style="display: block; margin: auto;" />

```r
ggsave(filename = paste0(results_outpath, 'microglia_population_dynamics.tiff'), plot = pop.dynamics, device = 'tiff', height = 3.75, width = 8)
```


```r
tmp <- FetchData(object = microglia, vars = c('seurat_clusters', 'orig.ident', VariableFeatures(microglia)), slot = 'scale.data')
tmp_annotation <- tmp[c('seurat_clusters', 'orig.ident')]
tmp_dat <- tmp[sapply(tmp, is.numeric)]
pheatmap::pheatmap(tmp_dat)
microglia.dist <- dist(x = tmp_dat)
microglia.tree <- hclust(d = microglia.dist, method = 'ward.D2')
microglia.dend <- dendsort::dendsort(as.dendrogram(microglia.tree, hang = 0.1))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
tmp_order <- order.dendrogram(microglia.dend)
microglia.cols <- gg_color_hue(n = length(levels(microglia$seurat_clusters)))
tmp_clust <- plyr::mapvalues(x = rownames(tmp_annotation)[tmp_order],
                             from = rownames(tmp_annotation)[tmp_order],
                             to = as.character(tmp_annotation$seurat_clusters[tmp_order]))
dendextend::labels_colors(microglia.dend) <- microglia.cols[tmp_clust]
dendextend::labels_colors(microglia.dend) <- microglia.cols
dendextend::labels_cex(microglia.dend) <- 0.75
plot(microglia.tree)
```


```r
sessionInfo()
```

```
## R version 3.6.1 (2019-07-05)
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
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] ggplot2_3.3.0.9000          dplyr_0.8.5                
##  [3] Seurat_3.1.2                SingleCellExperiment_1.6.0 
##  [5] SummarizedExperiment_1.14.1 DelayedArray_0.10.0        
##  [7] BiocParallel_1.18.1         matrixStats_0.56.0         
##  [9] Biobase_2.44.0              GenomicRanges_1.36.1       
## [11] GenomeInfoDb_1.20.0         IRanges_2.18.3             
## [13] S4Vectors_0.22.1            BiocGenerics_0.30.0        
## 
## loaded via a namespace (and not attached):
##   [1] TH.data_1.0-10         Rtsne_0.15             colorspace_1.4-1      
##   [4] ellipsis_0.3.0         ggridges_0.5.2         dendsort_0.3.3        
##   [7] XVector_0.24.0         farver_2.0.3           leiden_0.3.3          
##  [10] listenv_0.8.0          npsurv_0.4-0           DT_0.13               
##  [13] ggrepel_0.8.2          RSpectra_0.16-0        fansi_0.4.1           
##  [16] mvtnorm_1.1-0          codetools_0.2-16       splines_3.6.1         
##  [19] R.methodsS3_1.8.0      mnormt_1.5-6           lsei_1.2-0            
##  [22] knitr_1.28             TFisher_0.2.0          jsonlite_1.6.1        
##  [25] ica_1.0-2              cluster_2.1.0          png_0.1-7             
##  [28] R.oo_1.23.0            uwot_0.1.8             sctransform_0.2.1     
##  [31] compiler_3.6.1         httr_1.4.1             lazyeval_0.2.2        
##  [34] assertthat_0.2.1       Matrix_1.2-17          cli_2.0.2             
##  [37] htmltools_0.4.0        tools_3.6.1            rsvd_1.0.3            
##  [40] igraph_1.2.5           gtable_0.3.0           glue_1.4.0            
##  [43] GenomeInfoDbData_1.2.1 reshape2_1.4.4         RANN_2.6.1            
##  [46] rappdirs_0.3.1         Rcpp_1.0.4.6           vctrs_0.2.4           
##  [49] multtest_2.40.0        gdata_2.18.0           ape_5.3               
##  [52] nlme_3.1-140           crosstalk_1.1.0.1      gbRd_0.4-11           
##  [55] lmtest_0.9-37          xfun_0.13              stringr_1.4.0         
##  [58] globals_0.12.5         lifecycle_0.2.0        irlba_2.3.3           
##  [61] gtools_3.8.2           dendextend_1.13.4      future_1.16.0         
##  [64] zlibbioc_1.30.0        MASS_7.3-51.4          zoo_1.8-7             
##  [67] scales_1.1.0           sandwich_2.5-1         RColorBrewer_1.1-2    
##  [70] yaml_2.2.1             gridExtra_2.3          reticulate_1.15       
##  [73] pbapply_1.4-2          stringi_1.4.6          mutoss_0.1-12         
##  [76] plotrix_3.7-7          caTools_1.18.0         bibtex_0.4.2.2        
##  [79] Rdpack_0.11-1          SDMTools_1.1-221.1     rlang_0.4.5           
##  [82] pkgconfig_2.0.3        bitops_1.0-6           evaluate_0.14         
##  [85] lattice_0.20-38        ROCR_1.0-7             purrr_0.3.3           
##  [88] labeling_0.3           htmlwidgets_1.5.1      cowplot_1.0.0         
##  [91] tidyselect_1.0.0       RcppAnnoy_0.0.16       plyr_1.8.6            
##  [94] magrittr_1.5           R6_2.4.1               gplots_3.0.3          
##  [97] multcomp_1.4-13        withr_2.1.2            pillar_1.4.3          
## [100] sn_1.6-1               fitdistrplus_1.0-14    survival_3.1-12       
## [103] RCurl_1.98-1.1         tsne_0.1-3             tibble_3.0.0          
## [106] future.apply_1.4.0     crayon_1.3.4           KernSmooth_2.23-15    
## [109] plotly_4.9.2.1         rmarkdown_2.1          viridis_0.5.1         
## [112] grid_3.6.1             data.table_1.12.8      metap_1.3             
## [115] digest_0.6.25          tidyr_1.0.2            numDeriv_2016.8-1.1   
## [118] R.utils_2.9.2          munsell_0.5.0          viridisLite_0.3.0
```
