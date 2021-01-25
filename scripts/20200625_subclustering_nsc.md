---
title: "Subclustering NSCs"
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

<img src="20200625_subclustering_nsc_files/figure-html/load_data-1.png" style="display: block; margin: auto;" />


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

<img src="20200625_subclustering_nsc_files/figure-html/set_identities-1.png" style="display: block; margin: auto;" />

We extract the NSCs and perform the standard single-cell clustering procedure as implemented in Seurat.  


```r
nsc <- svz[,svz$celltype == 'NSC']
DefaultAssay(nsc) <- 'RNA'
```


### Normalization and PCA


```r
nsc <- nsc %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = c('CC.difference')) %>%
  RunPCA(npcs = 50)
```


```r
elbow <- ElbowPlot(nsc, ndims = 50)
elbow <- cowplot::plot_grid(NULL, elbow, NULL, rel_widths = c(0.3, 1, 0.3), ncol = 3)
nsc <- JackStraw(nsc, dims = 30)
nsc <- ScoreJackStraw(nsc, dims = 1:30)
jackstraw <- JackStrawPlot(nsc, dims = 1:30)
jackstraw <- cowplot::plot_grid(NULL, jackstraw, NULL, rel_widths = c(0.3, 1, 0.3), ncol = 3)

npcs <- which(nsc@reductions$pca@jackstraw$overall.p.values[,2] < 1e-5)
# npcs <- 1:10
total_var <- sum(matrixStats::rowVars(x = nsc[['RNA']]@scale.data))
pc_ev <- nsc[['pca']]@stdev^2
explained_var <- round(pc_ev/total_var * 100, digits = 1)

color_by <- 'orig.ident'
tmp <- cbind(nsc[['pca']]@cell.embeddings, nsc@meta.data)
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

<img src="20200625_subclustering_nsc_files/figure-html/pca_dim_selection-1.png" style="display: block; margin: auto;" />

```r
# ggsave(filename = paste0(results_outpath, 'PCA_biplot.tiff'), plot = tmp_plots, height = 14, width = 14.5)
```


```r
VizDimLoadings(nsc, dims = npcs, ncol = round(length(npcs)/2))
```

<img src="20200625_subclustering_nsc_files/figure-html/vizdimloadings-1.png" style="display: block; margin: auto;" />


### Dimensional Reduction and Clustering  



```r
# distance matrix in PCA-space
nsc_dist <- dist(x = nsc[['pca']]@cell.embeddings[,npcs])

# Assign values of k.param to iterate
k_range <- seq(from = 5, to = 35, by = 5) # sqrt(ncol(nsc)) == 25.6907
names(k_range) <- paste('k', k_range, sep ='_')
k_results <- vector(mode = 'list', length = length(k_range))
names(k_results) <- names(k_range)

# Cluster data through k.params using a set louvain resolution
for(k in 1:length(k_range)) {
  tmp <- FindNeighbors(object = nsc, k.param = k_range[k], dims = npcs)
  tmp <- FindClusters(tmp, resolution = 0.8, algorithm = 3)
  k_results[[names(k_range)[k]]] <- tmp@active.ident
}
k_results <- do.call(cbind, k_results)

# Compute silhouette coefficients for clustering results
silhouette_val <- vector(mode = 'list', length = ncol(k_results))
names(silhouette_val) <- colnames(k_results)
for(k in 1:ncol(k_results)) {
  silhouette_val[[colnames(k_results)[k]]] <- cluster::silhouette(x = k_results[,k], dist = nsc_dist)
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
nsc <- FindNeighbors(object = nsc, dims = npcs, k.param = 20)
```


```r
# distance matrix in PCA-space
nsc_dist <- dist(x = nsc[['pca']]@cell.embeddings[,npcs])

res_range <- seq(from = 0.2, to = 1.0, by = 0.1) # Exclude 0.1 because results in single cluster
names(res_range) <- paste('res', res_range, sep ='_')
res_results <- vector(mode = 'list', length = length(res_range))
names(res_results) <- names(res_range)

# Cluster data through k.params using a set louvain resolution
for(r in 1:length(res_range)) {
  tmp <- FindClusters(nsc, resolution = res_range[r])
  res_results[[names(res_range)[r]]] <- tmp@active.ident
}
res_results <- do.call(cbind, res_results)

# Compute silhouette coefficients for clustering results
silhouette_val <- vector(mode = 'list', length = ncol(res_results))
names(silhouette_val) <- colnames(res_results)
for(r in 1:ncol(res_results)) {
  silhouette_val[[colnames(res_results)[r]]] <- cluster::silhouette(x = res_results[,r], dist = nsc_dist)
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
nsc <- FindNeighbors(object = nsc, dims = npcs, k.param = 20)
nsc <- FindClusters(object = nsc, resolution = 0.8)
```


```r
nsc <- RunUMAP(nsc, dims = npcs)
plot1 <- DimPlot(nsc, group.by = 'orig.ident', reduction = 'umap', pt.size = 2) + 
  theme(axis.line = element_line(size = 0.5),
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        axis.title = element_text(size = 10), 
        panel.border = element_rect(color = 'black', fill = NA)) +
  labs(title = 'By time')
plot2 <- DimPlot(nsc, group.by = 'seurat_clusters', label = TRUE, label.size = 6, reduction = 'umap', pt.size = 2) + 
  theme(axis.line = element_line(size = 0.5),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_text(size = 10), 
        panel.border = element_rect(color = 'black', fill = NA)) +
  labs(title = 'By NSC subcluster')
plot3 <- DimPlot(nsc, group.by = 'full_clusters', label = TRUE, label.size = 6, reduction = 'umap', pt.size = 2) +
  theme(axis.line = element_line(size = 0.5),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_text(size = 10), 
        panel.border = element_rect(color = 'black', fill = NA)) +
  labs(title = 'By full dataset cluster')
plot <- cowplot::plot_grid(plot1, plot2, plot3, ncol = 3, rel_widths = c(1,0.925, 0.925))
plot
```

<img src="20200625_subclustering_nsc_files/figure-html/dimensional_reduction-1.png" style="display: block; margin: auto;" />

```r
ggsave(filename = paste0(results_outpath, 'nsc_umap.tiff'), plot = plot, device = 'tiff', height = 4.5, width = 16)
```

--------------------------------------------------------------------------------

### NSC sub-cluster DE genes


```r
nsc_markers <- FindAllMarkers(object = nsc, assay = 'RNA', slot = 'data', test.use = 'wilcox', only.pos = TRUE)
saveRDS(object = nsc_markers, file = paste0(results_outpath, 'nsc_de_markers.rds'))
write.csv(file = paste0(results_outpath, 'nsc_de_markers.csv'), x = nsc_markers[nsc_markers$p_val_adj < 0.001,], quote = FALSE)
```


```r
markers <- readRDS(file = paste0(results_outpath, 'nsc_de_markers.rds'))
markers <- markers[markers$p_val_adj < 0.001, c(6,7,2,5,1,3,4)]
markers[c('p_val_adj','avg_logFC','p_val','pct.1','pct.2')] <- signif(markers[c('p_val_adj','avg_logFC','p_val','pct.1','pct.2')]
, digits = 3)
DT::datatable(data = markers, rownames = FALSE, extensions = c('Scroller','Buttons'), options = list(scroller = TRUE, deferRender = TRUE, scrollY = 350, fixedColumns = TRUE, autoWidth = TRUE, pageLength = 10, dom = 'Bfrtip', buttons = list(extend = 'collection', buttons = c('csv', 'excel', 'pdf'))))
```

<!--html_preserve--><div id="htmlwidget-f843a680a1fcd62a0dd1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f843a680a1fcd62a0dd1">{"x":{"filter":"none","extensions":["Scroller","Buttons"],"data":[["0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","0","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","2","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3","3"],["Stmn2","Nrep","Tubb3","Dcx","Basp1","Pfn2","Rtn1","Tiam2","Dlx6os1","Map1b","Igfbpl1","Gad2","Gad1","Cdk5r1","Cd24a","Tmsb10","Tubb2a","Nsg2","Ncam1","Sp9","Mir124a-1hg","Btg1","Sox4","Crmp1","Dpysl3","Mllt11","Sox11","Meis2","Aplp1","Celf4","Zfp704","Gm9844","Elavl3","Apc2","Dlx1","Fscn1","Dpysl5","Dlgap4","Ttc3","Nfib","Pbx1","Gng3","Nrxn3","Nnat","Islr2","D430019H16Rik","Rtn4","Nav1","Marcksl1","Gng2","Trp53i11","Shtn1","Stmn3","Mpped2","Uchl1","Atcay","Ntrk2","Prdx6","Gja1","Ttyh1","Mlc1","Mt1","Clu","Slc1a3","Apoe","Cpe","Pbxip1","Fam107a","Prnp","Atp1a2","Plpp3","Aldoc","Htra1","Acsl3","Timp3","Gstm1","Dbi","Atp1b2","Itm2b","Sparc","Zfp36l1","Mt2","Itm2c","Glul","Slc1a2","Id3","F3","Cst3","Pla2g7","Sdc4","S1pr1","Bcan","Scd2","Fxyd1","Sparcl1","Dtna","Tsc22d4","Gpr37l1","Aldoa","Gabrb1","Ndrg2","Fam213a","Tppp3","Ddah1","Cldn10","Dkk3","Gpm6b","Sox9","Prex2","Ednrb","Mgll","Cd63","Hepacam","Tspan3","Scara3","Glud1","Serpine2","S100a6","Tmbim6","Nr1d1","Aqp4","Cnn3","Id2","Add3","Arhgef4","Fgfr1","Tspan7","Ntsr2","Klf9","Cspg5","Chst2","Ifitm3","Pdpn","Cd81","Tmem47","Klf15","Lsamp","Sorbs1","Asrgl1","Anxa5","Pttg1ip","Slc9a3r1","Fbxo2","C4a","Ramp1","Ngef","Cbs","Usp2","S100a1","Prelp","Mgst1","Ppp1r3c","Spsb1","Fgfr3","Tprkb","Thbs4","Foxo1","Ptprf","Tnfrsf19","Hmgcs2","Chchd10","Abhd4","Ptn","Ggta1","Cyp2d22","Acsbg1","Megf10","Arhgap5","Prodh","Anxa2","Enkur","Fnbp1","Pdlim4","Itih5","Papss2","Ccdc80","Naaa","Sfxn5","Cav2","Slc7a2","Aldh1l1","Oat","Gbp2","Ccdc74a","Farp1","Eps15","Pla2g16","Nrbp2","Sdc2","Rarres2","Crim1","Gstm5","Camk2n1","Reep5","Il6st","Lgi4","Grina","Plcd4","Psap","Crb2","Cd9","Aldh2","Hspa2","Id4","Gfap","Rorb","Flnc","Epas1","Entpd2","Sat1","Ctsl","Kcnj10","Ckb","Fth1","Lhfp","Acot1","Tsc22d3","Pltp","Osmr","Tspan18","Timp2","Gna13","Plagl1","Cd82","Mmd2","Fhad1","Qk","Degs1","Rgcc","Mical2","Slco1c1","Fjx1","Enpp5","Csdc2","Scg3","Neat1","Acaa2","Plekhb1","Abca1","Mfsd2a","Ppp2r2b","Nkain4","Oaf","Vcam1","Notch2","Ptrf","Acsl6","Sepp1","Abhd3","Gm973","Rabac1","Vim","Ass1","Kcnk1","Rora","Sik3","Serpinh1","Fam183b","Slc4a4","Dlgap1","Slc38a3","Rras","Lmbrd1","Mid1ip1","Cyr61","Gpld1","Fas","Gramd3","Rsph1","Gjb6","Wwtr1","Cfap126","Gldc","Laptm4a","Ntm","Saraf","Aldh6a1","Sash1","Arap2","Ppp1r1b","Tcf7l2","Tnfrsf1a","Nme5","Cyp4f15","Gfra1","Mro","Slc14a1","Iqca","Fam174a","Dab1","Cadm1","Phgdh","Pon2","Trps1","Mfap3l","Phkg1","Bcar3","Daam2","Slc16a2","Nrxn1","Tob2","S100a13","Cmtm5","Foxj1","Smpdl3a","Cfap54","Acot11","Bag3","Mir22hg","Psph","2410004P03Rik","Slc22a23","Lap3","Adhfe1","Sfrp1","Gas6","Slc2a1","Cfap100","2810417H13Rik","Hmgn2","Cdca7","Rrm2","Rfc5","Insm1","Mki67","H2afz","Fbxo5","Stmn1","H2afx","Tmpo","Ube2s","Tubb5","Nasp","Gm10282","Ptma","Smc2","Hnrnpab","Dck","Lmo3","Hmgb2","Ube2ql1","Sox11","Dhfr","Ezh2","Pcna","Kif11","Ccnd2","Smc4","Whsc1","Elavl4","Arx","Prox1","Lmnb1","Pcna-ps2","Stil","Dlx2","Dek","Ticrr","Nup62","Top2a","Atad2","Mdk","Anp32b","Cdk1","H2afv","Cenpm","Mcm5","Ckap2l","Cit","Tra2a","Ska2","Hmgn1","Hat1","Fam64a","Rrm1","Hjurp","Pbk","Foxm1","Hes6","Ascl1","Egfr","Thrsp","Uhrf1","Gas1","Lima1","Rpsa","Rplp1","Tipin","Ndufc2","Ckb","Hmgb1","Etv5","Naa50","Mcm3","Btg2","Lix1","Dut","Nptx2","Odc1","Ppia","Rpl18","Rpl10a","Gm8225","Rad51ap1","Rpl3","Eif5a","Dek","Dll4","Ranbp1","Tyro3","Hmgb2","Rpl8","mt-Co2","Ptma","Nbl1","Mfge8","Mcm7","Incenp","Dll3","Rps3a1","Rbbp7","Rpl41","Cenph","2700094K13Rik","Nup210","Cdk6","Jag1","mt-Co1","Rpl4","Nme1","Rbm3","Eef1b2","Ran","Aurkb","Rbl1","Mcm2","Cbx5","Rcc1","Impdh2","Efhd2","Car3","Hdgf","Dctpp1","Ddias","Ptbp1","Nup93","Gm8730","Nap1l1","Cad","Fabp7","mt-Co3","Pcna-ps2","Pcna","Rps8","Rad51","Tk1","Pabpc1","Hells","Sox8","Neil3","mt-Nd1","Anp32b","Pmf1","Ncald","Gm14207","Rps15","Sema5b","Gapdh","Nfic","Lzts1","Cdk2ap1","Rps11","Rps16","Fn1","Rps20","D030055H07Rik","Ddx39b","Dtl","E2f8","AI506816","Eri2","Psmg2","Iqgap3","Rpl36a","Rpl13","Ccnd1","Snrpb","Gm10073","Rps10","Nr2f6","Dll1","Mapkapk2","Rpl12","mt-Atp6","Wdr76","Wdr18","Ybx1","Gm10076","Rrm1","Usp39","Nop56","Rplp0","Mcm6","Abce1","Pold2","Cenpk","Shmt1","Wdr89","Mrps31","Rps13","2810417H13Rik","Rpa2","H2afz","Tbc1d31","Mtap","Gm6793","Rps3","Tmem97","Chchd3","Slbp","Nup205","Atf1","Lsm14a","Rps24","Tubb6","Usp1","Gm3839","Mrpl52","Mcm5","Fkbp3","Cdt1","Mybl2","Kars"],[1.56,1.48,1.26,1.36,1.17,1.21,0.984,1.23,1.2,1.23,0.977,0.972,0.991,0.875,0.947,0.908,0.732,0.917,0.759,0.974,0.803,0.875,0.843,0.821,0.814,0.818,0.673,0.847,0.818,0.784,1.21,0.812,0.717,0.803,0.749,0.816,0.764,0.751,0.544,0.695,0.641,0.73,0.905,0.332,0.651,0.651,0.61,0.809,0.494,0.787,0.655,0.688,0.67,0.987,0.811,0.418,2.31,2.04,1.75,1.88,1.54,2.57,3.16,1.96,1.85,2.03,1.22,2.27,1.45,1.98,2.17,2.54,1.39,1.66,1.27,1.66,1.55,1.64,1.36,1.68,1.21,2.08,1.46,1.84,2.26,1.88,1.16,1.91,1.13,1.22,1.65,1.43,1.09,1.66,2.21,1.32,1.21,1.38,1.11,1.18,1.34,0.945,1.08,1.5,0.922,1.23,1.04,1.15,1.14,1.13,1.15,0.985,0.982,0.913,1.05,1.01,1.35,1.48,0.833,0.778,1.37,1.05,1.69,1,0.591,1.09,1.4,1.07,1.3,1.51,1.08,1.37,0.703,0.939,1.1,0.718,0.812,1.04,0.835,1.01,0.805,0.732,0.814,0.995,0.612,0.755,0.786,0.86,1.08,0.961,0.764,0.736,0.631,0.802,0.876,1.46,0.876,0.79,0.837,0.83,0.916,0.886,1.39,0.613,0.79,0.775,0.464,0.833,0.646,1.15,0.555,0.767,0.739,0.425,0.531,0.571,0.666,0.646,0.582,0.837,0.597,0.956,0.794,0.582,0.705,0.675,0.697,0.825,0.576,0.984,0.612,0.816,0.837,0.702,0.655,0.705,0.791,0.468,0.801,0.428,0.722,0.84,0.587,0.963,1.91,0.905,0.591,0.74,0.587,0.756,0.742,0.756,0.614,0.614,0.574,0.885,0.962,0.558,0.484,0.497,0.59,0.712,0.735,0.417,0.878,0.668,0.793,0.534,0.497,0.406,0.842,0.691,0.609,0.47,0.75,1.13,0.585,0.735,0.817,0.364,0.576,0.584,0.647,0.448,0.789,0.739,0.929,0.856,0.386,0.541,0.65,1,0.433,0.748,0.67,0.585,0.647,0.552,0.71,0.642,0.589,0.539,0.622,0.839,0.807,0.564,0.46,0.494,0.446,0.439,0.803,0.54,0.594,0.646,0.46,0.523,0.666,0.668,0.583,0.546,0.676,0.524,0.42,0.436,0.476,0.386,0.663,0.412,0.48,0.582,0.693,0.839,0.621,0.573,0.291,0.382,0.467,0.491,0.714,0.98,0.598,0.374,0.538,0.742,0.454,0.515,0.452,0.395,0.323,0.638,0.525,0.485,0.632,0.576,1.05,0.624,0.581,0.351,1.29,0.867,0.877,0.984,0.53,0.872,1.06,0.817,0.686,0.736,0.908,0.835,0.822,0.718,0.767,0.772,0.545,0.81,0.631,0.63,0.439,0.86,0.598,0.727,0.665,0.728,0.714,0.808,0.771,0.84,0.692,0.62,0.699,0.638,0.563,0.722,0.449,0.689,0.619,0.414,0.46,1.06,0.643,0.66,0.661,0.924,0.648,0.477,0.614,0.694,0.461,0.481,0.572,0.467,0.557,0.463,0.581,0.548,0.629,0.467,1.17,1.48,1.16,1.23,0.806,1.35,1.09,0.773,0.742,0.686,1,0.993,0.744,0.681,0.722,0.72,0.987,0.672,0.777,1.05,0.884,0.647,0.681,0.563,0.664,0.441,0.65,0.672,0.958,0.486,0.728,0.417,0.958,0.689,0.367,0.532,0.54,0.917,0.791,0.664,0.587,0.589,0.675,0.594,0.332,0.645,0.29,0.894,0.47,0.377,0.553,0.667,0.898,0.56,0.758,0.485,0.479,0.579,0.657,0.472,0.542,0.718,0.41,0.684,0.551,0.26,0.591,0.362,0.507,0.659,0.429,0.717,0.353,0.82,0.872,0.489,0.447,0.599,0.538,0.519,0.623,0.533,0.373,0.656,0.641,0.497,0.263,0.503,0.542,0.589,0.758,0.653,0.67,0.526,0.466,0.492,0.535,0.319,0.559,0.519,0.354,0.336,0.284,0.402,0.455,0.562,0.554,0.325,0.623,0.653,0.562,0.264,0.585,0.443,0.574,0.279,0.314,0.369,0.535,0.604,0.647,0.283,0.631,0.472,0.576,0.415,0.386,0.3,0.488,0.44,0.288,0.463,0.607,0.495,0.734,0.295,0.3,0.563,0.528,0.266,0.439,0.612,0.292,0.284,0.443,0.522,0.328,0.625,0.541,0.452,0.487,0.514,0.4,0.364,0.479],[9.76e-022,1.31e-019,5.17e-019,1.56e-018,2.23e-017,1.12e-016,1.25e-015,3.17e-015,1.1e-014,6.9e-013,9.95e-013,2.33e-011,1.16e-010,3.57e-010,5.79e-010,2.23e-009,3.23e-009,3.85e-009,4.17e-009,6.3e-009,6.74e-009,1.09e-008,1.62e-008,1.71e-008,7.06e-008,9.57e-008,8.34e-007,1.34e-006,1.46e-006,1.48e-006,2.25e-006,2.84e-006,2.97e-006,3.32e-006,4.23e-006,4.37e-006,6.06e-006,7.45e-006,1.41e-005,2.08e-005,2.29e-005,2.41e-005,4.88e-005,6.31e-005,6.69e-005,0.000105,0.000175,0.000192,0.000358,0.000516,0.000556,0.000591,0.000791,0.000923,0.00099,0.000995,1.62e-027,1.68e-027,3.19e-027,1.05e-024,1.55e-024,1.63e-024,1.81e-024,1.44e-023,2.02e-023,2.35e-023,2.44e-023,2.7e-023,8.09e-023,5.67e-022,6.68e-022,1.68e-021,2.29e-021,6.88e-021,1.3e-020,1.48e-020,1.9e-019,4.43e-019,1.86e-018,2.21e-018,2.44e-018,2.52e-018,7.14e-018,7.33e-018,1.68e-017,2.54e-017,2.67e-017,3.15e-017,4.74e-017,6.38e-017,3.86e-016,5.4e-016,7.64e-016,1.05e-015,1.13e-015,1.72e-015,2.79e-015,2.88e-015,4.5e-015,6.59e-015,7.6e-015,2.28e-014,2.33e-014,5.43e-014,5.86e-014,5.87e-014,6.06e-014,1.11e-013,2.43e-013,2.59e-013,2.74e-013,2.81e-013,3.68e-013,4.17e-013,4.92e-013,1.1e-012,2.16e-012,3.88e-012,3.95e-012,5.5e-012,9.94e-012,1.27e-011,1.42e-011,1.52e-011,1.65e-011,1.66e-011,2.03e-011,4.78e-011,5.13e-011,5.71e-011,7.49e-011,1.05e-010,1.28e-010,1.78e-010,2.58e-010,3.65e-010,4.59e-010,5.64e-010,7.41e-010,8.04e-010,9.63e-010,1.07e-009,1.18e-009,1.79e-009,1.79e-009,1.83e-009,2.26e-009,2.37e-009,2.42e-009,2.45e-009,2.56e-009,3.28e-009,3.52e-009,3.62e-009,3.78e-009,3.95e-009,4.31e-009,5.87e-009,6.67e-009,6.92e-009,6.97e-009,8.32e-009,8.9e-009,1.72e-008,1.76e-008,1.84e-008,1.92e-008,2.29e-008,2.29e-008,2.84e-008,2.92e-008,3.4e-008,3.42e-008,3.56e-008,3.58e-008,3.79e-008,3.85e-008,4.7e-008,4.84e-008,5.01e-008,5.27e-008,6.02e-008,6.02e-008,1.07e-007,1.1e-007,1.17e-007,1.62e-007,1.65e-007,1.79e-007,1.92e-007,1.97e-007,2.04e-007,2.15e-007,2.24e-007,2.49e-007,2.63e-007,3.35e-007,3.4e-007,3.75e-007,4.5e-007,4.91e-007,5.37e-007,6.55e-007,6.77e-007,6.78e-007,7.23e-007,8.68e-007,9.81e-007,1.06e-006,1.23e-006,1.38e-006,1.39e-006,1.45e-006,1.49e-006,1.53e-006,1.61e-006,1.72e-006,2.14e-006,2.14e-006,2.89e-006,3.23e-006,4.33e-006,4.91e-006,4.96e-006,5.3e-006,5.43e-006,5.85e-006,6.26e-006,6.92e-006,7.4e-006,1.01e-005,1.08e-005,1.24e-005,1.28e-005,1.28e-005,1.43e-005,1.46e-005,1.48e-005,1.64e-005,1.72e-005,1.86e-005,1.93e-005,2.05e-005,2.14e-005,2.19e-005,2.36e-005,2.72e-005,2.77e-005,2.86e-005,3.09e-005,3.22e-005,3.73e-005,3.85e-005,4.27e-005,4.51e-005,4.64e-005,4.66e-005,4.69e-005,4.75e-005,4.75e-005,5.52e-005,5.64e-005,6.06e-005,6.55e-005,7.47e-005,7.52e-005,7.87e-005,7.91e-005,8.39e-005,8.39e-005,9.25e-005,9.51e-005,0.000107,0.000109,0.00011,0.00011,0.000118,0.000119,0.00012,0.000132,0.00014,0.000144,0.000153,0.000181,0.000195,0.000207,0.000216,0.000216,0.00025,0.000272,0.0003,0.000332,0.000338,0.000357,0.00038,0.000392,0.000394,0.000403,0.000425,0.000431,0.000444,0.000489,0.000497,0.0005,0.000548,0.000564,0.00058,0.000594,0.000598,0.000604,0.000634,0.000675,0.000703,0.000716,0.000719,0.000917,0.000964,0.000985,4.44e-011,3.31e-010,5.19e-010,8.88e-010,2.08e-009,2.66e-009,7.32e-009,1.43e-008,2.32e-008,2.95e-008,4.94e-008,7.9e-008,8.73e-008,2.04e-007,2.19e-007,4.18e-007,8.18e-007,9.81e-007,1.27e-006,1.74e-006,2.99e-006,3.91e-006,4.36e-006,5.32e-006,5.81e-006,5.97e-006,6.55e-006,7.6e-006,7.78e-006,9.31e-006,1.05e-005,2.4e-005,2.45e-005,2.57e-005,3.1e-005,4e-005,5.69e-005,6.02e-005,6.11e-005,6.26e-005,6.49e-005,9.07e-005,0.000124,0.000142,0.000211,0.000221,0.000288,0.000304,0.000323,0.000325,0.000347,0.000352,0.000354,0.000445,0.000478,0.000598,0.000616,0.000694,0.000763,0.000971,1.07e-015,3.89e-015,6.03e-015,1.02e-014,3.21e-013,5.43e-012,5.68e-012,8.51e-011,3.72e-010,7.71e-010,7.73e-010,1.45e-009,1.74e-009,2.21e-009,3.77e-009,5.49e-009,5.62e-009,7.49e-009,9.94e-009,1.08e-008,1.34e-008,1.69e-008,1.73e-008,2.19e-008,2.74e-008,3.93e-008,4.03e-008,6.62e-008,7.94e-008,8.68e-008,8.87e-008,1.46e-007,1.66e-007,1.67e-007,1.74e-007,2.07e-007,3.57e-007,4.86e-007,5.05e-007,5.52e-007,6.62e-007,7.01e-007,7.28e-007,7.81e-007,7.82e-007,7.86e-007,8.23e-007,8.64e-007,9.41e-007,1e-006,1.5e-006,1.8e-006,2.01e-006,2.2e-006,2.28e-006,2.34e-006,2.96e-006,3.07e-006,3.19e-006,4.14e-006,4.2e-006,4.5e-006,4.71e-006,4.89e-006,6.59e-006,8.8e-006,9.57e-006,1.06e-005,1.08e-005,1.11e-005,1.36e-005,1.5e-005,1.67e-005,1.9e-005,1.91e-005,1.95e-005,1.96e-005,2.04e-005,2.32e-005,2.42e-005,2.45e-005,2.47e-005,2.54e-005,2.69e-005,2.87e-005,2.89e-005,2.92e-005,3.07e-005,3.11e-005,3.12e-005,3.28e-005,3.38e-005,3.47e-005,4.03e-005,5.28e-005,5.39e-005,5.87e-005,6.31e-005,6.73e-005,6.74e-005,7.8e-005,7.93e-005,8.41e-005,8.6e-005,9.33e-005,9.98e-005,0.000106,0.000106,0.000113,0.000116,0.000119,0.000119,0.000131,0.000134,0.000137,0.000143,0.000145,0.000148,0.000149,0.000153,0.000169,0.000174,0.000182,0.000198,0.000199,0.000209,0.000217,0.000218,0.000245,0.000307,0.000322,0.000325,0.000326,0.000348,0.000417,0.000514,0.000527,0.000531,0.000548,0.000549,0.000595,0.000598,0.0006,0.000659,0.000738,0.000742,0.000746,0.000751,0.000759,0.000769,0.000806,0.000905,0.000936,0.000967,0.000974],[2.85e-026,3.82e-024,1.51e-023,4.56e-023,6.53e-022,3.28e-021,3.66e-020,9.27e-020,3.21e-019,2.02e-017,2.91e-017,6.81e-016,3.39e-015,1.04e-014,1.69e-014,6.51e-014,9.45e-014,1.12e-013,1.22e-013,1.84e-013,1.97e-013,3.18e-013,4.72e-013,5.01e-013,2.07e-012,2.8e-012,2.44e-011,3.92e-011,4.26e-011,4.32e-011,6.57e-011,8.31e-011,8.7e-011,9.69e-011,1.24e-010,1.28e-010,1.77e-010,2.18e-010,4.11e-010,6.09e-010,6.69e-010,7.04e-010,1.43e-009,1.84e-009,1.95e-009,3.08e-009,5.13e-009,5.62e-009,1.05e-008,1.51e-008,1.63e-008,1.73e-008,2.31e-008,2.7e-008,2.9e-008,2.91e-008,4.73e-032,4.91e-032,9.32e-032,3.06e-029,4.52e-029,4.75e-029,5.3e-029,4.2e-028,5.91e-028,6.87e-028,7.15e-028,7.88e-028,2.37e-027,1.66e-026,1.95e-026,4.9e-026,6.69e-026,2.01e-025,3.79e-025,4.33e-025,5.55e-024,1.3e-023,5.45e-023,6.45e-023,7.13e-023,7.37e-023,2.09e-022,2.14e-022,4.91e-022,7.44e-022,7.79e-022,9.21e-022,1.39e-021,1.87e-021,1.13e-020,1.58e-020,2.23e-020,3.06e-020,3.31e-020,5.04e-020,8.15e-020,8.42e-020,1.32e-019,1.93e-019,2.22e-019,6.66e-019,6.82e-019,1.59e-018,1.71e-018,1.72e-018,1.77e-018,3.24e-018,7.1e-018,7.58e-018,8.01e-018,8.23e-018,1.08e-017,1.22e-017,1.44e-017,3.22e-017,6.33e-017,1.13e-016,1.15e-016,1.61e-016,2.91e-016,3.7e-016,4.15e-016,4.44e-016,4.83e-016,4.85e-016,5.92e-016,1.4e-015,1.5e-015,1.67e-015,2.19e-015,3.06e-015,3.74e-015,5.21e-015,7.56e-015,1.07e-014,1.34e-014,1.65e-014,2.17e-014,2.35e-014,2.82e-014,3.13e-014,3.45e-014,5.25e-014,5.25e-014,5.34e-014,6.62e-014,6.91e-014,7.07e-014,7.16e-014,7.49e-014,9.59e-014,1.03e-013,1.06e-013,1.11e-013,1.16e-013,1.26e-013,1.72e-013,1.95e-013,2.02e-013,2.04e-013,2.43e-013,2.6e-013,5.04e-013,5.16e-013,5.39e-013,5.62e-013,6.7e-013,6.7e-013,8.31e-013,8.55e-013,9.94e-013,1e-012,1.04e-012,1.05e-012,1.11e-012,1.13e-012,1.37e-012,1.42e-012,1.47e-012,1.54e-012,1.76e-012,1.76e-012,3.12e-012,3.22e-012,3.41e-012,4.75e-012,4.82e-012,5.24e-012,5.61e-012,5.76e-012,5.98e-012,6.27e-012,6.56e-012,7.27e-012,7.68e-012,9.8e-012,9.95e-012,1.1e-011,1.32e-011,1.43e-011,1.57e-011,1.92e-011,1.98e-011,1.98e-011,2.11e-011,2.54e-011,2.87e-011,3.1e-011,3.59e-011,4.04e-011,4.07e-011,4.23e-011,4.34e-011,4.47e-011,4.7e-011,5.03e-011,6.26e-011,6.26e-011,8.44e-011,9.45e-011,1.27e-010,1.44e-010,1.45e-010,1.55e-010,1.59e-010,1.71e-010,1.83e-010,2.02e-010,2.16e-010,2.97e-010,3.17e-010,3.63e-010,3.74e-010,3.75e-010,4.17e-010,4.28e-010,4.33e-010,4.81e-010,5.02e-010,5.44e-010,5.66e-010,6e-010,6.26e-010,6.39e-010,6.91e-010,7.95e-010,8.08e-010,8.37e-010,9.05e-010,9.4e-010,1.09e-009,1.13e-009,1.25e-009,1.32e-009,1.36e-009,1.36e-009,1.37e-009,1.39e-009,1.39e-009,1.61e-009,1.65e-009,1.77e-009,1.91e-009,2.19e-009,2.2e-009,2.3e-009,2.31e-009,2.45e-009,2.45e-009,2.7e-009,2.78e-009,3.13e-009,3.2e-009,3.2e-009,3.22e-009,3.44e-009,3.49e-009,3.5e-009,3.85e-009,4.09e-009,4.21e-009,4.46e-009,5.3e-009,5.71e-009,6.05e-009,6.31e-009,6.31e-009,7.3e-009,7.95e-009,8.78e-009,9.7e-009,9.87e-009,1.04e-008,1.11e-008,1.15e-008,1.15e-008,1.18e-008,1.24e-008,1.26e-008,1.3e-008,1.43e-008,1.45e-008,1.46e-008,1.6e-008,1.65e-008,1.69e-008,1.74e-008,1.75e-008,1.77e-008,1.85e-008,1.97e-008,2.06e-008,2.09e-008,2.1e-008,2.68e-008,2.82e-008,2.88e-008,1.3e-015,9.68e-015,1.52e-014,2.6e-014,6.08e-014,7.77e-014,2.14e-013,4.17e-013,6.78e-013,8.61e-013,1.45e-012,2.31e-012,2.55e-012,5.96e-012,6.41e-012,1.22e-011,2.39e-011,2.87e-011,3.7e-011,5.07e-011,8.74e-011,1.14e-010,1.27e-010,1.55e-010,1.7e-010,1.75e-010,1.92e-010,2.22e-010,2.27e-010,2.72e-010,3.08e-010,7.02e-010,7.16e-010,7.52e-010,9.07e-010,1.17e-009,1.66e-009,1.76e-009,1.79e-009,1.83e-009,1.9e-009,2.65e-009,3.61e-009,4.15e-009,6.18e-009,6.45e-009,8.43e-009,8.88e-009,9.43e-009,9.51e-009,1.01e-008,1.03e-008,1.04e-008,1.3e-008,1.4e-008,1.75e-008,1.8e-008,2.03e-008,2.23e-008,2.84e-008,3.13e-020,1.14e-019,1.76e-019,2.97e-019,9.37e-018,1.59e-016,1.66e-016,2.49e-015,1.09e-014,2.25e-014,2.26e-014,4.24e-014,5.1e-014,6.47e-014,1.1e-013,1.6e-013,1.64e-013,2.19e-013,2.91e-013,3.16e-013,3.91e-013,4.95e-013,5.06e-013,6.41e-013,8e-013,1.15e-012,1.18e-012,1.94e-012,2.32e-012,2.54e-012,2.59e-012,4.27e-012,4.85e-012,4.88e-012,5.1e-012,6.06e-012,1.04e-011,1.42e-011,1.48e-011,1.61e-011,1.93e-011,2.05e-011,2.13e-011,2.28e-011,2.28e-011,2.3e-011,2.41e-011,2.53e-011,2.75e-011,2.93e-011,4.4e-011,5.26e-011,5.88e-011,6.44e-011,6.66e-011,6.83e-011,8.65e-011,8.97e-011,9.33e-011,1.21e-010,1.23e-010,1.32e-010,1.38e-010,1.43e-010,1.93e-010,2.57e-010,2.8e-010,3.09e-010,3.17e-010,3.24e-010,3.99e-010,4.39e-010,4.87e-010,5.56e-010,5.57e-010,5.69e-010,5.72e-010,5.96e-010,6.77e-010,7.08e-010,7.15e-010,7.22e-010,7.42e-010,7.87e-010,8.39e-010,8.44e-010,8.54e-010,8.97e-010,9.1e-010,9.12e-010,9.58e-010,9.87e-010,1.01e-009,1.18e-009,1.55e-009,1.58e-009,1.72e-009,1.84e-009,1.97e-009,1.97e-009,2.28e-009,2.32e-009,2.46e-009,2.52e-009,2.73e-009,2.92e-009,3.09e-009,3.09e-009,3.29e-009,3.39e-009,3.48e-009,3.49e-009,3.82e-009,3.91e-009,4e-009,4.17e-009,4.25e-009,4.33e-009,4.36e-009,4.46e-009,4.94e-009,5.1e-009,5.33e-009,5.79e-009,5.82e-009,6.11e-009,6.35e-009,6.36e-009,7.15e-009,8.97e-009,9.4e-009,9.51e-009,9.54e-009,1.02e-008,1.22e-008,1.5e-008,1.54e-008,1.55e-008,1.6e-008,1.6e-008,1.74e-008,1.75e-008,1.76e-008,1.93e-008,2.16e-008,2.17e-008,2.18e-008,2.19e-008,2.22e-008,2.25e-008,2.36e-008,2.65e-008,2.74e-008,2.83e-008,2.85e-008],[0.958,0.917,1,0.944,1,0.986,0.986,0.861,0.819,0.889,0.806,0.764,0.681,0.847,0.917,0.931,0.986,0.875,0.931,0.764,0.431,0.889,0.944,0.833,0.903,0.611,1,0.861,0.722,0.694,0.736,0.847,0.875,0.528,0.931,0.833,0.625,0.694,0.931,0.958,0.972,0.542,0.694,0.986,0.333,0.431,0.903,0.681,1,0.653,0.653,0.431,0.722,0.611,0.667,0.306,0.984,1,0.841,0.921,0.905,0.984,0.905,0.984,0.984,1,0.794,0.841,0.984,0.937,0.937,0.873,0.73,0.968,0.762,0.937,0.984,0.937,0.984,0.937,0.921,0.984,0.937,0.984,0.984,0.651,0.619,1,0.524,0.667,0.778,0.667,0.984,0.778,0.81,0.857,0.873,0.714,0.857,0.492,0.937,0.651,0.587,0.746,0.54,0.778,1,0.968,0.492,0.619,0.54,0.921,0.683,0.937,0.619,0.889,0.619,0.825,0.921,0.635,0.524,0.905,0.714,0.651,0.413,0.683,0.889,0.54,0.714,0.714,0.603,0.54,0.492,0.937,0.778,0.476,0.571,0.571,0.698,0.571,0.762,0.508,0.429,0.397,0.349,0.429,0.619,0.524,0.683,0.444,0.492,0.397,0.429,0.556,0.635,0.571,0.667,0.619,0.698,0.413,0.444,0.714,0.667,0.413,0.492,0.46,0.286,0.587,0.381,0.571,0.302,0.667,0.54,0.302,0.317,0.381,0.349,0.333,0.365,0.476,0.365,0.683,0.381,0.413,0.587,0.619,0.508,0.556,0.524,0.302,0.444,0.81,0.73,0.778,0.508,0.381,0.667,0.317,0.81,0.254,0.857,0.603,0.381,0.857,0.492,0.683,0.27,0.429,0.302,0.778,0.667,0.524,1,0.921,0.413,0.667,0.556,0.238,0.238,0.381,0.524,0.825,0.524,0.27,0.587,0.27,0.889,0.429,0.27,0.302,0.413,0.492,0.381,0.381,0.857,0.429,0.429,0.619,0.746,0.238,0.381,0.349,0.317,0.302,0.476,0.286,0.54,0.54,0.254,0.302,0.683,0.905,0.27,0.413,0.429,0.556,0.476,0.206,0.365,0.492,0.349,0.365,0.556,0.714,0.333,0.302,0.222,0.302,0.222,0.222,0.46,0.222,0.397,0.841,0.286,0.698,0.476,0.54,0.46,0.27,0.46,0.349,0.238,0.27,0.381,0.27,0.19,0.19,0.413,0.46,0.746,0.635,0.524,0.603,0.206,0.206,0.349,0.286,0.46,0.651,0.444,0.333,0.349,0.238,0.397,0.222,0.286,0.286,0.222,0.524,0.222,0.302,0.508,0.397,0.524,0.302,0.365,0.175,0.902,0.98,0.98,0.824,0.51,0.882,0.824,1,0.627,1,0.686,0.941,0.902,1,0.824,0.863,1,0.804,0.98,0.667,0.608,0.98,0.725,0.961,0.706,0.784,0.863,0.588,0.98,0.843,0.961,0.922,0.941,0.784,0.745,0.784,0.51,0.941,0.98,0.373,0.647,0.745,0.608,0.98,0.922,0.667,0.922,0.588,0.647,0.549,0.51,0.765,0.667,1,0.706,0.49,0.765,0.725,0.608,0.529,0.935,0.935,0.87,0.717,0.826,0.826,0.848,1,1,0.717,0.935,1,1,0.696,0.804,0.739,0.739,0.783,0.891,0.652,0.957,1,1,1,1,0.478,0.978,0.978,0.935,0.304,1,0.565,0.978,0.978,1,1,0.652,0.87,0.891,0.609,0.435,1,0.891,1,0.587,0.891,0.326,0.717,0.391,1,1,0.891,0.957,1,0.978,0.565,0.609,0.696,0.913,0.565,0.761,0.891,0.304,0.957,0.63,0.283,0.826,0.674,1,0.957,0.478,0.739,1,0.783,0.826,1,0.5,0.609,0.978,0.63,0.565,0.457,1,0.978,0.739,0.435,0.196,1,0.674,1,0.891,0.478,0.957,0.978,1,0.391,1,0.239,1,0.587,0.457,0.522,0.413,0.5,0.543,1,0.978,0.413,0.913,0.891,0.957,0.348,0.478,0.717,1,1,0.326,0.587,0.978,0.978,0.783,0.478,0.783,1,0.783,0.717,0.565,0.457,0.457,1,0.457,1,0.848,0.674,0.978,0.413,0.391,0.891,0.978,0.348,0.652,0.783,0.5,0.478,0.696,0.978,0.217,0.696,1,0.761,0.696,0.826,0.5,0.391,0.717],[0.356,0.375,0.538,0.356,0.638,0.8,0.838,0.325,0.256,0.619,0.288,0.269,0.181,0.438,0.619,0.606,0.85,0.538,0.706,0.281,0.05,0.625,0.738,0.381,0.562,0.206,0.769,0.706,0.388,0.275,0.5,0.544,0.55,0.175,0.619,0.538,0.25,0.375,0.825,0.906,0.95,0.194,0.325,0.881,0.05,0.112,0.769,0.375,0.844,0.406,0.312,0.112,0.444,0.288,0.388,0.05,0.396,0.45,0.101,0.231,0.213,0.473,0.29,0.402,0.97,0.787,0.089,0.16,0.444,0.355,0.396,0.266,0.077,0.42,0.107,0.55,0.592,0.467,0.769,0.556,0.308,0.538,0.527,0.722,0.669,0.083,0.065,0.935,0.018,0.095,0.219,0.107,0.822,0.243,0.308,0.331,0.343,0.148,0.32,0.024,0.621,0.112,0.077,0.237,0.053,0.231,0.828,0.68,0.036,0.101,0.059,0.467,0.148,0.592,0.118,0.444,0.136,0.367,0.675,0.107,0.071,0.55,0.254,0.172,0.018,0.213,0.592,0.089,0.26,0.284,0.118,0.095,0.059,0.68,0.355,0.059,0.118,0.13,0.225,0.124,0.296,0.083,0.041,0.03,0.012,0.047,0.148,0.095,0.237,0.053,0.077,0.03,0.041,0.118,0.195,0.142,0.225,0.166,0.183,0.041,0.059,0.272,0.266,0.047,0.083,0.071,0,0.166,0.03,0.142,0.006,0.243,0.112,0.006,0.012,0.036,0.024,0.018,0.03,0.089,0.03,0.278,0.036,0.053,0.172,0.195,0.107,0.142,0.112,0.012,0.071,0.414,0.272,0.349,0.118,0.047,0.243,0.018,0.574,0,0.538,0.195,0.047,0.556,0.118,0.302,0.006,0.077,0.018,0.385,0.266,0.136,0.734,0.864,0.065,0.284,0.183,0,0,0.053,0.136,0.521,0.154,0.012,0.219,0.012,0.698,0.083,0.012,0.024,0.077,0.124,0.065,0.053,0.592,0.089,0.095,0.243,0.379,0.006,0.065,0.047,0.036,0.03,0.13,0.024,0.195,0.189,0.012,0.03,0.296,0.722,0.018,0.095,0.095,0.178,0.13,0,0.065,0.13,0.053,0.059,0.183,0.385,0.047,0.036,0.006,0.036,0.006,0.006,0.13,0.006,0.089,0.609,0.03,0.343,0.142,0.195,0.13,0.024,0.136,0.059,0.012,0.024,0.077,0.024,0,0,0.095,0.13,0.373,0.296,0.195,0.231,0.006,0.006,0.065,0.036,0.148,0.379,0.124,0.053,0.065,0.018,0.095,0.012,0.036,0.036,0.012,0.195,0.012,0.041,0.207,0.095,0.219,0.047,0.083,0,0.37,0.773,0.525,0.282,0.072,0.414,0.293,0.713,0.144,0.873,0.204,0.613,0.448,0.967,0.354,0.403,1,0.298,0.807,0.188,0.144,0.691,0.238,0.807,0.249,0.331,0.376,0.182,0.928,0.409,0.547,0.481,0.514,0.309,0.298,0.32,0.122,0.558,0.508,0.061,0.193,0.365,0.199,0.702,0.547,0.282,0.68,0.188,0.243,0.171,0.144,0.326,0.276,0.768,0.293,0.122,0.315,0.309,0.215,0.16,0.263,0.366,0.188,0.102,0.177,0.258,0.274,0.935,0.823,0.167,0.667,0.758,0.919,0.151,0.28,0.188,0.183,0.188,0.312,0.151,0.565,0.968,0.855,0.935,0.898,0.059,0.844,0.78,0.532,0.011,0.575,0.097,0.699,0.828,1,1,0.14,0.376,0.414,0.145,0.059,0.919,0.435,0.892,0.108,0.387,0.022,0.247,0.043,1,0.984,0.392,0.672,0.79,0.715,0.118,0.14,0.188,0.527,0.124,0.247,0.349,0.022,0.672,0.161,0.016,0.317,0.183,0.812,0.64,0.091,0.263,1,0.333,0.398,0.952,0.097,0.172,0.876,0.167,0.129,0.086,1,0.543,0.269,0.075,0,0.849,0.204,0.903,0.468,0.097,0.591,0.876,0.935,0.065,0.892,0.011,0.828,0.145,0.086,0.108,0.065,0.102,0.124,0.677,0.898,0.065,0.532,0.527,0.806,0.043,0.102,0.237,0.817,1,0.038,0.145,0.909,0.656,0.323,0.091,0.392,0.801,0.323,0.22,0.145,0.086,0.097,0.925,0.091,0.892,0.398,0.226,0.726,0.075,0.07,0.484,0.833,0.048,0.204,0.323,0.108,0.102,0.247,0.882,0.011,0.285,0.903,0.28,0.242,0.425,0.124,0.075,0.253]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>cluster<\/th>\n      <th>gene<\/th>\n      <th>avg_logFC<\/th>\n      <th>p_val_adj<\/th>\n      <th>p_val<\/th>\n      <th>pct.1<\/th>\n      <th>pct.2<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"scroller":true,"deferRender":true,"scrollY":350,"fixedColumns":true,"autoWidth":true,"pageLength":10,"dom":"Bfrtip","buttons":{"extend":"collection","buttons":["csv","excel","pdf"]},"columnDefs":[{"className":"dt-right","targets":[2,3,4,5,6]}],"order":[],"orderClasses":false}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->


### Cluster Dendrograms


```r
Idents(nsc) <- 'seurat_clusters' 
avg <- AverageExpression(object = nsc, assays = 'RNA', slot = 'scale.data')[[1]]
nsc.dist <- dist(x = t(avg))
nsc.tree <- hclust(d = nsc.dist, method = 'ward.D2')
nsc.dend <- dendsort::dendsort(as.dendrogram(nsc.tree, hang = 0.1))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
nsc.cols <- gg_color_hue(n = length(levels(nsc$seurat_clusters)))
dendextend::labels_colors(nsc.dend) <- nsc.cols[nsc.tree$labels][order.dendrogram(nsc.dend)]
dendextend::labels_colors(nsc.dend) <- nsc.cols
dendextend::labels_cex(nsc.dend) <- 0.75
tiff(filename = paste0(results_outpath, 'Cluster_dendrogram.tiff'), height = 630, width = 630)
{par(mar = c(2,2,2,2), cex = 3, lwd = 2); dendextend::plot_horiz.dendrogram(nsc.dend, side = FALSE, main = 'Cluster dendrogram')}
dev.off()
{par(mar = c(2,2,2,2), cex = 2, lwd = 2, cex.main = 1); dendextend::plot_horiz.dendrogram(nsc.dend, side = FALSE, main = 'Cluster dendrogram')}
```

<img src="20200625_subclustering_nsc_files/figure-html/dendrogram_subcluster-1.png" style="display: block; margin: auto;" />

```r
rm(avg, nsc.dist, nsc.tree, nsc.dend)
```

```
## png 
##   2
```



```r
counts <- data.frame(table(nsc$seurat_clusters, nsc$orig.ident))
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

<img src="20200625_subclustering_nsc_files/figure-html/nsc_population_dynamics-1.png" style="display: block; margin: auto;" />

```r
ggsave(filename = paste0(results_outpath, 'nsc_population_dynamics.tiff'), plot = pop.dynamics, device = 'tiff', height = 3.75, width = 8)
```




```r
tmp <- FetchData(object = nsc, vars = c('seurat_clusters', 'orig.ident', VariableFeatures(nsc)), slot = 'scale.data')
tmp_annotation <- tmp[c('seurat_clusters', 'orig.ident')]
tmp_dat <- tmp[sapply(tmp, is.numeric)]
pheatmap::pheatmap(tmp_dat)
nsc.dist <- dist(x = tmp_dat)
nsc.tree <- hclust(d = nsc.dist, method = 'ward.D2')
nsc.dend <- dendsort::dendsort(as.dendrogram(nsc.tree, hang = 0.1))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
tmp_order <- order.dendrogram(nsc.dend)
nsc.cols <- gg_color_hue(n = length(levels(nsc$seurat_clusters)))
tmp_clust <- plyr::mapvalues(x = rownames(tmp_annotation)[tmp_order],
                             from = rownames(tmp_annotation)[tmp_order],
                             to = as.character(tmp_annotation$seurat_clusters[tmp_order]))
dendextend::labels_colors(nsc.dend) <- nsc.cols[tmp_clust]
dendextend::labels_colors(nsc.dend) <- nsc.cols
dendextend::labels_cex(nsc.dend) <- 0.75
plot(nsc.tree)
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
