

######### Data preparation for GEO submission ########
######### 
######### 

# Format and write-out data -----------------------------------------------

library('Matrix')
library('Seurat')
library('dplyr')

dir.create(path = 'data/GEO_submission')

svz = readRDS(file = 'data/svz.rds')
svz_mat = Matrix::Matrix(data = svz@assays$RNA@counts, sparse = FALSE)
svz_mat = Matrix(data = svz_mat)
genes_tsv = rownames(svz_mat)
barcodes_tsv = colnames(svz_mat)
metadata_tsv = svz@meta.data
metadata_tsv$barcode = rownames(svz@meta.data)
metadata_tsv = metadata_tsv[c(ncol(metadata_tsv), 1:(ncol(metadata_tsv)-1))]

writeMM(svz_mat, file = 'data/GEO_submission/svz_mat.mtx')
write(x = genes_tsv, file = 'data/GEO_submission/genes.tsv')
write(x = barcodes_tsv, file = 'data/GEO_submission/barcodes.tsv')
write.table(x = metadata_tsv, file = 'data/GEO_submission/metadata.tsv', sep = '\t', row.names = FALSE)



# Code chunk for verifying successful write -------------------------------

test = Matrix::readMM(file = 'data/GEO_submission/svz_mat.mtx')
barcodes = read.table(file = 'data/GEO_submission/barcodes.tsv', sep = '\t')
genes = read.table(file = 'data/GEO_submission/genes.tsv', sep = '\t')
metadata = read.table(file = 'data/GEO_submission/metadata.tsv', sep = '\t')

colnames(test) = barcodes$V1
rownames(test) = genes$V1

test = CreateSeuratObject(counts = test)
test = test %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:10) %>% 
  FindNeighbors(dims = 1:10) %>% 
  FindClusters()
FeaturePlot(test, features = c('Ascl1', 'Cx3cr1'), order = TRUE)
