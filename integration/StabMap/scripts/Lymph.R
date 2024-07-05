library(Seurat)
# library(ArchR)
library(StabMap)
library(glue)
library(scran)
library(Matrix)
library(Signac)

set.seed(2021)


out_dir = '/disco_500t/xuhua/gitrepo/BridgeNorm/stabmap/outputs/'

data_dir = '/disco_500t/xuhua/data/spatial_multi_omics/lymp_node/'
ln1_dir = paste0(data_dir, 'A1_LN-Jinmiao-sl6/outs/filtered_feature_bc_matrix/')
ln2_dir = paste0(data_dir, 'D1_LN-Jinmiao-sl5/outs/filtered_feature_bc_matrix/')

featn1 = read.table(gzfile(paste0(ln1_dir, "features.tsv.gz")), sep = "\t", header = F)
barcd1 = read.table(gzfile(paste0(ln1_dir, "barcodes.tsv.gz")), sep = "\t", header = F)
mat1 = readMM(gzfile(paste0(ln1_dir, 'matrix.mtx.gz')))
rownames(mat1) = featn1$V1
colnames(mat1) = barcd1$V1

featn2 = read.table(gzfile(paste0(ln2_dir, "features.tsv.gz")), sep = "\t", header = F)
barcd2 = read.table(gzfile(paste0(ln2_dir, "barcodes.tsv.gz")), sep = "\t", header = F)
mat2 = readMM(gzfile(paste0(ln2_dir, 'matrix.mtx.gz')))
rownames(mat2) = featn2$V1
colnames(mat2) = barcd2$V1


rna_mask <- which(sapply(featn1$V3, function(x) grepl("^Gene", x)))
adt_mask <- which(sapply(featn1$V3, function(x) grepl("^Antibody", x)))

mult_rna_count = mat1[rna_mask,]
mult_adt_count = mat1[adt_mask,]
single_rna_count = mat2[rna_mask,]
single_adt_count = mat2[adt_mask,]
colnames(single_rna_count) = paste0('rna_', colnames(single_rna_count))
colnames(single_adt_count) = paste0('adt_', colnames(single_adt_count))
n_genes = dim(mult_rna_count)[1]
n_peaks = dim(mult_adt_count)[1]

# select features
data_rna <- CreateSeuratObject(counts = cbind(mult_rna_count, single_rna_count), project = "rna")
data_adt <- CreateSeuratObject(counts = cbind(mult_adt_count, single_adt_count), project = "adt")

data_rna <- NormalizeData(data_rna) %>% FindVariableFeatures()
variable_rna_features = VariableFeatures(data_rna)

data_adt <- NormalizeData(data_adt, normalization.method = 'CLR', margin = 2)
variable_adt_features = rownames(data_adt)

n_mult = dim(mult_rna_count)[2]
n_single_rna = dim(single_rna_count)[2]
n_single_adt = dim(single_adt_count)[2]

rna_data = GetAssayData(data_rna, assay = "RNA", layer = "data")
adt_data = GetAssayData(data_adt, assay = "RNA", layer = "data")

mult_lognorm = rbind(
  rna_data[variable_rna_features, 1:n_mult], 
  adt_data[variable_adt_features, 1:n_mult]  # (n_mult+1):(n_mult+n_single)
)  
dim(mult_lognorm)

# stabmap !!!!!!!!!!!!!!
single_rna_lognorm = as(rna_data[variable_rna_features, (n_mult+1):(n_mult+n_single_rna)], 'dgCMatrix')
single_adt_lognorm = as(adt_data[variable_adt_features, (n_mult+1):(n_mult+n_single_adt)], 'dgCMatrix') 
colnames(single_rna_lognorm) = paste0('rna_', colnames(single_rna_lognorm))
colnames(single_adt_lognorm) = paste0('adt_', colnames(single_adt_lognorm))
assay_list_indirect = list(
  RNA = single_rna_lognorm,    # dgCMatrix
  Multiome = as(mult_lognorm, 'dgCMatrix'), # dgCMatrix
  ADT = single_adt_lognorm  # dgCMatrix
)

# lapply(assay_list_indirect, dim)

# jpeg(file=paste0(out_dir, "E15-18-18-overlapFeature.jpg"),  width = 4, height = 4, units = 'in', res = 300)
# mosaicDataUpSet(assay_list_indirect, plot = T)
# dev.off()

# Running !!!!!!!!!!
# return matrix array
stab_indirect = stabMap(assay_list_indirect,
                        reference_list = c("Multiome"),
                        maxFeatures=30000,
                        plot = FALSE)
write.csv(as.data.frame(stab_indirect), 
  file=paste0(out_dir, 'lymph.csv'),
  quote=F, row.names=T)

