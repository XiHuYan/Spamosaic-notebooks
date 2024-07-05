library(Seurat)
# library(ArchR)
library(StabMap)
library(glue)
library(scran)
library(Matrix)
library(Signac)

set.seed(2021)
pj <- file.path

read_data <- function(dir, name){
    mat = readMM(pj(dir, glue('{name}_mat.mtx')))
    cname = read.table(pj(dir, glue('{name}_cname.csv')), sep=',', header=T, row.names=1)
    feat_name = read.table(pj(dir, glue('{name}_feat.csv')), sep=',', header=T)
    rownames(mat) = cname[, 1]
    colnames(mat) = feat_name[, 2]

    return (t(mat))
}

out_dir = '/disco_500t/xuhua/gitrepo/BridgeNorm/stabmap/outputs/'

data_dir = '/disco_500t/xuhua/data/spatial_multi_omics/lymp_tonsil_ramen/R_data_modalMatch'
name = 'bridge_adt'
cname = read.table(pj(data_dir, glue('{name}_cname.csv')), sep=',', header=T, row.names=1)
feat_name = read.table(pj(data_dir, glue('{name}_feat.csv')), sep=',', header=T)

mult_rna_count = read_data(data_dir, 'bridge_rna')
mult_adt_count = read_data(data_dir, 'bridge_adt')
single_rna_count = read_data(data_dir, 'test_rna')
single_adt_count = read_data(data_dir, 'test_adt')

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
# colnames(single_rna_lognorm) = paste0('rna_', colnames(single_rna_lognorm))
# colnames(single_adt_lognorm) = paste0('adt_', colnames(single_adt_lognorm))
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
  file=paste0(out_dir, 'lymph_modalMatch.csv'),
  quote=F, row.names=T)

