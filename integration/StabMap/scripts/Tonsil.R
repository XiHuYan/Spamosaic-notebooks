library(Seurat)
# library(ArchR)
library(StabMap)
library(glue)
library(scran)
library(Matrix)
library(Signac)

set.seed(2021)

pj <- file.path

load_data <- function(path, prefix){
  data = readMM(pj(path, glue('{prefix}_data.mtx')))
  cell_name = read.table(pj(path, glue('{prefix}_cell_names.csv')), sep=',', header=T, row.names=1)
  feat_name = read.table(pj(path, glue('{prefix}_feat_names.csv')), sep=',', header=T, row.names=1)
  rownames(data) = cell_name$X0
  colnames(data) = feat_name$X0

  return (t(data))
}

out_dir = '/disco_500t/xuhua/gitrepo/BridgeNorm/stabmap/outputs/'

a1_dir = '/disco_500t/xuhua/data/spatial_multi_omics/lymp_tonsil_ramen/tonsil_A1/R_data'
d1_dir = '/disco_500t/xuhua/data/spatial_multi_omics/lymp_tonsil_ramen/tonsil_D1/R_data'

mult_rna_count = load_data(a1_dir, 'rna')
mult_adt_count = load_data(a1_dir, 'adt')
single_rna_count = load_data(d1_dir, 'rna')
single_adt_count = load_data(d1_dir, 'adt')

n_genes = dim(mult_rna_count)[1]
n_peaks = dim(mult_adt_count)[1]

head(colnames(single_rna_count))
head(colnames(single_adt_count))

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
  file=paste0(out_dir, 'tonsil.csv'),
  quote=F, row.names=T)

