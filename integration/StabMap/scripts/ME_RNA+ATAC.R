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
    colnames(mat) = feat_name[, 1]

    return (t(mat))
}

out_dir = '/disco_500t/xuhua/gitrepo/BridgeNorm/stabmap/outputs/'
data_dir = '/disco_500t/xuhua/data/real_mosaic_cases/mouse_embryo_rna+atac/R_data/'

# name = 'bridge_rna'
# cname = read.table(pj(data_dir, glue('Tsample_{name}_cname.csv')), sep=',', header=T, row.names=1)
# feat_name = read.table(pj(data_dir, glue('Tsample_{name}_feat.csv')), sep=',', header=T, row.names=1)
# head(cname)
# head(feat_name)

mult_rna_count = read_data(data_dir, 'bridge_rna')
mult_atac_count = read_data(data_dir, 'bridge_atac')
single_atac_count1 = read_data(data_dir, 'atac11')
single_atac_count2 = read_data(data_dir, 'atac13')
single_atac_count3 = read_data(data_dir, 'atac13_2')

dim(mult_rna_count)
dim(mult_atac_count)
dim(single_atac_count1)
dim(single_atac_count2)
dim(single_atac_count3)
head(colnames(mult_rna_count))
head(colnames(single_atac_count1))

n_genes = dim(mult_rna_count)[1]
n_peaks = dim(mult_atac_count)[1]

# select features
data_rna <- CreateSeuratObject(counts = mult_rna_count, project = "rna")
data_atac <- CreateSeuratObject(counts = cbind(mult_atac_count, single_atac_count1, single_atac_count2, single_atac_count3), project = "atac")

data_rna <- NormalizeData(data_rna) %>% FindVariableFeatures()
variable_rna_features = VariableFeatures(data_rna)

data_atac <- RunTFIDF(data_atac)
data_atac <- FindTopFeatures(data_atac, min.cutoff = 'q50')
variable_atac_features = VariableFeatures(data_atac)
length(variable_atac_features)

n_mult = dim(mult_rna_count)[2]
n_single_atac1 = dim(single_atac_count1)[2]
n_single_atac2 = dim(single_atac_count2)[2]
n_single_atac3 = dim(single_atac_count3)[2]

rna_data = GetAssayData(data_rna, assay = "RNA", layer = "data")
atac_data = GetAssayData(data_atac, assay = "RNA", layer = "data")

mult_lognorm = rbind(
  rna_data[variable_rna_features, 1:n_mult], 
  atac_data[variable_atac_features, 1:n_mult]  # (n_mult+1):(n_mult+n_single)
)  
dim(mult_lognorm)

# stabmap !!!!!!!!!!!!!!
single_atac_lognorm1 = as(atac_data[variable_atac_features, (n_mult+1):(n_mult+n_single_atac1)], 'dgCMatrix') 
single_atac_lognorm2 = as(atac_data[variable_atac_features, (n_mult+n_single_atac1+1):(n_mult+n_single_atac1+n_single_atac2)], 'dgCMatrix') 
single_atac_lognorm3 = as(atac_data[variable_atac_features, (n_mult+n_single_atac1+n_single_atac2+1):(n_mult+n_single_atac1+n_single_atac2+n_single_atac3)], 'dgCMatrix') 
# colnames(single_rna_lognorm) = paste0('rna_', colnames(single_rna_lognorm))
# colnames(single_atac_lognorm) = paste0('atac_', colnames(single_atac_lognorm))
assay_list_indirect = list(
  Multiome = as(mult_lognorm, 'dgCMatrix'), # dgCMatrix
  ATAC1 = single_atac_lognorm1,  # dgCMatrix,
  ATAC2 = single_atac_lognorm2,
  ATAC3 = single_atac_lognorm3
)

# lapply(assay_list_indirect, dim)

# jpeg(file=paste0(out_dir, "E15-18-18-overlapFeature.jpg"),  width = 4, height = 4, units = 'in', res = 300)
# mosaicDataUpSet(assay_list_indirect, plot = T)
# dev.off()

# Running !!!!!!!!!!
# return matrix array
stab_indirect = stabMap(assay_list_indirect,
                        reference_list = c("Multiome"),
                        maxFeatures=80000,
                        plot = FALSE)
write.csv(as.data.frame(stab_indirect), 
  file=paste0(out_dir, 'ME_RNA+ATAC.csv'),
  quote=F, row.names=T)

