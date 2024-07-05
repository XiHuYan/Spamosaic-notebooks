library(Seurat)
library(ArchR)
library(StabMap)
library(glue)
library(scran)
library(Matrix)
library(Signac)

set.seed(2021)

data_dir = '/home/cb213/local/cache3/yxh/Data/spatial_dataset/MISAR_seq/section1/'
setwd(data_dir)
rds1 = readRDS(paste0(data_dir, 'ArchR_Section1.rds'))
out_dir = '/home/cb213/local/cache3/yxh/Data/exp_cache/stabmap/misar_seq/'

# loading 
meta_data = (rds1@cellColData)

pkm = getMatrixFromProject(
  ArchRProj = rds1,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  # logFile = createLogFile("getMatrixFromProject")
)

# write GEM matrix
gem = getMatrixFromProject(
  ArchRProj = rds1,
  useMatrix = "GeneExpressionMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  # logFile = createLogFile("getMatrixFromProject")
)

# E15
E15_names <- which(sapply(colnames(pkm), function(x) grepl("^E15", x)))
E13_names <- which(sapply(colnames(pkm), function(x) grepl("^E13", x)))
E18_names <- which(sapply(colnames(pkm), function(x) grepl("^E18", x)))

mult_rna_count = assay(gem[, E15_names])
mult_atac_count = assay(pkm[, E15_names])
single_rna_count = assay(gem[, E13_names])
single_atac_count = assay(pkm[, E18_names])
n_genes = dim(mult_rna_count)[1]
n_peaks = dim(mult_atac_count)[1]
rownames(mult_rna_count) = paste0('gene', c(1:n_genes))
rownames(mult_atac_count) = paste0('peak', c(1:n_peaks))
rownames(single_rna_count) = paste0('gene', c(1:n_genes))
rownames(single_atac_count) = paste0('peak', c(1:n_peaks))

# select features
data_rna <- CreateSeuratObject(counts = cbind(mult_rna_count, single_rna_count), project = "rna")
data_atac <- CreateSeuratObject(counts = cbind(mult_atac_count, single_atac_count), project = "atac")

data_rna <- NormalizeData(data_rna) %>% FindVariableFeatures()
variable_rna_features = VariableFeatures(data_rna)

data_atac <- RunTFIDF(data_atac)
data_atac <- FindTopFeatures(data_atac, min.cutoff = 'q90')
variable_atac_features = VariableFeatures(data_atac)
length(variable_atac_features)

n_mult = dim(mult_rna_count)[2]
n_single_rna = dim(single_rna_count)[2]
n_single_atac = dim(single_atac_count)[2]
mult_lognorm = rbind(
  data_rna[['RNA']]@data[variable_rna_features, 1:n_mult], 
  data_atac[['RNA']]@data[variable_atac_features, 1:n_mult]  # (n_mult+1):(n_mult+n_single)
)  
dim(mult_lognorm)

# stabmap !!!!!!!!!!!!!!
single_rna_lognorm = as(data_rna[['RNA']]@data[variable_rna_features, (n_mult+1):(n_mult+n_single_rna)], 'dgCMatrix')
single_atac_lognorm = as(data_atac[['RNA']]@data[variable_atac_features, (n_mult+1):(n_mult+n_single_atac)], 'dgCMatrix') 
colnames(single_rna_lognorm) = paste0('rna_', colnames(single_rna_lognorm))
colnames(single_atac_lognorm) = paste0('atac_', colnames(single_atac_lognorm))
assay_list_indirect = list(
  RNA = single_rna_lognorm,    # dgCMatrix
  Multiome = as(mult_lognorm, 'dgCMatrix'), # dgCMatrix
  ATAC = single_atac_lognorm  # dgCMatrix
)

lapply(assay_list_indirect, dim)

jpeg(file=paste0(out_dir, "E15-13-18-overlapFeature.jpg"),  width = 4, height = 4, units = 'in', res = 300)
mosaicDataUpSet(assay_list_indirect, plot = T)
dev.off()

# Running !!!!!!!!!!
# return matrix array
stab_indirect = stabMap(assay_list_indirect,
                        reference_list = c("Multiome"),
                        maxFeatures=30000,
                        plot = FALSE)
write.csv(as.data.frame(stab_indirect), 
  file=paste0(out_dir, 'embed_Npeak=20000.csv'),
  quote=F, row.names=T)