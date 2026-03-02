suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(Seurat)
  library(mclust)
})
pdf(NULL)

args <- commandArgs(trailingOnly=T)
dataset <- args[1]
sample <- args[2]
n_domains <- as.integer(args[3])

################################ read data ################################
print('Reading data')

input_folder = file.path('resources', dataset)
output_folder = file.path('results', dataset, 'Seurat', sample)
if(!dir.exists(output_folder)){
  dir.create(output_folder, recursive=T)
}
counts = data.table::fread(file.path(input_folder, paste(sample, 'UMImatrix.csv', sep='_')), header=F, showProgress=T)
coords = data.table::fread(file.path(input_folder, paste(sample, 'coords.csv', sep='_')), header=F)

data <- t(as.matrix(counts))
rownames(data) = 1:dim(data)[1]
colnames(data) = 1:dim(data)[2]

################################ run method ###############################
print('\nRunning method\n')

sce <- CreateSeuratObject(counts=data)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce)
sce <- FindNeighbors(sce)

lo <- -20.
hi <- 20.
while(hi-lo > 1e-3){
  mi <- (lo+hi)/2
  sce <- FindClusters(sce, resolution=mi)
  n_pred <- length(levels(Idents(sce)))
  if(n_pred==n_domains){
    lo <- mi
    break
  } else if(n_pred<n_domains){
    lo <- mi
  } else {
    hi <- mi
  }
}

############################### save output ###############################

labels = Idents(sce)
assertthat::assert_that(length(labels)==dim(coords)[1])
write.table(labels, file=file.path(output_folder, 'clustering.tsv'), 
            sep='\t', quote=F, row.names=F, col.names=c('region'))
print('Saved domain labels to clustering.tsv')
