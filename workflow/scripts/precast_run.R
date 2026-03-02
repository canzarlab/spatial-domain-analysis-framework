suppressPackageStartupMessages({
  library(Seurat)
  library(PRECAST)
  library(scater)
  library(patchwork)
})
pdf(NULL)

args <- commandArgs(trailingOnly=T)
dataset <- args[1]
sample <- args[2]
n_domains <- as.integer(args[3])
if (length(args)==4) {
  tech <- args[4]
} else {
  tech <- stringr::str_split(dataset, '_')[[1]][1]
}
if (tech == 'ST' | tech == 'sciSpace') {
  platform = 'ST'
} else {
  platform = 'Visium'
}

################################ read data ################################
print('Reading data')

input_folder <- file.path('resources', dataset)
output_folder <- file.path('results', dataset, 'PRECAST', sample)
if(!dir.exists(output_folder)){
  dir.create(output_folder, recursive=T)
}
counts <- data.table::fread(file.path(input_folder, 
                                      paste(sample, 'UMImatrix.csv', sep='_')), 
                            header=F, showProgress=T)
coords <- data.table::fread(file.path(input_folder, 
                                      paste(sample, 'coords.csv', sep='_')), 
                            header=F)

################################ run method ###############################
print('Running method')

sobj <- CreateSeuratObject(counts= t(counts), assay='RNA')
colnames(coords) <- c('x','y')
sobj@meta.data$row <- coords[,1]
sobj@meta.data$col <- coords[,2]

preobj<- CreatePRECASTObject(seuList = list(sobj))

PRECASTObj<- AddAdjList(preobj, platform = platform)

PRECASTObj<- AddParSetting(PRECASTObj)

PRECASTObj <- PRECAST(PRECASTObj, K = n_domains)

resList<- PRECASTObj@resList
PRECASTObj <- SelectModel(PRECASTObj)

############################### save output ###############################

labels <- PRECASTObj@resList$cluster[[1]]
tryCatch({
  print('number of labels equals number of coordinates:')
  assertthat::assert_that(length(labels)==dim(coords)[1])
  write.table(labels, file=file.path(output_folder, 'clustering.tsv'), 
              sep='\t', quote=F, row.names=F, col.names=c('region'))
  print('Saved domain labels to clustering.tsv')
}, error = function(e){
  print('reinserting missing spots')
  row.names(labels) <- row.names(PRECASTObj@seulist[[1]]@meta.data)
  missing_spots <- setdiff(paste0('Cell_', 1:dim(coords)[1]), row.names(labels))
  labels_df <- as.data.frame(labels)
  labels_df[missing_spots,] <- NA
  row.names(labels_df) <- sub('^Cell_', '', row.names(labels_df))
  labels <- labels_df[order(as.numeric(row.names(labels_df))),]
  print('number of labels after reinserting:')
  print(length(labels))
  assertthat::assert_that(dim(labels_df)[1]==dim(coords)[1])
  write.table(labels, file=file.path(output_folder, 'clustering.tsv'),
              sep='\t', quote=F, row.names=F, col.names=c('region'))
  print('Saved domain labels to clustering.tsv')
})
