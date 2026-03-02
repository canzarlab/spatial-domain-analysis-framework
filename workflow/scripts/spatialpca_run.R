suppressPackageStartupMessages({
  library(SpatialPCA)
  library(bluster)
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
  shape <- 'square'
} else {
  shape <- 'hexagon'
}

################################ read data ################################
print('Reading data')

input_folder <- file.path('resources', dataset)
output_folder <- file.path('results', dataset, 'SpatialPCA', sample)
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

counts <- t(counts)
colnames(counts) <- paste0('cell', 1:ncol(counts))
coords <- as.matrix(coords)
row.names(coords) <- colnames(counts)

LIBD <- CreateSpatialPCAObject(counts=counts, 
                              location=coords)

LIBD <- SpatialPCA_buildKernel(LIBD)

LIBD <- SpatialPCA_EstimateLoading(LIBD) 

LIBD <- SpatialPCA_SpatialPCs(LIBD)

clusterlabel <- walktrap_clustering(clusternum=n_domains,
                                  latent_dat=LIBD@SpatialPCs) 

clusterlabel_refine <- refine_cluster_10x(clusterlabels=clusterlabel,
                                          location=LIBD@location,
                                          shape=shape)

############################### save output ###############################

labels <- as.integer(clusterlabel_refine)
tryCatch({
  print('number of labels equals number of coordinates:')
  assertthat::assert_that(length(labels)==dim(coords)[1])
  write.table(labels, file=file.path(output_folder, 'clustering.tsv'), 
              sep='\t', quote=F, row.names=F, col.names=c('region'))
  print('Saved domain labels to clustering.tsv')
}, error = function(e){
  print('reinserting missing cells')
  labels_df <- as.data.frame(labels)
  row.names(labels_df) <- colnames(LIBD@SpatialPCs)
  missing_cells <- setdiff(colnames(counts), row.names(labels_df))
  labels_df[missing_cells,] <- NA
  row.names(labels_df) <- sub('^cell', '', row.names(labels_df))
  labels <- labels_df[order(as.numeric(row.names(labels_df))),]
  print('number of labels after reinserting:')
  print(length(labels))
  assertthat::assert_that(length(labels)==dim(coords)[1])
  write.table(labels, file=file.path(output_folder, 'clustering.tsv'),
              sep='\t', quote=F, row.names=F, col.names=c('region'))
  print('Saved domain labels to clustering.tsv')
})