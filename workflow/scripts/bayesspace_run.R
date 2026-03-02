suppressPackageStartupMessages({
  library(BayesSpace)
  library(Matrix)
})
pdf(NULL)
set.seed(104)

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

input_folder = file.path('resources', dataset)
output_folder = file.path('results', dataset, 'BayesSpace', sample)
if(!dir.exists(output_folder)){
  dir.create(output_folder, recursive=T)
}
counts <- data.table::fread(file.path(input_folder, 
                                      paste(sample, 'UMImatrix.csv', sep='_')), 
                            header=F, showProgress=T)
coordsfile <- file.path(input_folder, 
                        paste(sample, 'index_coords.csv', sep='_'))
if(file.exists(coordsfile)){
  coords <- data.table::fread(coordsfile, header=F)
} else {
  print('No index_coords available, processing data')
  orig <- data.table::fread(file.path(input_folder,
                                      paste(sample, 'coords.csv', sep='_')),
                            header=F)
  bins_index <- data.table::fread(file.path(input_folder, paste(sample, 'nbins.csv', sep='_')))
  data.table::setkey(bins_index, sample)
  nbins_row <- bins_index[sample, 'row']
  nbins_col <- bins_index[sample, 'col']
  indexrow <- cut(sapply(orig[, 1], as.numeric), 
                  nbins_row, 
                  include.lowest=T, 
                  labels=F)
  indexcol <- cut(sapply(orig[, 2], as.numeric), 
                  nbins_col, 
                  include.lowest=T, 
                  labels=F)
  coords <- data.table::data.table(indexrow, indexcol)
  grouped_indices <- coords[, .(indices = list(.I)), by = .(indexrow, indexcol)]
  grouped_indices <- grouped_indices[sapply(indices, length) > 1]
  to_drop <- meanrows <- rowindices <- c()
  rownames_orig <- row.names(orig)
  if(length(grouped_indices$indices)>0){
    pb = txtProgressBar(max=length(grouped_indices$indices), style=3)
    ii=0
    for(inds in grouped_indices$indices){
      meanrows <- append(meanrows, as.list(colMeans(counts[inds,])))
      rowindices <- append(rowindices, inds[1])
      to_drop <- append(to_drop, inds[2:length(inds)])
      ii = ii+1
      setTxtProgressBar(pb, ii)
    }
    counts[rowindices,] <- data.table::data.table(matrix(meanrows, ncol=dim(counts)[2], byrow=TRUE))
    counts <- counts[-to_drop,]
    coords <- coords[-to_drop,]
    rownames_drop <- rownames_orig[-to_drop]
    close(pb)
  }
  print(paste('number of spots reduced from', dim(orig)[1], 
              'to', dim(coords)[1]))
}

################################ run method ###############################
print('Running method')

names(coords) <- c('row', 'col')
dlpfc <- SingleCellExperiment(assays=list(counts=as(t(counts), 'dgCMatrix')),
                              colData=coords)

dlpfc <- spatialPreprocess(dlpfc, platform=platform)

dlpfc <- spatialCluster(dlpfc, 
                        q=n_domains, 
                        platform=platform)

############################### save output ###############################

labels <- colData(dlpfc)$spatial.cluster
if(!file.exists(coordsfile)){
  if(length(grouped_indices$indices)>0){
    labels_df <- as.data.frame(labels)
    row.names(labels_df) <- rownames_drop
    for(inds in grouped_indices$indices){
      labels_df[as.character(inds[2:length(inds)]),] <- labels_df[as.character(inds[1]),]
    }
    labels <- labels_df[order(as.numeric(row.names(labels_df))),]
  }
  assertthat::assert_that(length(labels)==dim(orig)[1])
} else {
  assertthat::assert_that(length(labels)==dim(coords)[1])
}
write.table(labels, file=file.path(output_folder, 'clustering.tsv'), 
            sep='\t', quote=F, row.names=F, col.names=c('region'))
print('Saved domain labels to clustering.tsv')
