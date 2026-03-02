suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(vesalius)
})
pdf(NULL)
set.seed(1)

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
output_folder <- file.path('results', dataset, 'Vesalius', sample)
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

colnames(coords) <- c('x','y')
coords[, 'barcodes'] <- 1:nrow(coords)
counts = t(counts)
colnames(counts) <- 1:ncol(counts)

ves <- build_vesalius_assay(coords, counts)
ves <- generate_embeddings(ves)
ves <- regularise_image(ves)
ves <- smooth_image(ves, sigma=5, iter=10)
ves <- equalize_image(ves, sleft=5, sright=5)

print('segmenting image and isolating territories')
for(n_seg in rep(n_domains:1, each=5)){
  ves <- segment_image(ves, 
                       col_resolution=n_seg,
                       verbose=F)
  ves <- isolate_territories(ves,
                             verbose=F)
  n_pred <- dim(unique(ves@territories[ncol(ves@territories)]))[1]
  if(n_pred==n_domains){
    print(paste('got correct n_domains with', n_seg, 'segments'))
    break
  }
}

if(n_pred!=n_domains){
  print('varying capture radius.')
  seg_cols <- grep('Segment', colnames(ves@territories))
  for(rad in rep(seq(from=0.05, to=0.5, by=0.05), each=5)){
    for(ii in 1:length(seg_cols)){
      segname <- colnames(ves@territories[seg_cols[ii]])
      ves <- isolate_territories(ves,
                                 trial=segname,
                                 capture_radius=rad,
                                 verbose=F)
      n_pred <- dim(unique(ves@territories[ncol(ves@territories)]))[1]
      if(n_pred==n_domains){ 
        print(paste('got correct n_domains with', 
                    dim(unique(ves@territories[seg_cols[ii]]))[1],
                    'segments and capture radius', rad))
        break 
      }
    }
    if(n_pred==n_domains){ break }
  }
}
assertthat::assert_that(n_pred==n_domains)

############################### save output ###############################

labels <- ves@territories[ncol(ves@territories)]
tryCatch({
  print('number of labels equals number of coordinates:')
  assertthat::assert_that(dim(labels)[1]==dim(coords)[1])
  write.table(labels, file=file.path(output_folder, 'clustering.tsv'), 
              sep='\t', quote=F, row.names=F, col.names=c('region'))
  print('Saved domain labels to clustering.tsv')
}, error = function(e){
  print('reinserting missing spots')
  missing_spots <- setdiff(row.names(coords), row.names(labels))
  labels_df <- as.data.frame(labels)
  labels_df[missing_spots,] <- NA
  labels <- labels_df[order(as.numeric(row.names(labels_df))),]
  print('number of labels after reinserting:')
  print(length(labels))
  assertthat::assert_that(length(labels)==dim(coords)[1])
  write.table(labels, file=file.path(output_folder, 'clustering.tsv'),
              sep='\t', quote=F, row.names=F, col.names=c('region'))
  print('Saved domain labels to clustering.tsv')
})