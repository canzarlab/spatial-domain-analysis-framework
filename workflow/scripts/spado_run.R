suppressPackageStartupMessages({
  library(SpaDo)
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

config_file <- read.csv('config/technology_meta.csv')

if (tech == 'MERFISH') {
  tech <- 'merFISH'
} else {
  tech <- tech
}

################################ read data ################################
print('Reading data')

input_folder <- file.path('resources', dataset)
output_folder <- file.path('results', dataset, 'SpaDo', sample)
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

expr_norm <- SpatialNormalize(expression_profile=counts,
                              ST_method=tech)
initial_result <- InitialClustering(expression_profile=expr_norm,
                                    user_offered=F)
ct_distribution <- SpatialCellTypeDistribution(sample_information_coordinate=coords,
                                               # sequence_resolution=seq_res,
                                               # with sequence_resolution=='spot', need ct deconvolution
                                               sample_information_cellType=initial_result$sample_information)
if (dim(coords)[1] < 30000){
  DD <- DistributionDistance(ct_distribution,
                             distance="jensen-shannon")
} else {
  DD <- DistributionDistance(ct_distribution,
                             distance="manhattan")
}
domain_hclust <- DomainHclust(distribution_distance=DD,
                              autoselection=F,
                              domain_num=n_domains)

############################### save output ###############################

labels <- domain_hclust$hclust_result_df[paste0('Domain_level_', n_domains)]
labels = labels[order(as.numeric(sub('^cell', '', row.names(labels)))), ]
assertthat::assert_that(length(labels)==dim(coords)[1])
write.table(labels, file=file.path(output_folder, 'clustering.tsv'), 
            sep='\t', quote=F, row.names=F, col.names=c('region'))
print('Saved domain labels to clustering.tsv')
