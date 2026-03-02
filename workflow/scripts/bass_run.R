suppressPackageStartupMessages({
  library(Matrix)
  library(BASS)
})
pdf(NULL)
set.seed(0)

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[1]
sample <- args[2]
n_domains <- as.integer(args[3])

################################ read data ################################
print('Reading data')

input_folder = file.path('resources', dataset)
output_folder = file.path('results', dataset, 'BASS', sample)
if(!dir.exists(output_folder)){
  dir.create(output_folder, recursive=TRUE)
}
counts = data.table::fread(file.path(input_folder, 
                                     paste(sample, 'UMImatrix.csv', sep='_')), 
                           header=F, showProgress=T)
coords = data.table::fread(file.path(input_folder, 
                                     paste(sample, 'coords.csv', sep='_')), 
                           header=F)

################################ run method ###############################
print('Running method')

BASS <- createBASSObject(list(t(counts)), 
                         list(coords), 
                         C = 20, 
                         R = n_domains)

BASS <- BASS.preprocess(BASS)

BASS <- BASS.run(BASS)

BASS <- BASS.postprocess(BASS)

############################### save output ###############################

labels = BASS@results$z[[1]]
assertthat::assert_that(length(labels)==dim(coords)[1])
write.table(labels, file=file.path(output_folder, 'clustering.tsv'), 
            sep='\t', quote=F, row.names=F, col.names=c('region'))
print('Saved domain labels to clustering.tsv')
