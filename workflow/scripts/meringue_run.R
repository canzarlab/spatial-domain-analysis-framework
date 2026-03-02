suppressPackageStartupMessages({
  library(MERINGUE)
})
pdf(NULL)

# this function from MERINGUE, added bisection for number of domains
my_getSpatiallyInformedClusters <- function(pcs, W, k, num, 
                                           alpha=1, beta=1,
                                           method = igraph::cluster_louvain,
                                           verbose=FALSE,
                                           details=FALSE) {
  
  if (verbose) {
    print("finding approximate nearest neighbors ...")
  }
  # nearest neighbors in PC space
  nn = RANN::nn2(pcs, k = k) ## KNN
  names(nn) <- c('idx', 'dists')
  
  if (verbose) {
    print("using spatial weights ...")
  }
  # create weights from binary adjacency matrix W
  nw.simple = igraph::graph_from_adjacency_matrix(W)
  nw.simple = igraph::simplify(nw.simple)
  pweight <- do.call(rbind, lapply(1:nrow(nn$idx), function(i) {
    d <- unlist(lapply(nn$idx[i,], function(j) {
      if(i==j) { return(0) }
      else {
        return(igraph::distances(nw.simple, rownames(W)[i], colnames(W)[j]))
        #return(igraph::min_cut(nw.simple, rownames(W)[i], colnames(W)[j]))
      }
    }))
    return(d)
  }))
  
  if (verbose) {
    print(paste0("calculating weights with alpha:", alpha, " and beta:", beta))
  }
  weight <- 1/(alpha + as.vector(pweight)) + beta
  
  if (verbose) {
    print("calculating clustering ...")
  }
  # graph-based clustering with spatial weights
  nn.df = data.frame(from = rep(1:nrow(nn$idx), k),
                     to = as.vector(nn$idx),
                     weight = weight
  )
  g <- igraph::graph_from_data_frame(nn.df, directed = FALSE)
  g <- igraph::simplify(g)
  
  ################### modified from original ###################
  lo <- 0.
  hi <- 20.
  while(hi-lo > 1e-3){
    mi <- (lo+hi)/2
    print(mi)
    km <- method(g, resolution=mi)
    n_pred <- length(unique(km$membership))
    print(paste0('number of clusters: ', n_pred))
    if(n_pred==num){
      lo <- mi
      break
    } else if(n_pred<num){
      lo <- mi
    } else {
      hi <- mi
    }
  }
  ######################## end modified ########################
  
  if (verbose) {
    mod <- igraph::modularity(km)
    if (mod < 0.3) {
      print("WARNING")
    }
    print(paste0("graph modularity: ", mod))
  }
  com <- km$membership
  names(com) <- rownames(pcs)
  com <- factor(com)
  if (verbose) {
    print("identifying cluster membership ...")
    print(table(com))
  }
  if (details) {
    return(list(com = com, mod = mod, g = g))
  }
  else {
    return(com)
  }
}


args <- commandArgs(trailingOnly=T)
dataset <- args[1]
sample <- args[2]
n_domains <- as.integer(args[3])

################################ read data ################################
print('Reading data')

input_folder <- file.path('resources', dataset)
output_folder <- file.path('results', dataset, 'MERINGUE', sample)
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

counts = Matrix::Matrix(t(counts), sparse=TRUE)
dimnames(counts) = list(paste0('gene', 1:nrow(counts)), paste0('cell', 1:ncol(counts)))
counts <- cleanCounts(counts=counts, min.reads=100, min.lib.size=100, plot=F)

counts <- normalizeCounts(counts=counts)

filterdist <- max(sapply(coords, max)-sapply(coords, min)) / 20
w <- getSpatialNeighbors(coords)

com2 <- my_getSpatiallyInformedClusters(pcs=t(as.matrix(counts)), 
                                        W=w, 
                                        k=min(50, dim(coords)[1]),
                                        num=n_domains,
                                        verbose=T)

############################### save output ###############################

labels <- com2
tryCatch({
  print('number of labels equals number of coordinates:')
  assertthat::assert_that(length(labels)==dim(coords)[1])
  write.table(labels, file=file.path(output_folder, 'clustering.tsv'), 
              sep='\t', quote=F, row.names=F, col.names=c('region'))
  print('Saved domain labels to clustering.tsv')
}, error = function(e){
  print('reinserting missing spots')
  missing_spots <- setdiff(paste0('cell', 1:dim(coords)[1]), names(labels))
  labels_df <- as.data.frame(labels)
  labels_df[missing_spots,] <- NA
  row.names(labels_df) <- sub('^cell', '', row.names(labels_df))
  labels <- labels_df[order(as.numeric(row.names(labels_df))),]
  print('number of labels after reinserting:')
  print(length(labels))
  assertthat::assert_that(dim(labels_df)[1]==dim(coords)[1])
  write.table(labels, file=file.path(output_folder, 'clustering.tsv'),
              sep='\t', quote=F, row.names=F, col.names=c('region'))
  print('Saved domain labels to clustering.tsv')
})
