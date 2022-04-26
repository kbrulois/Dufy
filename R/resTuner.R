#' Generate a specific number of clusters by tuning the resolution parameter of graph-based clustering algorithms
#' 
#' Existing functions for graph-based clustering allow users to obtain a larger or smaller number of clusters by varying parameters, in particular the resolution parameter. This function allows users to directly control the number of clusters by iteratively updating the resolution parameter such that a certain number of clusters are generated. The cluster_number parameter can be an exact number of clusters or a range. Clustering can be repeated to obtain variant clusterings with the desired number of clusters.
#' 
#' @param input_data A matrix where features correspond to rows and cells correspond to columns.
#' @param cluster_number Either a positive integer specifying the exact number of clusters desired or a vector of length two specifying range.
#' @param cluster_repeats The number of rounds of clusterings to perform. Default: 5.
#' @param seed_resolution The initial value of resolution.
#' @param max_iter The maximum number of iterations per clustering. Default: 10.
#' @param label_prefix A character string to be prepended to the colnames of the returned matrix.
#' @param cluster_function A resTuner-compatible function that performs nearest-neighbor graph generation and graph-based clustering. Must contain these three named arguments: input_data, resolution, and k. Must return a vector of length \code{ncol(input_data)} containing cluster assignments. Default: seurat_cluster_wrapper
#' @param ... Additional arguments passed to cluster_function
#' @return A \code{ncol(input_data)} by \code{cluster_repeats} \code{matrix} containing cluster assignments.
#'
#' @author Kevin Brulois
#' @export

resTuner <- function(input_data,
                     cluster_number,
                     cluster_repeats = 5, 
                     seed_resolution = 0.005, 
                     max_iter = 10,
                     cluster_function = Dufy::seurat_cluster_wrapper,
                     label_prefix = "seurat_",
                     ...) {
  
  assertthat::assert_that(sum(names(formals(cluster_function)) %in% c("input_data", ".resolution", "k")) == 3)
  
  x <- 1
  y <- 1
  resolution <- seed_resolution
  if(length(cluster_number) == 1) {cluster_number <- rep(cluster_number, 2)}
  clusts <- list()

  while(x <= cluster_repeats & y <= max_iter) {
    
    message("\nclustering round ", x, " iteration ", y, "\ntrying resolution ", round(resolution, 5))
    
    clusters <- cluster_function(input_data = input_data, .resolution = resolution, ...)
    
    clust_num <- length(unique(as.character(clusters)))
    
    message("found ", clust_num, " clusters")
    
    if(clust_num < cluster_number[1]) {
      message("too few")
      resolution <- resolution * cluster_number[2] / clust_num
      y <- y + 1
    }
    
    if(clust_num > cluster_number[2]) {
      message("too many")
      resolution <- resolution * cluster_number[1] / clust_num
      y <- y + 1
    }
    
    if(clust_num >= cluster_number[1] & clust_num <= cluster_number[2]) {
      message("within range")
      clusts[[paste0(label_prefix, "res", round(resolution, 5), "_n", clust_num)]] <- clusters
      x <- x + 1
      y <- 1
      resolution <- resolution * runif(1, 0.75, 1.5)
      
    }
  }
  if(length(clusts) > 1) {
  final_result <- do.call(cbind, clusts)
  ord <- order(sapply(1:ncol(final_result), function(x) length(unique(final_result[,x]))))
  final_result <- final_result[,ord]
  } else {
    final_result <- clusts[[1]]
  }
  final_result
}
    




#' A resTuner-compatible wrapper function for leiden clustering using monocle3
#' 
#' @param input_data A matrix where features correspond to rows and cells correspond to columns.
#' @param resolution Either a positive integer specifying the exact number of clusters desired or a vector of length two specifying range.
#' @param k The number of rounds of clusterings to perform. Default: 5.
#' @param ... Arguments passed to monocle3:::leiden_clustering
#' @return A vector of length \code{ncol(input_data)} containing cluster assignments.
#'
#' @author Kevin Brulois
#' @export



monocle3_leiden_wrapper <- function(input_data, resolution, k, ...) {
  
fake_meta_data <- data.frame(subsets = rep(1, nrow(input_data)), 
                             row.names = rownames(input_data))
  
igraph::membership(monocle3:::leiden_clustering(input_data, 
                                               pd = fake_meta_data, 
                                               resolution_parameter=resolution,
                                               k = k,
                                               random_seed = 12,
                                               weight = FALSE,
                                               ...)[["optim_res"]])
}

#' A resTuner-compatible wrapper function for clustering using Seurat
#' 
#' @param input_data A matrix where features correspond to rows and cells correspond to columns.
#' @param resolution Either a positive integer specifying the exact number of clusters desired or a vector of length two specifying range.
#' @param k The number of rounds of clusterings to perform. Default: 5.
#' @param ... Arguments passed to Seurat::FindClusters
#' @return A vector of length \code{ncol(input_data)} containing cluster assignments.
#'
#' @author Kevin Brulois
#' @export

seurat_cluster_wrapper <- function(input_data, .resolution = 0.05, k = 5, ...) {
  
  nn <- Seurat::FindNeighbors(input_data, k.param = k)
  Seurat::FindClusters(nn[["snn"]], resolution = .resolution, ...)[,1]
  
}







