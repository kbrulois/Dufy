
#' Density-based k-means clustering
#' 
#' This algorithm performs k-means clustering and sub-clustering such that the size of each cluster is proportional to the local density. This is done iteratively and the initial size of the kmeans clusters is updated until the desired number of total clusters is reached.
#' 
#' @param data A numeric matrix of data where cells correspond to rows and genes correspond to columns
#' @param clusters The number of clusters
#' @param d A numeric vector of density values corresponding to rows of data. 
#' @param iter_max The maximum number of iterations allowed.
#' @param label (Optional) A label to be prepended to cluster names and used in messages
#' @return A \code{list} that includes the following elements:
#' \describe{
#'   \item{cluster}{A \code{vector} of length \code{nrow(data)} containing cluster assignments.}
#'   \item{d}{A \code{vector} of length \code{nrow(data)} containing density values.}
#' }
#'
#' @author Kevin Brulois
#' @export

db_kmeans <- function(data = dat,
                      clusters = 400,
                      d = NULL,
                      iter_max = 100,
                      label = "") {
  
  message("calculating clusters ", label)
  
  error_margin <- max(2, clusters * 0.01)
  
  target_cell_range <- c(clusters - error_margin, clusters + error_margin)
  
  og_cell_num <- nrow(data)
  
  if(is.null(d)) d <- optiDensity(data)
  
  max_flag <- FALSE
  
  cell_num <- 0
  
  iters <- 0
  
  min_clust_size <- 1
  
  init_clust_size_factor <- 2
  
  while(cell_num < target_cell_range[1] | cell_num > target_cell_range[2]) {
    
    iters <- iters + 1
    
    print(init_clust_size_factor)
    
    init_clust_size <- init_clust_size_factor * og_cell_num / clusters
    
    clust_num <- round(og_cell_num / init_clust_size)
    
    if(clust_num < 24) {
      max_flag <- TRUE
      clust_num <- max(clust_num, 12)
      min_clust_size <- min_clust_size + 1  
      main_clust <- as.character(Dufy::resTuner(input_data = data, cluster_number = c(clust_num -3, clust_num +3), max_iter = 40, cluster_repeats = 1))
      #main_clust <- unname(kmeans(data, data[sample(1:og_cell_num, size = clust_num, prob = d/max(d)), ], iter.max = 40)[["cluster"]])
      
    } else {
    message("clust_num is ", clust_num)
    
    main_clust <- unname(kmeans(data, data[sample(1:og_cell_num, size = clust_num, prob = d/max(d)), ], iter.max = 40)[["cluster"]])
    }
    main_clust_og <- main_clust
    for(x in unique(main_clust)) {
      clust <- main_clust == x
      clust_size <- sum(clust)
      if(clust_size > init_clust_size) {
        sub_clust_num <- ceiling(clust_size / init_clust_size)
        if(sub_clust_num < clust_size) {
          sub_clust <- unname(kmeans(data[clust,], sub_clust_num)[["cluster"]])
          sub_clust <- paste(x, sub_clust, sep = "_")
          main_clust[clust] <- sub_clust
        }
      }
    }
    
    uni_main_clust <- unique(main_clust)
    
    message("uni_main_clust is ", length(unique(uni_main_clust)))
    
    dpc <- do.call(c, lapply(uni_main_clust, function(x) {
      mean(d[main_clust == x])
    }))
    
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    dpc_scaled <- range01(dpc)
    clust_size <- table(main_clust)[uni_main_clust]
    mean_target_clust_size <- mean(clust_size)
    target_clust_size <- (mean_target_clust_size - min_clust_size) * dpc_scaled + min_clust_size
    names(target_clust_size) <- uni_main_clust
    

    for(x in uni_main_clust) {
      clust <- main_clust == x
      cl_size <- sum(clust)
      targ_size <- target_clust_size[x]
      sub_cl_num <- round(cl_size / targ_size) - 1
     
      if(sub_cl_num > 1) {
        sub_clust <- unname(kmeans(data[clust,], sub_cl_num)[["cluster"]])
        sub_clust <- paste(x, sub_clust, sep = "_")
        main_clust[clust] <- sub_clust
      }
    }
    
    cell_num <- length(unique(main_clust))
    
    message("cell_num again is ", cell_num)
    init_clust_size_factor <- init_clust_size_factor * (cell_num / clusters)
    
    if(iters == iter_max) {
      warning("max iterations reached (", label, "). Continuing with ", cell_num, " target cells")
      break
    }
  }
  
  main_clust <- paste(label, main_clust, sep = "_")
  
  return(list(cluster = main_clust,
              initial_cluster = main_clust_og,
              density = d))
  
}





#' Aggregate numeric data within clusters
#' 
#' 
#' 
#' @param data A numeric matrix of data where cells correspond to rows and genes correspond to columns
#' @param centering_function A function used to aggregate data for a given gene within a given cluster
#' @return A \code{list} that includes the following elements:
#' \describe{
#'   \item{cluster}{A \code{vector} of length \code{nrow(data)} containing cluster assignments.}
#'   \item{d}{A \code{vector} of length \code{nrow(data)} containing density values.}
#' }
#'
#' @author Kevin Brulois
#' @export


compute_centroids_n <- function(data = dat,
                              cluster = NULL,
                              centering_function = median) {
  
  do.call(rbind, lapply(unique(cluster), function(x) {
    clust <- cluster == x
    if(sum(clust) > 1) to.return <- apply(data[clust,], 2, centering_function)
    else to.return <- data[clust, ]
    return(to.return)
  }))
  
}



#' Aggregate categorical data within clusters
#' 
#' 
#' 
#' @param data A numeric matrix of data where cells correspond to rows and genes correspond to columns
#' @param ambiguous_cutoff Proportion of the most represented category in a cluster. Clusters falling below this cutoff will be considered ambiguous.
#' @return A \code{list} that includes the following elements:
#' \describe{
#'   \item{cluster}{A \code{vector} of length \code{nrow(data)} containing cluster assignments.}
#'   \item{d}{A \code{vector} of length \code{nrow(data)} containing density values.}
#' }
#'
#' @author Kevin Brulois
#' @export


compute_centroids_c <- function(data = dat,
                                cluster = NULL,
                                ambiguous_cutoff = .5) {
  
  do.call(rbind, lapply(unique(cluster), function(x) {
    clust <- cluster == x
    clust_size <- sum(clust)
    data <- as.data.frame(data)
    apply(data, 2, function(y) {
    if(clust_size > 1) {
      tab_dat <- table(y[clust])
      max_ind <- which.max(tab_dat)
      if(tab_dat[max_ind]/clust_size >= ambiguous_cutoff) {
        to.return <- names(tab_dat[max_ind])
      } else {
        to.return <- "ambiguous"
      }
    }
    else to.return <- y[clust]
    return(to.return)
  })
  }))
  
}



nz_median <- function(x) {nz <- x != 0
if(sum(nz) > 0) to.return <- median(x[nz])
else to.return <- 0
return(to.return)}




