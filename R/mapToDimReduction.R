
#' Use nearest neighbor relationships to approximate the positions of new data points within pre-computed low-dimensional representations of a reference dataset
#' 
#' 
#' 
#' @param refData A matrix where cells correspond to rows and features correspond to columns.
#' @param newData A matrix where cells correspond to rows and features correspond to columns.
#' @param int_dat A matrix of batch corrected data. If left unset, batch correction will be performed on refData and newData using the batchelor::fastMNN function.
#' @param refReducedDims A list where each element is a dimensionality reduction of refData.
#' @param nn_k For a given new data point, the number of nearest neighbors in the reference data to use.
#' @param min_k The minimum number of nearest neighbors out of nn_k nearest neighbors from which to select the "best" combination of nearest neighbors. If set, mapping will be performed using the combination of l (where min_k =< l <= nn_k) nearest neighbors in the reference data whose centroid is the closest to its corresponding new data point in batch-corrected space. Recommended for use with nn_k < 12 as the number of combinations (and computation time) increases exponentially with the size of nn_k.
#' @param nn_span_max For a given refReducedDim, nn_span is the ratio between the diagnal spanning the nearest neighbors and the diagnal spanning all data points. If set and if nn_span of the reference nearest neighbors for a given new data point is greater than \code{nn_span_max}, then nearest neighbors will be iteratively removed by selecting the most distant neighbor (defined as the nn with the largest distance to the average the other nn) and recomputing nn_span until falls below \code{nn_span_max}. 
#' @return A \code{list} of the same length as refReducedDims
#'
#' @author Kevin Brulois
#' @export

mapToDimReduction <- function(refData,
                              newData,
                              int_dat = NULL,
                              refReducedDims, 
                              nn_k = 10, 
                              min_k = 4, 
                              nn_span = 0.1) {
  


if(is.null(int_dat)) {
  int_dat <- SingleCellExperiment::reducedDim(batchelor::fastMNN(refData, newData, prop.k = 0.01), "corrected")
}

euc_dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2)) 

int_dat_dist <- euc_dist(apply(int_dat, 2, min), apply(int_dat, 2, max))

batch <- factor(c(rep("refData", ncol(refData)), rep("newData", ncol(newData))))

nn_res <- RANN::nn2(int_dat[batch == "refData", ], query = int_dat[batch == "newData", ], k = nn_k)

int_dat_dist2 <- lapply(1:ncol(newData), function(y) {
  x1 <- int_dat[batch == "newData", ][y, ]
  euc_dist(x1,apply(int_dat[nn_res[["nn.idx"]][y, ], ,drop = FALSE ], 2, mean))
})

if(!is.null(min_k)) {
  start <- Sys.time()
  best_combs <- lapply(1:nrow(nn_res[["nn.idx"]]), function(x)  {
    x1 <- int_dat[batch == "newData", ][x, ]
    tmp <- do.call(c, lapply(min_k:nn_k, function(y) combn(nn_res[["nn.idx"]][x, ], y, simplify = FALSE)))
    tmp2 <- lapply(tmp, function(y) apply(int_dat[y, , drop = FALSE], 2, mean))
    tmp3 <- lapply(tmp2, function(y) euc_dist(x1, y))
    ind <- which.min(do.call(c,tmp3))
    c(tmp[ind], tmp3[ind][[1]]/int_dat_dist)
  })
  end <- Sys.time()
  end - start
}

int_rd_dat <- lapply(refReducedDims, function(rd_dat) {
  
  dims <- ncol(rd_dat)
  rd_dist <- euc_dist(apply(rd_dat, 2, min), apply(rd_dat, 2, max))
  
  new_dat <- do.call(rbind, lapply(1:ncol(newData), function(y) {
    if(!is.null(min_k)) {
      nn_rd_dat <- rd_dat[best_combs[[y]][[1]],]
      int_dat_dist2 <- best_combs[[y]][[2]]
    } else {
      nn_rd_dat <- rd_dat[nn_res[["nn.idx"]][y, ],]
      int_dat_dist2 <- int_dat_dist2[[y]]
    }
    nn_rd_dist <- euc_dist(apply(nn_rd_dat, 2, min), apply(nn_rd_dat, 2, max))
    rd_dat_dist <- nn_rd_dist/rd_dist
    to.return <- c(unname(apply(nn_rd_dat, 2, mean)), int_dat_dist2, rd_dat_dist, rd_dat_dist/int_dat_dist2, nrow(nn_rd_dat))
    
    if(!is.null(nn_span)) {
      if(rd_dat_dist > nn_span) {
        while(rd_dat_dist > nn_span & nrow(nn_rd_dat) > 2) {
          to.remove <- which.max(do.call(c, lapply(1:nrow(nn_rd_dat), function(x) {
            euc_dist(nn_rd_dat[x,], apply(nn_rd_dat[-x, ], 2, mean))
          })))
          nn_rd_dat <- nn_rd_dat[-to.remove,]
          nn_rd_dist <- euc_dist(apply(nn_rd_dat, 2, min), apply(nn_rd_dat, 2, max))
          rd_dat_dist <- nn_rd_dist/rd_dist
        }
      }
      to.return <- c(unname(apply(nn_rd_dat, 2, mean)), to.return, rd_dat_dist, rd_dat_dist/int_dat_dist2, nrow(nn_rd_dat))
    }
    return(to.return)
  }))
  new_dat <- as.data.frame(new_dat)
  colnames(new_dat)[1:dims] <- colnames(rd_dat)
  col_labels <- c("er", "ier", "k")
  col_labels1 <- c("ir", paste0(col_labels, ".1"))
  col_labels2 <- paste0(col_labels, ".2")
  tmp <- plyr::rbind.fill(list(as.data.frame(rd_dat), as.data.frame(new_dat)))
  
  colnames(tmp)[1] <- paste0(colnames(rd_dat)[1], "|", "nn_k", nn_k, "_", "min_k", min_k, "_", "nn_span", nn_span, "_", paste0(c("sel1", "sel2")[c(!is.null(min_k),  !is.null(nn_span))], collapse = "_"))
  if(!is.null(nn_span)) {
    colnames(tmp)[dims+1] <- paste0(colnames(rd_dat)[1], "|", "nn_k", nn_k, "_", "min_k", min_k, "_", "nn_span", nn_span, "_", c("sel1"))
    colnames(tmp)[(dims+2):(2*dims)] <- colnames(rd_dat)[2:dims]
    colnames(tmp)[(2*dims +1):(2*dims + 7)] <- c(col_labels1, col_labels2)
  } else {
    colnames(tmp)[(2*dims +1):(2*dims + 7)] <- c(col_labels1)
  }
  tmp
})

}

