
 
#' Graphically-based Downsampling
#' 
#' This method starts by constructing a lattice encompassing the data and identifies the points in the data that are closest to one of the lattice points. 
#' This is done iteratively and the distance between lattice points is adjusted until the number of closest data points is approximately equal to a user-specified target_number.
#' @param data A numeric matrix of data where cells correspond to rows and genes correspond to columns
#' @param target_number The number of cells to downsample to
#' @param iter_max The maximum number of iterations allowed.
#' @param label (Optional) A label to be used in messages
#' @return A \code{list} that includes the following elements:
#' \describe{
#'   \item{indices}{A \code{vector} of indices of \code{data}}
#'   \item{lattice}{A \code{data.frame} of lattice coordinates used to select data points}
#' }
#'
#' @author Kevin Brulois
#' @export


graphicalDownsample <- function(data = dat,
                                target_number = 400,
                                iter_max = 100,
                                label = "") {
  
  error_margin <- target_number * 0.01

  target_cell_range <- c(target_number - error_margin, target_number + error_margin)
  
  n <- target_number ^ (1/ncol(data))

  cell_num <- 0
  
  iters <- 0
  
  while(cell_num < target_cell_range[1] | cell_num > target_cell_range[2]) {

    iters <- iters + 1
  
    if(iters == iter_max & ncol(data) > 2) {
      rm("upper_n", "lower_n")
      iters <- 1
      data <- data[,-ncol(data)]
      n <- target_number ^ (1/ncol(data))
      warning(label, " failed to converge using ", ncol(data) + 1, " dimensional data. Retrying with ", ncol(data), " dimensional data.")
    }
  
    bins <- apply(data, 2, function(x) seq(min(x), max(x), by =  (max(x) - min(x))/n))

    lat_mat <- expand.grid(as.data.frame(bins))

    nn_res <- RANN::nn2(data, query = lat_mat, k = 1)[["nn.idx"]][ ,1]

    indices <- unique(nn_res)
  
    cell_num <- length(indices)
  
    if(cell_num < target_cell_range[1]) lower_n <- n
  
    if(cell_num > target_cell_range[2]) upper_n <- n
  
    if(!exists("closest")) closest <- indices
  
    if(abs(length(closest) - target_number) > abs(cell_num - target_number)) {
      closest <- indices
      if(exists("upper_n") & exists("lower_n")) {
        n <- mean(c(upper_n, lower_n))
      } else {
        n <- n * target_number/cell_num
        }
      } else {
      n <- n * target_number/cell_num
    }
  
    if(iters == iter_max) {
      warning("max iterations reached for ", label, ". Obtained ", cell_num, " target cells")
      indices <- closest
      break
    }
    
  }

list(indices = indices,
     lattice = lat_mat)
  
}








