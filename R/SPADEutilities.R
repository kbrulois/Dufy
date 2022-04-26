



#' Optimized version of SPADE's local density computation
#' 
#' This function computes local density using the spade::SPADE.density function and a range of values for the kernel_mult and apprx_mult parameters. It then selects and returns the density with the greatest variance.
#' 
#' @param data A numeric matrix of data where cells correspond to rows
#' @return A \code{vector} of density values corresponding to each row of \code{data}.
#'
#' @author Kevin Brulois
#' @export


optiDensity <- function(data, selection_method = c("ks.test", "variance")) {
  
  kernal_param <- seq(1,6, 0.5)
  dens_ities <- lapply(kernal_param, function(x) {
    spade_d_tsp <- spade::SPADE.density(tbl = data, kernel_mult = x, apprx_mult = x*(1.5/5))
    if(selection_method == "ks.test") {
    res <- ks.test(spade_d_tsp, "punif", min(spade_d_tsp), max(spade_d_tsp))[["statistic"]]
    } 
    if(selection_method == "variance") {
    res <- var(spade_d_tsp)
    }
    list(stat = res,
         d = spade_d_tsp)
  })

  if(selection_method == "ks.test") {
  d_stat <- which.min(sapply(dens_ities, `[[`,1))
  }
  if(selection_method == "variance") {
  d_stat <- which.max(sapply(dens_ities, `[[`,1))
  }
  f_dens <- dens_ities[d_stat][[1]][[2]]
  
  print(paste0("optimal kernel_mult = ", kernal_param[d_stat]))
  
  f_dens
}


#' A flexible version of SPADE's density-based downsampling function
#' 
#' This is a version of SPADE's density-based downsampling function, adapted for use with standard R matrices instead of FCS files.
#' 
#' @param data A numeric matrix of data where cells correspond to rows and genes correspond to columns
#' @return A \code{vector} of indices corresponding to rows of \code{data}
#'
#' @author Kevin Brulois
#' @export


dbSpadeDownsample <- function(data,
                             d = NULL,
                             target_number=400,
                             target_pctile=NULL,
                             target_percent=NULL
) {
  
  # boundary[1]: exclusion, boundary[2]: potential target
  boundary <- quantile(d,target_pctile,names=FALSE)
  
  if (!is.null(target_percent)) {
    target_number = round(target_percent * nrow(data))
  }
  if(is.null(d)) {
    density <- Dufy::optiDensity(data, selection_method = "ks.test")
  }
  density <- d
  if (is.null(target_number)) {
    data <- subset(data,boundary/density > runif(nrow(data)))
  } else if (target_number < nrow(data)) {
    # Need to find target density such there are approximately target_number
    # remaining after downsampling. To do so we solve for the density such that
    # the sum of samples below that density plus the expected value of
    # samples retained above that density equals approximately the desired
    # number of samples
    density_s <- sort(density)
    cdf       <- rev(cumsum(1.0/rev(density_s)))
    
    # Default solution if target density smaller than any present
    boundary <- target_number/cdf[1] 
    if (boundary > density_s[1]) {  # Boundary actually falls amongst densities present
      targets <- (target_number-1:length(density_s)) / cdf 
      boundary <- targets[which.min(targets-density_s > 0)]
    }
    res  <- which(boundary/density > runif(length(density)))
  } else if (target_number > nrow(out_data)) {
    stop("More events requested than present in data")
  }
  
  return(list(indices = res,
              d = density))
}








