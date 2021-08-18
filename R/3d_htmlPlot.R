
#' Generate a 3D scatter plot html widget with color changer
#' 
#' This function generates an interactive 3d scatter plot. Includes a color changer and
#' multi-panel wrapping. 
#' 
#' @param coordinates A m x 3 matrix of numeric xyz coordinates.
#' @param color A m x n data.frame of data to plot as color. Categorical data (factors and character vectors) will be mapped to a discrete color scheme. Numerical data will be mapped to a continuous color scale.
#' @param discrete_colors_default A character vector of hexidecimal color codes to use as discrete colors. Default: readRDS(system.file("data/jmp_discrete_colors.rds", package = "Dufy")).
#' @param discrete_colors_custom Optional. A custom set of discrete color schemes specified as a named list of named vectors of hexidecimal color codes. Only list elements whose name matches a column name of color will be used. The names attribute of each list element must correspond to each unique element within the matched color column. Default: NULL.
#' @param continuous_color_scale A character vector of hexidecimal color codes to use as a continuous color scale. Defalut: NULL.
#' @param convert_integer_to_string Whether to convert integer columns of color to character. Default: TRUE
#' @param texts Optional. A string corresponding to a column name of color containing a categorical variable. If specified,  
#' @param arrows Optional. A p x q data.frame containing arrow coordinates. The base of the arrows must be specified as xyz coodinates in columns 1:3. The tip of the arrows must be provided as xyz coordinates in column 4:6. Two additional optional columns can be included: one called "color" containing hexidecimal color codes for each arrow and one with the same name as wrap containing categorical data. 
#' @param wrap Optional. A string corresponding to a column name of color (and optionally a column name of arrows) containing a categorical variable.
#' @param include_all_var Whether to include an extra panel with all data points. Ignored if wrap is NULL.
#' @param save_widget Whether to save a html widget. Default: TRUE.
#' @param capture_2D Whether to capture a 2D orthogonal projection of the current 3d view. Default: FALSE.
#' @param out_dir A character string specifying the directory where the plot should be saved. Default: current working directory.
#' @param file_name A character string specifying the file name basename. The ".html" file extension will be appended to file_name automatically.
#' @param selfcontained Whether to save the plot as a self contained html file. See htmlwidgets::saveWidget. Default: FALSE.
#' @return html_3dPlot is called for the side effect of saving an html file. If capture_2D is TRUE, a \code{matrix} containing a 2D orthogonal projection of the current view of the data is invisibly returned.
#'
#' @author Kevin Brulois
#' @export




html_3dPlot <- function(coordinates = NULL,
                        color = NULL,
                        discrete_colors_default = readRDS(system.file("data/jmp_discrete_colors.rds", package = "Dufy")),
                        discrete_colors_custom = NULL,
                        continuous_color_scale = viridis::viridis(100, option = "inferno"),
                        convert_integer_to_string = TRUE,
                        texts = NULL,
                        arrows = NULL,
                        wrap = NULL,
                        include_all_var = TRUE,
                        save_widget = TRUE,
                        capture_2D = FALSE,
                        out_dir = getwd(),
                        file_name = "plot",
                        selfcontained = FALSE) {
  
  color <- as.data.frame(color)
  if(convert_integer_to_string) {
  color[, sapply(color, is.integer)] <- as.character(color[, sapply(color, is.integer)])
  }
  coordinates <- as.matrix(coordinates)
  cell.num <- nrow(coordinates)
  assertthat::assert_that(nrow(color) == cell.num)
  
  color_features <- colnames(color)
  isnum.dat <- do.call(c, lapply(color, is.numeric))
  
  num.dat <- color[,isnum.dat, drop = FALSE]
  cat.dat <- color[,!isnum.dat, drop = FALSE]
  
  if(dim(cat.dat)[2] != 0) {
  discrete_colors <- Map(function(x) {
    if(x %in% names(discrete_colors_custom)) {
      to.return <- discrete_colors_custom[[x]]
    } else {
      col_num <- length(unique(cat.dat[[x]]))
      while(length(discrete_colors_default) < col_num) {
        discrete_colors_default <- c(discrete_colors_default, 
                                     discrete_colors_default)
      }
      
      disc_colors <- discrete_colors_default[1:col_num]
      cat_dat_names <- unique(cat.dat[[x]])
      names(disc_colors) <- cat_dat_names[gtools::mixedorder(cat_dat_names)]
      to.return <- disc_colors
      
    }
    return(to.return)
  }, colnames(cat.dat))
  
  
  hex.dat1 <- do.call(cbind, lapply(colnames(cat.dat), function(x) {
    discrete_colors[[x]][cat.dat[[x]]]
  }))
  }
  
  if(dim(num.dat)[2] != 0 ) {
  map2colorScale <- function(x) continuous_color_scale[round(scales::rescale(x, to = c(1, length(continuous_color_scale))), 0)]
  
  hex.dat2 <- apply(num.dat, 2, map2colorScale)
  
  continuous_colors <- Map(function(x) {
    params <- c("min", "25p", "mean", "median", "75p", "nz_25p", "nz_mean", "nz_median", "nz_75p", "max")
    dat <- num.dat[[x]]
    nz_dat <- dat[dat != 0]
    inds <- c(which.min(dat), 
              which.min(abs(dat - quantile(dat, 0.25))),
              which.min(abs(dat - mean(dat))),
              which.min(abs(dat - median(dat))),
              which.min(abs(dat - quantile(dat, 0.75))),
              which.min(abs(dat - quantile(nz_dat, 0.25))),
              which.min(abs(dat - mean(nz_dat))),
              which.min(abs(dat - median(nz_dat))),
              which.min(abs(dat - quantile(nz_dat, 0.75))),
              which.max(dat))
    setNames(object = hex.dat2[,x][inds],
             nm = paste0(params, "_", round(do.call(c, lapply(inds, function(x) dat[x])), 2)))
  }, colnames(num.dat))
  }
  
  if(dim(cat.dat)[2] != 0 & dim(num.dat)[2] != 0) {
  leg_colors <- c(discrete_colors, continuous_colors)
  hex.dat <- cbind(hex.dat1, hex.dat2)
  }
  if(dim(cat.dat)[2] != 0 & dim(num.dat)[2] == 0) {
    leg_colors <- discrete_colors
    hex.dat <- hex.dat1
  }
  if(dim(cat.dat)[2] == 0 & dim(num.dat)[2] != 0) {
    leg_colors <- continuous_colors
    hex.dat <- hex.dat2
  }
  
  rgb.dat <- t(apply(hex.dat, 2, col2rgb, alpha = TRUE)) / 255
  
  rgb.dat <- cbind(rgb.dat[, seq(1, ncol(rgb.dat), 4)],
                   rgb.dat[, seq(2, ncol(rgb.dat), 4)],
                   rgb.dat[, seq(3, ncol(rgb.dat), 4)])
  
  rownames(rgb.dat) <- c(colnames(cat.dat), colnames(num.dat))
  
  rgb.dat <- rgb.dat[color_features, ]
  
  leg_colors <- leg_colors[color_features]
  
  if(is.null(wrap)) {
    wrap_vars <- "All"
  } else {
  wrap_vars <- unique(color[[wrap]])
  wrap_vars <- wrap_vars[!is.na(wrap_vars)]
  wrap_vars <- wrap_vars[gtools::mixedorder(wrap_vars)]
  if(include_all_var) wrap_vars <- c("All", wrap_vars)
  }
  
  rap <- Map(function(x) {
    if(x == "All") {to.return <- rep(TRUE, cell.num)}
    else {to.return <- color[[wrap]] == x}
    to.return[is.na(to.return)] <- FALSE
    return(to.return)
  }, wrap_vars)
  
  num_rows <- function(x) {
    if(1 <= x & x <= 2) y <- 1
    if(3 <= x & x <= 5) y <- 2
    if(6 <= x & x <= 14) y <- 3
    if(15 <= x & x <= 23) y <- 4
    if(24 <= x & x <= 39) y <- 5
    y
  }
  
  data_panels <- 1 + length(rap)
  row_num <- num_rows(length(rap))
  data_panels <- data_panels + (row_num - data_panels %% row_num)

  layout_mat <- matrix(data = c(1:data_panels),
                       nrow = row_num,
                       byrow = TRUE)
  
  rgl::clear3d()
  
  rgl::rglFonts(sans = c(system.file("data/FreeSans.ttf", package = "Dufy"), 
                         system.file("data/FreeSansBold.ttf", package = "Dufy"),
                         system.file("data/FreeSansOblique.ttf", package = "Dufy"), 
                         system.file("data/FreeSansBoldOblique.ttf", package = "Dufy")))
  
  rgl::par3d(windowRect = c(100, 100, 1224, 612), 
             family = "sans", 
             font = 4)
  parent <- rgl::currentSubscene3d()
  rgl::mfrow3d(dim(layout_mat)[1], dim(layout_mat)[2], sharedMouse = TRUE)
  rgl::layout3d(mat = layout_mat, sharedMouse = TRUE)
  mplot <- list()
  sub <- list()
  cones <- list()
  
  for(x in names(rap)) {

      mplot[[x]] <- rgl::plot3d(x = coordinates[rap[[x]], ],
                                col = hex.dat[,1][rap[[x]]],
                                size = 6, xlab = "", ylab = "", zlab = "", 
                                axes = FALSE, ann = FALSE,
                                point_antialias = TRUE)
    
  
  
  if(!is.null(arrows)) {
    # dists <- do.call(c, lapply(1:nrow(arrows), function(x) {sqrt((arrows[x,1] - arrows[x,4])^2 +
    #                                                               (arrows[x,2] - arrows[x,5])^2 +
    #                                                               (arrows[x,3] - arrows[x,6])^2)}))
    
    if(!is.null(wrap) & x != "All") {
      if(!is.null(arrows[[wrap]])) {
    arrows_sub <- arrows[arrows[[wrap]] == x, ]
      }
    } else {
      arrows_sub <- arrows
    }
    for(y in 1:nrow(arrows_sub)) {
    cones[[paste0(x,"_", y)]] <- Dufy:::cone3d(base = unname(unlist(arrows_sub[y,1:3])), tip = unname(unlist(arrows_sub[y,4:6])), rad = 0.2,
                 color = if(is.null(arrows_sub[["color"]])) {"red"} else {arrows_sub[["color"]][y]})
    
    }
  }
      if(!is.null(texts)) {
        label.centroids <- do.call(rbind, Map(function(x) {
          apply(coordinates[color[[texts]] == x, , drop = FALSE], 2, median, na.rm = TRUE)
        }, unique(color[[texts]])))
        
        rgl::text3d(x = label.centroids, 
                    texts = rownames(label.centroids), 
                    pos = 4)
        
      }
      
    rgl::bgplot3d({
      plot.new()
      title(main = if(length(rap) == 1) {""} else {x}, line = 3)
    })
  
  rgl::aspect3d(1,1,1)
  sub[[x]] <- rgl::subsceneInfo()$id
  rgl::next3d()
  }
  
#   legends <- list()
#   
#   sub2 <- rgl::subsceneInfo()$id
#   
#   
#       for(x in names(leg_colors)) {
#         legends[[x]] <- rgl::legend3d("center", names(leg_colors[[x]]), pch = 16, col = unname(leg_colors[[x]]), title = x)
# }
#     
#   leg_ids <- do.call(c, Map(function(x) legends[[x]][1], names(leg_colors)))
#   
  
  legends <- list()
  leg_text <- list()
  leg_title <- list()
  
  for(x in names(leg_colors)) {
    num <- length(leg_colors[[x]])
    legends[[x]] <- rgl::points3d(0, 0, c(1,0, -1*1:num), 
                                  size = 12, 
                                  point_antialias = TRUE, 
                                  col = c("#FFFFFF00", 
                                          "#FFFFFF00", 
                                          unname(leg_colors[[x]])), 
                                  alpha = c(0,0, rep(1, num)))
    leg_title[[x]] <- rgl::text3d(0, 0, 1, texts = x)
    leg_text[[x]] <- rgl::text3d(0, 0, 
                                 c(0, -1*1:num), 
                                 texts = c("", names(leg_colors[[x]])), 
                                 pos = 4)
  }
  
  sub2 <- rgl::subsceneInfo()$id
  
  
  if(save_widget) {
    message("saving html widget")
  widget <- magrittr::`%>%`(rgl::rglwidget(elementId = file_name),
                            
  rgl::asRow(
    "Feature Changer:",
    rgl::playwidget(
      sceneId = rgl::getWidgetId(.),
      controls = c(Map(function(x) {
      Dufy:::vertexControlMod(values = rgb.dat[,rep(rap[[x]], 3)],
                              attributes = c(rep("red", sum(rap[[x]])),
                                             rep("green", sum(rap[[x]])),
                                             rep("blue", sum(rap[[x]]))),
                              vertices = 1:sum(rap[[x]]), 
                              objid = mplot[[x]]["data"], 
                              interp = FALSE,
                              labels = color_features)}, names(rap)),
                   list(rgl::subsetControl(value = 1, 
                                           subscenes = sub2, 
                                           subsets = legends),
                        rgl::subsetControl(value = 1, 
                                           subscenes = sub2, 
                                           subsets = leg_title),
                        rgl::subsetControl(value = 1, 
                                           subscenes = sub2, 
                                           subsets = leg_text))),
      step = 1, 
      interval = 1, 
      start = 0, 
      labels = color_features,
      stop = length(color_features) - 1),
    "Rotator:",
    rgl::playwidget(
      sceneId = rgl::getWidgetId(.), 
      controls = Map(function(x) { 
        rgl::par3dinterpControl(rgl::spin3d(), 
                                0, 
                                12, 
                                steps=40, 
                                subscene = x)}, 
                     sub),
      step = 0.01, 
      loop = TRUE, 
      rate = 0.5, 
      components = c("Reverse", "Play", "Slower", "Faster", "Reset", "Slider")),
    last = 4, 
    colsize = c(1,2,1,2)
    )
  )
  if(selfcontained) {
  message("coverting to selfcontained html")
  }
  data_size <- dim(color)
  data_size <- data_size[1]*data_size[2]
  if(data_size > 500000) {
    message("warning: this is a large dataset. Converting to a self-contained html file will take a while... Consider setting selfcontained to FALSE.")
  }
  htmlwidgets::saveWidget(widget, 
                          file = path.expand(paste0(out_dir, "/", file_name, ".html")), 
                          libdir = path.expand(paste0(out_dir, "/", paste0(file_name, "_files"))), 
                          title = file_name,
                          selfcontained = selfcontained)
  
  }
  
  if(capture_2D) {
    
    
    readline(prompt="Orient plot and press [enter] to save a 2d projection") 
    
    rot.mat <- t(rgl::par3d()[["userMatrix"]][1:3, 1:3])
    x_axis <- c(1,0,0)
    y_axis <- c(0,1,0)
    ref.plane <- c(0,0,1)
    
    scalar1 <- function(x) {x / sqrt(sum(x^2))}
    doRotation <- function(x) c(rot.mat %*% x)
    x_axis.rot <- scalar1(doRotation(x_axis))
    y_axis.rot <- scalar1(doRotation(y_axis))
    ref.plane.rot <- scalar1(doRotation(ref.plane))
    
    projected_data <- t(apply(coordinates, 1, function(x) {
      c(x_axis.rot %*% (x - ref.plane.rot), y_axis.rot %*% (x - ref.plane.rot))
    }))
    
    colnames(projected_data) <- paste0("2d_proj_of_", colnames(coordinates)[1:2])
    
    projected_data <- as.data.frame(projected_data)
    plot(x = projected_data[[1]], y = projected_data[[2]], pch = 16, col = hex.dat[,1])
    
    
    return(invisible(projected_data))
  }
  rgl::close3d()
  
}







bc_stripper <- function(x) {
  stringi::stri_extract(x, regex = '[ATGC]{2}[ATGC]{1,20}[ATGC]{2}')
}

###rgl source code modification

vertexControlMod <- function (value = 0, values = NULL, vertices = 1, attributes, 
                              objid, param = seq_len(NROW(values)) - 1, interp = TRUE, labels = NULL) 
{
  attributes <- match.arg(attributes, choices = c("x", "y", 
                                                  "z", "red", "green", "blue", "alpha", "radii", "nx", 
                                                  "ny", "nz", "ox", "oy", "oz", "ts", "tt", "offset"), 
                          several.ok = TRUE)
  if (!is.null(values)) {
    ncol <- max(length(vertices), length(attributes))
    if (is.matrix(values)) 
      stopifnot(ncol == ncol(values))
    else {
      stopifnot(ncol == 1)
      values <- matrix(values, ncol = 1)
    }
    param <- c(-Inf, param, Inf)
    values <- rbind(values[1, ], values, values[nrow(values), 
                                                ])
  }
  structure(list(type = "vertexSetter", value = value, values = values, 
                 vertices = vertices - 1, attributes = attributes, objid = as.integer(objid), 
                 param = param, interp = interp, labels = labels), class = "rglControl")
}


cone3d <- function(base=c(0,0,0),tip=c(0,0,1),rad=1,n=30,draw.base=TRUE,qmesh=FALSE,
                   trans = par3d("userMatrix"), ...) {
  ax <- tip-base
  if (missing(trans) && !rgl::cur3d()) trans <- diag(4)
  ### is there a better way?
  if (ax[1]!=0) {
    p1 <- c(-ax[2]/ax[1],1,0)
    p1 <- p1/sqrt(sum(p1^2))
    if (p1[1]!=0) {
      p2 <- c(-p1[2]/p1[1],1,0)
      p2[3] <- -sum(p2*ax)
      p2 <- p2/sqrt(sum(p2^2))
    } else {
      p2 <- c(0,0,1)
    }
  } else if (ax[2]!=0) {
    p1 <- c(0,-ax[3]/ax[2],1)
    p1 <- p1/sqrt(sum(p1^2))
    if (p1[1]!=0) {
      p2 <- c(0,-p1[3]/p1[2],1)
      p2[3] <- -sum(p2*ax)
      p2 <- p2/sqrt(sum(p2^2))
    } else {
      p2 <- c(1,0,0)
    }
  } else {
    p1 <- c(0,1,0); p2 <- c(1,0,0)
  }
  degvec <- seq(0,2*pi,length=n+1)[-1]
  ecoord2 <- function(theta) {
    base+rad*(cos(theta)*p1+sin(theta)*p2)
  }
  i <- rbind(1:n,c(2:n,1),rep(n+1,n))
  v <- cbind(sapply(degvec,ecoord2),tip)
  if (qmesh) 
    ## minor kluge for quads -- draw tip twice
    i <- rbind(i,rep(n+1,n))
  if (draw.base) {
    v <- cbind(v,base)
    i.x <- rbind(c(2:n,1),1:n,rep(n+2,n))
    if (qmesh)  ## add base twice
      i.x <-  rbind(i.x,rep(n+2,n))
    i <- cbind(i,i.x)
  }
  if (qmesh) v <- rbind(v,rep(1,ncol(v))) ## homogeneous
  if (!qmesh)
    rgl::triangles3d(v[1,i],v[2,i],v[3,i],...)
  else
    return(rgl::rotate3d(rgl::qmesh3d(v,i,material=list(...)), matrix=trans))
}     



html_3dPlot_subetter_slow <- function(coordinates = NULL,
                        text = "subsets",
                        subsetter = "sample",
                        color = NULL,
                        discrete_colors_default = readRDS(system.file("data/discrete_colors.rds", package = "Dufy")),
                        discrete_colors_custom = NULL,
                        continuous_color_scale = viridis::viridis(100, option = "inferno"),
                        out.dir = ".",
                        file.name = "plot") {
  
  color <- as.data.frame(color)
  color[, sapply(color, is.integer)] <- as.character(color[, sapply(color, is.integer)])
  coordinates <- as.matrix(coordinates)
  cell.num <- nrow(coordinates)
  assertthat::assert_that(nrow(color) == cell.num)
  color_features <- colnames(color)
  isnum.dat <- do.call(c, lapply(color, is.numeric))
  
  num.dat <- color[,isnum.dat]
  cat.dat <- color[,!isnum.dat]
  
  discrete_colors <- Map(function(x) {
    if(x %in% names(discrete_colors_custom)) {
      to.return <- discrete_colors_custom[[x]]
    } else {
      col_num <- length(unique(cat.dat[[x]]))
      while(length(discrete_colors_default) < col_num) {
        discrete_colors_default <- c(discrete_colors_default, 
                                     discrete_colors_default)
      }
      
      disc_colors <- discrete_colors_default[1:col_num]
      names(disc_colors) <- unique(cat.dat[[x]])
      to.return <- disc_colors
      
    }
    return(to.return)
  }, colnames(cat.dat))
  
  
  hex.dat1 <- do.call(cbind, lapply(colnames(cat.dat), function(x) {
    discrete_colors[[x]][cat.dat[[x]]]
  }))
  
  map2colorScale <- function(x) continuous_color_scale[round(scales::rescale(x, to = c(1, length(continuous_color_scale))), 0)]
  
  hex.dat2 <- apply(num.dat, 2, map2colorScale)
  
  continuous_colors <- Map(function(x) {
    params <- c("min", "25p", "mean", "median", "75p", "nz_25p", "nz_mean", "nz_median", "nz_75p", "max")
    dat <- num.dat[[x]]
    nz_dat <- dat[dat != 0]
    inds <- c(which.min(dat), 
              which.min(abs(dat - quantile(dat, 0.25))),
              which.min(abs(dat - mean(dat))),
              which.min(abs(dat - median(dat))),
              which.min(abs(dat - quantile(dat, 0.75))),
              which.min(abs(dat - quantile(nz_dat, 0.25))),
              which.min(abs(dat - mean(nz_dat))),
              which.min(abs(dat - median(nz_dat))),
              which.min(abs(dat - quantile(nz_dat, 0.75))),
              which.max(dat))
    setNames(object = hex.dat2[,x][inds],
             nm = paste0(params, "_", round(do.call(c, lapply(inds, function(x) dat[x])), 2)))
  }, colnames(num.dat))
  
  leg_colors <- c(discrete_colors, continuous_colors)
  
  hex.dat <- cbind(hex.dat1, hex.dat2)
  colnames(hex.dat) <- c(colnames(cat.dat), colnames(num.dat))
  hex.dat <- hex.dat[,color_features]
  
  rgb.dat <- apply(hex.dat, 2, col2rgb, alpha = TRUE) / 255
  
  rgb.dat <- list(red = rgb.dat[seq(1, nrow(rgb.dat), 4), ],
                   green = rgb.dat[seq(2, nrow(rgb.dat), 4),],
                   blue = rgb.dat[seq(3, nrow(rgb.dat), 4), ])
  
  rgb.dat <- Map(function(x) {
    rownames(x) <- rownames(coordinates)
    colnames(x) <- c(colnames(cat.dat), colnames(num.dat))
    x <- x[,color_features]
    x
  }, rgb.dat)
  
  leg_colors <- leg_colors[color_features]
  
  to.plot <- as.data.frame(cbind(coordinates, hex.dat))
  
  sd <- crosstalk::SharedData$new(to.plot, key = rownames(to.plot), group = "html_3dplot")
  
  sd_red <- crosstalk::SharedData$new(rgb.dat$red, group = "html_3dplot")
  sd_green <- crosstalk::SharedData$new(rgb.dat$green, group = "html_3dplot")
  sd_blue <- crosstalk::SharedData$new(rgb.dat$blue, group = "html_3dplot")
  
  rgl::clear3d()
  
  rgl::par3d(windowRect = c(100, 100, 1224, 612))
  rgl::mfrow3d(1, 2)
  rgl::layout3d(mat = matrix(c(1, 2), nrow = 1), widths = c(1.25, 0.75))
  parent <- rgl::currentSubscene3d()
  
  mplot <- list()
  share <- list()
  
  for(x in color_features) {
  rgl::useSubscene3d(subscene = parent)
  mplot[[x]] <- rgl::plot3d(x = sd$origData()[,1],
                            y = sd$origData()[,2],
                            z = sd$origData()[,3],
                       col = sd$origData()[,x],
                       size = 6, 
                       point_antialias = TRUE, add = TRUE)
  share[[x]] <- rgl::rglShared(mplot[[x]]["data"], key = sd$key(), group = sd$groupName(), deselectedFade = 0)
  
  }
  

  if(!is.null(text)) {
    label.centroids <- do.call(rbind, Map(function(x) {
      apply(coordinates[color[[text]] == x, , drop = FALSE], 2, median, na.rm = TRUE)
    }, unique(color[[text]])))
    
    rgl::text3d(x = label.centroids, texts = rownames(label.centroids), pos = 4)
    
  }
  
  rgl::aspect3d(1,1,1)
  sub1 <- rgl::subsceneInfo()$id
  rgl::next3d()
  
  legends <- list()
  leg_text <- list()
  leg_title <- list()
  
  for(x in names(leg_colors)) {
    num <- length(leg_colors[[x]])
    legends[[x]] <- rgl::points3d(0, 0, c(1,0, -1*1:num), 
                                  size = 12, 
                                  point_antialias = TRUE, 
                                  col = c("#FFFFFF00", "#FFFFFF00", unname(leg_colors[[x]])), alpha = c(0,0, rep(1, num)))
    leg_title[[x]] <- rgl::text3d(0, 0, 1, texts = x)
    leg_text[[x]] <- rgl::text3d(0, 0, c(0, -1*1:num), texts = c("", names(leg_colors[[x]])), pos = 4)
  }
  
  sub2 <- rgl::subsceneInfo()$id
  
  widget <- magrittr::`%>%`(rgl::rglwidget(shared = share),
                            
                            rgl::asRow(
                              "Subset:",
                              crosstalk::filter_checkbox("subsetselector22", 
                                                         "", sd, ~ sample, inline = TRUE),
                              "Feature Changer:",
                              rgl::playwidget(sceneId = rgl::getWidgetId(.),
                                              controls = list(rgl::subsetControl(value = 1, subscenes = sub1, subsets = mplot),
                                                              rgl::subsetControl(value = 1, subscenes = sub2, subsets = legends),
                                                              rgl::subsetControl(value = 1, subscenes = sub2, subsets = leg_title),
                                                              rgl::subsetControl(value = 1, subscenes = sub2, subsets = leg_text)),
                                              
                                              step = 1, 
                                              interval = 1, 
                                              start = 0, 
                                              labels = color_features,
                                              stop = length(color_features) - 1),
                              "Rotator:",
                              rgl::playwidget(sceneId = rgl::getWidgetId(.), 
                                              controls = rgl::par3dinterpControl(rgl::spin3d(), 0, 12, steps=40, subscene = sub1),
                                              step = 0.01, 
                                              loop = TRUE, 
                                              rate = 0.5, 
                                              components = c("Reverse", "Play", "Slower", "Faster", "Reset", "Slider")),
  last = 1, colsize = c(1,2,1,2,1,2))
  )
  
  htmlwidgets::saveWidget(widget, file = paste0(out.dir, "/", file.name, ".html"))
  
  rgl::close3d()
  
}

