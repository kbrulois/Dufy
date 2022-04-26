
#' Download GSE140348 and format as SingleCellExperiment
#' 
#' 
#' @param local.dir A directory where data should be downloaded. Default: tempdir().
#' @return A \code{SingleCellExperiment}
#' 
#' @author Kevin Brulois
#' @export




extractButcherBEC <- function(local.dir = tempdir()) {
  dir.create(local.dir)
  setwd(local.dir)
  system("curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140348/suppl/GSE140348_RAW.tar")
  system("curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140348/suppl/GSE140348_cell_meta_data.csv.gz")
  system("curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140348/suppl/GSE140348_features.tsv.gz")
  untar(grep(".tar", list.files(), value = TRUE))
  
  sample_base_names <- c("PLN1", "PLN2", "PLN3")
  file.names <- list.files()
  getFile <- function(samp_name, data_type) {
    paste0(local.dir, "/", file.names[grepl(paste(samp_name, data_type, sep = "_"), file.names)])}
  
  dat.comb <- list()
  
  for(x in sample_base_names) {
    
    dat <- as(Matrix::readMM(file = getFile(x, "matrix")), "dgCMatrix")
    feature.names = read.delim(getFile(x, "features"), 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    barcode.names = read.delim(getFile(x, "barcodes"), 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    
    colnames(dat) = sub("-1", "", barcode.names$V1)
    colnames(dat) = paste0(colnames(dat), "_", x)
    rownames(dat) = feature.names$V1
    
    dat.comb[[x]] <- dat
    
  }
  
  
  samps <- do.call(c, lapply(sample_base_names, function(x) rep(x, ncol(dat.comb[[x]]))))
  dat.comb <- Map(function(x) {
    message(paste0("normalizing ", x))
    
    sce.temp <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=dat.comb[[x]]), 
                                                           rowData = rownames(dat.comb[[x]]),
                                                           colData = colnames(dat.comb[[x]]))
    
    
    nGene <- apply(SingleCellExperiment::counts(sce.temp), 2, function(y) sum(y > 0))
    sce.temp <- sce.temp[ ,nGene > 100]
    nUMI = Matrix::colSums(SingleCellExperiment::counts(sce.temp))
    sce.temp <- sce.temp[, nUMI > 500]
    
    clusters <- scran::quickCluster(sce.temp,
                                    min.size=min(70, ncol(sce.temp)))
    
    scran::computeSumFactors(sce.temp, 
                             clusters=clusters, 
                             positive = T)
    
  }, names(dat.comb))
  
  dat.comb <- batchelor::multiBatchNorm(dat.comb, 
                                        #subset.row = norm.genes, 
                                        min.mean = 1, normalize.all = TRUE)
  
  dat.comb <- Map(function(x) {
    sce.mb.temp <- dat.comb[[x]]
    SummarizedExperiment::colData(sce.mb.temp)$sample <- factor(rep(x, ncol(sce.mb.temp)))
    return(sce.mb.temp)
  }, names(dat.comb))
  
  
  dat.comb <- do.call(SingleCellExperiment::cbind, dat.comb)
  
  dat.comb <- scater::logNormCounts(dat.comb)
  
  meta_data <- S4Vectors::DataFrame(data.table::fread(paste0(local.dir, "/GSE140348_cell_meta_data.csv.gz")))
  unlink(local.dir, recursive = TRUE)
  meta_data <- S4Vectors::DataFrame(meta_data, sizeFactor = SingleCellExperiment::colData(dat.comb)[["sizeFactor"]], row.names = meta_data$barcodes)
  SummarizedExperiment::colData(dat.comb) <- meta_data
  dat.comb
}






#' Download a subset of GSE140348 and format as SingleCellExperiment
#' 
#' 
#' @param local.dir A directory where data should be downloaded. Default: tempdir().
#' @return A \code{SingleCellExperiment}
#' 
#' @author Kevin Brulois
#' @export




extractPLN1 <- function(local.dir = tempdir()) {
  dir.create(local.dir)
  og.dir <- getwd()
  on.exit({setwd(og.dir)})
  setwd(local.dir)
  system("curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140348/suppl/GSE140348_RAW.tar")
  system("curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140348/suppl/GSE140348_cell_meta_data.csv.gz")
  system("curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140348/suppl/GSE140348_features.tsv.gz")
  untar(grep(".tar", list.files(), value = TRUE))
  
  x <- "PLN1"
  file.names <- list.files()
  getFile <- function(samp_name, data_type) {
    paste0(local.dir, "/", file.names[grepl(paste(samp_name, data_type, sep = "_"), file.names)])}
    
    dat <- as(Matrix::readMM(file = getFile(x, "matrix")), "dgCMatrix")
    feature.names = read.delim(getFile(x, "features"), 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    barcode.names = read.delim(getFile(x, "barcodes"), 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    
    colnames(dat) = sub("-1", "", barcode.names$V1)
    colnames(dat) = paste0(colnames(dat), "_", x)
    rownames(dat) = feature.names$V1
    
    sce.temp <- SingleCellExperiment::SingleCellExperiment(assays = list(counts=dat), 
                                                           rowData = rownames(dat),
                                                           colData = colnames(dat))
    clusters <- scran::quickCluster(sce.temp,
                                    min.size=min(70, ncol(sce.temp)))
    
    sce.temp <- scran::computeSumFactors(sce.temp, 
                                         clusters=clusters, 
                                         positive = T)
  
    sce.temp <- scater::logNormCounts(sce.temp)
  
  SingleCellExperiment::counts(sce.temp) <- NULL
  
  meta_data <- S4Vectors::DataFrame(data.table::fread(paste0(local.dir, "/GSE140348_cell_meta_data.csv.gz")))
  unlink(local.dir, recursive = TRUE)
  color.key <- meta_data[!is.na(meta_data[["color.scheme"]]),c("color.scheme", "color.scheme.key")]
  color.key <- setNames(color.key[["color.scheme"]], nm = color.key[["color.scheme.key"]])
  meta_data <- meta_data[meta_data$barcodes %in% colnames(sce.temp), !colnames(meta_data) %in% c("color.scheme", "color.scheme.key")]
  meta_data <- S4Vectors::DataFrame(meta_data, sizeFactor = SingleCellExperiment::colData(sce.temp)[["sizeFactor"]], row.names = meta_data$barcodes)
  SummarizedExperiment::colData(sce.temp) <- meta_data
  SummarizedExperiment::metadata(sce.temp) <- list(color.key = color.key)
  sce.temp[ ,!is.na(SingleCellExperiment::colData(sce.temp)[["fig1_tSpace1"]])]
}

