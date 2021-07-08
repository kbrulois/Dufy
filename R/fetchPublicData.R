
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
  system("curl -O ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140348/suppl/GSE140348_RAW.tar")
  system("curl -O ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140348/suppl/GSE140348_cell_meta_data.csv.gz")
  system("curl -O ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140348/suppl/GSE140348_features.tsv.gz")
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
  
  SingleCellExperiment::counts(dat.comb) <- NULL
  
  meta_data <- S4Vectors::DataFrame(data.table::fread("~/Desktop/ButcherBEC/GSE140348_cell_meta_data.csv.gz"))
  unlink(local.dir, recursive = TRUE)
  meta_data <- S4Vectors::DataFrame(meta_data, sizeFactor = SingleCellExperiment::colData(dat.comb)[["sizeFactor"]], row.names = meta_data$barcodes)
  SummarizedExperiment::colData(dat.comb) <- meta_data
  dat.comb
}

