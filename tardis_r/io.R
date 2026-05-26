#' Load BGI GEM file into Seurat object
#'
#' Port of Python \code{load_bin}. Returns a Seurat object with counts and
#' spatial coordinates in \code{misc$spatial}.
#'
#' @param gem_file Path to GEM TSV.
#' @param bin_size Bin size for aggregation.
#' @param assay Assay name.
#' @return Seurat object.
#' @export
load_bin <- function(gem_file, bin_size, assay = "RNA") {
  dat <- utils::read.table(
    gem_file, sep = "\t", header = TRUE, comment.char = "#",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (names(dat)[1L] != "geneID") {
    names(dat)[1L] <- "geneID"
  }

  dat$x <- dat$x - min(dat$x)
  dat$y <- dat$y - min(dat$y)
  dat$xp <- (dat$x %/% bin_size) * bin_size
  dat$yp <- (dat$y %/% bin_size) * bin_size
  dat$xb <- floor(dat$xp / bin_size + 1)
  dat$yb <- floor(dat$yp / bin_size + 1)
  dat$bin_ID <- max(dat$xb) * (dat$yb - 1) + dat$xb

  trans_x_xb <- aggregate(x ~ xb, data = unique(dat[, c("x", "xb")]), FUN = function(v) floor(mean(v)))
  trans_y_yb <- aggregate(y ~ yb, data = unique(dat[, c("y", "yb")]), FUN = function(v) floor(mean(v)))
  trans_matrix <- merge(
    expand.grid(xb = trans_x_xb$xb, yb = trans_y_yb$yb),
    trans_x_xb, by = "xb"
  )
  trans_matrix <- merge(trans_matrix, trans_y_yb, by = "yb")
  trans_matrix$bin_ID <- max(trans_matrix$xb) * (trans_matrix$yb - 1) + trans_matrix$xb
  trans_matrix$in_tissue <- 1L

  tissue_positions <- data.frame(
    barcodes = as.character(trans_matrix$bin_ID),
    in_tissue = trans_matrix$in_tissue,
    array_row = trans_matrix$yb,
    array_col = trans_matrix$xb,
    pxl_row_in_fullres = trans_matrix$y,
    pxl_col_in_fullres = trans_matrix$x,
    row.names = as.character(trans_matrix$bin_ID),
    stringsAsFactors = FALSE
  )

  count_col <- if ("MIDCount" %in% names(dat)) "MIDCount" else "MIDCounts"
  dat <- stats::aggregate(
    dat[[count_col]],
    by = list(geneID = dat$geneID, xb = dat$xb, yb = dat$yb),
    FUN = sum
  )
  names(dat)[4L] <- count_col
  dat$bin_ID <- max(dat$xb) * (dat$yb - 1) + dat$xb

  genes <- unique(dat$geneID)
  barcodes <- unique(dat$bin_ID)
  gene_hash <- stats::setNames(seq_along(genes) - 1L, genes)
  barcode_hash <- stats::setNames(seq_along(barcodes) - 1L, as.character(barcodes))

  i <- barcode_hash[as.character(dat$bin_ID)] + 1L
  j <- gene_hash[dat$geneID] + 1L
  x <- dat[[count_col]]

  mat <- Matrix::sparseMatrix(
    i = i, j = j, x = x,
    dims = c(length(barcodes), length(genes))
  )
  rownames(mat) <- names(barcode_hash)
  colnames(mat) <- names(gene_hash)

  obj <- Seurat::CreateSeuratObject(counts = mat, assay = assay)
  obj <- Seurat::AddMetaData(obj, metadata = tissue_positions[colnames(obj), , drop = FALSE])
  coords <- as.matrix(tissue_positions[colnames(obj), c("pxl_row_in_fullres", "pxl_col_in_fullres")])
  .set_spatial(obj, coords)
  obj
}
