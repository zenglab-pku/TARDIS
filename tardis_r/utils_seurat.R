#' Internal helpers for Seurat guide objects (Python AnnData parity)
#' @keywords internal
NULL

#' @keywords internal
.default_assay <- function(object, assay = NULL) {
  if (!is.null(assay)) return(assay)
  Seurat::DefaultAssay(object)
}

#' @keywords internal
.get_counts <- function(object, assay = NULL) {
  assay <- .default_assay(object, assay)
  mat <- Seurat::GetAssayData(object, assay = assay, slot = "counts")
  if (!inherits(mat, "dgCMatrix")) {
    mat <- Matrix::Matrix(mat, sparse = TRUE)
  }
  mat
}

#' @keywords internal
.set_counts <- function(object, mat, assay = NULL) {
  assay <- .default_assay(object, assay)
  object[[assay]] <- Seurat::SetAssayData(object[[assay]], slot = "counts", new.data = mat)
  object
}

#' @keywords internal
.feature_names <- function(object, assay = NULL) {
  rownames(.get_counts(object, assay))
}

#' @keywords internal
.cell_names <- function(object, assay = NULL) {
  colnames(.get_counts(object, assay))
}

#' @keywords internal
.get_meta_features <- function(object, assay = NULL) {
  assay <- .default_assay(object, assay)
  mf <- object[[assay]]@meta.features
  if (nrow(mf) == 0L || is.null(rownames(mf))) {
    mf <- data.frame(row.names = .feature_names(object, assay))
  }
  mf
}

#' @keywords internal
.set_meta_feature <- function(object, name, values, assay = NULL) {
  assay <- .default_assay(object, assay)
  mf <- .get_meta_features(object, assay)
  if (length(values) == nrow(mf)) {
    mf[[name]] <- values
  } else {
    mf[[name]] <- NA_real_
    idx <- match(names(values), rownames(mf))
    mf[[name]][idx] <- values
  }
  object[[assay]]@meta.features <- mf
  object
}

#' @keywords internal
.get_meta_feature <- function(object, name, assay = NULL) {
  mf <- .get_meta_features(object, assay)
  if (!name %in% colnames(mf)) return(rep(NA_real_, nrow(mf)))
  mf[[name]]
}

#' @keywords internal
.set_meta_feature_col <- function(object, col, assay = NULL) {
  assay <- .default_assay(object, assay)
  mf <- .get_meta_features(object, assay)
  for (nm in names(col)) mf[[nm]] <- col[[nm]]
  object[[assay]]@meta.features <- mf
  object
}

#' @keywords internal
.get_obs <- function(object) {
  object@meta.data
}

#' @keywords internal
.set_obs_col <- function(object, col_name, values) {
  object@meta.data[[col_name]] <- values
  object
}

#' @keywords internal
.get_spatial <- function(object) {
  if (!is.null(object@misc$spatial)) return(object@misc$spatial)
  if ("spatial" %in% names(object@reductions)) {
    return(Seurat::Embeddings(object, "spatial"))
  }
  obs <- .get_obs(object)
  if (all(c("pxl_row_in_fullres", "pxl_col_in_fullres") %in% colnames(obs))) {
    return(as.matrix(obs[, c("pxl_row_in_fullres", "pxl_col_in_fullres")]))
  }
  if (all(c("array_row", "array_col") %in% colnames(obs))) {
    return(as.matrix(obs[, c("array_row", "array_col")]))
  }
  stop(
    "Cannot infer spatial coordinates. Set object@misc$spatial, ",
    "a 'spatial' reduction, or obs columns pxl_row/array coordinates."
  )
}

#' @keywords internal
.set_spatial <- function(object, coords) {
  object@misc$spatial <- coords
  object
}

#' @keywords internal
.kl_div <- function(pk, qk) {
  pk <- pk / sum(pk)
  qk <- qk / sum(qk)
  idx <- pk > 0
  sum(pk[idx] * log(pk[idx] / qk[idx]))
}

#' @keywords internal
.euclidean <- function(x, y) {
  sqrt(sum((x - y)^2))
}

#' @keywords internal
.progress <- function(i, n, label = "", show = TRUE) {
  if (!show || n <= 1L) return(invisible(NULL))
  if (i %% max(1L, floor(n / 20L)) == 0L || i == n) {
    cat(sprintf("\r%s %d/%d", label, i, n))
    if (i == n) cat("\n")
  }
  invisible(NULL)
}

#' Read 10x HDF5 or h5ad guide matrix into Seurat
#'
#' Supports 10x HDF5 (\code{Seurat::Read10X_h5}) and AnnData h5ad (CSR \code{X} group)
#' as used in the TARDIS tutorial (\code{filtered_guide_bc_matrix.h5}).
#'
#' @param path Path to .h5 / .h5ad file.
#' @param assay Assay name (default "RNA").
#' @export
read_guide_h5 <- function(path, assay = "RNA") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required.")
  }
  fmt <- .detect_h5_format(path)
  if (fmt == "tenx") {
    counts <- Seurat::Read10X_h5(path)
    obj <- Seurat::CreateSeuratObject(counts = counts, assay = assay)
    return(obj)
  }
  if (fmt == "h5ad") {
    return(.read_h5ad_seurat(path, assay = assay))
  }
  stop("Unrecognized HDF5 format: ", path)
}

#' @keywords internal
.detect_h5_format <- function(path) {
  if (requireNamespace("hdf5r", quietly = TRUE)) {
    f <- hdf5r::H5File$new(path, mode = "r")
    on.exit(f$close_all())
    if ("matrix" %in% names(f)) return("tenx")
    if ("X" %in% names(f)) return("h5ad")
  }
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    top <- rhdf5::h5ls(path, recursive = FALSE)$name
    if ("matrix" %in% top) return("tenx")
    if ("X" %in% top) return("h5ad")
  }
  # Heuristic: tutorial files are h5ad despite .h5 extension
  if (grepl("bc_matrix", path)) return("h5ad")
  "tenx"
}

#' @keywords internal
.read_h5ad_seurat <- function(path, assay = "RNA") {
  if (requireNamespace("zellkonverter", quietly = TRUE)) {
    sce <- zellkonverter::readH5AD(path)
    return(Seurat::as.Seurat(sce, data = NULL))
  }
  if (requireNamespace("reticulate", quietly = TRUE)) {
    ok <- tryCatch({
      reticulate::import("anndata", delay_load = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if (ok) return(.read_h5ad_reticulate(path, assay = assay))
  }
  if (!requireNamespace("hdf5r", quietly = TRUE) &&
      !requireNamespace("rhdf5", quietly = TRUE)) {
    stop(
      "Reading h5ad requires 'hdf5r', 'rhdf5', 'zellkonverter', or 'reticulate' (+ Python anndata)."
    )
  }

  spatial <- NULL
  obs <- NULL

  if (requireNamespace("hdf5r", quietly = TRUE)) {
    f <- hdf5r::H5File$new(path, mode = "r")
    on.exit(f$close_all())
    data <- f[["X/data"]][]
    indices <- f[["X/indices"]][] + 1L
    indptr <- f[["X/indptr"]][]
    barcodes <- as.character(f[["obs/_index"]][])
    features <- as.character(f[["var/_index"]][])
    if ("obsm" %in% names(f) && "spatial" %in% names(f[["obsm"]])) {
      spatial <- f[["obsm/spatial"]][]
    }
    obs <- data.frame(row.names = barcodes)
    for (nm in names(f[["obs"]])) {
      if (nm == "_index") next
      obs[[nm]] <- f[[paste0("obs/", nm)]][]
    }
  } else if (requireNamespace("rhdf5", quietly = TRUE)) {
    data <- rhdf5::h5read(path, "X/data")
    indices <- rhdf5::h5read(path, "X/indices") + 1L
    indptr <- rhdf5::h5read(path, "X/indptr")
    barcodes <- as.character(rhdf5::h5read(path, "obs/_index"))
    features <- as.character(rhdf5::h5read(path, "var/_index"))
    spatial <- tryCatch(rhdf5::h5read(path, "obsm/spatial"), error = function(e) NULL)
    obs <- data.frame(row.names = barcodes)
  } else {
    return(.read_h5ad_reticulate(path, assay = assay))
  }

  n_cells <- length(indptr) - 1L
  row_idx <- rep(seq_len(n_cells), diff(indptr))
  cell_gene <- Matrix::sparseMatrix(
    i = row_idx,
    j = indices,
    x = data,
    dims = c(n_cells, length(features)),
    index1 = TRUE
  )
  counts <- Matrix::t(cell_gene)
  rownames(counts) <- features
  colnames(counts) <- barcodes

  obj <- Seurat::CreateSeuratObject(counts = counts, assay = assay, meta.data = obs)
  if (!is.null(spatial)) obj <- .set_spatial(obj, spatial)
  obj
}

#' @keywords internal
.read_h5ad_reticulate <- function(path, assay = "RNA") {
  tmp <- tempfile(fileext = ".mtx")
  obs_path <- tempfile(fileext = ".csv")
  on.exit(unlink(c(tmp, obs_path)), add = TRUE)
  reticulate::py_run_string(sprintf(
    paste(
      "import anndata as ad, scipy.io as sio, pandas as pd, numpy as np",
      "a = ad.read_h5ad(%s)",
      "sio.mmwrite(%s, a.X.T.tocoo())",
      "a.obs.to_csv(%s)",
      "barcodes = list(a.obs_names)",
      "features = list(a.var_names)",
      "spatial = a.obsm['spatial'] if 'spatial' in a.obsm else None",
      sep = "\n"
    ),
    shQuote(path), shQuote(tmp), shQuote(obs_path)
  ))
  counts <- Matrix::readMM(tmp)
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  barcodes <- as.character(reticulate::py$barcodes)
  features <- as.character(reticulate::py$features)
  rownames(counts) <- features
  colnames(counts) <- barcodes
  obs <- utils::read.csv(obs_path, row.names = 1, check.names = FALSE)
  obj <- Seurat::CreateSeuratObject(counts = counts, assay = assay, meta.data = obs)
  Seurat::RowNames(obj[[assay]]) <- features
  sp <- tryCatch(reticulate::py_to_r(reticulate::py$spatial), error = function(e) NULL)
  if (!is.null(sp)) obj <- .set_spatial(obj, sp)
  obj
}

#' Subset Seurat object (features and cells)
#' @keywords internal
.subset_object <- function(object, cells = NULL, features = NULL, assay = NULL) {
  assay <- .default_assay(object, assay)
  if (!is.null(features)) {
    features <- intersect(features, .feature_names(object, assay))
    object <- object[features, ]
  }
  if (!is.null(cells)) {
    cells <- intersect(cells, .cell_names(object, assay))
    object <- object[, cells]
  }
  object
}
