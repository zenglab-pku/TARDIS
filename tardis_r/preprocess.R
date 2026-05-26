#' @importFrom stats aggregate
#' @importFrom utils read.table write.table
NULL

#' QC barplot for guide GEM duplicates
#' @param gem_path GEM file path.
#' @param guide_prefix Optional guide prefix filter.
#' @param fig_path Optional output figure path.
#' @export
qc_guide_bins <- function(gem_path, guide_prefix = NULL, fig_path = NULL) {
  df <- utils::read.table(gem_path, sep = "\t", header = TRUE, comment.char = "#",
                          stringsAsFactors = FALSE, check.names = FALSE)
  if (names(df)[1L] != "geneID") names(df)[1L] <- "geneID"
  if (!is.null(guide_prefix)) {
    df <- df[startsWith(df$geneID, guide_prefix), , drop = FALSE]
  }
  dup_df <- df[duplicated(df[, c("x", "y")]) | duplicated(df[, c("x", "y")], fromLast = TRUE), ]
  dedup_df <- dup_df[!duplicated(dup_df[, c("x", "y")]), , drop = FALSE]
  single_df <- df[!duplicated(df[, c("x", "y")]) & !duplicated(df[, c("x", "y")], fromLast = TRUE), ]

  x <- c("Total CID", "Singlet CID", "Doublet CID", "Dedup CID", "Combined CID")
  y <- c(nrow(df), nrow(single_df), nrow(dup_df), nrow(dedup_df), nrow(dedup_df) + nrow(single_df))

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mar = c(5, 4, 2, 1))
  bp <- barplot(y, names.arg = x, col = grDevices::adjustcolor(grDevices::tab.colors(5), 0.5),
                border = NA, ylab = "Count", las = 2)
  text(bp, y + max(y) / 200, labels = y, pos = 3, cex = 0.8)
  if (!is.null(fig_path)) {
    grDevices::png(fig_path, width = 1200, height = 900, res = 300)
    on.exit(dev.off(), add = TRUE)
    barplot(y, names.arg = x, col = grDevices::adjustcolor(grDevices::tab.colors(5), 0.5),
            border = NA, ylab = "Count", las = 2)
    text(bp, y + max(y) / 200, labels = y, pos = 3, cex = 0.8)
  }
  invisible(NULL)
}

#' @keywords internal
.filter_guide_long_df <- function(
    df,
    guide_prefix = NULL,
    binarilize = FALSE,
    assign_pattern = "max",
    filter_threshold = NULL
) {
  if (nrow(df) == 0L) return(df)
  if (!is.null(guide_prefix)) df <- df[startsWith(df$geneID, guide_prefix), , drop = FALSE]
  if (!is.null(filter_threshold)) df <- df[df$MIDCount > filter_threshold, , drop = FALSE]
  if (nrow(df) == 0L) return(df)

  single_df <- df[!duplicated(df[, c("x", "y")]) & !duplicated(df[, c("x", "y")], fromLast = TRUE), , drop = FALSE]
  dup_df <- df[duplicated(df[, c("x", "y")]) | duplicated(df[, c("x", "y")], fromLast = TRUE), , drop = FALSE]

  if (assign_pattern == "max") {
    dup_df$max_count <- ave(dup_df$MIDCount, dup_df$x, dup_df$y, FUN = max)
    dedup_df <- dup_df[dup_df$MIDCount == dup_df$max_count, , drop = FALSE]
    dedup_df$max_count <- NULL
  } else if (assign_pattern == "drop") {
    dedup_df <- dup_df[0L, , drop = FALSE]
  } else if (assign_pattern == "all") {
    dedup_df <- dup_df
  } else {
    stop("Invalid assign_pattern. Use 'max', 'drop', or 'all'.")
  }

  out <- rbind(single_df, dedup_df)
  if (binarilize) out$MIDCount <- 1L
  out
}

#' @keywords internal
.get_spatial_xy_df <- function(object) {
  coords <- .get_spatial(object)
  if (ncol(coords) >= 3L) {
    y <- coords[, 2L]
    x <- coords[, 3L]
  } else {
    y <- coords[, 1L]
    x <- coords[, 2L]
  }
  data.frame(barcode = .cell_names(object), x = x, y = y, stringsAsFactors = FALSE)
}

#' @keywords internal
.adata_to_guide_long <- function(object, assay = NULL) {
  mat <- .get_counts(object, assay)
  mat <- Matrix::summary(mat)
  if (nrow(mat) == 0L) {
    return(data.frame(geneID = character(), barcode = character(), MIDCount = numeric(),
                      x = numeric(), y = numeric(), stringsAsFactors = FALSE))
  }
  df <- data.frame(
    geneID = rownames(mat)[mat$j],
    barcode = colnames(mat)[mat$i],
    MIDCount = mat$x,
    stringsAsFactors = FALSE
  )
  spatial <- .get_spatial_xy_df(object)
  merge(df, spatial, by = "barcode", all.x = TRUE)
}

#' @keywords internal
.long_to_seurat <- function(df, template, guide_names, assay = NULL) {
  assay <- .default_assay(template, assay)
  if (nrow(df) == 0L) {
    out <- template[guide_names, ]
    out <- .set_counts(out, Matrix::sparseMatrix(
      i = integer(), j = integer(), x = numeric(),
      dims = c(ncol(out), nrow(out))
    ), assay)
    return(out)
  }

  barcodes <- unique(df$barcode)
  genes <- guide_names
  gi <- match(df$geneID, genes)
  bi <- match(df$barcode, barcodes)
  mat <- Matrix::sparseMatrix(
    i = bi, j = gi, x = df$MIDCount,
    dims = c(length(barcodes), length(genes))
  )
  colnames(mat) <- genes
  rownames(mat) <- barcodes

  cells <- intersect(.cell_names(template, assay), barcodes)
  mat <- mat[match(cells, barcodes), match(genes, colnames(mat)), drop = FALSE]
  rownames(mat) <- cells
  colnames(mat) <- genes
  out <- Seurat::CreateSeuratObject(counts = Matrix::t(mat), assay = assay,
                                    meta.data = template@meta.data[cells, , drop = FALSE])

  coords <- template@misc$spatial
  if (!is.null(coords)) {
    idx <- match(cells, .cell_names(template, assay))
    out <- .set_spatial(out, coords[idx, , drop = FALSE])
  }
  out
}

#' Filter guide reads from GEM file
#' @export
filter_guide_reads <- function(
    gem_path,
    guide_prefix = NULL,
    output_path = NULL,
    binarilize = FALSE,
    assign_pattern = "max",
    filter_threshold = NULL
) {
  df <- utils::read.table(gem_path, sep = "\t", header = TRUE, comment.char = "#",
                          stringsAsFactors = FALSE, check.names = FALSE)
  if (names(df)[1L] != "geneID") names(df)[1L] <- "geneID"
  out <- .filter_guide_long_df(df, guide_prefix, binarilize, assign_pattern, filter_threshold)
  if (is.null(output_path)) return(out)
  utils::write.table(out, output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(NULL)
}

#' Filter guide reads on Seurat / 10x h5 object
#' @param object_or_path Seurat object or path to h5/h5ad (10x h5 via Seurat).
#' @param copy Return modified copy (default TRUE).
#' @export
filter_guide_reads_h5 <- function(
    object_or_path,
    guide_prefix = NULL,
    output_path = NULL,
    binarilize = FALSE,
    assign_pattern = "max",
    filter_threshold = NULL,
    copy = TRUE,
    assay = "RNA"
) {
  if (is.character(object_or_path)) {
    object <- read_guide_h5(object_or_path, assay = assay)
  } else if (copy) {
    mat <- .get_counts(object_or_path, assay)
    object <- Seurat::CreateSeuratObject(counts = mat, assay = assay, meta.data = .get_obs(object_or_path))
    if (!is.null(object_or_path@misc$spatial)) {
      object <- .set_spatial(object, object_or_path@misc$spatial)
    }
  } else {
    object <- object_or_path
  }

  long_df <- .adata_to_guide_long(object, assay)
  filtered <- .filter_guide_long_df(long_df, guide_prefix, binarilize, assign_pattern, filter_threshold)
  genes <- if (!is.null(guide_prefix)) {
    .feature_names(object, assay)[startsWith(.feature_names(object, assay), guide_prefix)]
  } else {
    .feature_names(object, assay)
  }
  out <- .long_to_seurat(filtered, object, genes, assay)
  if (is.null(output_path)) return(out)
  tryCatch(
    Seurat::SaveH5Seurat(out, filename = output_path, overwrite = TRUE),
    error = function(e) saveRDS(out, file = sub("\\.h5[^.]*$", ".rds", output_path))
  )
  invisible(NULL)
}

#' Remove mito, ribo, housekeeping, and lncRNA-like genes
#' @export
remove_mito_ribo_hk_lnc_genes <- function(object, housekeeping_list = "He2020Nature_mouseHK.txt") {
  genes <- .feature_names(object)
  drop <- grepl("^Mt", genes) | grepl("^mt-", genes) | grepl("^Gm", genes) |
    grepl("^Rp", genes) | grepl("Rik", genes)
  if (file.exists(housekeeping_list)) {
    hk <- scan(housekeeping_list, what = "", sep = "\t", quiet = TRUE)
    drop <- drop | genes %in% hk
  }
  object[!(drop), ]
}

#' Combine guide replicates (sum columns sharing prefix before '_')
#' @export
combine_guide_replicates <- function(object, assay = NULL) {
  assay <- .default_assay(object, assay)
  mat <- as.matrix(.get_counts(object, assay))
  prefixes <- sub("_.*$", "", rownames(mat))
  combined <- rowsum(mat, group = prefixes)
  out <- Seurat::CreateSeuratObject(counts = combined, assay = assay, meta.data = .get_obs(object))
  if (!is.null(object@misc$spatial)) out <- .set_spatial(out, object@misc$spatial)
  out
}

#' NMF clustering on highly variable genes
#' @export
nmf_clustering <- function(
    object,
    n_components = 10L,
    random_state = 42L,
    max_iter = 1000L,
    verbose = 0L,
    n_top_genes = 2000L,
    assay = NULL
) {
  assay <- .default_assay(object, assay)
  obj <- object
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = n_top_genes, verbose = FALSE)
  hvg <- Seurat::VariableFeatures(obj)
  mat <- as.matrix(.get_counts(obj, assay)[hvg, ])
  mat <- t(mat)

  if (requireNamespace("RcppML", quietly = TRUE)) {
    model <- RcppML::nmf(mat, k = n_components, maxit = max_iter, seed = random_state, verbose = verbose > 0L)
    W <- model$w
    H <- model$h
  } else if (requireNamespace("NMF", quietly = TRUE)) {
    fit <- NMF::nmf(mat, rank = n_components, method = "brunet", seed = random_state, maxIter = max_iter)
    W <- NMF::basis(fit)
    H <- NMF::coef(fit)
  } else {
    stop("Install 'RcppML' or 'NMF' for NMF clustering.")
  }

  obj@misc$X_nmf <- W
  obj@misc$X_nmf_components <- H
  obj@misc$nmf_features <- hvg
  obj
}

#' NMF consensus clustering
#' @export
nmf_consensus <- function(
    object,
    min_clusters = 4L,
    max_clusters = 10L,
    n_resamples = 100L,
    resample_frac = 0.8,
    random_state = 42L,
    n_cluster_genes = 50L
) {
  if (is.null(object@misc$X_nmf)) stop("Run nmf_clustering() first.")
  W <- object@misc$X_nmf
  H <- object@misc$X_nmf_components
  hvg <- object@misc$nmf_features
  n <- ncol(W)
  corr <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in i:n) {
      corr[i, j] <- stats::cor(W[, i], W[, j])
      corr[j, i] <- corr[i, j]
    }
  }

  set.seed(random_state)
  # Simple hierarchical clustering fallback (consensusclustering is Python-only)
  d <- as.dist(1 - corr)
  hc <- stats::hclust(d, method = "average")
  k <- min(max_clusters, max(min_clusters, 2L))
  labels <- stats::cutree(hc, k = k)

  for (i in seq_len(k)) {
    idx <- which(labels == i)
    top_genes <- hvg[order(rowMeans(H[idx, , drop = FALSE]), decreasing = TRUE)][seq_len(min(n_cluster_genes, length(hvg)))]
    object <- Seurat::AddModuleScore(object, features = list(top_genes), name = paste0("nmf_cluster_", i - 1L))
    sc <- object@meta.data[[paste0("nmf_cluster_", i - 1L, "1")]]
    object@meta.data[[paste0("nmf_cluster_", i - 1L)]] <- (sc - mean(sc)) / stats::sd(sc)
    object@meta.data[[paste0("nmf_cluster_", i - 1L, "1")]] <- NULL
  }
  object
}

#' DBSCAN density region for a guide
#' @export
dbscan_density_region <- function(
    object,
    guide,
    eps = 5,
    min_samples = 20,
    label = "0",
    mode = c("most", "label"),
    density_level = 7L,
    step = 1,
    bandwidth = 2,
    region_label_columns = c("pxl_row_in_fullres", "pxl_col_in_fullres"),
    assay = NULL
) {
  mode <- match.arg(mode)
  assay <- .default_assay(object, assay)
  mat <- .get_counts(object, assay)
  guide_idx <- which(rownames(mat) == guide)
  if (length(guide_idx) != 1L) stop("Guide not found: ", guide)
  cells <- which(mat[guide_idx, ] > 0)
  if (length(cells) == 0L) stop("No cells with guide ", guide)

  sub <- object[, colnames(mat)[cells]]
  coords <- .get_spatial(sub)
  if (ncol(coords) >= 3L) guide_dots <- coords[, 2:3, drop = FALSE] else guide_dots <- coords

  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required for dbscan_density_region().")
  }
  cl <- dbscan::dbscan(guide_dots, eps = eps, minPts = min_samples)$cluster
  n_clust <- length(setdiff(unique(cl), -1L))
  message("Number of Clusters (excluding noise): ", n_clust)

  sub@meta.data$dbscan_cluster <- ifelse(cl == -1L, "Noise", as.character(cl))
  labels <- sub@meta.data$dbscan_cluster

  if (mode == "most") {
    tab <- table(labels[labels != "-1" & labels != "Noise"])
    if (length(tab) == 0L) stop("No clusters (other than noise) found.")
    region_label <- names(which.max(tab))
  } else {
    region_label <- as.character(label)
  }

  sel <- labels == region_label
  if (!any(sel)) stop("No points found for label ", region_label)

  cluster_points <- guide_dots[sel, , drop = FALSE]
  x_min <- min(cluster_points[, 1]) - 1
  x_max <- max(cluster_points[, 1]) + 1
  y_min <- min(cluster_points[, 2]) - 1
  y_max <- max(cluster_points[, 2]) + 1
  x_range <- seq(x_min, x_max, by = step)
  y_range <- seq(y_min, y_max, by = step)
  grid <- expand.grid(x = x_range, y = y_range)

  if (requireNamespace("ks", quietly = TRUE)) {
    kde_fit <- ks::kde(cluster_points, bandwidth = bandwidth)
    z <- exp(ks::kde(grid, H = kde_fit$H, binned = FALSE)$eval)
  } else {
    d <- as.matrix(dist(rbind(cluster_points, grid)))
    d <- d[seq_len(nrow(cluster_points)), -seq_len(nrow(cluster_points)), drop = FALSE]
    w <- exp(-rowMeans(d) / bandwidth)
    z <- colSums(w) / nrow(cluster_points)
  }

  levels <- seq(min(z), max(z), length.out = density_level)
  thresh <- levels[min(3L, length(levels))]
  core <- grid[z >= thresh, , drop = FALSE]
  colnames(core) <- region_label_columns
  list(guide_adata = sub, core_points_df = core)
}

#' Plot guide gene summary (spots vs expression rank)
#' @export
plot_guide_gene_summary <- function(object, assay = NULL) {
  mat <- .get_counts(object, assay)
  expr_count <- Matrix::colSums(mat > 0)
  total_expr <- Matrix::colSums(mat)
  mean_expr <- Matrix::colMeans(mat)
  gene_rank <- length(total_expr) - rank(total_expr) + 1

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  cols <- if (any(mean_expr > 0)) {
    grDevices::colorRampPalette(c("#440154", "#21908C", "#FDE725"))(100)
  } else rep("gray", 100)
  brks <- seq(min(mean_expr[mean_expr > 0], na.rm = TRUE), max(mean_expr), length.out = 101)
  if (!any(is.finite(brks))) brks <- c(0, 1)

  plot(expr_count, gene_rank, pch = 16, cex = 0.4, col = cols[findInterval(mean_expr, brks)],
       log = "x", xlab = "Number of spots detected", ylab = "Gene total expression rank",
       main = "")
  axis(2, labels = FALSE)
  legend_image <- grDevices::colorRampPalette(c("#440154", "#21908C", "#FDE725"))(50)
  # simple colorbar
  invisible(list(expr_count = expr_count, gene_rank = gene_rank, mean_expr = mean_expr))
}
