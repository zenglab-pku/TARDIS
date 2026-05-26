#' Local (cluster-composition) statistics
#' @importFrom stats dist
NULL

#' @keywords internal
.clr_transform <- function(x) {
  lx <- log10(x + 1)
  10^(lx - sum(lx) / length(lx))
}

#' @keywords internal
.aitchison_perm_swap <- function(p_guide, p_ntc, p_swap) {
  for (k in seq_along(p_guide)) {
    if (stats::runif(1) < p_swap) {
      diff <- as.integer(p_ntc[k] - p_guide[k])
      if (diff > 0L) {
        t <- sample.int(diff, 1L)
        p_guide[k] <- p_guide[k] + t
        p_ntc[k] <- p_ntc[k] - t
      } else if (diff < 0L) {
        t <- sample.int(-diff, 1L)
        p_guide[k] <- p_guide[k] - t
        p_ntc[k] <- p_ntc[k] + t
      }
    }
  }
  list(guide = p_guide, ntc = p_ntc)
}

#' Aitchison distance between guides and reference across clusters
#' @export
aitchison_distance <- function(
    object,
    cluster_field,
    result_field = "aitchison_dist",
    reference_guide = "sgNon-targeting",
    library_key = NULL,
    n_permutations = 1000L,
    show_progress = TRUE,
    p_swap = 0.1,
    copy = FALSE,
    assay = NULL
) {
  if (copy) {
    mat <- .get_counts(object, assay)
    object <- Seurat::CreateSeuratObject(counts = mat, assay = .default_assay(object, assay),
                                         meta.data = .get_obs(object))
    if (!is.null(object@misc$spatial)) object <- .set_spatial(object, object@misc$spatial)
  }

  obs <- .get_obs(object)
  if (!cluster_field %in% colnames(obs)) stop("cluster_field not in meta.data: ", cluster_field)
  obs[[cluster_field]] <- as.character(obs[[cluster_field]])
  object@meta.data <- obs

  mat <- as.matrix(.get_counts(object, assay))
  features <- rownames(mat)
  cells <- colnames(mat)
  all_clusters <- unique(obs[[cluster_field]])

  .count_df <- function(idx_cells) {
    sub <- mat[, idx_cells, drop = FALSE]
    cl <- obs[[cluster_field]][match(colnames(sub), rownames(obs))]
    out <- matrix(0, nrow = nrow(sub), ncol = length(all_clusters))
    rownames(out) <- rownames(sub)
    colnames(out) <- all_clusters
    for (g in seq_len(nrow(sub))) {
      for (c in seq_along(all_clusters)) {
        out[g, c] <- sum(sub[g, cl == all_clusters[c]])
      }
    }
    out
  }

  if (is.null(library_key)) {
    count_df <- .count_df(seq_along(cells))
    n_clusters <- ncol(count_df)
    transform_df <- t(apply(count_df, 1, .clr_transform))
    if (!(reference_guide %in% rownames(transform_df))) {
      stop("reference_guide '", reference_guide, "' not found.")
    }
    ref <- transform_df[reference_guide, ]
    dists <- apply(transform_df, 1, function(x) .euclidean(x, ref))
    object <- .set_meta_feature(object, result_field, dists[match(features, names(dists))])

    if (!is.null(n_permutations)) {
      pvals <- rep(0, length(features))
      names(pvals) <- features
      guides <- setdiff(features, reference_guide)
      for (guide in guides) {
        obs_dist <- dists[guide]
        greater <- 0L
        for (p in seq_len(n_permutations)) {
          .progress(p, n_permutations, guide, show_progress)
          pg <- count_df[guide, ]
          pn <- count_df[reference_guide, ]
          sw <- .aitchison_perm_swap(pg, pn, p_swap)
          tg <- .clr_transform(sw$guide)
          tn <- .clr_transform(sw$ntc)
          if (.euclidean(tg, tn) > obs_dist) greater <- greater + 1L
        }
        pvals[guide] <- greater / n_permutations
      }
      object <- .set_meta_feature(object, paste0(result_field, ".p_value"), pvals)
    }
  } else {
    libs <- unique(obs[[library_key]])
    for (sample in libs) {
      idx <- which(obs[[library_key]] == sample)
      count_df <- .count_df(idx)
      transform_df <- t(apply(count_df, 1, .clr_transform))
      ref <- transform_df[reference_guide, ]
      dists <- apply(transform_df, 1, function(x) .euclidean(x, ref))
      col <- paste0(sample, result_field)
      object <- .set_meta_feature(object, col, dists[match(features, names(dists))])
    }
    if (!is.null(n_permutations)) {
      pvals <- rep(0, length(features))
      names(pvals) <- features
      for (p in seq_len(n_permutations)) {
        .progress(p, n_permutations, "Permutations", show_progress)
        perm <- list()
        for (sample in libs) {
          idx <- which(obs[[library_key]] == sample)
          count_df <- .count_df(idx)
          for (guide in setdiff(features, reference_guide)) {
            sw <- .aitchison_perm_swap(count_df[guide, ], count_df[reference_guide, ], p_swap)
            tg <- .clr_transform(sw$guide)
            tn <- .clr_transform(sw$ntc)
            perm[[sample]][[guide]] <- .euclidean(tg, tn)
          }
        }
        for (guide in setdiff(features, reference_guide)) {
          obs_vals <- vapply(libs, function(s) {
            mf <- .get_meta_features(object)
            v <- mf[guide, paste0(s, result_field)]
            if (is.na(v)) NA_real_ else v
          }, numeric(1))
          if (any(!is.na(obs_vals))) {
            obs_dist <- mean(obs_vals, na.rm = TRUE)
            sim_vals <- vapply(libs, function(s) perm[[s]][[guide]], numeric(1))
            if (mean(sim_vals, na.rm = TRUE) > obs_dist) pvals[guide] <- pvals[guide] + 1L
          }
        }
      }
      pvals <- pvals / n_permutations
      object <- .set_meta_feature(object, paste0(result_field, ".p_value"), pvals)
    }
  }

  if (copy) object else invisible(NULL)
}

#' @keywords internal
.nan_safe_dist <- function(mat) {
  mat <- as.matrix(mat)
  if (nrow(mat) != ncol(mat)) return(NULL)
  diag(mat) <- 0
  mat[is.na(mat)] <- 0
  (mat + t(mat)) / 2
}

#' @keywords internal
.kde_or_zeros <- function(vals, x_grid) {
  vals <- vals[!is.na(vals)]
  if (length(vals) > 1L) {
    d <- stats::density(vals, from = min(x_grid), to = max(x_grid), n = length(x_grid))
    y <- approx(d$x, d$y, xout = x_grid, rule = 2)$y
    y / sum(y)
  } else {
    rep(0, length(x_grid))
  }
}

#' @keywords internal
.fast_pseudo_f <- function(distance_mat, labels) {
  groups <- unique(labels)
  idx0 <- which(labels == groups[1L])
  idx1 <- which(labels == groups[2L])
  n0 <- length(idx0)
  n1 <- length(idx1)
  ss_between <- mean(distance_mat[idx0, idx1])
  ss_within_0 <- if (n0 > 1L) mean(distance_mat[idx0, idx0]) else 0
  ss_within_1 <- if (n1 > 1L) mean(distance_mat[idx1, idx1]) else 0
  n <- n0 + n1
  ss_within <- (ss_within_0 * n0 * (n0 - 1L) / 2 + ss_within_1 * n1 * (n1 - 1L) / 2) / (n * (n - 1L) / 2)
  ss_between / (ss_within + 1e-10)
}

#' @keywords internal
.braycurtis_dist <- function(x, y) {
  sum(abs(x - y)) / sum(abs(x + y))
}

#' @keywords internal
.pdist_braycurtis <- function(mat) {
  n <- nrow(mat)
  d <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      d[i, j] <- .braycurtis_dist(mat[i, ], mat[j, ])
    }
  }
  d
}

#' PERMANOVA-style pseudo-F on cluster density profiles
#' @export
permanova <- function(
    object,
    cluster_field,
    result_field = "permanova_f_value",
    reference_guide = "sgNon-targeting",
    library_key = NULL,
    count_bins = 10L,
    n_permutations = 1000L,
    show_progress = TRUE,
    copy = FALSE,
    assay = NULL
) {
  if (copy) {
    mat <- .get_counts(object, assay)
    object <- Seurat::CreateSeuratObject(counts = mat, assay = .default_assay(object, assay),
                                         meta.data = .get_obs(object))
  }

  obs <- .get_obs(object)
  obs[[cluster_field]] <- as.character(obs[[cluster_field]])
  object@meta.data <- obs
  mat <- as.matrix(.get_counts(object, assay))
  features <- rownames(mat)
  cluster_tags <- unique(obs[[cluster_field]])

  .perm_profiles <- function(idx_cells) {
    sub_mat <- mat[, idx_cells, drop = FALSE]
    cl <- obs[[cluster_field]][match(colnames(sub_mat), rownames(obs))]
    out <- list()
    for (guide in setdiff(features, reference_guide)) {
      g_counts <- sub_mat[guide, ]
      r_counts <- sub_mat[reference_guide, ]
      guide_df <- tapply(g_counts, cl, function(v) v[v > 0])
      ref_df <- tapply(r_counts, cl, function(v) v[v > 0])
      gx <- max(unlist(guide_df), 0)
      rx <- max(unlist(ref_df), 0)
      x_grid_g <- seq(0, gx + 0.5, length.out = count_bins)
      x_grid_r <- seq(0, rx + 0.5, length.out = count_bins)
      guide_cnts <- t(vapply(cluster_tags, function(ct) {
        vals <- if (ct %in% names(guide_df)) guide_df[[ct]] else numeric()
        .kde_or_zeros(vals, x_grid_g)
      }, numeric(count_bins)))
      ref_cnts <- t(vapply(cluster_tags, function(ct) {
        vals <- if (ct %in% names(ref_df)) ref_df[[ct]] else numeric()
        .kde_or_zeros(vals, x_grid_r)
      }, numeric(count_bins)))
      m <- max(ncol(guide_cnts), ncol(ref_cnts))
      if (ncol(guide_cnts) < m) guide_cnts <- cbind(guide_cnts, matrix(0, nrow(guide_cnts), m - ncol(guide_cnts)))
      if (ncol(ref_cnts) < m) ref_cnts <- cbind(ref_cnts, matrix(0, nrow(ref_cnts), m - ncol(ref_cnts)))
      perm_data <- rbind(guide_cnts, ref_cnts)
      out[[guide]] <- perm_data
    }
    out
  }

  if (is.null(library_key)) {
    perm_data_dict <- .perm_profiles(seq_len(ncol(mat)))
    pseudo_f <- vapply(names(perm_data_dict), function(guide) {
      pd <- perm_data_dict[[guide]]
      labels <- c(rep("guide", nrow(pd) / 2), rep("ref", nrow(pd) / 2))
      dm <- .nan_safe_dist(.pdist_braycurtis(pd))
      if (is.null(dm)) NA_real_ else .fast_pseudo_f(dm, labels)
    }, numeric(1))
    object <- .set_meta_feature(object, result_field, pseudo_f[match(features, names(pseudo_f))])

    if (!is.null(n_permutations)) {
      guides <- names(perm_data_dict)
      num_clusters <- nrow(perm_data_dict[[guides[1L]]]) / 2L
      perm_swaps <- matrix(stats::runif(n_permutations * num_clusters) < 0.5,
                           nrow = n_permutations, ncol = num_clusters)
      pvals <- rep(NA_real_, length(features))
      names(pvals) <- features
      for (guide in guides) {
        arr <- perm_data_dict[[guide]]
        obs_f <- pseudo_f[guide]
        higher <- 0L
        gi <- seq_len(num_clusters)
        ri <- gi + num_clusters
        for (p in seq_len(n_permutations)) {
          .progress(p, n_permutations, guide, show_progress)
          a <- arr
          idx <- which(perm_swaps[p, ])
          if (length(idx) > 0L) {
            tmp <- a[gi[idx], , drop = FALSE]
            a[gi[idx], ] <- a[ri[idx], , drop = FALSE]
            a[ri[idx], ] <- tmp
          }
          dm <- .nan_safe_dist(.pdist_braycurtis(a))
          if (!is.null(dm)) {
            f <- .fast_pseudo_f(dm, c(rep("guide", num_clusters), rep("ref", num_clusters)))
            if (!is.na(f) && !is.na(obs_f) && f > obs_f) higher <- higher + 1L
          }
        }
        pvals[guide] <- higher / n_permutations
      }
      object <- .set_meta_feature(object, paste0(result_field, ".p_value"), pvals)
    }
  } else {
    libs <- unique(obs[[library_key]])
    all_pseudo_f <- list()
    perm_data_dict <- list()
    for (sample in libs) {
      idx <- which(obs[[library_key]] == sample)
      perm_data_dict[[sample]] <- .perm_profiles(idx)
      for (guide in names(perm_data_dict[[sample]])) {
        pd <- perm_data_dict[[sample]][[guide]]
        labels <- c(rep("guide", nrow(pd) / 2), rep("ref", nrow(pd) / 2))
        dm <- .nan_safe_dist(.pdist_braycurtis(pd))
        pf <- if (is.null(dm)) NA_real_ else .fast_pseudo_f(dm, labels)
        all_pseudo_f[[guide]][[sample]] <- pf
        object <- .set_meta_feature(object, paste0(sample, result_field),
                                    stats::setNames(pf, guide))
      }
    }
    if (!is.null(n_permutations)) {
      guides <- setdiff(features, reference_guide)
      num_clusters <- nrow(perm_data_dict[[libs[1L]]][[guides[1L]]]) / 2L
      perm_swaps <- matrix(stats::runif(n_permutations * num_clusters) < 0.5,
                           nrow = n_permutations, ncol = num_clusters)
      pvals <- rep(0, length(guides))
      names(pvals) <- guides
      gi <- seq_len(num_clusters)
      ri <- gi + num_clusters
      for (guide in guides) {
        obs_f <- mean(unlist(all_pseudo_f[[guide]]), na.rm = TRUE)
        higher <- 0L
        for (p in seq_len(n_permutations)) {
          fvals <- numeric()
          for (sample in libs) {
            if (!guide %in% names(perm_data_dict[[sample]])) next
            a <- perm_data_dict[[sample]][[guide]]
            idx <- which(perm_swaps[p, ])
            if (length(idx) > 0L) {
              tmp <- a[gi[idx], , drop = FALSE]
              a[gi[idx], ] <- a[ri[idx], , drop = FALSE]
              a[ri[idx], ] <- tmp
            }
            dm <- .nan_safe_dist(.pdist_braycurtis(a))
            if (!is.null(dm)) fvals <- c(fvals, .fast_pseudo_f(dm, c(rep("guide", num_clusters), rep("ref", num_clusters))))
          }
          if (length(fvals) > 0L && mean(fvals) > obs_f) higher <- higher + 1L
        }
        pvals[guide] <- higher / n_permutations
      }
      object <- .set_meta_feature(object, paste0(result_field, ".p_value"),
                                  pvals[match(features, names(pvals))])
    }
  }

  if (copy) object else invisible(NULL)
}
