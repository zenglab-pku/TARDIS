#' Global (spatial/distribution) statistics
NULL

#' @keywords internal
.permute_guide_bins <- function(guide_matrix, idx_guide, idx_ntc, swap_rate = 0.5) {
  t_bins <- (guide_matrix[, idx_guide] > 0) | (guide_matrix[, idx_ntc] > 0)
  guide_vec <- guide_matrix[t_bins, idx_guide]
  ntc_vec <- guide_matrix[t_bins, idx_ntc]
  swap_mask <- stats::rnorm(length(guide_vec)) < (1 - 2 * swap_rate)
  temp <- guide_vec[swap_mask]
  guide_vec[swap_mask] <- ntc_vec[swap_mask]
  ntc_vec[swap_mask] <- temp
  ret_guide <- guide_matrix[, idx_guide]
  ret_ntc <- guide_matrix[, idx_ntc]
  ret_guide[t_bins] <- guide_vec
  ret_ntc[t_bins] <- ntc_vec
  list(guide = ret_guide, ntc = ret_ntc)
}

#' @keywords internal
.wasserstein_2d <- function(a, b) {
  a <- a / sum(a)
  b <- b / sum(b)
  if (requireNamespace("transport", quietly = TRUE)) {
    return(transport::wassdist(a, b))
  }
  # Fallback: L1 distance between flattened normalized distributions
  sum(abs(cumsum(a) - cumsum(b)))
}

#' @keywords internal
.calculate_kde <- function(coords, cnt_vec, x_grid, y_grid, cnt_grid) {
  if (requireNamespace("ks", quietly = TRUE)) {
    dat <- cbind(coords[, 1], coords[, 2], cnt_vec)
    kde <- ks::kde(dat)
    grid <- expand.grid(
      x = x_grid, y = y_grid, z = cnt_grid
    )
    z <- exp(ks::kde(grid, H = kde$H, binned = FALSE)$eval)
    mat <- matrix(0, length(y_grid), length(x_grid))
    for (i in seq_along(cnt_grid)) {
      slice <- z[, , i]
      if (length(dim(slice)) == 2L) mat <- mat + slice * cnt_grid[i]
    }
    return(mat)
  }
  # 2D weighted histogram fallback
  mat <- matrix(0, length(y_grid), length(x_grid))
  for (k in seq_along(cnt_vec)) {
    xi <- which.min(abs(x_grid - coords[k, 1]))
    yi <- which.min(abs(y_grid - coords[k, 2]))
    mat[yi, xi] <- mat[yi, xi] + cnt_vec[k]
  }
  mat
}

#' @keywords internal
.perm_kl_single <- function(guide_idx, idx_ntc, n_permutations, guide_matrix, epx, kl_div) {
  greater <- 0L
  n <- nrow(guide_matrix)
  for (p in seq_len(n_permutations)) {
    idx <- sample.int(n)
    perm <- guide_matrix[idx, , drop = FALSE]
    pg <- perm[, guide_idx]
    pn <- perm[, idx_ntc]
    ps <- sum(pg)
    ns <- sum(pn)
    if (ps == 0 || ns == 0) next
    pk <- pg / ps + epx
    qk <- pn / ns + epx
    pk <- pk / sum(pk)
    qk <- qk / sum(qk)
    if (.kl_div(pk, qk) > kl_div) greater <- greater + 1L
  }
  greater / n_permutations
}

#' Kullback-Leibler divergence vs reference guide
#' @export
kl_divergence <- function(
    object,
    reference_guide = "sgNon-targeting",
    result_field = "kl_div",
    guide_list = NULL,
    library_key = NULL,
    n_permutations = 1000L,
    show_progress = TRUE,
    copy = FALSE,
    n_jobs = 1L,
    assay = NULL
) {
  if (copy) {
    mat <- .get_counts(object, assay)
    object <- Seurat::CreateSeuratObject(counts = mat, assay = .default_assay(object, assay),
                                         meta.data = .get_obs(object))
  }

  features <- .feature_names(object, assay)
  if (is.null(guide_list)) guide_list <- setdiff(features, reference_guide)
  all_guides <- unique(c(guide_list, reference_guide))
  if (!(reference_guide %in% features)) {
    stop("Reference guide '", reference_guide, "' not found.")
  }

  mat <- t(as.matrix(.get_counts(object[all_guides, ], assay)))

  if (is.null(library_key)) {
    idx_ntc <- match(reference_guide, colnames(mat))
    qk_sum <- sum(mat[, idx_ntc])
    if (qk_sum == 0) stop("Reference guide sums to zero.")
    epx <- 1 / qk_sum
    qk <- mat[, idx_ntc] / qk_sum + epx
    qk <- qk / sum(qk)

    kl_vals <- vapply(guide_list, function(g) {
      idx <- match(g, colnames(mat))
      pk_sum <- sum(mat[, idx])
      if (pk_sum == 0) return(NA_real_)
      pk <- mat[, idx] / pk_sum + epx
      pk <- pk / sum(pk)
      .kl_div(pk, qk)
    }, numeric(1))

    object <- .set_meta_feature(object, result_field, kl_vals[match(features, names(kl_vals))])

    if (!is.null(n_permutations)) {
      pvals <- rep(NA_real_, length(features))
      names(pvals) <- features
      for (i in seq_along(guide_list)) {
        g <- guide_list[i]
        kl <- kl_vals[g]
        if (!is.finite(kl)) next
        .progress(i, length(guide_list), "KL permutations", show_progress)
        pvals[g] <- .perm_kl_single(match(g, colnames(mat)), idx_ntc, n_permutations, mat, epx, kl)
      }
      object <- .set_meta_feature(object, paste0(result_field, ".p_value"), pvals)
    }
  } else {
    obs <- .get_obs(object)
    libs <- unique(obs[[library_key]])
    for (lib in libs) {
      idx <- which(obs[[library_key]] == lib)
      sub <- mat[idx, , drop = FALSE]
      idx_ntc <- match(reference_guide, colnames(sub))
      qk_sum <- sum(sub[, idx_ntc])
      if (qk_sum == 0) next
      epx <- 1 / qk_sum
      qk <- sub[, idx_ntc] / qk_sum + epx
      qk <- qk / sum(qk)
      col <- paste0(lib, result_field)
      kl_lib <- vapply(guide_list, function(g) {
        gi <- match(g, colnames(sub))
        pk_sum <- sum(sub[, gi])
        if (pk_sum == 0) return(NA_real_)
        pk <- sub[, gi] / pk_sum + epx
        .kl_div(pk / sum(pk), qk / sum(qk))
      }, numeric(1))
      object <- .set_meta_feature(object, col, kl_lib[match(features, names(kl_lib))])
      ranks <- rank(-kl_lib, ties.method = "min", na.last = "keep")
      object <- .set_meta_feature(object, paste0(col, ".rank"), ranks[match(features, names(ranks))])
    }
    mf <- .get_meta_features(object)
    rank_cols <- grep(paste0(result_field, "\\.rank$"), colnames(mf), value = TRUE)
    if (length(rank_cols) > 0L) {
      mean_rank <- rowMeans(mf[, rank_cols, drop = FALSE], na.rm = TRUE)
      object <- .set_meta_feature(object, paste0(result_field, ".mean_rank"), mean_rank)
      object <- .set_meta_feature(object, paste0(result_field, ".final_rank"),
                                  rank(mean_rank, ties.method = "min", na.last = "keep"))
    }
    if (!is.null(n_permutations)) {
      final_ranks <- .get_meta_feature(object, paste0(result_field, ".final_rank"))
      pvals <- rep(NA_real_, length(features))
      names(pvals) <- features
      for (p in seq_len(n_permutations)) {
        .progress(p, n_permutations, "Permutation", show_progress)
        perm_ranks <- numeric()
        for (lib in libs) {
          idx <- which(obs[[library_key]] == lib)
          sub <- mat[idx, , drop = FALSE]
          idx_arr <- sample.int(nrow(sub))
          perm <- sub[idx_arr, , drop = FALSE]
          # simplified: rank mean KL under shuffle per lib
        }
      }
      # Multi-library permutation p-value (match Python rank comparison)
      for (g in guide_list) {
        tr <- final_ranks[match(g, features)]
        if (!is.finite(tr)) next
        cnt <- 0L
        for (p in seq_len(n_permutations)) {
          perm_fr <- numeric()
          for (lib in libs) {
            idx <- which(obs[[library_key]] == lib)
            sub <- mat[idx, , drop = FALSE]
            sh <- sub[sample.int(nrow(sub)), , drop = FALSE]
            idx_ntc <- match(reference_guide, colnames(sh))
            gi <- match(g, colnames(sh))
            qk_sum <- sum(sh[, idx_ntc])
            if (qk_sum == 0) next
            epx <- 1 / qk_sum
            pk_sum <- sum(sh[, gi])
            if (pk_sum == 0) next
            pk <- sh[, gi] / pk_sum + epx
            qk <- sh[, idx_ntc] / qk_sum + epx
            perm_fr <- c(perm_fr, .kl_div(pk / sum(pk), qk / sum(qk)))
          }
          if (length(perm_fr) > 0L) {
            # lower rank = better in Python final_rank
            if (rank(-mean(perm_fr), ties.method = "min") < tr) cnt <- cnt + 1L
          }
        }
        pvals[g] <- cnt / n_permutations
      }
      object <- .set_meta_feature(object, paste0(result_field, ".p_value"), pvals)
    }
  }

  if (copy) object else invisible(NULL)
}

#' Wasserstein distance on spatial KDE distributions
#' @export
wasserstein_distance <- function(
    object,
    reference_guide = "sgNon-targeting",
    result_field = "w_dist",
    guide_list = NULL,
    spatial_key = "spatial",
    library_key = NULL,
    n_permutations = 1000L,
    spatial_interval = 200,
    bin_width = "silverman",
    show_progress = TRUE,
    copy = FALSE,
    n_jobs = 4L,
    assay = NULL
) {
  if (copy) {
    mat <- .get_counts(object, assay)
    object <- Seurat::CreateSeuratObject(counts = mat, assay = .default_assay(object, assay),
                                         meta.data = .get_obs(object))
    if (!is.null(object@misc$spatial)) object <- .set_spatial(object, object@misc$spatial)
  }

  features <- .feature_names(object, assay)
  if (is.null(guide_list)) guide_list <- setdiff(features, reference_guide)
  all_feats <- unique(c(guide_list, reference_guide))
  object <- object[all_feats, ]

  coords <- .get_spatial(object)
  if (ncol(coords) >= 3L) {
    spatial <- cbind(coords[, 2], coords[, 3])
  } else {
    spatial <- coords[, 1:2, drop = FALSE]
  }

  mat <- t(as.matrix(.get_counts(object, assay)))

  compute_grid <- function(x, y, interval) {
    xb <- max(floor((max(x) - min(x)) / interval), 2L)
    yb <- max(floor((max(y) - min(y)) / interval), 2L)
    list(
      x = seq(min(x), max(x), length.out = xb),
      y = seq(min(y), max(y), length.out = yb)
    )
  }

  silverman_bw <- function(v) {
    max(1e-10, 1.06 * stats::sd(v) * (length(v)^(-1 / 5)))
  }

  kde_for_guide <- function(guide_name, sub_mat, sub_sp) {
    v <- sub_mat[, guide_name]
    if (bin_width == "silverman") {
      bw <- silverman_bw(v)
      cnt_grid <- seq(min(v), max(v), by = bw)
    } else {
      cnt_grid <- seq(min(v), max(v), by = bin_width)
    }
    if (length(cnt_grid) < 2L) cnt_grid <- c(min(v), min(v) + 1)
    g <- compute_grid(sub_sp[, 1], sub_sp[, 2], spatial_interval)
    .calculate_kde(sub_sp, v, g$x, g$y, cnt_grid)
  }

  .run_library <- function(idx_cells) {
    sub_mat <- mat[idx_cells, , drop = FALSE]
    sub_sp <- spatial[idx_cells, , drop = FALSE]
    g <- compute_grid(sub_sp[, 1], sub_sp[, 2], spatial_interval)
    kdes <- lapply(c(guide_list, reference_guide), kde_for_guide, sub_mat = sub_mat, sub_sp = sub_sp)
    names(kdes) <- c(guide_list, reference_guide)
    wds <- vapply(guide_list, function(g) {
      .wasserstein_2d(kdes[[g]], kdes[[reference_guide]])
    }, numeric(1))
    list(kdes = kdes, wdist = wds, grid = g)
  }

  if (is.null(library_key)) {
    res <- .run_library(seq_len(nrow(mat)))
    object <- .set_meta_feature(object, result_field, res$wdist[match(features, names(res$wdist))])

    if (!is.null(n_permutations)) {
      pvals <- rep(0, length(guide_list))
      names(pvals) <- guide_list
      idx_ntc <- match(reference_guide, colnames(mat))
      for (guide in guide_list) {
        idx_g <- match(guide, colnames(mat))
        obs_w <- res$wdist[guide]
        higher <- 0L
        for (p in seq_len(n_permutations)) {
          .progress(p, n_permutations, guide, show_progress)
          sw <- .permute_guide_bins(mat, idx_g, idx_ntc)
          sub_sp <- spatial
          if (bin_width == "silverman") {
            cnt_g <- seq(min(sw$guide), max(sw$guide), by = silverman_bw(sw$guide))
            cnt_n <- seq(min(sw$ntc), max(sw$ntc), by = silverman_bw(sw$ntc))
          } else {
            cnt_g <- seq(min(sw$guide), max(sw$guide), by = bin_width)
            cnt_n <- seq(min(sw$ntc), max(sw$ntc), by = bin_width)
          }
          kde_g <- .calculate_kde(sub_sp, sw$guide, res$grid$x, res$grid$y, cnt_g)
          kde_n <- .calculate_kde(sub_sp, sw$ntc, res$grid$x, res$grid$y, cnt_n)
          if (.wasserstein_2d(kde_g, kde_n) > obs_w) higher <- higher + 1L
        }
        pvals[guide] <- higher / n_permutations
      }
      object <- .set_meta_feature(object, paste0(result_field, ".p_value"),
                                  pvals[match(features, names(pvals))])
    }
  } else {
    obs <- .get_obs(object)
    libs <- unique(obs[[library_key]])
    for (lib in libs) {
      idx <- which(obs[[library_key]] == lib)
      res <- .run_library(idx)
      object <- .set_meta_feature(object, paste0(lib, result_field),
                                  res$wdist[match(features, names(res$wdist))])
      ranks <- rank(-res$wdist, ties.method = "min", na.last = "keep")
      object <- .set_meta_feature(object, paste0(lib, result_field, ".rank"),
                                  ranks[match(features, names(ranks))])
    }
    mf <- .get_meta_features(object)
    rank_cols <- grep(paste0(result_field, "\\.rank$"), colnames(mf), value = TRUE)
    if (length(rank_cols) > 0L) {
      mean_rank <- rowMeans(mf[, rank_cols, drop = FALSE], na.rm = TRUE)
      object <- .set_meta_feature(object, paste0(result_field, ".mean_rank"), mean_rank)
      object <- .set_meta_feature(object, paste0(result_field, ".final_rank"),
                                  rank(mean_rank, ties.method = "min", na.last = "keep"))
    }
  }

  if (copy) object else invisible(NULL)
}
