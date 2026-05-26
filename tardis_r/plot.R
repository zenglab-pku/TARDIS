#' Plotting utilities (ggplot2, Python/matplotlib style)
#' @import ggplot2
NULL

#' @keywords internal
.theme_tardis <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.3),
      legend.position = "right"
    )
}

#' KDE plot of ranking metric with NTC and top guides highlighted
#' @export
plot_top_kde <- function(
    object,
    result_field,
    show_sgnt = TRUE,
    top_n = 2L,
    violin = TRUE,
    sgnt_label = "sgNon-targeting",
    figsize = c(5, 2.7),
    min_count = 0L,
    assay = NULL
) {
  mf <- .get_meta_features(object, assay)
  if (!(result_field %in% colnames(mf))) {
    stop(result_field, " not found in feature metadata.")
  }
  var <- mf
  if ("TotalCount" %in% colnames(var)) {
    var <- var[var$TotalCount > min_count, , drop = FALSE]
  }
  vals <- var[[result_field]]
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) stop("No valid values for ", result_field)

  dens <- stats::density(vals)
  ddf <- data.frame(x = dens$x, y = dens$y)

  p <- ggplot2::ggplot(ddf, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_area(fill = "gray70", alpha = 1, colour = NA) +
    ggplot2::geom_line(colour = "gray40") +
    ggplot2::labs(x = "", y = "Density") +
    .theme_tardis() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())

  if (show_sgnt && sgnt_label %in% rownames(var)) {
    sg <- var[sgnt_label, result_field]
    yv <- approx(dens$x, dens$y, xout = sg)$y
    p <- p + ggplot2::geom_vline(xintercept = sg, colour = "red", linewidth = 1) +
      ggplot2::annotate("text", x = sg, y = yv, label = sgnt_label, colour = "red",
                        hjust = 0, vjust = 0, size = 3, fontface = "bold")
  }

  top_genes <- setdiff(rownames(var), sgnt_label)
  top_genes <- top_genes[order(-var[top_genes, result_field], na.last = NA)][seq_len(min(top_n, length(top_genes)))]
  cols <- grDevices::hcl.colors(length(top_genes), palette = "Set3")
  for (i in seq_along(top_genes)) {
    g <- top_genes[i]
    val <- var[g, result_field]
    yv <- approx(dens$x, dens$y, xout = val)$y
    p <- p + ggplot2::geom_vline(xintercept = val, colour = cols[i], linewidth = 0.8, linetype = "dashed") +
      ggplot2::annotate("text", x = val, y = yv, label = g, colour = cols[i], hjust = 0, vjust = 0, size = 2.8)
  }

  if (violin) {
    vdf <- data.frame(val = var[[result_field]])
    p2 <- ggplot2::ggplot(vdf, ggplot2::aes(x = .data$val, y = 1)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = 0, ymax = 1), alpha = 0.3, linewidth = 1, colour = "gray50") +
      ggplot2::labs(x = gsub("_", " ", result_field, fixed = TRUE), y = "") +
      ggplot2::theme_void() +
      ggplot2::theme(panel.background = ggplot2::element_rect(fill = "whitesmoke", colour = NA))
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- patchwork::wrap_plots(p, p2, ncol = 1, heights = c(3, 1))
    }
  }

  p <- p + ggplot2::ggtitle(NULL)
  print(p)
  invisible(p)
}

#' Spatial scatter of dominant guide per bin
#' @export
plot_spatial_guides <- function(
    object,
    scale_factor,
    image = NULL,
    figsize = c(10, 10),
    palette = NULL,
    s = 3,
    edgecolor = NA,
    legend = FALSE,
    alpha = 1,
    assay = NULL
) {
  mat <- t(as.matrix(.get_counts(object, assay)))
  rs <- rowSums(mat)
  keep <- rs > 0
  mat <- mat[keep, , drop = FALSE]
  coords <- .get_spatial(object)[keep, , drop = FALSE]
  if (ncol(coords) >= 3L) {
    px <- coords[, 3] * scale_factor
    py <- coords[, 2] * scale_factor
  } else {
    px <- coords[, 2] * scale_factor
    py <- coords[, 1] * scale_factor
  }
  guide <- colnames(mat)[max.col(mat, ties.method = "first")]

  plot_df <- data.frame(x = px, y = py, guide = guide)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data$y, colour = .data$guide)) +
    ggplot2::geom_point(size = s, alpha = alpha, stroke = 0) +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed() +
    .theme_tardis() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())

  if (!is.null(palette)) {
    if (is.character(palette) && length(palette) == 1L) {
      n <- length(unique(guide))
      pal <- grDevices::hcl.colors(n, palette)
      p <- p + ggplot2::scale_colour_manual(values = pal)
    } else {
      p <- p + ggplot2::scale_colour_manual(values = palette)
    }
  }
  if (!legend) p <- p + ggplot2::theme(legend.position = "none")

  print(p)
  invisible(p)
}

#' Rank vs distance scatter coloured by -log10(p)
#' @export
plot_ranking_scatter <- function(object, NTC, result_field = "dist", assay = NULL) {
  mf <- .get_meta_features(object, assay)
  pcol <- paste0(result_field, ".p_value")
  if (!(result_field %in% colnames(mf)) || !(pcol %in% colnames(mf))) {
    stop("Need ", result_field, " and ", pcol, " in feature metadata.")
  }
  top_clones <- setdiff(rownames(mf), NTC)
  plot_data <- mf[top_clones, c(result_field, pcol), drop = FALSE]
  plot_data <- plot_data[complete.cases(plot_data), , drop = FALSE]
  plot_data$rank <- rank(-plot_data[[result_field]], ties.method = "min")
  plot_data$neg_log10_p <- -log10(plot_data[[pcol]] + 1e-10)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data[[result_field]],
    y = .data$rank,
    colour = .data$neg_log10_p
  )) +
    ggplot2::geom_point(size = 1.5, alpha = 0.8) +
    ggplot2::scale_colour_viridis_c(option = "viridis", trans = "log10") +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(x = result_field, y = paste0(result_field, ".rank"), colour = "-log10(p)") +
    .theme_tardis()

  print(p)
  invisible(p)
}
