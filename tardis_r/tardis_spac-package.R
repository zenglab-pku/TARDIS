#' @keywords internal
"_PACKAGE"

#' tardisSpac: TARDIS spatial perturbation analysis (R port)
#'
#' R port of the Python \code{tardis_spac} package for spatial CRISPR screens.
#' Uses Seurat objects (read 10x HDF5 via \code{\link{read_guide_h5}}) as the
#' primary data structure, matching AnnData layout in Python:
#' \itemize{
#'   \item Cell/bin metadata -> \code{object@meta.data}
#'   \item Guide-level results -> \code{object[[assay]]@meta.features}
#'   \item Spatial coordinates -> \code{object@misc$spatial}
#' }
#'
#' @section Statistics:
#' \code{\link{kl_divergence}}, \code{\link{wasserstein_distance}},
#' \code{\link{aitchison_distance}}, \code{\link{permanova}}
#'
#' @section Preprocessing:
#' \code{\link{filter_guide_reads}}, \code{\link{filter_guide_reads_h5}},
#' \code{\link{combine_guide_replicates}}, \code{\link{load_bin}}, etc.
#'
#' @name tardis_spac-package
#' @aliases tardis_spac tardisSpac
NULL

#' @export
#' @rdname tardis_spac-package
stats <- function() {
  list(
    kl_divergence = kl_divergence,
    wasserstein_distance = wasserstein_distance,
    aitchison_distance = aitchison_distance,
    permanova = permanova
  )
}

#' @export
#' @rdname tardis_spac-package
utils <- function() {
  list(
    qc_guide_bins = qc_guide_bins,
    remove_mito_ribo_hk_lnc_genes = remove_mito_ribo_hk_lnc_genes,
    filter_guide_reads = filter_guide_reads,
    filter_guide_reads_h5 = filter_guide_reads_h5,
    combine_guide_replicates = combine_guide_replicates,
    nmf_clustering = nmf_clustering,
    nmf_consensus = nmf_consensus,
    dbscan_density_region = dbscan_density_region,
    plot_guide_gene_summary = plot_guide_gene_summary,
    load_bin = load_bin,
    read_guide_h5 = read_guide_h5,
    plot_spatial_guides = plot_spatial_guides,
    plot_top_kde = plot_top_kde,
    plot_ranking_scatter = plot_ranking_scatter
  )
}
