#' Step 4: Heatmap of top ligand–receptor (LR) pairs using pheatmap
#'
#' Selects top LR pairs from statistical results and plots a pair (rows) × sample (columns)
#' heatmap of LR scores.
#'
#' Pair selection logic:
#'   1) If any pairs meet the FDR threshold, take the top \code{n_top} by ascending \code{p.adj}
#'      (ties broken by |t|).
#'   2) Else if any pairs have finite \code{p.adj}, take the top \code{n_top} by ascending \code{p.adj}.
#'   3) Else fall back to top \code{n_top} by descending |delta|.
#'
#' @param lr_scores Data frame from \code{get_ChimeLRscore()}.
#' @param stat_results Data frame from \code{get_ChimeStatistics()} (must contain
#'   \code{pair}, \code{p.adj}, \code{t}, and \code{delta}).
#' @param n_top Integer. Number of top pairs to plot (default 20).
#' @param fdr_thr Numeric. FDR threshold used to prefer significant pairs (default 0.05).
#' @param pair_column Which LR-pair column in \code{lr_scores} to use: \code{"pair"} (display)
#'   or \code{"pair_human"} (both sides in human symbols).
#' @param sample_order Character vector. Desired column order in the heatmap. Columns not
#'   listed will be appended after the specified order. Default:
#'   \code{c(paste0("Control_", 1:3), paste0("Treatment_", 1:3))}.
#' @param cluster_rows Logical. Cluster rows? (default TRUE)
#' @param cluster_cols Logical. Cluster columns? (default FALSE)
#' @param scale_rows Logical. If TRUE, row-scale scores (i.e., \code{scale="row"}) (default TRUE).
#' @param palette Brewer palette name (default \code{"RdBu"}).
#' @param palette_n Integer. Number of colors to interpolate (default 100).
#' @param palette_rev Logical. Reverse palette? (default TRUE).
#' @param breaks Optional numeric vector of breaks passed to \code{pheatmap}. If NULL,
#'   symmetric breaks are computed; when \code{scale_rows=TRUE}, a centered scheme in \code{[-2, 2]}
#'   is used by default.
#' @param fontsize_row Numeric. Row font size (default 12).
#' @param fontsize_col Numeric. Column font size (default 12).
#' @param border_color Border color (default \code{"white"}).
#' @param show_rownames Logical. Show row names? (default TRUE).
#' @param ... Additional arguments passed to \code{pheatmap::pheatmap()}.
#'
#' @return A list with three elements:
#'   \itemize{
#'     \item \code{pairs} — Character vector of plotted pairs (in the order used).
#'     \item \code{matrix} — Numeric matrix (pairs × samples) used for plotting.
#'     \item \code{plot} — The \code{pheatmap} object (invisibly returned by pheatmap).
#'   }
#' @export
run_ChimeHeatmap <- function(lr_scores,
                             stat_results,
                             n_top = 20,
                             fdr_thr = 0.05,
                             pair_column = c("pair", "pair_human"),
                             sample_order = c(paste0("Control_", 1:3),
                                              paste0("Treatment_", 1:3)),
                             cluster_rows = TRUE,
                             cluster_cols = FALSE,
                             scale_rows   = TRUE,
                             palette      = "RdBu",
                             palette_n    = 100,
                             palette_rev  = TRUE,
                             breaks = NULL,
                             fontsize_row = 12,
                             fontsize_col = 12,
                             border_color = "white",
                             show_rownames = TRUE,
                             ...) {

  if (!requireNamespace("pheatmap", quietly = TRUE))   stop("Package 'pheatmap' is required.")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Package 'RColorBrewer' is required.")

  pair_column <- match.arg(pair_column)

  if (!all(c("pair","p.adj","t","delta") %in% colnames(stat_results))) {
    stop("stat_results must contain columns: pair, p.adj, t, delta (use get_ChimeStatistics()).")
  }
  sr <- stat_results

  if (any(is.finite(sr$p.adj) & sr$p.adj < fdr_thr, na.rm = TRUE)) {
    sig_pairs <- sr |>
      dplyr::filter(is.finite(.data$p.adj), .data$p.adj < fdr_thr) |>
      dplyr::arrange(.data$p.adj, dplyr::desc(abs(.data$t))) |>
      dplyr::slice_head(n = n_top) |>
      dplyr::pull(.data$pair) |>
      unique()
  } else if (any(is.finite(sr$p.adj), na.rm = TRUE)) {
    sig_pairs <- sr |>
      dplyr::filter(is.finite(.data$p.adj)) |>
      dplyr::arrange(.data$p.adj, dplyr::desc(abs(.data$t))) |>
      dplyr::slice_head(n = n_top) |>
      dplyr::pull(.data$pair) |>
      unique()
  } else {
    sig_pairs <- sr |>
      dplyr::arrange(dplyr::desc(abs(.data$delta))) |>
      dplyr::slice_head(n = n_top) |>
      dplyr::pull(.data$pair) |>
      unique()
  }

  if (length(sig_pairs) == 0L) stop("No pairs selected for heatmap.")

  if (!all(c("sample","score") %in% colnames(lr_scores))) {
    stop("lr_scores must contain 'sample' and 'score' columns.")
  }
  if (!pair_column %in% colnames(lr_scores)) {
    stop(sprintf("Column '%s' not found in lr_scores.", pair_column))
  }

  mat_df <- lr_scores |>
    dplyr::filter(!is.na(.data[[pair_column]]),
                  .data[[pair_column]] %in% sig_pairs) |>
    dplyr::select(pair = dplyr::all_of(pair_column),
                  sample, score) |>
    dplyr::group_by(.data$pair, .data$sample) |>
    dplyr::summarise(score = mean(.data$score, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = .data$sample, values_from = .data$score)

  if (!"pair" %in% colnames(mat_df)) stop("Failed to build pair x sample matrix.")

  rownames_mat <- mat_df$pair
  sig_scores_mat <- as.data.frame(mat_df[, setdiff(colnames(mat_df), "pair"), drop = FALSE])
  rownames(sig_scores_mat) <- rownames_mat
  sig_scores_mat <- as.matrix(sig_scores_mat)

  keep_cols  <- intersect(sample_order, colnames(sig_scores_mat))
  other_cols <- setdiff(colnames(sig_scores_mat), keep_cols)
  sig_scores_mat <- sig_scores_mat[, c(keep_cols, other_cols), drop = FALSE]

  cols <- RColorBrewer::brewer.pal(11, palette)
  if (palette_rev) cols <- rev(cols)
  colors.use <- grDevices::colorRampPalette(cols)(palette_n)

  if (is.null(breaks)) {
    if (scale_rows) {
      n_half_hi <- floor(palette_n / 2)
      n_half_lo <- ceiling(palette_n / 2) + 1
      breaks <- c(seq(-2, 0, length.out = n_half_lo),
                  seq(0.05, 2, length.out = n_half_hi))
    } else {
      vmax <- max(abs(sig_scores_mat), na.rm = TRUE)
      if (!is.finite(vmax) || vmax <= 0) vmax <- 1
      breaks <- seq(-vmax, vmax, length.out = palette_n + 1)
    }
  }

  p <- pheatmap::pheatmap(
    sig_scores_mat,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    scale        = if (scale_rows) "row" else "none",
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    color        = colors.use,
    border_color = border_color,
    breaks       = breaks,
    show_rownames = show_rownames,
    ...
  )

  invisible(list(
    pairs   = sig_pairs,
    matrix  = sig_scores_mat,
    plot    = p
  ))
}
