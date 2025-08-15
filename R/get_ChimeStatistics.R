#' Step 3: Perform a per-pair t-test (Treatment vs Control) with multiple-testing correction
#'
#' @param lr_scores A data frame from get_ChimeLRscore(), containing at least:
#'   \itemize{
#'     \item \code{pair}: Ligand–receptor pair identifier.
#'     \item \code{score}: Numeric score for each sample.
#'     \item \code{condition}: Condition label ("Control" or "Treatment").
#'     \item \code{direction} (optional): Interaction direction, if applicable.
#'   }
#' @param adjust_method Method for \code{p.adjust} (default "fdr").
#' @param by_direction If TRUE and \code{lr_scores} contains \code{direction},
#'   adjust p-values separately within each direction.
#'
#' @return A tibble with one row per pair (and per direction if present), including:
#'   \itemize{
#'     \item \code{direction}: Interaction direction (if present).
#'     \item \code{pair}: Ligand–receptor pair identifier.
#'     \item \code{p.value}: Raw p-value from the t-test.
#'     \item \code{t}: t-statistic from the t-test.
#'     \item \code{mean_Treatment}: Mean score in Treatment samples.
#'     \item \code{mean_Control}: Mean score in Control samples.
#'     \item \code{delta}: Difference \code{mean_Treatment - mean_Control}.
#'     \item \code{n_Treatment}: Number of Treatment samples.
#'     \item \code{n_Control}: Number of Control samples.
#'     \item \code{p.adj}: Adjusted p-value after multiple-testing correction.
#'   }
#' @export
get_ChimeStatistics <- function(lr_scores,
                                adjust_method = "fdr",
                                by_direction = TRUE) {
  if (!requireNamespace("dplyr", quietly = TRUE))  stop("Package 'dplyr' is required.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package 'tibble' is required.")

  need_cols <- c("pair", "score", "condition")
  miss <- setdiff(need_cols, colnames(lr_scores))
  if (length(miss)) stop("lr_scores is missing required columns: ", paste(miss, collapse = ", "))

  cond_levels <- unique(lr_scores$condition)
  if (!all(c("Control","Treatment") %in% cond_levels)) {
    stop("lr_scores$condition must contain both 'Control' and 'Treatment'.")
  }

  grp_keys <- if (by_direction && "direction" %in% colnames(lr_scores)) c("direction", "pair") else "pair"

  stats <- dplyr::group_by(lr_scores, dplyr::across(dplyr::all_of(grp_keys)))
  stats <- dplyr::summarise(
    stats,
    {
      x_t <- score[condition == "Treatment"]
      x_c <- score[condition == "Control"]
      n_Treatment <- sum(is.finite(x_t))
      n_Control   <- sum(is.finite(x_c))
      tt <- tryCatch(stats::t.test(x_t, x_c), error = function(e) NULL)
      p.value <- ifelse(is.null(tt), NA_real_, tt$p.value)
      t_stat  <- ifelse(is.null(tt), NA_real_, unname(tt$statistic))
      mean_Treatment <- mean(x_t, na.rm = TRUE)
      mean_Control   <- mean(x_c, na.rm = TRUE)
      delta <- mean_Treatment - mean_Control
      tibble::tibble(p.value, t = t_stat,
                     mean_Treatment, mean_Control, delta,
                     n_Treatment, n_Control)
    },
    .groups = "drop_last"
  )
  stats <- dplyr::ungroup(stats)

  if ("direction" %in% colnames(stats) && by_direction) {
    stats <- dplyr::group_by(stats, .data$direction)
    stats <- dplyr::mutate(stats, p.adj = stats::p.adjust(.data$p.value, method = adjust_method))
    stats <- dplyr::ungroup(stats)
  } else {
    stats <- dplyr::mutate(stats, p.adj = stats::p.adjust(.data$p.value, method = adjust_method))
  }

  stats <- dplyr::arrange(stats, .data$p.adj, dplyr::desc(abs(.data$t)))
  stats
}
