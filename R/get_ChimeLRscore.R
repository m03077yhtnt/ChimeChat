#' Step 2: Build a ligand–receptor (LR) table and compute LR scores
#'
#' @description
#' Select LR interactions from CellChatDB and compute LR scores as the sum of
#' VST-transformed ligand and receptor values for matched samples.
#'
#' @param combined A numeric matrix (rows = human gene symbols) returned by \code{run_ChimeDB()}.
#' @param direction Direction of scoring: \code{"mouse_ligand_human_receptor"} or
#'   \code{"human_ligand_mouse_receptor"}.
#' @param db_scope Character vector of any of \code{c("secreted","contact","ecm")} to select
#'   CellChatDB categories.
#' @param rename_to_mouse Logical; if \code{TRUE}, convert the mouse-side gene names for display.
#' @param ortholog_map A list or data.frame used to map human-to-mouse names (typically
#'   \code{run_ChimeDB()$ortholog_map}).
#' @param rename_mode Display mode when multiple orthologs exist: \code{"one"} (first only) or
#'   \code{"compact"} (first plus “+n more”).
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{direction}{Scoring direction.}
#'   \item{ligand}{Ligand gene symbol.}
#'   \item{receptor}{Receptor gene symbol.}
#'   \item{score}{Sum of VST-transformed ligand and receptor values.}
#'   \item{condition}{Condition label (Control/Treatment).}
#'   \item{sample}{Sample identifier.}
#'   \item{pair_human}{Pair shown with human symbols on both sides.}
#'   \item{pair}{Display pair, with one side renamed to mouse symbols if requested.}
#' }
#' @export
get_ChimeLRscore <- function(combined,
                             direction = c("mouse_ligand_human_receptor",
                                           "human_ligand_mouse_receptor"),
                             db_scope = c("secreted", "contact", "ecm"),
                             rename_to_mouse = FALSE,
                             ortholog_map = NULL,
                             rename_mode = c("compact", "one")) {
  if (!requireNamespace("CellChat", quietly = TRUE)) stop("Package 'CellChat' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE))    stop("Package 'dplyr' is required.")
  if (!requireNamespace("tibble", quietly = TRUE))   stop("Package 'tibble' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE))    stop("Package 'tidyr' is required.")

  direction   <- match.arg(direction)
  rename_mode <- match.arg(rename_mode)

  CellChatDB <- CellChat::CellChatDB.human
  db_list <- list()
  if ("secreted" %in% db_scope) db_list$secreted <- CellChat::subsetDB(CellChatDB, search = "Secreted Signaling")
  if ("contact"  %in% db_scope) db_list$contact  <- CellChat::subsetDB(CellChatDB, search = "Cell-Cell Contact")
  if ("ecm"      %in% db_scope) db_list$ecm      <- CellChat::subsetDB(CellChatDB, search = "ECM-Receptor")

  LRdb <- dplyr::bind_rows(lapply(db_list, function(db) {
    as.data.frame(db$interaction[, c("ligand", "receptor")])
  }))
  LRdb <- dplyr::distinct(LRdb)

  mouse_cols <- c("Control_Mouse_1","Control_Mouse_2","Control_Mouse_3",
                  "Treatment_Mouse_1","Treatment_Mouse_2","Treatment_Mouse_3")
  human_cols <- c("Control_Human_1","Control_Human_2","Control_Human_3",
                  "Treatment_Human_1","Treatment_Human_2","Treatment_Human_3")

  if (direction == "mouse_ligand_human_receptor") {
    ligand_cols   <- mouse_cols
    receptor_cols <- human_cols
  } else {
    ligand_cols   <- human_cols
    receptor_cols <- mouse_cols
  }
  conditions <- c(rep("Control", 3), rep("Treatment", 3))

  has_gene <- function(g) !is.na(g) && nzchar(g) && g %in% rownames(combined)

  out_list <- vector("list", length(ligand_cols))
  for (i in seq_along(ligand_cols)) {
    lig_col <- ligand_cols[i]
    rec_col <- receptor_cols[i]
    cond    <- conditions[i]
    cond_ix <- sum(conditions[1:i] == cond)

    scores <- vapply(seq_len(nrow(LRdb)), function(j) {
      lig <- LRdb$ligand[j]; rec <- LRdb$receptor[j]
      if (has_gene(lig) && has_gene(rec)) {
        vl <- combined[lig, lig_col]; vr <- combined[rec, rec_col]
        if (is.na(vl) || is.na(vr)) NA_real_ else (vl + vr)
      } else NA_real_
    }, numeric(1))

    out_list[[i]] <- tibble::tibble(
      ligand    = LRdb$ligand,
      receptor  = LRdb$receptor,
      score     = scores,
      condition = cond,
      sample    = paste0(cond, "_", cond_ix)
    )
  }

  df <- dplyr::bind_rows(out_list)
  df <- tidyr::drop_na(df, score)
  df <- dplyr::mutate(df, pair_human = paste(.data$ligand, .data$receptor, sep = "-"))

  make_mouse_label <- function(hsym_vec, h2m_all, mode = c("compact", "one")) {
    mode <- match.arg(mode)
    sapply(hsym_vec, function(h) {
      ms <- unique(h2m_all[[h]])
      ms <- ms[!is.na(ms) & nzchar(ms)]
      if (length(ms) == 0) return(NA_character_)
      if (mode == "one") return(ms[1])
      if (length(ms) == 1) return(ms[1])
      paste0(ms[1], " +", length(ms) - 1, " more")
    })
  }

  get_h2m_all <- function(omap) {
    if (is.list(omap) && !is.null(omap$h2m_all)) return(omap$h2m_all)
    if (is.data.frame(omap) &&
        all(c("mouse_symbol","human_symbol") %in% colnames(omap))) {
      return(split(omap$mouse_symbol, omap$human_symbol))
    }
    NULL
  }
  h2m_all <- get_h2m_all(ortholog_map)

  if (isTRUE(rename_to_mouse) && !is.null(h2m_all)) {
    if (direction == "mouse_ligand_human_receptor") {
      lab <- make_mouse_label(df$ligand, h2m_all, rename_mode)
      df  <- dplyr::mutate(df,
                           ligand = ifelse(!is.na(lab), lab, .data$ligand),
                           pair   = paste(.data$ligand, .data$receptor, sep = "-"))
    } else {
      lab <- make_mouse_label(df$receptor, h2m_all, rename_mode)
      df  <- dplyr::mutate(df,
                           receptor = ifelse(!is.na(lab), lab, .data$receptor),
                           pair     = paste(.data$ligand, .data$receptor, sep = "-"))
    }
  } else {
    df <- dplyr::mutate(df, pair = .data$pair_human)
  }

  df <- dplyr::mutate(df, direction = direction, .before = 1)
  df
}
