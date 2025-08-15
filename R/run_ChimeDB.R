#' Step 1: Map mouse -> human at raw-count level, VST per species, then merge
#'
#' Reads two raw count tables (columns: Control_1..3, Treatment_1..3, symbol).
#' First, mouse gene symbols are mapped to representative human symbols
#' (one2one preferred, otherwise the highest %id) and mouse counts are
#' aggregated per human symbol **at the raw-count level**. Then DESeq2 size
#' factors and VST are computed separately for the (mapped) mouse matrix and
#' the human matrix. Finally, the union of genes is taken and matrices are merged.
#'
#' @param human_file Path to human raw counts (must include columns Control_1..3, Treatment_1..3, symbol).
#' @param mouse_file Path to mouse raw counts (same columns as human).
#' @param host Ensembl archive host for biomaRt.
#' @return list(combined, ortholog_map) where `combined` is a numeric matrix with columns
#'   Control_Mouse_1..3, Treatment_Mouse_1..3, Control_Human_1..3, Treatment_Human_1..3;
#'   `ortholog_map` is a list with elements: m2h_full (data.frame), h2m_all (list),
#'   h2m_rep (named character), m2h_rep (named character).
#' @examples
#' # ChimeDB <- run_ChimeDB("Human_CountData.xlsx","Mouse_CountData.xlsx")
#' @export
run_ChimeDB <- function(human_file,
                        mouse_file,
                        host = "https://dec2021.archive.ensembl.org/") {

  for (p in c("DESeq2","biomaRt","readxl","dplyr","tibble","SummarizedExperiment")) {
    if (!requireNamespace(p, quietly = TRUE)) stop(sprintf("Package '%s' is required.", p))
  }

  .read_counts <- function(path) {
    df  <- readxl::read_xlsx(path)
    req <- c("Control_1","Control_2","Control_3",
             "Treatment_1","Treatment_2","Treatment_3","symbol")
    if (!all(req %in% colnames(df))) {
      stop("Input must contain columns: Control_1..3, Treatment_1..3, symbol.")
    }
    df <- dplyr::filter(df, !is.na(.data$symbol) & nzchar(.data$symbol))

    mat <- as.matrix(df[, req[1:6]])
    rownames(mat) <- df$symbol

    if (any(duplicated(rownames(mat)))) {
      tmp <- tibble::rownames_to_column(as.data.frame(mat), "symbol")
      tmp <- dplyr::group_by(tmp, .data$symbol)
      tmp <- dplyr::summarise(
        tmp,
        dplyr::across(c("Control_1","Control_2","Control_3",
                        "Treatment_1","Treatment_2","Treatment_3"),
                      ~ sum(.x, na.rm = TRUE)),
        .groups = "drop"
      )
      mat <- as.matrix(tibble::column_to_rownames(tmp, "symbol"))
    }
    mat <- mat[rowSums(mat) > 0.5, , drop = FALSE]
    mat
  }

  .reindex_to <- function(mat, genes) {
    out <- matrix(NA_real_, nrow = length(genes), ncol = ncol(mat),
                  dimnames = list(genes, colnames(mat)))
    common <- intersect(genes, rownames(mat))
    if (length(common)) out[common, ] <- mat[common, , drop = FALSE]
    out
  }

  human_raw <- .read_counts(human_file)
  mouse_raw <- .read_counts(mouse_file)

  human_ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host = host)
  mouse_ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)

  mm_syms <- unique(rownames(mouse_raw))

  fetch_m2h <- function(mm_syms) {
    res <- tryCatch({
      df <- biomaRt::getBM(
        attributes = c("mgi_symbol",
                       "hsapiens_homolog_associated_gene_name",
                       "hsapiens_homolog_ensembl_gene",
                       "hsapiens_homolog_orthology_type",
                       "hsapiens_homolog_perc_id",
                       "hsapiens_homolog_perc_id_r1",
                       "hsapiens_homolog_orthology_confidence"),
        filters    = "mgi_symbol",
        values     = mm_syms,
        mart       = mouse_ensembl
      )
      colnames(df) <- c("mouse_symbol","human_symbol","human_ensembl",
                        "orthology_type","perc_id","perc_id_r1","orthology_confidence")
      df
    }, error = function(e) NULL)

    if (!is.null(res)) return(res)

    basic <- biomaRt::getLDS(
      attributes  = c("mgi_symbol"),
      filters     = "mgi_symbol",
      values      = mm_syms,
      mart        = mouse_ensembl,
      attributesL = c("hgnc_symbol"),
      martL       = human_ensembl,
      uniqueRows  = TRUE
    )
    colnames(basic) <- c("mouse_symbol","human_symbol")
    basic$human_ensembl         <- NA_character_
    basic$orthology_type        <- NA_character_
    basic$perc_id               <- NA_real_
    basic$perc_id_r1            <- NA_real_
    basic$orthology_confidence  <- NA_real_
    basic
  }

  m2h_full <- fetch_m2h(mm_syms)
  m2h_full <- dplyr::filter(
    m2h_full,
    !is.na(.data$human_symbol) & nzchar(.data$human_symbol),
    !is.na(.data$mouse_symbol) & nzchar(.data$mouse_symbol)
  )
  m2h_full <- dplyr::distinct(m2h_full)

  type_rank <- function(x) {
    pri <- c("ortholog_one2one","ortholog_one2many","ortholog_many2many")
    r   <- match(x, pri); r[is.na(r)] <- length(pri) + 1; r
  }

  m2h_ranked <- dplyr::mutate(
    m2h_full,
    .type_rank = type_rank(.data$orthology_type),
    .pid       = dplyr::coalesce(.data$perc_id, .data$perc_id_r1, 0)
  )

  h2m_rep_tbl <- m2h_ranked |>
    dplyr::group_by(.data$human_symbol) |>
    dplyr::arrange(.data$.type_rank, dplyr::desc(.data$.pid), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(human_symbol, mouse_symbol)
  h2m_rep <- stats::setNames(h2m_rep_tbl$mouse_symbol, h2m_rep_tbl$human_symbol)

  h2m_all <- split(m2h_ranked$mouse_symbol, m2h_ranked$human_symbol)

  m2h_rep_tbl <- m2h_ranked |>
    dplyr::group_by(.data$mouse_symbol) |>
    dplyr::arrange(.data$.type_rank, dplyr::desc(.data$.pid), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(mouse_symbol, human_symbol)
  m2h_rep <- stats::setNames(m2h_rep_tbl$human_symbol, m2h_rep_tbl$mouse_symbol)

  mouse_df <- tibble::rownames_to_column(as.data.frame(mouse_raw), "mouse_symbol")
  mouse_df <- dplyr::inner_join(mouse_df, m2h_rep_tbl, by = "mouse_symbol")
  mouse_df <- dplyr::group_by(mouse_df, .data$human_symbol)
  mouse_df <- dplyr::summarise(
    mouse_df,
    dplyr::across(c("Control_1","Control_2","Control_3",
                    "Treatment_1","Treatment_2","Treatment_3"),
                  ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )
  mouse_df <- tibble::column_to_rownames(mouse_df, "human_symbol")
  mouse_df <- as.matrix(mouse_df)

  cond <- factor(rep(c("Control","Treatment"), each = 3), levels = c("Control","Treatment"))

  human_dds <- DESeq2::DESeqDataSetFromMatrix(round(human_raw),
                                              colData = data.frame(condition = cond),
                                              design  = ~ condition)
  human_dds <- DESeq2::estimateSizeFactors(human_dds)
  human_vsd <- SummarizedExperiment::assay(DESeq2::vst(human_dds, blind = FALSE))

  mouse_dds <- DESeq2::DESeqDataSetFromMatrix(round(mouse_df),
                                              colData = data.frame(condition = cond),
                                              design  = ~ condition)
  mouse_dds <- DESeq2::estimateSizeFactors(mouse_dds)
  mouse_vsd <- SummarizedExperiment::assay(DESeq2::vst(mouse_dds, blind = FALSE))

  all_genes  <- sort(unique(c(rownames(human_vsd), rownames(mouse_vsd))))
  all_genes  <- all_genes[!is.na(all_genes) & nzchar(all_genes)]
  human_vsd2 <- .reindex_to(human_vsd, all_genes)
  mouse_vsd2 <- .reindex_to(mouse_vsd, all_genes)

  colnames(mouse_vsd2) <- c("Control_Mouse_1","Control_Mouse_2","Control_Mouse_3",
                            "Treatment_Mouse_1","Treatment_Mouse_2","Treatment_Mouse_3")
  colnames(human_vsd2) <- c("Control_Human_1","Control_Human_2","Control_Human_3",
                            "Treatment_Human_1","Treatment_Human_2","Treatment_Human_3")

  combined <- cbind(mouse_vsd2, human_vsd2)

  ortholog_map <- list(
    m2h_full = m2h_full,
    h2m_all  = h2m_all,
    h2m_rep  = h2m_rep,
    m2h_rep  = m2h_rep
  )

  invisible(list(
    combined     = combined,
    ortholog_map = ortholog_map
  ))
}
