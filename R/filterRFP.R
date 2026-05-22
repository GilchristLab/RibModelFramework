#' Filter and Subsample RFP Data
#'
#' Filters a tidy long-format ribosome footprint (RFP) data frame and a
#' matching wide-format expression (phi) data frame by gene length, read
#' counts, pause-site outliers, position trimming, gene-ID exclusion,
#' and random gene subsampling. Filters are applied in a fixed order
#' (see Details). The accompanying phi table is subset to the surviving
#' gene set after RFP filtering.
#'
#' @param rfp A data frame with at least the columns \code{GeneID},
#'   \code{Position}, \code{Codon}, and \code{RFPCount}. One row per
#'   position per gene (tidy long format).
#' @param phi A data frame with a \code{GeneID} column plus one or more
#'   phi/expression value columns.
#' @param min_length_gene Optional integer. Genes with fewer than this
#'   many positions are removed.
#' @param last_position_to_include Optional integer position. Codons
#'   beyond this position are dropped (inclusive cutoff). Applied before
#'   \code{trim_5}/\code{trim_3}; a warning is emitted if used together
#'   with \code{trim_5}.
#' @param trim_5 Optional integer. Drop the first \code{trim_5} codons
#'   of every gene and renumber positions to start at 1.
#' @param trim_3 Optional integer. Drop the last \code{trim_3} codons of
#'   every gene.
#' @param min_read_counts Optional integer. Minimum total
#'   \code{RFPCount} per gene. Mutually exclusive with
#'   \code{min_read_counts_per_pos}.
#' @param min_read_counts_per_pos Optional numeric. Minimum mean
#'   \code{RFPCount}/position per gene. Mutually exclusive with
#'   \code{min_read_counts}.
#' @param apply_filter_to_halves Logical. If \code{TRUE}, the
#'   \code{min_read_counts} or \code{min_read_counts_per_pos} threshold
#'   is evaluated independently on the 5' and 3' halves of each gene;
#'   both halves must pass. (Renamed from the original CLI flag
#'   \code{--apply_filter_to_front_back}.)
#' @param remove_pause_sites Optional positive numeric. Z-score
#'   threshold for pause-site removal: genes containing any position
#'   whose per-gene Z-score of \code{RFPCount} exceeds the threshold are
#'   removed in full.
#' @param exclude_genes Optional character vector of \code{GeneID}s to
#'   exclude before any other filter.
#' @param num_genes_to_include Optional integer. After all other
#'   filters, draw a uniform random subsample of
#'   \code{num_genes_to_include} genes. Has no effect if fewer genes
#'   remain after filtering.
#' @param seed Optional integer passed to \code{\link[base]{set.seed}}
#'   immediately before \code{num_genes_to_include} sampling, for
#'   reproducibility. Note that this modifies the global RNG state as a
#'   side effect.
#'
#' @return A list with two data frames: \code{rfp} (filtered RFP table)
#'   and \code{phi} (rows subset to the surviving gene set, in the same
#'   order as \code{unique(rfp$GeneID)}).
#'
#' @details
#' Filters are applied in this fixed order:
#' \enumerate{
#'   \item Exclude by gene ID (\code{exclude_genes})
#'   \item Minimum gene length (\code{min_length_gene})
#'   \item Minimum total read count, global or per-half
#'     (\code{min_read_counts})
#'   \item Minimum read count per position, global or per-half
#'     (\code{min_read_counts_per_pos})
#'   \item Pause-site removal (\code{remove_pause_sites})
#'   \item Position ceiling (\code{last_position_to_include})
#'   \item 5' trim (\code{trim_5})
#'   \item 3' trim (\code{trim_3})
#'   \item Random subsample (\code{num_genes_to_include})
#' }
#' After RFP filtering completes, the phi table is subset and reordered
#' to match \code{unique(rfp$GeneID)}.
#'
#' @examples
#' \dontrun{
#' rfp <- read.csv(system.file("extdata", "test_rfp.csv",
#'                             package = "AnaCoDa"))
#' phi <- read.csv(system.file("extdata", "test_phi.csv",
#'                             package = "AnaCoDa"))
#' out <- filterRFPData(rfp, phi,
#'                      min_length_gene = 100,
#'                      min_read_counts = 200)
#' nrow(out$rfp); nrow(out$phi)
#' }
#'
#' @importFrom dplyr filter group_by summarize select mutate left_join n pull
#' @importFrom tidyr pivot_wider
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
filterRFPData <- function(rfp, phi,
                          min_length_gene          = NULL,
                          last_position_to_include = NULL,
                          trim_5                   = NULL,
                          trim_3                   = NULL,
                          min_read_counts          = NULL,
                          min_read_counts_per_pos  = NULL,
                          apply_filter_to_halves   = FALSE,
                          remove_pause_sites       = NULL,
                          exclude_genes            = NULL,
                          num_genes_to_include     = NULL,
                          seed                     = NULL) {

  required    <- c("GeneID", "Position", "Codon", "RFPCount")
  missing_cols <- setdiff(required, colnames(rfp))
  if (length(missing_cols) > 0) {
    stop("rfp is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (!"GeneID" %in% colnames(phi)) {
    stop("phi must contain a GeneID column")
  }
  if (!is.null(min_read_counts) && !is.null(min_read_counts_per_pos)) {
    stop("only one of min_read_counts and min_read_counts_per_pos ",
         "should be specified")
  }
  if (!is.null(trim_5) && !is.null(last_position_to_include)) {
    warning("both last_position_to_include and trim_5 specified. ",
            "last_position_to_include will be applied first, followed ",
            "by additional trimming from the 5'-end.")
  }

  if (!is.null(exclude_genes)) {
    rfp <- dplyr::filter(rfp, !.data$GeneID %in% exclude_genes)
  }
  if (!is.null(min_length_gene)) {
    rfp <- .filterRFPByGeneLength(rfp, min_length_gene)
  }
  if (!is.null(min_read_counts)) {
    rfp <- .filterRFPByReadCount(rfp, min_read_counts,
                                 per_position = FALSE,
                                 halves       = apply_filter_to_halves)
  }
  if (!is.null(min_read_counts_per_pos)) {
    rfp <- .filterRFPByReadCount(rfp, min_read_counts_per_pos,
                                 per_position = TRUE,
                                 halves       = apply_filter_to_halves)
  }
  if (is.numeric(remove_pause_sites) && remove_pause_sites > 0) {
    rfp <- .filterRFPByPauseSites(rfp, remove_pause_sites)
  }
  rfp <- .trimRFPPositions(rfp,
                           trim_5           = trim_5,
                           trim_3           = trim_3,
                           position_ceiling = last_position_to_include)
  if (!is.null(num_genes_to_include)) {
    gene_names <- unique(rfp$GeneID)
    if (length(gene_names) > num_genes_to_include) {
      if (!is.null(seed)) set.seed(seed)
      keep <- sample(gene_names, size = num_genes_to_include,
                     replace = FALSE)
      rfp  <- dplyr::filter(rfp, .data$GeneID %in% keep)
    }
  }

  phi <- .coFilterPhi(phi, rfp)
  list(rfp = rfp, phi = phi)
}


# -- internal helpers (not exported) ------------------------------------

.filterRFPByGeneLength <- function(rfp, min_length) {
  short <- rfp %>%
    dplyr::group_by(.data$GeneID) %>%
    dplyr::summarize(Length = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(.data$Length < min_length)
  dplyr::filter(rfp, !.data$GeneID %in% short$GeneID)
}

.filterRFPByReadCount <- function(rfp, threshold, per_position, halves) {
  if (!halves) {
    summary_df <- rfp %>%
      dplyr::group_by(.data$GeneID) %>%
      dplyr::summarize(Total  = sum(.data$RFPCount),
                       Length = dplyr::n(),
                       .groups = "drop") %>%
      dplyr::mutate(Average = .data$Total / .data$Length)
    keep <- if (per_position) {
      dplyr::filter(summary_df, .data$Average >= threshold)$GeneID
    } else {
      dplyr::filter(summary_df, .data$Total >= threshold)$GeneID
    }
  } else {
    halved <- rfp %>%
      dplyr::group_by(.data$GeneID) %>%
      dplyr::mutate(Loc = ifelse(.data$Position / dplyr::n() <= 0.5,
                                 "Front", "Back")) %>%
      dplyr::group_by(.data$GeneID, .data$Loc) %>%
      dplyr::summarize(Length = dplyr::n(),
                       Total  = sum(.data$RFPCount),
                       .groups = "drop") %>%
      dplyr::mutate(PerPos = .data$Total / .data$Length)
    value_col <- if (per_position) "PerPos" else "Total"
    wide <- tidyr::pivot_wider(halved,
                               id_cols    = "GeneID",
                               names_from = "Loc",
                               values_from = tidyr::all_of(value_col))
    keep <- dplyr::filter(wide,
                          .data$Front >= threshold,
                          .data$Back  >= threshold)$GeneID
  }
  dplyr::filter(rfp, .data$GeneID %in% keep)
}

.filterRFPByPauseSites <- function(rfp, z_threshold) {
  extreme <- rfp %>%
    dplyr::group_by(.data$GeneID) %>%
    dplyr::mutate(Z = (.data$RFPCount - mean(.data$RFPCount)) /
                     stats::sd(.data$RFPCount)) %>%
    dplyr::filter(.data$Z > z_threshold) %>%
    dplyr::pull("GeneID") %>%
    unique()
  dplyr::filter(rfp, !.data$GeneID %in% extreme)
}

.trimRFPPositions <- function(rfp,
                              trim_5           = NULL,
                              trim_3           = NULL,
                              position_ceiling = NULL) {
  if (!is.null(position_ceiling)) {
    rfp <- dplyr::filter(rfp, .data$Position <= position_ceiling)
  }
  if (!is.null(trim_5)) {
    rfp <- rfp %>%
      dplyr::filter(.data$Position > trim_5) %>%
      dplyr::mutate(Position = .data$Position - trim_5)
  }
  if (!is.null(trim_3)) {
    lengths <- rfp %>%
      dplyr::group_by(.data$GeneID) %>%
      dplyr::summarize(Length = dplyr::n(), .groups = "drop")
    rfp <- rfp %>%
      dplyr::left_join(lengths, by = "GeneID") %>%
      dplyr::mutate(RelTo3End = .data$Length - .data$Position) %>%
      dplyr::filter(.data$RelTo3End > trim_3) %>%
      dplyr::select(-"Length", -"RelTo3End")
  }
  rfp
}

.coFilterPhi <- function(phi, rfp) {
  surviving <- unique(rfp$GeneID)
  idx <- match(surviving, phi$GeneID)
  missing_in_phi <- surviving[is.na(idx)]
  if (length(missing_in_phi) > 0) {
    warning(length(missing_in_phi),
            " gene(s) present in rfp but not in phi; ",
            "their phi rows will be NA. First few: ",
            paste(utils::head(missing_in_phi, 5), collapse = ", "))
  }
  out <- phi[idx, , drop = FALSE]
  # preserve GeneID for missing-in-phi rows (match returns NA so phi[NA,]
  # produces all-NA row including GeneID; reassign from surviving)
  out$GeneID <- surviving
  rownames(out) <- NULL
  out
}
