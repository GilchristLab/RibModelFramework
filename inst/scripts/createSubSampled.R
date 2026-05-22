#!/usr/bin/env Rscript
# createSubSampled.R
#
# CLI wrapper around AnaCoDa::filterRFPData(). Reads tidy long-format
# RFP CSV and wide-format phi CSV, applies a configurable set of
# filters, and writes the filtered tables to an output directory.
#
# Replaces the standalone script of the same name in
# Nonsense_error_rates/R_scripts/createSubSampled.R. The library
# function does the work; this wrapper just handles I/O, legacy
# column-name fallbacks (gene -> GeneID, ORF -> GeneID), and argument
# parsing.
#
# Usage:
#   Rscript createSubSampled.R \
#     --input rfp.csv --phi phi.csv \
#     --output_dir out/ --output_rfp filtered_rfp.csv \
#     --output_phi filtered_phi.csv \
#     [--min_length_gene N] [--min_read_counts N | --min_read_counts_per_position X] \
#     [--apply_filter_to_front_back] [--remove_pause_sites Z] \
#     [--last_position_to_include N] [--trim_5 N] [--trim_3 N] \
#     [--exclude_these_genes file] [--num_genes_to_include N] [--seed N]

if (!requireNamespace("argparse", quietly = TRUE)) {
  stop("This script requires the 'argparse' package. ",
       "Install with: install.packages('argparse')")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  stop("This script requires the 'readr' package. ",
       "Install with: install.packages('readr')")
}
if (!requireNamespace("AnaCoDa", quietly = TRUE)) {
  stop("This script requires the 'AnaCoDa' package.")
}

parser <- argparse::ArgumentParser()
parser$add_argument("-i", "--input", type = "character",
                    help = "Input RFP CSV (tidy long format)")
parser$add_argument("--output_dir", type = "character",
                    help = "Output directory")
parser$add_argument("--output_rfp", type = "character",
                    help = "Output RFP filename (written under output_dir)")
parser$add_argument("--phi", type = "character",
                    help = "Input phi CSV (first column = GeneID or ORF)")
parser$add_argument("--output_phi", type = "character",
                    help = "Output phi filename (written under output_dir)")
parser$add_argument("--min_length_gene", type = "integer",
                    help = "Minimum gene length")
parser$add_argument("--last_position_to_include", type = "integer",
                    default = NULL,
                    help = "Inclusive position cutoff; codons beyond are dropped")
parser$add_argument("--trim_5", type = "integer", default = NULL,
                    help = "Drop the first N codons; renumber positions")
parser$add_argument("--trim_3", type = "integer", default = NULL,
                    help = "Drop the last N codons")
parser$add_argument("--min_read_counts", type = "integer", default = NULL,
                    help = "Minimum total RFP counts per gene")
parser$add_argument("--min_read_counts_per_position", type = "double",
                    default = NULL,
                    help = "Minimum mean RFP counts per position per gene")
parser$add_argument("--apply_filter_to_front_back",
                    help = paste("Apply read-count threshold to 5'",
                                 "and 3' halves independently"),
                    action = "store_true")
parser$add_argument("--remove_pause_sites", type = "double", default = NULL,
                    help = "Z-score threshold; genes with any site exceeding Z are dropped")
parser$add_argument("--exclude_these_genes", type = "character",
                    default = NULL,
                    help = "Plaintext file of GeneIDs to exclude (one per line)")
parser$add_argument("--num_genes_to_include", type = "integer",
                    default = NULL,
                    help = "Random subsample to N genes after filtering")
parser$add_argument("--seed", type = "integer", default = NULL,
                    help = "Optional RNG seed for reproducible subsampling")

args <- parser$parse_args()

# -- read inputs --------------------------------------------------------

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

rfp <- readr::read_csv(args$input, show_col_types = FALSE)
phi <- readr::read_csv(args$phi,   show_col_types = FALSE)

# Legacy column-name handling: the library expects GeneID in both tables.
if (!"GeneID" %in% colnames(rfp)) {
  if ("gene" %in% colnames(rfp)) {
    rfp <- dplyr::rename(rfp, GeneID = "gene")
  } else {
    stop("rfp input must contain a GeneID (or legacy 'gene') column")
  }
}
if (!"GeneID" %in% colnames(phi)) {
  if ("ORF" %in% colnames(phi)) {
    phi <- dplyr::rename(phi, GeneID = "ORF")
  } else {
    stop("phi input must contain a GeneID (or legacy 'ORF') column")
  }
}

exclude_genes <- if (is.character(args$exclude_these_genes)) {
  readLines(args$exclude_these_genes)
} else {
  NULL
}

# -- call library function ----------------------------------------------

out <- AnaCoDa::filterRFPData(
  rfp                      = rfp,
  phi                      = phi,
  min_length_gene          = args$min_length_gene,
  last_position_to_include = args$last_position_to_include,
  trim_5                   = args$trim_5,
  trim_3                   = args$trim_3,
  min_read_counts          = args$min_read_counts,
  min_read_counts_per_pos  = args$min_read_counts_per_position,
  apply_filter_to_halves   = args$apply_filter_to_front_back,
  remove_pause_sites       = args$remove_pause_sites,
  exclude_genes            = exclude_genes,
  num_genes_to_include     = args$num_genes_to_include,
  seed                     = args$seed
)

cat("Remaining number of genes: ", nrow(out$phi), "\n", sep = "")

readr::write_csv(out$rfp, file.path(args$output_dir, args$output_rfp))
readr::write_csv(out$phi, file.path(args$output_dir, args$output_phi))
