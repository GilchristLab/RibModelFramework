library(testthat)
library(AnaCoDa)

context("filterRFPData")

# --------------------------------------------------------------------
# Fixture overview (inst/extdata/test_rfp.csv, inst/extdata/test_phi.csv):
#
#   GENE_A: 8 positions, RFPCount=50 each       (passes everything)
#   GENE_B: 4 positions, RFPCount=50 each       (fails min_length>=5)
#   GENE_C: 10 positions, RFPCount=1 each       (fails read-count filters)
#   GENE_D: 10 positions, front 5 @50, back 5 @1
#                                                (passes global,
#                                                 fails halves filter)
#   GENE_E: 10 positions, mostly 50, position 9 = 500 (pause-site spike)
#   GENE_X: 10 positions, RFPCount=40 each      (kept; used for exclusion test)
#
#   phi has an extra row GENE_MISSING_IN_RFP that should be dropped
#   by the phi co-filter (not present in rfp).
# --------------------------------------------------------------------

read_fixture <- function() {
  rfp_path <- system.file("extdata", "test_rfp.csv", package = "AnaCoDa")
  phi_path <- system.file("extdata", "test_phi.csv", package = "AnaCoDa")
  list(rfp = read.csv(rfp_path, stringsAsFactors = FALSE),
       phi = read.csv(phi_path, stringsAsFactors = FALSE))
}

test_that("input validation: missing rfp columns errors", {
  f <- read_fixture()
  bad_rfp <- f$rfp[, c("GeneID", "Position", "RFPCount")]  # missing Codon
  expect_error(filterRFPData(bad_rfp, f$phi),
               regexp = "missing required columns")
})

test_that("input validation: mutually exclusive read-count thresholds", {
  f <- read_fixture()
  expect_error(filterRFPData(f$rfp, f$phi,
                             min_read_counts = 100,
                             min_read_counts_per_pos = 10),
               regexp = "only one of min_read_counts")
})

test_that("warning when both last_position_to_include and trim_5 set", {
  f <- read_fixture()
  expect_warning(filterRFPData(f$rfp, f$phi,
                               last_position_to_include = 5,
                               trim_5 = 2),
                 regexp = "last_position_to_include will be applied first")
})

test_that("min_length_gene removes short genes only", {
  f <- read_fixture()
  out <- filterRFPData(f$rfp, f$phi, min_length_gene = 5)
  surviving <- unique(out$rfp$GeneID)
  expect_false("GENE_B" %in% surviving)
  expect_true(all(c("GENE_A", "GENE_C", "GENE_D", "GENE_E", "GENE_X")
                  %in% surviving))
})

test_that("min_read_counts (global) drops low-read genes", {
  f <- read_fixture()
  out <- filterRFPData(f$rfp, f$phi, min_read_counts = 100)
  surviving <- unique(out$rfp$GeneID)
  # GENE_C has total 10, fails. Others all have totals > 100.
  expect_false("GENE_C" %in% surviving)
  expect_true("GENE_A" %in% surviving)
  expect_true("GENE_D" %in% surviving)  # 255 total, passes
})

test_that("min_read_counts_per_pos (global) drops low-density genes", {
  f <- read_fixture()
  out <- filterRFPData(f$rfp, f$phi, min_read_counts_per_pos = 10)
  surviving <- unique(out$rfp$GeneID)
  # GENE_C has avg=1, fails. GENE_D has avg=25.5, passes.
  expect_false("GENE_C" %in% surviving)
  expect_true("GENE_D" %in% surviving)
})

test_that("halves mode (count): GENE_D fails because back half is low", {
  f <- read_fixture()
  out <- filterRFPData(f$rfp, f$phi,
                       min_read_counts = 100,
                       apply_filter_to_halves = TRUE)
  surviving <- unique(out$rfp$GeneID)
  # GENE_D global passes (255>100) but halves fails (back=5 < 100).
  expect_false("GENE_D" %in% surviving)
  expect_true("GENE_A" %in% surviving)
})

test_that("halves mode (per-pos): GENE_D fails because back half avg is low", {
  f <- read_fixture()
  out <- filterRFPData(f$rfp, f$phi,
                       min_read_counts_per_pos = 10,
                       apply_filter_to_halves = TRUE)
  surviving <- unique(out$rfp$GeneID)
  # GENE_D back-half avg = 1, fails threshold 10.
  expect_false("GENE_D" %in% surviving)
})

test_that("remove_pause_sites drops GENE_E (planted spike)", {
  f <- read_fixture()
  out <- filterRFPData(f$rfp, f$phi, remove_pause_sites = 2)
  surviving <- unique(out$rfp$GeneID)
  expect_false("GENE_E" %in% surviving)
  expect_true("GENE_A" %in% surviving)  # no spike
})

test_that("last_position_to_include drops codons beyond cutoff", {
  f <- read_fixture()
  out <- filterRFPData(f$rfp, f$phi, last_position_to_include = 5)
  expect_true(all(out$rfp$Position <= 5))
})

test_that("trim_5 removes first N codons and renumbers positions", {
  f <- read_fixture()
  out <- filterRFPData(f$rfp, f$phi, trim_5 = 2)
  # No remaining position should be <= 0 after renumbering.
  expect_true(all(out$rfp$Position >= 1))
  # GENE_A original had positions 1..8; after trim_5=2, expect 1..6.
  gene_a <- subset(out$rfp, GeneID == "GENE_A")
  expect_equal(sort(gene_a$Position), 1:6)
})

test_that("trim_3 removes last N codons", {
  f <- read_fixture()
  out <- filterRFPData(f$rfp, f$phi, trim_3 = 2)
  # GENE_A original had positions 1..8; after trim_3=2, expect 1..6
  # (drops last 2 positions; uses strict inequality Rel.to.3.end > trim_3,
  # so positions where Length-Position <= trim_3 are removed: 7, 8).
  gene_a <- subset(out$rfp, GeneID == "GENE_A")
  expect_equal(sort(gene_a$Position), 1:6)
})

test_that("exclude_genes removes the listed gene IDs", {
  f <- read_fixture()
  out <- filterRFPData(f$rfp, f$phi,
                       exclude_genes = c("GENE_X"))
  expect_false("GENE_X" %in% out$rfp$GeneID)
})

test_that("num_genes_to_include subsamples to N and is reproducible", {
  f <- read_fixture()
  out1 <- filterRFPData(f$rfp, f$phi, num_genes_to_include = 3, seed = 42)
  out2 <- filterRFPData(f$rfp, f$phi, num_genes_to_include = 3, seed = 42)
  expect_equal(length(unique(out1$rfp$GeneID)), 3)
  expect_equal(unique(out1$rfp$GeneID), unique(out2$rfp$GeneID))
})

test_that("num_genes_to_include is no-op when N >= surviving genes", {
  f <- read_fixture()
  # 6 genes in fixture; ask for 10.
  out <- filterRFPData(f$rfp, f$phi, num_genes_to_include = 10, seed = 42)
  expect_equal(length(unique(out$rfp$GeneID)), 6)
})

test_that("phi co-filter drops genes not in rfp, preserves order", {
  f <- read_fixture()
  out <- filterRFPData(f$rfp, f$phi)
  # GENE_MISSING_IN_RFP is in phi but not rfp; should be gone.
  expect_false("GENE_MISSING_IN_RFP" %in% out$phi$GeneID)
  # phi row order matches unique(out$rfp$GeneID).
  expect_equal(out$phi$GeneID, unique(out$rfp$GeneID))
})

test_that("phi co-filter warns when rfp has genes missing in phi", {
  f <- read_fixture()
  # Drop GENE_X from phi.
  phi_short <- subset(f$phi, GeneID != "GENE_X")
  expect_warning(
    out <- filterRFPData(f$rfp, phi_short),
    regexp = "present in rfp but not in phi"
  )
  # GENE_X should appear in the output phi with NA value but GeneID preserved.
  gene_x_row <- subset(out$phi, GeneID == "GENE_X")
  expect_equal(nrow(gene_x_row), 1)
  expect_true(is.na(gene_x_row$phi))
})

test_that("end-to-end: combined filters reproduce expected gene set", {
  f <- read_fixture()
  out <- filterRFPData(
    f$rfp, f$phi,
    min_length_gene   = 5,        # drops GENE_B
    min_read_counts   = 100,      # drops GENE_C
    remove_pause_sites = 2,       # drops GENE_E
    exclude_genes     = "GENE_X"  # drops GENE_X
  )
  surviving <- sort(unique(out$rfp$GeneID))
  expect_equal(surviving, c("GENE_A", "GENE_D"))
})
