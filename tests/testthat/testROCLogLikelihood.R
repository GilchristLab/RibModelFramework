library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Tests for the ROC LOG-SPACE likelihood path.
#
# testROCNumerical.R already verifies the NON-log codon-probability
# softmax (CalculateProbabilitiesForCodons) against a hand oracle at
# 1e-12.  The log-space variant used inside the MCMC --
# ROCModel::calculateLogCodonProbabilityVector, reached via
# ROCModel::calculateLogLikelihood(genome) -- had no unit test.
#
# These tests close that gap using small, fully-controlled genomes
# containing a SINGLE multi-codon amino acid, so the per-AA aggregation
# is unambiguous (no serine-split or codon-ordering bookkeeping):
#
#   calculateLogLikelihood(genome)
#       == sum_genes sum_codons  count_codon * log(p_codon)
#
# where p_codon is the ROC codon probability for that amino acid.
# A single deterministic evaluation (no RNG, no MCMC) -> stable across
# platforms, so we use a tight 1e-8 tolerance.
# ======================================================================

context("ROC calculateLogLikelihood (log-space path)")

mutFile <- file.path("UnitTestingData", "testMCMCROCFiles", "mutation_1.csv")
selFile <- file.path("UnitTestingData", "testMCMCROCFiles", "selection_1.csv")

# Leucine codons, alphabetical; the reference codon (last alphabetical, TTG)
# is implicit and excluded from the parameter vectors.
leuCodons    <- c("CTA", "CTC", "CTG", "CTT", "TTA", "TTG")
leuNonRef    <- leuCodons[-length(leuCodons)]   # CTA,CTC,CTG,CTT,TTA
numLeuCodons <- length(leuCodons)               # 6

# Build a 2-gene, Leu-only genome. Each gene repeats all 6 Leu codons
# `reps` times, so every codon has an identical, known count per gene.
# (>=2 genes are required by parameter initialization.)
makeLeuGenome <- function(reps = c(2, 3)) {
  seqs <- vapply(reps, function(r)
    paste0("ATG", paste(rep(leuCodons, r), collapse = ""), "TAA"), character(1))
  ff <- tempfile(fileext = ".fasta")
  writeLines(as.vector(rbind(paste0(">g", seq_along(seqs)), seqs)), ff)
  ff
}

writeZeroCSV <- function(srcFile) {
  d <- read.csv(srcFile, stringsAsFactors = FALSE)
  d[, 3:ncol(d)] <- 0
  out <- tempfile(fileext = ".csv")
  write.csv(d, out, row.names = FALSE, quote = FALSE)
  out
}

# ----------------------------------------------------------------------
# Test 1: zero dM/dEta -> every codon equiprobable -> ll = N * log(1/k)
# ----------------------------------------------------------------------
test_that("calculateLogLikelihood is uniform under zero dM/dEta", {
  reps <- c(2, 3)
  ff   <- makeLeuGenome(reps)
  g    <- initializeGenomeObject(file = ff)
  p    <- initializeParameterObject(genome = g, sphi = 1, num.mixtures = 1,
                                    gene.assignment = rep(1, length(g)),
                                    mixture.definition = "allUnique")
  p$initMutationCategories(c(writeZeroCSV(mutFile)), 1, FALSE)
  p$initSelectionCategories(c(writeZeroCSV(selFile)), 1, FALSE)
  m  <- initializeModelObject(p, "ROC")

  totalLeu <- sum(reps) * numLeuCodons      # codons per gene = reps*6
  expected <- totalLeu * log(1 / numLeuCodons)
  expect_equal(m$calculateLogLikelihood(g), expected, tolerance = 1e-8)
})

# ----------------------------------------------------------------------
# Test 2: real dM/dEta -> ll must equal the codon-probability oracle
# (cross-validates the log-space path against CalculateProbabilitiesForCodons)
# ----------------------------------------------------------------------
test_that("calculateLogLikelihood matches codon-probability oracle (nonzero params)", {
  reps <- c(2, 3)
  ff   <- makeLeuGenome(reps)
  g    <- initializeGenomeObject(file = ff)
  p    <- initializeParameterObject(genome = g, sphi = 1, num.mixtures = 1,
                                    gene.assignment = rep(1, length(g)),
                                    mixture.definition = "allUnique")
  p$initMutationCategories(c(mutFile), 1, FALSE)
  p$initSelectionCategories(c(selFile), 1, FALSE)
  m  <- initializeModelObject(p, "ROC")

  mut <- read.csv(mutFile, stringsAsFactors = FALSE)
  sel <- read.csv(selFile, stringsAsFactors = FALSE)
  dM  <- mut$Posterior[match(leuNonRef, mut$Codon)]
  dE  <- sel$Posterior[match(leuNonRef, sel$Codon)]
  phi <- p$getSynthesisRate()[[1]]            # expression category 1, per gene

  llOracle <- 0
  for (gi in seq_along(reps)) {
    pr  <- m$CalculateProbabilitiesForCodons(dM, dE, phi[gi])  # [nonref..., ref]
    cnt <- rep(reps[gi], numLeuCodons)         # equal count per Leu codon in gene
    llOracle <- llOracle + sum(cnt * log(pr))
  }
  expect_equal(m$calculateLogLikelihood(g), llOracle, tolerance = 1e-8)
})

# ----------------------------------------------------------------------
# Test 3: calculateLogLikelihood is a pure function (deterministic, no leak
# side effects) -- guards the CalculateProbabilitiesForCodons buffer fix.
# ----------------------------------------------------------------------
test_that("calculateLogLikelihood is deterministic across repeated calls", {
  ff <- makeLeuGenome(c(2, 3))
  g  <- initializeGenomeObject(file = ff)
  p  <- initializeParameterObject(genome = g, sphi = 1, num.mixtures = 1,
                                  gene.assignment = rep(1, length(g)),
                                  mixture.definition = "allUnique")
  p$initMutationCategories(c(mutFile), 1, FALSE)
  p$initSelectionCategories(c(selFile), 1, FALSE)
  m  <- initializeModelObject(p, "ROC")

  expect_identical(m$calculateLogLikelihood(g), m$calculateLogLikelihood(g))
})
