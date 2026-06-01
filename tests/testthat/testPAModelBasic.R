library(testthat)
library(AnaCoDa)

# ======================================================================
# Active unit tests for the PA (pausing-time) model.
#
# Context: the legacy testPAModel.R is an intentional MANUAL vignette
# (a 20,000-sample parameter-recovery run, skip()-guarded) and is not a
# CI test.  Until now PA had NO active automated coverage.  These tests
# provide a fast, deterministic smoke + regression layer using a small,
# modern RFP fixture.
#
# The fixture mini_rfp.csv uses the current readRFPData format
# (GeneID,Position,Codon,Mixture,RFPCount with Mixture >= 1); the older
# headerless PA fixtures (pa_rfpdata.csv, testPAModel.csv) crash the
# current parser, which is why a new fixture was added.
#
# Reproducibility note: PA synthesis-rate acceptance flows through the
# shared MCMCAlgorithm::acceptRejectSynthesisRateLevelForAllGenes, so the
# 2026-06 RNG-determinism fix makes PA runMCMC(..., ncores = 1)
# bit-reproducible under set.seed (see testMCMCDeterminism.R for ROC).
# ======================================================================

context("PA model (basic)")

rfpFile <- file.path("UnitTestingData", "testPAFiles", "mini_rfp.csv")

test_that("PA RFP fixture exists and loads via readRFPData", {
  expect_true(file.exists(rfpFile))
  g <- initializeGenomeObject(file = rfpFile, fasta = FALSE)
  expect_equal(length(g), 2)
})

# Helper: build PA model + run a short deterministic chain, return the
# log-posterior trace.
runShortPA <- function(seed, samples = 5, thinning = 2) {
  set.seed(seed)
  g <- initializeGenomeObject(file = rfpFile, fasta = FALSE)
  p <- initializeParameterObject(g, sphi = c(2), num.mixtures = 1,
                                 gene.assignment = rep(1, length(g)),
                                 model = "PA", mixture.definition = "allUnique")
  m <- initializeModelObject(p, "PA")
  mcmc <- initializeMCMCObject(samples = samples, thinning = thinning,
                               adaptive.width = thinning,
                               est.expression = TRUE, est.csp = TRUE, est.hyper = TRUE)
  invisible(capture.output(suppressMessages(runMCMC(mcmc, g, m, ncores = 1))))
  mcmc$getLogPosteriorTrace()
}

test_that("PA MCMC produces a finite log-posterior trace of expected length", {
  samples <- 5
  tr <- runShortPA(seed = 42, samples = samples, thinning = 2)
  expect_equal(length(tr), samples + 1)
  expect_true(all(is.finite(tr)))
})

test_that("PA MCMC is bit-reproducible at ncores = 1 (same seed)", {
  a <- runShortPA(seed = 42)
  b <- runShortPA(seed = 42)
  expect_identical(a, b)
})

test_that("PA MCMC diverges for different seeds (sanity)", {
  a <- runShortPA(seed = 42)
  c <- runShortPA(seed = 7)
  expect_false(isTRUE(all.equal(a, c)))
})
