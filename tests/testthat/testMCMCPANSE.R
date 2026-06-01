library(testthat)
library(AnaCoDa)

# ======================================================================
# Tests for the PANSE MCMC code pathway.
#
# Scope: exercises runMCMC() end-to-end with model="PANSE" -- the
# acceptRejectSynthesisRateLevelForAllGenes / acceptRejectCodonSpecific-
# Parameter / acceptRejectForHyperParameters loop -- NOT the Stan port.
#
# Key PANSE differences vs PA in the MCMC path:
#   - Three CSP types: Alpha, LambdaPrime, NSERate (per codon)
#   - Survival probability log P(no nonsense through codon i) via
#     forward-Lentz continued fraction (exact, post-2026 fix); older code
#     used an approximation that produced quiet_NaN on extreme params.
#   - A validity guard floors log(pSuccess) to 0 on roundoff excursions
#     and counts genuine breaches (prob_success_breach_count).
#   - split.serine = FALSE (PANSEParameter convention)
#
# Reproducibility: the 2026-06 RNG-determinism fix seeds both the CSP
# scan-order shuffle and the per-gene Metropolis draws from R's RNG, so
# runMCMC(..., ncores = 1) is bit-reproducible under set.seed.
# ncores > 1 remains nondeterministic (OpenMP FP reduction order); all
# tests pin ncores = 1.
# ======================================================================

context("MCMC with PANSE")

rfpFile <- file.path("UnitTestingData", "testMCMCPANSEFiles", "mini_rfp.csv")

test_that("PANSE MCMC fixture exists and loads", {
  expect_true(file.exists(rfpFile))
  g <- initializeGenomeObject(file = rfpFile, fasta = FALSE)
  expect_equal(length(g), 2)
})

# Run a short PANSE MCMC chain; return the log-posterior trace.
runPANSEMCMC <- function(seed, samples = 5, thinning = 2) {
  set.seed(seed)
  g <- initializeGenomeObject(file = rfpFile, fasta = FALSE)
  p <- initializeParameterObject(genome = g, sphi = c(1), num.mixtures = 1,
                                 gene.assignment = rep(1, length(g)),
                                 model = "PANSE",
                                 split.serine = FALSE,
                                 mixture.definition = "allUnique")
  m <- initializeModelObject(p, "PANSE")
  mcmc <- initializeMCMCObject(samples = samples, thinning = thinning,
                               adaptive.width = thinning,
                               est.expression = TRUE, est.csp = TRUE, est.hyper = TRUE)
  invisible(capture.output(suppressMessages(runMCMC(mcmc, g, m, ncores = 1))))
  mcmc$getLogPosteriorTrace()
}

test_that("PANSE MCMC produces a finite log-posterior trace of expected length", {
  samples <- 5
  tr <- runPANSEMCMC(seed = 42, samples = samples, thinning = 2)
  expect_equal(length(tr), samples + 1)
  expect_true(all(is.finite(tr)))
})

test_that("PANSE MCMC is bit-reproducible at ncores = 1 (same seed)", {
  a <- runPANSEMCMC(seed = 42)
  b <- runPANSEMCMC(seed = 42)
  expect_identical(a, b)
})

test_that("PANSE MCMC diverges for different seeds (sanity)", {
  a <- runPANSEMCMC(seed = 42)
  c <- runPANSEMCMC(seed = 7)
  expect_false(isTRUE(all.equal(a, c)))
})
