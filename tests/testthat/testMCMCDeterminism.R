library(testthat)
library(AnaCoDa)

# ======================================================================
# Determinism regression test for the RNG-determinism fix (2026-06).
#
# Two sources of non-set.seed-controlled randomness were routed through
# R's RNG so that runMCMC(..., ncores = 1) is bit-reproducible:
#   1. MCMCAlgorithm::acceptRejectSynthesisRateLevelForAllGenes -- the
#      per-gene Metropolis threshold and mixture-assignment draws were
#      taken from a default-seeded std::default_random_engine that was
#      also shared across OpenMP threads (data race). They are now
#      pre-drawn from R's RNG outside the parallel region.
#   2. MCMCAlgorithm::acceptRejectCodonSpecificParameter -- the
#      random-scan std::shuffle engine was seeded from the wall clock.
#      It is now seeded from R's RNG.
#
# If either source regresses (clock seed, shared engine, RNG call inside
# an OpenMP region), two same-seed chains will diverge and this fails.
# ======================================================================

context("MCMC determinism (ncores = 1)")

fastaFile <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")

run_short_roc_chain <- function(seed) {
  set.seed(seed)
  genome    <- initializeGenomeObject(file = fastaFile)
  parameter <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1,
                                         gene.assignment = rep(1, length(genome)),
                                         mixture.definition = "allUnique")
  model <- initializeModelObject(parameter, "ROC")
  mcmc  <- initializeMCMCObject(samples = 10, thinning = 10, adaptive.width = 10,
                                est.expression = TRUE, est.csp = TRUE, est.hyper = TRUE)
  invisible(capture.output(suppressMessages(
    runMCMC(mcmc, genome, model, ncores = 1, divergence.iteration = 0))))
  mcmc$getLogPosteriorTrace()
}

test_that("two ncores=1 chains with the same seed are bit-identical", {
  a <- run_short_roc_chain(987654)
  b <- run_short_roc_chain(987654)
  expect_identical(a, b)
})

test_that("different seeds produce different chains (sanity)", {
  a <- run_short_roc_chain(987654)
  c <- run_short_roc_chain(123456)
  # Initial element is the pre-sampling 0.0 in both; the chain must diverge.
  expect_false(isTRUE(all.equal(a, c)))
})
