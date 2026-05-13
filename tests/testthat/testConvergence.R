library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))
context("Convergence diagnostics")

# Reuses the same simulated genome used by testMCMCROC.R. We want a real chain
# (not synthetic data) so as.mcmc() and convergence.test() are exercised end-to-end
# against actual Rcpp_MCMCAlgorithm / Rcpp_Trace objects. Sample count is kept
# small -- these tests check structure and the windowing contract, not statistical
# convergence of the chain itself.

fileName <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")

set.seed(20260513)
samples <- 30
thinning <- 1
adaptiveWidth <- 10
numMixtures <- 2
sphi_init <- c(1, 1)

mcmc <- initializeMCMCObject(samples = samples, thinning = thinning,
                             adaptive.width = adaptiveWidth,
                             est.expression = TRUE, est.csp = TRUE, est.hyper = TRUE)
genome <- initializeGenomeObject(file = fileName)
geneAssignment <- sample(c(1, 2), size = length(genome), replace = TRUE,
                         prob = c(0.3, 0.7))
parameter <- initializeParameterObject(genome, sphi_init, numMixtures,
                                       geneAssignment, split.serine = TRUE,
                                       mixture.definition = "allUnique")
model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)

sink(file.path("UnitTestingOut", "testConvergenceLog.txt"))
runMCMC(mcmc = mcmc, genome = genome, model = model, ncores = 1)
sink()

trace <- parameter$getTraceObject()
trace_len <- length(mcmc$getLogPosteriorTrace())


## ---- as.mcmc on Rcpp_MCMCAlgorithm -----------------------------------------

test_that("as.mcmc(mcmc) returns coda::mcmc of LogPosterior by default", {
  m <- as.mcmc(mcmc)
  expect_s3_class(m, "mcmc")
  expect_equal(length(m), trace_len)
  expect_equal(as.numeric(m), mcmc$getLogPosteriorTrace())
})

test_that("as.mcmc(mcmc, what='LogLikelihood') returns the LL trace", {
  m <- as.mcmc(mcmc, what = "LogLikelihood")
  expect_s3_class(m, "mcmc")
  expect_equal(as.numeric(m), mcmc$getLogLikelihoodTrace())
})

test_that("as.mcmc.Rcpp_MCMCAlgorithm rejects unknown `what`", {
  expect_error(as.mcmc(mcmc, what = "Bogus"))
})


## ---- as.mcmc on Rcpp_Trace --------------------------------------------------

test_that("as.mcmc(trace, what='Mutation') returns a samples x ncodons matrix", {
  m <- as.mcmc(trace, what = "Mutation")
  expect_s3_class(m, "mcmc")
  expect_equal(nrow(m), trace_len)
  expect_gt(ncol(m), 30) # ~40 non-stop, non-M/W codons
})

test_that("as.mcmc(trace, what='Selection') has same shape as Mutation", {
  mut <- as.mcmc(trace, what = "Mutation")
  sel <- as.mcmc(trace, what = "Selection")
  expect_equal(dim(sel), dim(mut))
})

test_that("as.mcmc(trace, what='Sphi') has one column per mixture", {
  m <- as.mcmc(trace, what = "Sphi")
  expect_s3_class(m, "mcmc")
  expect_equal(nrow(m), trace_len)
  expect_equal(ncol(m), numMixtures)
})

test_that("as.mcmc(trace, what='Mphi') equals -Sphi^2/2", {
  sphi <- as.mcmc(trace, what = "Sphi")
  mphi <- as.mcmc(trace, what = "Mphi")
  expect_equal(as.matrix(mphi), -(as.matrix(sphi)^2) / 2, tolerance = 1e-12)
})

test_that("as.mcmc(trace, what='ExpectedPhi') returns a per-sample vector", {
  m <- as.mcmc(trace, what = "ExpectedPhi")
  expect_s3_class(m, "mcmc")
  expect_equal(length(m), trace_len)
})

test_that("as.mcmc(trace, what='AcceptanceCSP') has one column per amino acid", {
  # AcceptanceCSP is recorded once per adaptive_width steps, so nrow is
  # samples / adaptiveWidth (= 3 here), not samples+1.
  m <- as.mcmc(trace, what = "AcceptanceCSP")
  expect_s3_class(m, "mcmc")
  expect_gte(nrow(m), 1)
  expect_gt(ncol(m), 10) # ~18 non-M/W/X amino acids
})

test_that("as.mcmc.Rcpp_Trace rejects unknown `what`", {
  expect_error(as.mcmc(trace, what = "Bogus"))
})

test_that("as.mcmc.Rcpp_Trace errors on not-yet-implemented `what`", {
  expect_error(as.mcmc(trace, what = "Aphi"), "not yet implemented")
  expect_error(as.mcmc(trace, what = "Sepsilon"), "not yet implemented")
  expect_error(as.mcmc(trace, what = "Expression"), "not yet implemented")
})


## ---- convergence.test windowing semantics ----------------------------------

test_that("convergence.test(mcmc) returns coda::geweke.diag", {
  g <- convergence.test(mcmc)
  expect_s3_class(g, "geweke.diag")
  expect_length(g$z, 1)
})

test_that("convergence.test(mcmc, samples=NULL) uses full trace", {
  g_null <- convergence.test(mcmc, samples = NULL)
  g_full <- convergence.test(mcmc, samples = trace_len)
  expect_equal(g_null$z, g_full$z)
})

test_that("convergence.test(mcmc, samples > trace) is clamped to full", {
  g_big  <- convergence.test(mcmc, samples = 1e6)
  g_null <- convergence.test(mcmc)
  expect_equal(g_big$z, g_null$z)
})

test_that("convergence.test(trace, what='Mutation') returns one z per codon", {
  g <- convergence.test(trace, what = "Mutation")
  expect_s3_class(g, "geweke.diag")
  expect_gt(length(g$z), 30)
})

test_that("convergence.test(trace, what='Sphi') returns one z per mixture", {
  g <- convergence.test(trace, what = "Sphi")
  expect_equal(length(g$z), numMixtures)
})


## ---- ESS via coda::effectiveSize on extracted traces -----------------------

test_that("coda::effectiveSize on LogPosterior is finite and positive", {
  ess <- coda::effectiveSize(as.mcmc(mcmc))
  expect_true(is.finite(ess))
  expect_gt(ess, 0)
})

test_that("coda::effectiveSize on Mutation trace returns one ESS per codon", {
  ess <- coda::effectiveSize(as.mcmc(trace, what = "Mutation"))
  expect_gt(length(ess), 30)
  expect_true(all(is.finite(ess[ess > 0])))
})


## ---- gelman.test (between-chain) -------------------------------------------

# A second independent chain with the same model structure but a different RNG
# stream. Required for Gelman-Rubin (needs >=2 chains). Kept small for speed --
# PSRF estimates won't be meaningful at this length, but the structural and
# error-handling contract is what these tests assert.

set.seed(20260514)
mcmc2 <- initializeMCMCObject(samples = samples, thinning = thinning,
                              adaptive.width = adaptiveWidth,
                              est.expression = TRUE, est.csp = TRUE, est.hyper = TRUE)
geneAssignment2 <- sample(c(1, 2), size = length(genome), replace = TRUE,
                          prob = c(0.3, 0.7))
parameter2 <- initializeParameterObject(genome, sphi_init, numMixtures,
                                        geneAssignment2, split.serine = TRUE,
                                        mixture.definition = "allUnique")
model2 <- initializeModelObject(parameter2, "ROC", with.phi = FALSE)

sink(file.path("UnitTestingOut", "testConvergenceLog2.txt"))
runMCMC(mcmc = mcmc2, genome = genome, model = model2, ncores = 1)
sink()

trace2 <- parameter2$getTraceObject()


test_that("gelman.test on two MCMC chains returns gelman.diag", {
  g <- gelman.test(list(mcmc, mcmc2))
  expect_s3_class(g, "gelman.diag")
  expect_true(is.matrix(g$psrf))
  expect_equal(nrow(g$psrf), 1) # logPosterior is a single series
  expect_true(all(is.finite(g$psrf[, "Point est."])))
})

test_that("gelman.test on two Trace, what='Sphi' returns PSRF per mixture", {
  g <- gelman.test(list(trace, trace2), what = "Sphi")
  expect_s3_class(g, "gelman.diag")
  expect_equal(nrow(g$psrf), numMixtures)
})

test_that("gelman.test windowing reduces chain length symmetrically", {
  expect_silent(gelman.test(list(mcmc, mcmc2), samples = 20))
})

test_that("gelman.test errors on a single chain", {
  expect_error(gelman.test(list(mcmc)),
               "at least 2 chain objects")
})

test_that("gelman.test errors on mixed MCMC + Trace types", {
  expect_error(gelman.test(list(mcmc, trace)),
               "same type")
})

test_that("gelman.test errors on a non-chain element", {
  expect_error(gelman.test(list(mcmc, "not a chain")),
               "Rcpp_MCMCAlgorithm or Rcpp_Trace")
})
