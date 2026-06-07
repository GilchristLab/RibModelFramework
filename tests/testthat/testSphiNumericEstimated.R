library(testthat)
library(AnaCoDa)

# ======================================================================
# Regression tests: sphi=<numeric> must be ESTIMATED, not frozen.
#
# Background: prior to this fix (PR #37 regression), passing a numeric
# sphi to initializeParameterObject() silently called fixSphi(), which
# set fix_stdDevSynthesis=TRUE and made proposeStdDevSynthesisRate() a
# no-op.  sphi was frozen even with est.hyper=TRUE.
#
# The fix: numeric sphi= is now the INITIAL VALUE only; estimation uses
# Uniform(0,10) prior by default.  To freeze sphi, use
# phi.sphi=fixed(value) or est.hyper=FALSE.
#
# Covers:
#  (a) sphi=<numeric> + est.hyper=TRUE -> sphi trace MOVES
#  (b) phi.sphi=fixed(1) explicit      -> sphi trace FROZEN
#  (c) sphi=NA                         -> sphi trace MOVES (baseline)
#  (d) init smoke: all four models accept numeric sphi without error
# ======================================================================

context("sphi=numeric estimated (regression)")

# ---- shared fixtures --------------------------------------------------

.roc_genome <- function() {
  initializeGenomeObject(
    file = system.file("extdata", "genome.fasta", package = "AnaCoDa")
  )
}

.pa_rfp_file  <- file.path("UnitTestingData", "testPAFiles",      "mini_rfp.csv")
.panse_rfp_file <- file.path("UnitTestingData", "testMCMCPANSEFiles", "mini_rfp.csv")

# Run a short ROC MCMC (50 samples) and return the sphi trace for mixture 1.
# Caller supplies the parameter object so it can be inspected after the run.
.run_roc_and_get_sphi <- function(parameter, genome, seed = 42) {
  model <- initializeModelObject(parameter, "ROC")
  mcmc  <- initializeMCMCObject(
    samples        = 50,
    thinning       = 1,
    adaptive.width = 25,
    est.expression = TRUE,
    est.csp        = TRUE,
    est.hyper      = TRUE
  )
  set.seed(seed)
  invisible(capture.output(
    suppressMessages(runMCMC(mcmc, genome, model, ncores = 1,
                             divergence.iteration = 0))
  ))
  parameter$getTraceObject()$getStdDevSynthesisRateTraces()[[1]]
}

# ======================================================================
# (a) sphi=<numeric> + est.hyper=TRUE -> sphi trace moves
# ======================================================================

test_that("sphi=1 (numeric) + est.hyper=TRUE: sphi trace moves (regression)", {
  g <- .roc_genome()
  p <- initializeParameterObject(g, sphi = 1, num.mixtures = 1,
                                 gene.assignment = rep(1L, length(g)))
  # Confirm parameter resolved to estimated (not fixed)
  # SPHI_PRIOR_UNIFORM = 2
  expect_equal(p$getSphiPriorType(), 2L,
               label = "sphi=1 should resolve to estimated(uniform) -> sphiPriorType==2")

  sphi_trace <- .run_roc_and_get_sphi(p, g)
  # trace[1] is a pre-MCMC sentinel (0); samples start at trace[2]
  samples_trace <- sphi_trace[-1]
  n_unique   <- length(unique(round(samples_trace, 8)))
  expect_gt(n_unique, 5L,
            label = paste("sphi trace should move; got n_unique =", n_unique))
  expect_gt(diff(range(samples_trace)), 1e-6,
            label = "sphi trace sample range should be non-zero")
})

# ======================================================================
# (b) phi.sphi=fixed(1) explicit -> sphi trace FROZEN
# ======================================================================

test_that("phi.sphi=fixed(1) explicit: sphi trace is frozen at 1", {
  g <- .roc_genome()
  p <- initializeParameterObject(g, sphi = 1, num.mixtures = 1,
                                 gene.assignment = rep(1L, length(g)),
                                 phi.sphi = fixed(value = 1))
  sphi_trace <- .run_roc_and_get_sphi(p, g)
  # trace[1] is a pre-MCMC sentinel (0); samples start at trace[2]
  samples_trace <- sphi_trace[-1]
  expect_lt(diff(range(samples_trace)), 1e-8,
            label = "explicit fixed() should freeze sphi: sample range must be ~0")
})

# ======================================================================
# (c) sphi=NA -> also estimated, trace moves (baseline)
# ======================================================================

test_that("sphi=NA: sphi trace moves (baseline)", {
  g <- .roc_genome()
  p <- initializeParameterObject(g, sphi = NA, num.mixtures = 1,
                                 gene.assignment = rep(1L, length(g)))
  sphi_trace <- .run_roc_and_get_sphi(p, g)
  # trace[1] is a pre-MCMC sentinel (0); samples start at trace[2]
  samples_trace <- sphi_trace[-1]
  n_unique   <- length(unique(round(samples_trace, 8)))
  expect_gt(n_unique, 5L,
            label = paste("sphi=NA trace should move; got n_unique =", n_unique))
})

# ======================================================================
# (d) init smoke: all four models accept numeric sphi without error
# ======================================================================

test_that("ROC initializeParameterObject accepts numeric sphi without error", {
  g <- .roc_genome()
  expect_error(
    initializeParameterObject(g, sphi = 1.5, num.mixtures = 1,
                              gene.assignment = rep(1L, length(g)),
                              model = "ROC"),
    NA
  )
})

test_that("FONSE initializeParameterObject accepts numeric sphi without error", {
  g <- .roc_genome()
  expect_error(
    initializeParameterObject(g, sphi = 1.5, num.mixtures = 1,
                              gene.assignment = rep(1L, length(g)),
                              model = "FONSE"),
    NA
  )
})

test_that("PA initializeParameterObject accepts numeric sphi without error", {
  skip_if_not(file.exists(.pa_rfp_file), "PA mini_rfp fixture not found")
  g <- initializeGenomeObject(file = .pa_rfp_file, fasta = FALSE)
  expect_error(
    initializeParameterObject(g, sphi = 1.5, num.mixtures = 1,
                              gene.assignment = rep(1L, length(g)),
                              model = "PA",
                              mixture.definition = "allUnique"),
    NA
  )
})

test_that("PANSE initializeParameterObject accepts numeric sphi without error", {
  skip_if_not(file.exists(.panse_rfp_file), "PANSE mini_rfp fixture not found")
  g <- initializeGenomeObject(file = .panse_rfp_file, fasta = FALSE)
  expect_error(
    initializeParameterObject(g, sphi = 1.5, num.mixtures = 1,
                              gene.assignment = rep(1L, length(g)),
                              model = "PANSE",
                              split.serine = FALSE,
                              mixture.definition = "allUnique"),
    NA
  )
})
