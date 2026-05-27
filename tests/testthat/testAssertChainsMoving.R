library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))
context("assertChainsMoving diagnostic")

# ---------------------------------------------------------------------------
# These tests verify assertChainsMoving() without running a full MCMC fit.
# Strategy: monkey-patch the trace accessors on a lightweight ROC parameter
# object (cheap to initialise) so we can inject controlled AR values.
#
# True-positive test (AR = 0, stuck chain):
#   Override getCodonSpecificAcceptanceRateTraceForAA to return 0 for every
#   AA.  assertChainsMoving() must fire a warning naming "Elongation".
#
# False-positive test (AR = 0.1, healthy chain):
#   Override to return a realistic value.  assertChainsMoving() must NOT
#   fire a warning and must return stuck = FALSE.
#
# Unrecognised class test:
#   Pass a plain R list (wrong class).  Must produce a warning about the
#   unrecognised class and return stuck = FALSE (non-breaking behaviour).
# ---------------------------------------------------------------------------

## ---- minimal ROC parameter (no genome needed for trace introspection) ----
fileName <- file.path("UnitTestingData", "testMCMCROCFiles",
                      "simulatedAllUniqueR.fasta")

skip_if_not(file.exists(fileName),
            "ROC unit-test genome not found; skipping assertChainsMoving tests")

genome    <- initializeGenomeObject(file = fileName)
parameter <- initializeParameterObject(
  genome,
  sphi           = 1,
  num.mixtures   = 1,
  gene.assignment = rep(1L, length(genome)),
  model          = "ROC"
)

## Run a minimal MCMC so the trace object is populated with AR entries.
mcmc  <- initializeMCMCObject(samples = 10L, thinning = 1L,
                              adaptive.width = 10L,
                              est.expression = FALSE,
                              est.csp        = TRUE,
                              est.hyper      = FALSE)
model <- initializeModelObject(parameter, "ROC")
suppressMessages(
  suppressWarnings(
    capture.output(runMCMC(mcmc, genome, model, ncores = 1L))
  )
)

## ---- Test 1: false-positive guard (healthy chain, no warning expected) ----
test_that("assertChainsMoving returns stuck=FALSE on a live ROC run", {
  ## After a real (if short) MCMC, at least some AAs should have AR > 0.
  ## Use a very tight threshold (0) to guarantee no false alarm.
  result <- assertChainsMoving(parameter, threshold = 0, action = "warn")
  expect_false(result$stuck)
  expect_length(result$stuck.types, 0L)
  expect_true("Elongation" %in% names(result$ar))
})

## ---- Test 2: true-positive (monkey-patch trace to return AR = 0) ----------
test_that("assertChainsMoving warns and names Elongation when AR=0", {
  ## Save real trace, then monkey-patch.
  trace.real <- parameter$getTraceObject()

  ## Build a fake trace environment with overridden accessor.
  fake.trace <- new.env(parent = emptyenv())
  fake.trace$getCodonSpecificAcceptanceRateTraceForAA <- function(aa) {
    return(c(0.0))   # AR = 0 for every AA
  }

  ## Patch parameter's getTraceObject to return fake.trace.
  ## We do this by wrapping parameter in a local reference object.
  fake.param <- new.env(parent = emptyenv())
  class(fake.param) <- "Rcpp_ROCParameter"
  fake.param$getTraceObject <- function() fake.trace

  expect_warning(
    result <- assertChainsMoving(fake.param, threshold = 0.001, action = "warn"),
    regexp = "STUCK CHAIN"
  )
  expect_true(result$stuck)
  expect_true("Elongation" %in% result$stuck.types)
})

## ---- Test 3: high threshold forces a warning even on a live run -----------
test_that("assertChainsMoving warns when threshold exceeds max observed AR", {
  ## Threshold = 1 guarantees AR < threshold for all groups.
  expect_warning(
    result <- assertChainsMoving(parameter, threshold = 1.0, action = "warn"),
    regexp = "STUCK CHAIN"
  )
  expect_true(result$stuck)
  expect_true("Elongation" %in% result$stuck.types)
})

## ---- Test 4: action="stop" raises an error --------------------------------
test_that("assertChainsMoving stops when action='stop' and chain is stuck", {
  expect_error(
    assertChainsMoving(parameter, threshold = 1.0, action = "stop"),
    regexp = "STUCK CHAIN"
  )
})

## ---- Test 5: context string appears in message ----------------------------
test_that("assertChainsMoving includes context string in message", {
  expect_warning(
    assertChainsMoving(parameter, threshold = 1.0, context = "chunk_test"),
    regexp = "chunk_test"
  )
})

## ---- Test 6: unrecognised class produces a non-breaking warning -----------
test_that("assertChainsMoving warns on unknown class and returns stuck=FALSE", {
  bad.param <- list()
  class(bad.param) <- "SomeUnknownClass"
  expect_warning(
    result <- assertChainsMoving(bad.param),
    regexp = "unrecognised parameter class"
  )
  expect_false(result$stuck)
})
