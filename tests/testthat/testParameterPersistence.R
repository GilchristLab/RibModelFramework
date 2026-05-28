library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Regression tests for parameter persistence -- writeParameterObject /
# loadParameterObject must preserve all C++ state that is consulted by
# downstream accessors. Specifically guards against the 2026-05-13
# regression in which numElongationMixtures was added to checkIndex()
# bounds in src/Parameter.cpp (getCodonSpecificPosteriorMean and friends)
# but was not persisted by extractBaseInfo/setBaseInfo, leaving every
# loaded parameter object with numElongationMixtures = 0 and every
# accessor failing with "Index 1 is out of bounds. Index must be
# between 1 & 0".
# ======================================================================

genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
genome <- initializeGenomeObject(file = genome_file)
gene_assignment <- rep(1, length(genome))

# Helper: build a parameter object with traces populated by a tiny MCMC.
# writeParameterObject / loadParameterObject indexes into the traces during
# I/O, so they require at least one MCMC sample.
buildSavedParameter <- function() {
    set.seed(42)
    parameter <- initializeParameterObject(genome = genome, sphi = 1,
                                            num.mixtures = 1,
                                            gene.assignment = gene_assignment)
    mcmc <- initializeMCMCObject(samples = 2, thinning = 2, adaptive.width = 2,
                                  est.expression = TRUE, est.csp = TRUE,
                                  est.hyper = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)
    sink(tempfile())  # swallow MCMC log output
    runMCMC(mcmc = mcmc, genome = genome, model = model, ncores = 1,
            divergence.iteration = 0)
    sink()
    parameter
}

tmpFile <- tempfile(fileext = ".Rdata")
on.exit(unlink(tmpFile), add = TRUE)


test_that("numElongationMixtures round-trips through write/load", {
    parameter <- buildSavedParameter()
    expect_equal(parameter$numMixtures, 1)
    expect_equal(parameter$numElongationMixtures, 1)

    writeParameterObject(parameter, tmpFile)
    loaded <- loadParameterObject(tmpFile)

    expect_equal(loaded$numMixtures, 1)
    expect_equal(loaded$numElongationMixtures, 1)
})


test_that("loadParameterObject falls back to numMixtures when paramBase lacks numElongMix", {
    # Simulate an "old" .Rdata file (pre-2026-05-13) that does not contain
    # numElongMix in paramBase. The loader should restore numElongationMixtures
    # to numMixtures rather than leaving it at the C++ default of 0.
    parameter <- buildSavedParameter()
    writeParameterObject(parameter, tmpFile)

    # Reach into the saved file and strip numElongMix from paramBase.
    e <- new.env()
    load(tmpFile, envir = e)
    expect_true("numElongMix" %in% names(e$paramBase))
    e$paramBase$numElongMix <- NULL
    save(list = ls(envir = e), envir = e, file = tmpFile)

    loaded <- loadParameterObject(tmpFile)
    expect_equal(loaded$numElongationMixtures, loaded$numMixtures)
})


test_that("getCodonSpecificPosteriorMean works after loadParameterObject", {
    # The actual surface area broken by the 2026-05-13 regression. Without
    # numElongationMixtures persistence, this accessor would fail with
    # "Index 1 is out of bounds. Index must be between 1 & 0".
    parameter <- buildSavedParameter()
    writeParameterObject(parameter, tmpFile)
    loaded <- loadParameterObject(tmpFile)

    expect_silent({
        val <- loaded$getCodonSpecificPosteriorMean(
            mixtureElement = 1L, samples = 1L, codon = "GCA",
            paramType = 0L, withoutReference = TRUE, log_scale = FALSE
        )
    })
    expect_true(is.finite(val))
})


# ======================================================================
# Regression tests for multi-file loadParameterObject (acope3 #424,
# filed 2026-05-28).
#
# combineTwoDimensionalTrace / combineThreeDimensionalTrace had two bugs:
#   1. combineThreeDimensionalTrace did not return its result, so the
#      caller's binding was never updated.
#   2. setBaseInfo callers passed the trace length positionally as the
#      3rd arg (start), instead of as the named 4th arg (end), and
#      then discarded the return value.
# Net effect: loadParameterObject(c(file1, file2, ...)) bumped
# lastIteration correctly but silently dropped all but the first
# file's trace data.
# ======================================================================

# Helper: write two parameter Rdata files with traces produced by
# successive runMCMC calls on the same parameter object.  The parameter
# object's C++ trace is overwritten by each runMCMC call, so file2
# captures only the *second* run's trace, not the cumulative trace.
# Useful for verifying multi-file load concatenation behavior.
buildTwoSavedParameters <- function() {
    set.seed(43)
    parameter <- initializeParameterObject(genome = genome, sphi = 1,
                                            num.mixtures = 1,
                                            gene.assignment = gene_assignment)
    mcmc <- initializeMCMCObject(samples = 3, thinning = 2, adaptive.width = 2,
                                  est.expression = TRUE, est.csp = TRUE,
                                  est.hyper = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)
    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome, model = model, ncores = 1,
            divergence.iteration = 0)
    sink()
    f1 <- tempfile(fileext = ".Rdata")
    writeParameterObject(parameter, f1)

    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome, model = model, ncores = 1,
            divergence.iteration = 0)
    sink()
    f2 <- tempfile(fileext = ".Rdata")
    writeParameterObject(parameter, f2)

    list(f1 = f1, f2 = f2)
}


test_that("loadParameterObject(c(f1, f2)) concatenates sphi trace", {
    files <- buildTwoSavedParameters()
    on.exit(unlink(c(files$f1, files$f2)), add = TRUE)

    p1  <- loadParameterObject(files$f1)
    p2  <- loadParameterObject(files$f2)
    p12 <- loadParameterObject(c(files$f1, files$f2))

    sphi1  <- as.numeric(p1$getTraceObject()$getStdDevSynthesisRateTraces()[[1]])
    sphi2  <- as.numeric(p2$getTraceObject()$getStdDevSynthesisRateTraces()[[1]])
    sphi12 <- as.numeric(p12$getTraceObject()$getStdDevSynthesisRateTraces()[[1]])

    # Concatenated length = sum, minus 1 for the overlap-skip slot.
    expect_equal(length(sphi12), length(sphi1) + length(sphi2) - 1)
    # First value comes from file1; last from file2.
    expect_equal(sphi12[1],            sphi1[1])
    expect_equal(tail(sphi12, 1),      tail(sphi2, 1))
})


test_that("loadParameterObject(c(f1, f2)) concatenates CSP traces", {
    files <- buildTwoSavedParameters()
    on.exit(unlink(c(files$f1, files$f2)), add = TRUE)

    p1  <- loadParameterObject(files$f1)
    p2  <- loadParameterObject(files$f2)
    p12 <- loadParameterObject(c(files$f1, files$f2))

    # Mutation trace (paramType = 0), category 1, codon index 2.
    mut1  <- p1$getTraceObject()$getCodonSpecificParameterTrace(0)[[1]][[2]]
    mut2  <- p2$getTraceObject()$getCodonSpecificParameterTrace(0)[[1]][[2]]
    mut12 <- p12$getTraceObject()$getCodonSpecificParameterTrace(0)[[1]][[2]]

    expect_equal(length(mut12), length(mut1) + length(mut2) - 1)
    expect_equal(mut12[1],         mut1[1])
    expect_equal(tail(mut12, 1),   tail(mut2, 1))

    # Selection trace (paramType = 1), same pattern.
    sel12 <- p12$getTraceObject()$getCodonSpecificParameterTrace(1)[[1]][[2]]
    sel1  <- p1$getTraceObject()$getCodonSpecificParameterTrace(1)[[1]][[2]]
    sel2  <- p2$getTraceObject()$getCodonSpecificParameterTrace(1)[[1]][[2]]
    expect_equal(length(sel12), length(sel1) + length(sel2) - 1)
    expect_equal(sel12[1],         sel1[1])
    expect_equal(tail(sel12, 1),   tail(sel2, 1))
})


test_that("loadParameterObject(c(f1, f2)) sets lastIteration = sum", {
    files <- buildTwoSavedParameters()
    on.exit(unlink(c(files$f1, files$f2)), add = TRUE)

    p1  <- loadParameterObject(files$f1)
    p2  <- loadParameterObject(files$f2)
    p12 <- loadParameterObject(c(files$f1, files$f2))

    expect_equal(p12$getLastIteration(),
                 p1$getLastIteration() + p2$getLastIteration())
})


# ======================================================================
# Regression test for loadMCMCObject first-element drop (acope3 #388,
# filed 2022-06-23).
#
# Pre-fix: loadMCMCObject used `curLogPostTrace[2:max]` for all files,
# including the first.  Single-file load returned a trace with the
# initial-evaluation slot at index 1 silently removed.
# ======================================================================

test_that("loadMCMCObject preserves the first element on single-file load", {
    # Construct an MCMC save in the writeMCMCObject convention:
    # logPostTrace length = samples + 1, where index 1 is the initial 0.0
    # evaluation slot.
    logPostTrace  <- c(0, -100, -90, -85)
    logLikeTrace  <- c(0, -200, -180, -170)
    samples       <- 3L
    thinning      <- 1L
    adaptiveWidth <- 10L
    tf <- tempfile(fileext = ".Rda")
    on.exit(unlink(tf), add = TRUE)
    save(logPostTrace, logLikeTrace, samples, thinning, adaptiveWidth,
         file = tf)

    m  <- loadMCMCObject(tf)
    tr <- m$getLogPosteriorTrace()
    expect_equal(length(tr), length(logPostTrace))
    expect_equal(tr[1], logPostTrace[1])
    expect_equal(tail(tr, 1), tail(logPostTrace, 1))
})


test_that("loadMCMCObject(c(f1, f2)) concatenates with one overlap-skip", {
    logPostTrace  <- c(0, -100, -90, -85, -80, -78)
    logLikeTrace  <- c(0, -200, -180, -170, -160, -155)
    samples       <- 5L
    thinning      <- 1L
    adaptiveWidth <- 10L
    tf1 <- tempfile(fileext = ".Rda")
    save(logPostTrace, logLikeTrace, samples, thinning, adaptiveWidth,
         file = tf1)

    logPostTrace  <- c(-78, -75, -72, -70)   # first elem duplicates tf1's last
    logLikeTrace  <- c(-155, -150, -147, -145)
    samples       <- 3L
    tf2 <- tempfile(fileext = ".Rda")
    save(logPostTrace, logLikeTrace, samples, thinning, adaptiveWidth,
         file = tf2)
    on.exit(unlink(c(tf1, tf2)), add = TRUE)

    m  <- loadMCMCObject(c(tf1, tf2))
    tr <- m$getLogPosteriorTrace()
    # 6 + 4 - 1 overlap-skip = 9
    expect_equal(length(tr), 9L)
    expect_equal(tr[1],            0)
    expect_equal(tail(tr, 1),      -70)
    # No duplicated overlap value.
    expect_equal(tr[6], -78)
    expect_equal(tr[7], -75)
})
