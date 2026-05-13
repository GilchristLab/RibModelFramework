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
