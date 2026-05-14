library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Regression tests for the init.with.restart.file path of
# initializeParameterObject. Specifically guards against the 2026-05-13
# segfault in which:
#
#  (1) Parameter::initBaseValuesFromFile blindly reads 4 fields from the
#      ">categories:" line (delM, delEta, nse, phi). Legacy restart files
#      (pre-Dec 2022) wrote only 2 fields (delM, delEta). The trailing
#      reads failed silently, leaving K.nse and K.phi at the
#      mixtureDefinition default of -1. getSynthesisRateCategory() then
#      returned UINT_MAX, producing an out-of-bounds index on
#      currentSynthesisRateLevel and a segfault inside the MCMC loop.
#
#  (2) The restart path bypassed initPhiMixtureStorage(), leaving the
#      phi-mixture parameter vectors size-0. Defensive code paths gate
#      on phiPriorType == SINGLE_LN so this didn't crash today, but any
#      future use that doesn't gate would hit OOB.
#
# Both issues are tested below using a freshly-written restart file
# whose categories line is then rewritten to mimic the 2-field legacy
# format.
# ======================================================================

genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
genome      <- initializeGenomeObject(file = genome_file)
gene_assign <- rep(1, length(genome))

# Helper: run a tiny MCMC with restart settings so a .rst is written.
# Returns the path to the final restart file (write.multiple = FALSE
# overwrites a single file, easy to find).
buildRestartFile <- function() {
    set.seed(42)
    parameter <- initializeParameterObject(genome        = genome,
                                            sphi          = 1,
                                            num.mixtures  = 1,
                                            gene.assignment = gene_assign,
                                            model         = "ROC")
    mcmc  <- initializeMCMCObject(samples = 4, thinning = 2,
                                   adaptive.width = 2,
                                   est.expression = TRUE,
                                   est.csp        = TRUE,
                                   est.hyper      = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)

    rst_base <- tempfile(pattern = "restart_test_", fileext = ".rst")
    setRestartSettings(mcmc, filename = rst_base, samples = 2,
                       write.multiple = FALSE)
    sink(tempfile())  # swallow MCMC log output
    runMCMC(mcmc = mcmc, genome = genome, model = model, ncores = 1,
            divergence.iteration = 0)
    sink()
    rst_base
}


test_that("restart with current 4-field categories format works (sanity)", {
    rst <- buildRestartFile()
    expect_true(file.exists(rst))

    # Quick sanity: restart loads, MCMC runs without segfault.
    p2 <- initializeParameterObject(genome = genome,
                                     init.with.restart.file = rst,
                                     model = "ROC")
    m2 <- initializeModelObject(p2, "ROC", with.phi = FALSE)
    mcmc2 <- initializeMCMCObject(samples = 2, thinning = 2,
                                   adaptive.width = 2,
                                   est.expression = TRUE,
                                   est.csp = TRUE, est.hyper = FALSE)
    sink(tempfile())
    runMCMC(mcmc = mcmc2, genome = genome, model = m2,
            ncores = 1, divergence.iteration = 0)
    sink()
    expect_equal(length(mcmc2$getLogPosteriorTrace()), 3L)
    unlink(rst)
})


test_that("restart with legacy 2-field categories format does not segfault", {
    # Build a fresh .rst, then rewrite its ">categories:" block to mimic
    # the pre-Dec-2022 format that wrote only delM and delEta. Without
    # the fallback in initBaseValuesFromFile, K.phi stays at -1 (the
    # mixtureDefinition default), getSynthesisRateCategory returns
    # UINT_MAX, and the MCMC loop segfaults on the first
    # currentSynthesisRateLevel access.
    rst <- buildRestartFile()
    expect_true(file.exists(rst))

    lines <- readLines(rst)
    cat_header <- which(lines == ">categories:")
    expect_length(cat_header, 1L)

    # Find end of the categories block (next line starting with ">").
    next_section <- which(grepl("^>", lines[seq_len(length(lines)) >
                                              cat_header]))[1] + cat_header
    expect_true(!is.na(next_section) && next_section > cat_header)

    # Sanity-check that the current writer emits 4 fields per category.
    cat_lines <- lines[(cat_header + 1L):(next_section - 1L)]
    cat_lines <- cat_lines[nzchar(cat_lines)]
    expect_true(all(lengths(strsplit(trimws(cat_lines), "\\s+")) == 4L))

    # Rewrite each category line to keep only the first two fields.
    legacy_cat_lines <- vapply(cat_lines, function(l) {
        parts <- strsplit(trimws(l), "\\s+")[[1]]
        paste(parts[1:2], collapse = " ")
    }, character(1), USE.NAMES = FALSE)

    legacy <- c(lines[seq_len(cat_header)],
                legacy_cat_lines,
                lines[next_section:length(lines)])
    legacy_rst <- tempfile(pattern = "restart_test_legacy_", fileext = ".rst")
    on.exit(unlink(c(rst, legacy_rst)), add = TRUE)
    writeLines(legacy, legacy_rst)

    # Load + run MCMC. With the fix, this completes; without it, segfault.
    p3 <- initializeParameterObject(genome = genome,
                                     init.with.restart.file = legacy_rst,
                                     model = "ROC")

    # Fix #1 verification: K.phi was loaded (or back-filled) such that
    # getSynthesisRateCategoryForMixture returns a valid index. R-side
    # mixtures are 1-based; 0 means category 0 (the only one for a
    # 1-mixture run), and the buggy state would surface as either an
    # error or a value congruent to (UINT_MAX %% something).
    expect_equal(p3$getSynthesisRateCategoryForMixture(1), 1L)

    # Fix #2 verification: phi-mixture storage is allocated on the
    # restart path. Setters bounds-check against the vector size and
    # return NaN when the vector is size-0; if storage is allocated
    # they round-trip the default value seeded by initPhiMixtureStorage.
    expect_false(is.nan(p3$getPhiMixtureP(0L, FALSE)))
    expect_equal(p3$getPhiMixtureP(0L, FALSE), 0.9, tolerance = 1e-12)

    m3 <- initializeModelObject(p3, "ROC", with.phi = FALSE)
    mcmc3 <- initializeMCMCObject(samples = 2, thinning = 2,
                                   adaptive.width = 2,
                                   est.expression = TRUE,
                                   est.csp = TRUE, est.hyper = FALSE)
    sink(tempfile())
    runMCMC(mcmc = mcmc3, genome = genome, model = m3,
            ncores = 1, divergence.iteration = 0)
    sink()
    expect_equal(length(mcmc3$getLogPosteriorTrace()), 3L)
})
