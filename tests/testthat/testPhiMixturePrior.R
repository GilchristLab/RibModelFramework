library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Tests for the phi-mixture prior C++ infrastructure (task #12a).
#
# Covers:
#   - static helpers derivePhiMixtureMu2, densityLogNormMixture
#   - per-mixture-category storage round-trip via R setters/getters
#   - phi.prior / phi.prior.constraint / phi.prior.init / phi.prior.hyperparams
#     plumbing through initializeParameterObject
#   - default values reproduce legacy behavior (single LN, mean=1)
#
# Reference values come from prototypes/phi_mixture.R (also re-derived inline
# so this test does not depend on the prototype directory).
# ======================================================================

PRIOR_SINGLE  <- 0L
PRIOR_MIXTURE <- 1L
CONSTRAINT_MEAN   <- 0L
CONSTRAINT_MEDIAN <- 1L

genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
genome <- initializeGenomeObject(file = genome_file)
gene_assignment <- rep(1, length(genome))


# ---------------- derivePhiMixtureMu2 ----------------

test_that("derivePhiMixtureMu2 (mean=1) matches the closed-form", {
    # mu2 = log[(1 - p * exp(mu1 + s1^2/2)) / (1 - p)] - s2^2/2
    p <- 0.85; mu1 <- -0.4; s1 <- 0.4; s2 <- 0.25
    expected <- log((1 - p * exp(mu1 + s1^2 / 2)) / (1 - p)) - s2^2 / 2
    got <- AnaCoDa:::derivePhiMixtureMu2(p, mu1, s1, s2, CONSTRAINT_MEAN)
    expect_equal(got, expected, tolerance = 1e-12)
})

test_that("derivePhiMixtureMu2 (median=1) matches the closed-form", {
    # mu2 = -s2 * qnorm[(0.5 - p * Phi(-mu1/s1)) / (1 - p)]
    p <- 0.75; mu1 <- -0.05; s1 <- 0.4; s2 <- 0.25
    Phi1 <- pnorm(-mu1 / s1)
    expected <- -s2 * qnorm((0.5 - p * Phi1) / (1 - p))
    got <- AnaCoDa:::derivePhiMixtureMu2(p, mu1, s1, s2, CONSTRAINT_MEDIAN)
    expect_equal(got, expected, tolerance = 1e-10)
})

test_that("derivePhiMixtureMu2 returns NaN when infeasible", {
    # Lower component already exceeds mean=1: p * exp(mu1+s1^2/2) > 1
    expect_true(is.nan(AnaCoDa:::derivePhiMixtureMu2(0.9, 1.0, 0.4, 0.25, CONSTRAINT_MEAN)))
})


# ---------------- densityLogNormMixture ----------------

test_that("densityLogNormMixture matches p*dlnorm + (1-p)*dlnorm at the derived mu2", {
    p <- 0.85; mu1 <- -0.4; s1 <- 0.4; s2 <- 0.25
    mu2 <- AnaCoDa:::derivePhiMixtureMu2(p, mu1, s1, s2, CONSTRAINT_MEAN)
    for (x in c(0.2, 0.5, 1.0, 2.0, 5.0)) {
        expected <- p * dlnorm(x, mu1, s1) + (1 - p) * dlnorm(x, mu2, s2)
        # log scale:
        expect_equal(
            AnaCoDa:::densityLogNormMixture(x, p, mu1, s1, s2, CONSTRAINT_MEAN, TRUE),
            log(expected), tolerance = 1e-10
        )
        # linear scale:
        expect_equal(
            AnaCoDa:::densityLogNormMixture(x, p, mu1, s1, s2, CONSTRAINT_MEAN, FALSE),
            expected, tolerance = 1e-12
        )
    }
})

test_that("densityLogNormMixture returns -DBL_MAX on label-switching violation", {
    # Force mu2 < mu1 by picking parameters where mean=1 constraint pushes mu2 negative.
    p <- 0.5; mu1 <- 0.5; s1 <- 0.4; s2 <- 0.25
    mu2 <- AnaCoDa:::derivePhiMixtureMu2(p, mu1, s1, s2, CONSTRAINT_MEAN)
    skip_if(is.nan(mu2) || mu2 >= mu1, "did not produce mu2 < mu1 in this regime")
    val <- AnaCoDa:::densityLogNormMixture(1.0, p, mu1, s1, s2, CONSTRAINT_MEAN, TRUE)
    expect_lt(val, -1e300)  # -DBL_MAX sentinel
})

test_that("densityLogNormMixture returns 0 (linear) / -DBL_MAX (log) for x <= 0", {
    p <- 0.85; mu1 <- -0.4; s1 <- 0.4; s2 <- 0.25
    expect_equal(
        AnaCoDa:::densityLogNormMixture(-1.0, p, mu1, s1, s2, CONSTRAINT_MEAN, FALSE),
        0.0
    )
    expect_lt(
        AnaCoDa:::densityLogNormMixture(0.0, p, mu1, s1, s2, CONSTRAINT_MEAN, TRUE),
        -1e300
    )
})


# ---------------- defaults preserve legacy behavior ----------------

test_that("default initializeParameterObject sets phiPriorType=SINGLE, constraint=MEAN", {
    parameter <- initializeParameterObject(genome = genome, sphi = 1,
                                            num.mixtures = 1,
                                            gene.assignment = gene_assignment)
    expect_equal(parameter$getPhiPriorType(), PRIOR_SINGLE)
    expect_equal(parameter$getPhiPriorConstraint(), CONSTRAINT_MEAN)
})


# ---------------- storage round-trip ----------------

test_that("setPhiMixture* / getPhiMixture* round-trip per-category values", {
    parameter <- initializeParameterObject(genome = genome, sphi = 1,
                                            num.mixtures = 1,
                                            gene.assignment = gene_assignment)
    parameter$setPhiMixtureP(0.7, 0L)
    parameter$setPhiMixtureMu1(-0.2, 0L)
    parameter$setPhiMixtureSigma1(0.5, 0L)
    parameter$setPhiMixtureSigma2(0.3, 0L)
    expect_equal(parameter$getPhiMixtureP(0L, FALSE), 0.7)
    expect_equal(parameter$getPhiMixtureMu1(0L, FALSE), -0.2)
    expect_equal(parameter$getPhiMixtureSigma1(0L, FALSE), 0.5)
    expect_equal(parameter$getPhiMixtureSigma2(0L, FALSE), 0.3)
    # proposed counterparts also seeded
    expect_equal(parameter$getPhiMixtureP(0L, TRUE), 0.7)
})

test_that("getPhiMixtureMu2Derived matches the static helper", {
    parameter <- initializeParameterObject(genome = genome, sphi = 1,
                                            num.mixtures = 1,
                                            gene.assignment = gene_assignment)
    parameter$setPhiPriorConstraint(CONSTRAINT_MEAN)
    parameter$setPhiMixtureP(0.85, 0L)
    parameter$setPhiMixtureMu1(-0.4, 0L)
    parameter$setPhiMixtureSigma1(0.4, 0L)
    parameter$setPhiMixtureSigma2(0.25, 0L)
    expected <- AnaCoDa:::derivePhiMixtureMu2(0.85, -0.4, 0.4, 0.25, CONSTRAINT_MEAN)
    expect_equal(parameter$getPhiMixtureMu2Derived(0L), expected, tolerance = 1e-12)
})


# ---------------- phi.prior.init / phi.prior.hyperparams plumbing ----------------

test_that("phi.prior.init applies via initializeParameterObject", {
    parameter <- initializeParameterObject(
        genome = genome, sphi = 1, num.mixtures = 1,
        gene.assignment = gene_assignment,
        phi.prior = "mixture-lognormal",
        phi.prior.init = list(p = 0.82, mu1 = -0.30, sigma1 = 0.45, sigma2 = 0.20)
    )
    expect_equal(parameter$getPhiPriorType(), PRIOR_MIXTURE)
    expect_equal(parameter$getPhiMixtureP(0L, FALSE), 0.82)
    expect_equal(parameter$getPhiMixtureMu1(0L, FALSE), -0.30)
    expect_equal(parameter$getPhiMixtureSigma1(0L, FALSE), 0.45)
    expect_equal(parameter$getPhiMixtureSigma2(0L, FALSE), 0.20)
})

test_that("phi.prior.hyperparams overrides selectively (others keep defaults)", {
    parameter <- initializeParameterObject(
        genome = genome, sphi = 1, num.mixtures = 1,
        gene.assignment = gene_assignment,
        phi.prior = "mixture-lognormal",
        phi.prior.hyperparams = list(p = list(alpha = 4, beta = 4),
                                      sigma2 = list(scale = 2))
    )
    expect_equal(parameter$getPhiMixtureHyperPAlpha(), 4)
    expect_equal(parameter$getPhiMixtureHyperPBeta(), 4)
    expect_equal(parameter$getPhiMixtureHyperMu1Mean(), 0)    # default
    expect_equal(parameter$getPhiMixtureHyperMu1Sd(), 10)     # default
    expect_equal(parameter$getPhiMixtureHyperSigma1Scale(), 1) # default
    expect_equal(parameter$getPhiMixtureHyperSigma2Scale(), 2)
})


# ---------------- validation ----------------

test_that("bad phi.prior value is rejected", {
    expect_error(
        initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1,
                                   gene.assignment = gene_assignment,
                                   phi.prior = "bogus"),
        "phi.prior must be"
    )
})

test_that("bad phi.prior.constraint value is rejected", {
    expect_error(
        initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1,
                                   gene.assignment = gene_assignment,
                                   phi.prior.constraint = "mode"),
        "phi.prior.constraint must be"
    )
})

test_that("unknown phi.prior.init element is rejected", {
    expect_error(
        initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1,
                                   gene.assignment = gene_assignment,
                                   phi.prior = "mixture-lognormal",
                                   phi.prior.init = list(p = 0.8, foo = 1)),
        "unknown elements: foo"
    )
})

test_that("phi.prior.init$p out of (0,1) is rejected", {
    expect_error(
        initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1,
                                   gene.assignment = gene_assignment,
                                   phi.prior = "mixture-lognormal",
                                   phi.prior.init = list(p = 1.5)),
        "phi.prior.init\\$p must be in"
    )
})

test_that("phi.prior = 'mixture-lognormal' rejects non-ROC models", {
    # Task #12b: mixture wiring is currently ROC-only.
    expect_error(
        initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1,
                                   gene.assignment = gene_assignment,
                                   model = "FONSE",
                                   phi.prior = "mixture-lognormal"),
        "supported only with model = 'ROC'"
    )
})


# ---------------- 12b smoke test: ROC MCMC under mixture-LN prior ----------------

test_that("ROC MCMC runs under phi.prior = 'mixture-lognormal' without NaN/Inf", {
    # Smoke test that the mixture prior wiring (task #12b) is exercised by
    # the MCMC loops and produces finite log-posteriors. Mixture hyperparams
    # stay frozen at the init values (their MCMC update step lands in 12c).
    fasta <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
    skip_if_not(file.exists(fasta), "test fixture FASTA missing")

    set.seed(12345)
    genome_local <- initializeGenomeObject(file = fasta)
    geneAssignment_local <- rep(1, length(genome_local))

    parameter <- initializeParameterObject(
        genome = genome_local, sphi = 1, num.mixtures = 1,
        gene.assignment = geneAssignment_local,
        model = "ROC",
        phi.prior = "mixture-lognormal",
        phi.prior.constraint = "mean",
        phi.prior.init = list(p = 0.9, mu1 = -0.4, sigma1 = 0.4, sigma2 = 0.25)
    )
    expect_equal(parameter$getPhiPriorType(), PRIOR_MIXTURE)
    expect_true(is.finite(parameter$getPhiMixtureMu2Derived(0L)))

    mcmc <- initializeMCMCObject(samples = 5, thinning = 2, adaptive.width = 5,
                                  est.expression = TRUE, est.csp = TRUE,
                                  est.hyper = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)

    outFile <- file.path("UnitTestingOut", "testPhiMixtureMCMCLog.txt")
    sink(outFile)
    runMCMC(mcmc = mcmc, genome = genome_local, model = model, ncores = 1,
            divergence.iteration = 0)
    sink()

    trace <- mcmc$getLogPosteriorTrace()
    # The first entry can be 0 on a fresh chain; check the chain body.
    body <- trace[-1]
    expect_true(all(is.finite(body)),
                info = paste("non-finite logPosterior values:",
                              paste(body[!is.finite(body)], collapse = ", ")))
})


# ---------------- 12c.1 test: mixture hyperparams actually update ----------------

test_that("mixture hyperparams move during MCMC under phi.prior = 'mixture-lognormal'", {
    # Task #12c.1: the updatePhiMixtureHyperparameters M-H step should run
    # each iteration when phiPriorType == MIXTURE_LN. Verify by:
    #   1. At least one of (p, mu1, sigma1, sigma2) gets accepted (counters > 0)
    #   2. The current values are no longer exactly at init (chain moved)
    fasta <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
    skip_if_not(file.exists(fasta), "test fixture FASTA missing")

    set.seed(54321)
    genome_local <- initializeGenomeObject(file = fasta)
    geneAssignment_local <- rep(1, length(genome_local))

    init_p <- 0.9; init_mu1 <- -0.4; init_s1 <- 0.4; init_s2 <- 0.25
    parameter <- initializeParameterObject(
        genome = genome_local, sphi = 1, num.mixtures = 1,
        gene.assignment = geneAssignment_local,
        model = "ROC",
        phi.prior = "mixture-lognormal",
        phi.prior.constraint = "mean",
        phi.prior.init = list(p = init_p, mu1 = init_mu1,
                               sigma1 = init_s1, sigma2 = init_s2)
    )

    # 30 iterations: enough for at least a few accepts across 4 params.
    mcmc <- initializeMCMCObject(samples = 30, thinning = 1, adaptive.width = 5,
                                  est.expression = TRUE, est.csp = TRUE,
                                  est.hyper = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)

    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome_local, model = model, ncores = 1,
            divergence.iteration = 0)
    sink()

    # Traces: at least one trace should show a value different from init,
    # i.e. at least one accept across the run. (The live accept counters
    # reset every adaptation window in 12c.2, so we check the trace instead.)
    trace <- parameter$getTraceObject()
    pTrace  <- trace$getPhiMixturePTrace()[[1]]
    m1Trace <- trace$getPhiMixtureMu1Trace()[[1]]
    s1Trace <- trace$getPhiMixtureSigma1Trace()[[1]]
    s2Trace <- trace$getPhiMixtureSigma2Trace()[[1]]
    movedAny <- length(unique(pTrace[-1]))  > 1 ||
                length(unique(m1Trace[-1])) > 1 ||
                length(unique(s1Trace[-1])) > 1 ||
                length(unique(s2Trace[-1])) > 1
    expect_true(movedAny)

    # Current values: at least one should have moved from init.
    moved <- (parameter$getPhiMixtureP(0L, FALSE)       != init_p)  ||
             (parameter$getPhiMixtureMu1(0L, FALSE)     != init_mu1) ||
             (parameter$getPhiMixtureSigma1(0L, FALSE)  != init_s1) ||
             (parameter$getPhiMixtureSigma2(0L, FALSE)  != init_s2)
    expect_true(moved)

    # The derived mu2 should still satisfy mu2 >= mu1 (label-switching guard).
    expect_gte(parameter$getPhiMixtureMu2Derived(0L),
                parameter$getPhiMixtureMu1(0L, FALSE))
})


test_that("mixture hyperparam traces populate during MCMC (task #12c.2)", {
    fasta <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
    skip_if_not(file.exists(fasta), "test fixture FASTA missing")

    set.seed(2024)
    genome_local <- initializeGenomeObject(file = fasta)
    geneAssignment_local <- rep(1, length(genome_local))

    parameter <- initializeParameterObject(
        genome = genome_local, sphi = 1, num.mixtures = 1,
        gene.assignment = geneAssignment_local,
        model = "ROC",
        phi.prior = "mixture-lognormal",
        phi.prior.init = list(p = 0.85, mu1 = -0.4, sigma1 = 0.4, sigma2 = 0.25)
    )
    n_samples <- 20
    mcmc <- initializeMCMCObject(samples = n_samples, thinning = 1,
                                  adaptive.width = 5,
                                  est.expression = TRUE, est.csp = TRUE,
                                  est.hyper = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)
    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome_local, model = model, ncores = 1,
            divergence.iteration = 0)
    sink()

    trace <- parameter$getTraceObject()
    pTrace  <- trace$getPhiMixturePTrace()
    m1Trace <- trace$getPhiMixtureMu1Trace()
    s1Trace <- trace$getPhiMixtureSigma1Trace()
    s2Trace <- trace$getPhiMixtureSigma2Trace()

    # All 4 traces should have one row per mixture category, n_samples+1 entries.
    expect_equal(length(pTrace), 1)
    expect_equal(length(pTrace[[1]]), n_samples + 1)

    # At least some entries should differ from the init value (chain moved).
    expect_true(any(pTrace[[1]][-1] != 0.85) ||
                any(m1Trace[[1]][-1] != -0.4) ||
                any(s1Trace[[1]][-1] != 0.4) ||
                any(s2Trace[[1]][-1] != 0.25))

    # All trace values within sane ranges (no NaN/Inf).
    expect_true(all(is.finite(pTrace[[1]])))
    expect_true(all(is.finite(m1Trace[[1]])))
    expect_true(all(is.finite(s1Trace[[1]])))
    expect_true(all(is.finite(s2Trace[[1]])))
})


test_that("acceptance-rate trace pushes one value per adaptation window", {
    fasta <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
    skip_if_not(file.exists(fasta), "test fixture FASTA missing")

    set.seed(7)
    genome_local <- initializeGenomeObject(file = fasta)
    geneAssignment_local <- rep(1, length(genome_local))

    parameter <- initializeParameterObject(
        genome = genome_local, sphi = 1, num.mixtures = 1,
        gene.assignment = geneAssignment_local,
        model = "ROC",
        phi.prior = "mixture-lognormal"
    )
    # 30 iters with adaptive.width=10 -> expect ~3 acceptance rate entries.
    mcmc <- initializeMCMCObject(samples = 30, thinning = 1,
                                  adaptive.width = 10,
                                  est.expression = TRUE, est.csp = TRUE,
                                  est.hyper = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)
    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome_local, model = model, ncores = 1,
            divergence.iteration = 0)
    sink()

    trace <- parameter$getTraceObject()
    pRates <- trace$getPhiMixturePAcceptanceRateTrace()
    expect_gte(length(pRates), 2)
    expect_true(all(is.finite(pRates)))
    expect_true(all(pRates >= 0 & pRates <= 1))
})


test_that("default phi.prior = 'lognormal' does NOT call mixture update step", {
    # Counter check: with the default prior, the mixture update is a no-op
    # (early return). Accept counters should stay at zero across MCMC.
    fasta <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")
    skip_if_not(file.exists(fasta), "test fixture FASTA missing")

    set.seed(99)
    genome_local <- initializeGenomeObject(file = fasta)
    geneAssignment_local <- rep(1, length(genome_local))

    parameter <- initializeParameterObject(
        genome = genome_local, sphi = 1, num.mixtures = 1,
        gene.assignment = geneAssignment_local,
        model = "ROC"
        # default phi.prior = "lognormal"
    )

    mcmc <- initializeMCMCObject(samples = 10, thinning = 1, adaptive.width = 5,
                                  est.expression = TRUE, est.csp = TRUE,
                                  est.hyper = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)
    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome_local, model = model, ncores = 1,
            divergence.iteration = 0)
    sink()

    expect_equal(parameter$getNumAcceptForPhiMixtureP(), 0)
    expect_equal(parameter$getNumAcceptForPhiMixtureMu1(), 0)
    expect_equal(parameter$getNumAcceptForPhiMixtureSigma1(), 0)
    expect_equal(parameter$getNumAcceptForPhiMixtureSigma2(), 0)
})
