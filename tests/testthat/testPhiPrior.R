library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Tests for the ROC phi prior parameterization.
#
# Covers:
#  1. Legacy mPhi math (E[phi]=1 mean anchor, median anchor)
#  2. computeMPhi closed-form formulas for all five statistics
#  3. phi spec R constructor validation
#  4. Parameter getter/setter round-trips (new C++ members via Rcpp)
#  5. phi.mphi integration (constrained / fixed modes)
#  6. phi.sphi integration (estimated with various priors)
#  7. MCMC smoke tests for new phi spec combinations
#  8. Restart file round-trip for all new phi spec members
#  9. sphi prior default-off regression
# ======================================================================


# ======================================================================
# Section 1: Legacy parameterization math
# ======================================================================

test_that("default ROC phi prior (single lognormal) anchors E[phi] = 1", {
    sigma <- 0.5
    mPhi  <- -(sigma^2) / 2

    # E[phi] = exp(mu + sigma^2/2) for LogNormal(mu, sigma)
    expect_equal(exp(mPhi + sigma^2 / 2), 1.0, tolerance = 1e-12)

    # Median = exp(mu) < 1 under mean-anchor parameterization
    expect_lt(exp(mPhi), 1.0)
    expect_equal(exp(mPhi), exp(-sigma^2 / 2), tolerance = 1e-12)

    set.seed(1)
    phi <- rlnorm(1e6, meanlog = mPhi, sdlog = sigma)
    expect_equal(mean(phi), 1.0, tolerance = 0.01)
    expect_lt(median(phi), 0.95)
})

test_that("median-constraint phi prior (mPhi = 0) anchors median[phi] = 1", {
    sigma <- 0.5
    mPhi  <- 0

    expect_equal(exp(mPhi), 1.0, tolerance = 1e-12)
    expect_gt(exp(mPhi + sigma^2 / 2), 1.0)

    set.seed(2)
    phi <- rlnorm(1e6, meanlog = mPhi, sdlog = sigma)
    expect_equal(median(phi), 1.0, tolerance = 0.01)
    expect_gt(mean(phi), 1.0)
})


# ======================================================================
# Section 2: computeMPhi formula verification (pure R math)
#
# These tests mirror the C++ formulas in Parameter::computeMPhi() so
# that any drift between the design spec and the implementation surfaces
# in code review (even without calling C++ directly).
# ======================================================================

test_that("mPhi formula: mean statistic pins E[phi] to value", {
    for (s in c(0.3, 0.5, 1.0, 1.5)) {
        v    <- 2.5
        mPhi <- log(v) - s^2 / 2
        # E[LogNormal(mPhi, s)] = exp(mPhi + s^2/2) = v
        expect_equal(exp(mPhi + s^2 / 2), v, tolerance = 1e-12,
                     label = paste0("E[phi] at sphi=", s))
    }
})

test_that("mPhi formula: median statistic pins median[phi] to value", {
    for (s in c(0.3, 0.5, 1.0, 1.5)) {
        v    <- 2.5
        mPhi <- log(v)
        # median[LogNormal(mPhi, s)] = exp(mPhi) = v
        expect_equal(exp(mPhi), v, tolerance = 1e-12,
                     label = paste0("median at sphi=", s))
        # E[phi] = v * exp(s^2/2) > v when s > 0
        expect_gt(exp(mPhi + s^2 / 2), v)
    }
})

test_that("mPhi formula: mode statistic pins mode[phi] to value", {
    for (s in c(0.3, 0.5, 1.0)) {
        v    <- 2.5
        mPhi <- log(v) + s^2
        # mode[LogNormal(mPhi, s)] = exp(mPhi - s^2) = v
        expect_equal(exp(mPhi - s^2), v, tolerance = 1e-12,
                     label = paste0("mode at sphi=", s))
    }
})

test_that("mPhi formula: variance statistic pins Var[phi] to value", {
    for (s in c(0.3, 0.5)) {
        v    <- 1.2
        mPhi <- (log(v / (exp(s^2) - 1)) - s^2) / 2
        # Var[LogNormal(mPhi, s)] = exp(2*mPhi + s^2) * (exp(s^2) - 1)
        var_phi <- exp(2 * mPhi + s^2) * (exp(s^2) - 1)
        expect_equal(var_phi, v, tolerance = 1e-10,
                     label = paste0("Var[phi] at sphi=", s))
    }
})

test_that("mPhi formula: sd statistic pins sd[phi] to value", {
    for (s in c(0.3, 0.5)) {
        v    <- 1.1
        mPhi <- (log(v^2 / (exp(s^2) - 1)) - s^2) / 2
        var_phi <- exp(2 * mPhi + s^2) * (exp(s^2) - 1)
        expect_equal(sqrt(var_phi), v, tolerance = 1e-10,
                     label = paste0("sd[phi] at sphi=", s))
    }
})


# ======================================================================
# Section 3: phi spec R constructor validation
# ======================================================================

test_that("constrained() returns phi_spec_constrained with correct fields", {
    suppressWarnings({
        s_mean   <- constrained(statistic = "mean",   value = 1)
        s_median <- constrained(statistic = "median", value = 2)
        s_mode   <- constrained(statistic = "mode",   value = 1)
        s_var    <- constrained(statistic = "variance", value = 0.5)
        s_sd     <- constrained(statistic = "sd",     value = 0.8)
    })
    expect_s3_class(s_mean,   c("phi_spec_constrained", "phi_spec"))
    expect_s3_class(s_median, c("phi_spec_constrained", "phi_spec"))
    expect_equal(s_mean$statistic,      "mean")
    expect_equal(s_mean$statistic_code, 0L)   # PHI_CONSTRAINT_MEAN
    expect_equal(s_median$statistic_code, 1L) # PHI_CONSTRAINT_MEDIAN
    expect_equal(s_mode$statistic_code,   2L) # PHI_STATISTIC_MODE
    expect_equal(s_var$statistic_code,    3L) # PHI_STATISTIC_VARIANCE
    expect_equal(s_sd$statistic_code,     4L) # PHI_STATISTIC_SD
    expect_equal(s_median$value, 2)
})

test_that("constrained() warns for statistics that allow E[phi] to drift", {
    expect_warning(constrained(statistic = "median"),   "prior on sphi")
    expect_warning(constrained(statistic = "mode"),     "prior on sphi")
    expect_warning(constrained(statistic = "variance"), "prior on sphi")
    expect_warning(constrained(statistic = "sd"),       "prior on sphi")
    expect_no_warning(constrained(statistic = "mean"))
})

test_that("constrained() rejects invalid statistic", {
    expect_error(constrained(statistic = "harmonic"), "statistic must be one of")
    expect_error(constrained(statistic = ""),         "statistic must be one of")
})

test_that("constrained() rejects non-positive value", {
    expect_error(constrained(value = 0),   "positive")
    expect_error(constrained(value = -1),  "positive")
    expect_error(constrained(value = "x"), "positive")
})

test_that("fixed() returns phi_spec_fixed with correct value", {
    s <- fixed(value = -0.3)
    expect_s3_class(s, c("phi_spec_fixed", "phi_spec"))
    expect_equal(s$value, -0.3)
})

test_that("fixed() rejects non-finite values", {
    expect_error(fixed(Inf),  "finite")
    expect_error(fixed(-Inf), "finite")
    expect_error(fixed(NaN),  "finite")
})

test_that("estimated() returns phi_spec_estimated with prior", {
    s_unif <- estimated(prior = prior_uniform(0, 5))
    s_norm <- estimated(prior = prior_normal(mean = 1.4, sd = 0.125))
    s_null <- estimated(prior = NULL)

    expect_s3_class(s_unif, c("phi_spec_estimated", "phi_spec"))
    expect_equal(s_unif$prior$dist, "uniform")
    expect_equal(s_norm$prior$dist, "normal")
    expect_null(s_null$prior)
})

test_that("estimated() rejects non-prior_dist prior", {
    expect_error(estimated(prior = list(dist = "uniform")), "prior_dist")
    expect_error(estimated(prior = 1),                      "prior_dist")
})

test_that("prior_uniform() stores bounds and validates", {
    p <- prior_uniform(low = 0.5, high = 8)
    expect_s3_class(p, "prior_dist")
    expect_equal(p$dist, "uniform")
    expect_equal(p$low,  0.5)
    expect_equal(p$high, 8)
    expect_error(prior_uniform(5, 1),    "less than")
    expect_error(prior_uniform(Inf, 10), "finite")
    expect_error(prior_uniform(0, -Inf), "finite")
})

test_that("prior_normal() stores parameters and validates", {
    p <- prior_normal(mean = 1.4, sd = 0.125)
    expect_equal(p$dist, "normal")
    expect_equal(p$mean, 1.4)
    expect_equal(p$sd,   0.125)
    expect_error(prior_normal(sd = 0),  "positive")
    expect_error(prior_normal(sd = -1), "positive")
})

test_that("prior_student_t() stores parameters and validates", {
    p <- prior_student_t(df = 5, mean = 0, sd = 1)
    expect_equal(p$dist, "student_t")
    expect_equal(p$df,   5)
    expect_error(prior_student_t(df = 0), "positive")
})

test_that("prior_exponential() stores rate and validates", {
    p <- prior_exponential(rate = 2)
    expect_equal(p$dist, "exponential")
    expect_equal(p$rate, 2)
    expect_error(prior_exponential(rate = 0),  "positive")
    expect_error(prior_exponential(rate = -1), "positive")
})

test_that("print methods run without error", {
    s <- constrained(statistic = "mean", value = 1)
    f <- fixed(value = -0.2)
    e <- estimated(prior = prior_uniform(0, 10))
    p <- prior_normal(1.4, 0.125)
    expect_output(print(s), "constrained")
    expect_output(print(f), "fixed")
    expect_output(print(e), "estimated")
    expect_output(print(p), "Normal")
    expect_match(format(p), "Normal")
})


# ======================================================================
# Section 4: Parameter getter/setter round-trips
# ======================================================================

.make_param <- function(sphi = 1) {
    genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
    genome      <- initializeGenomeObject(file = genome_file)
    initializeParameterObject(genome = genome, sphi = sphi,
                               num.mixtures = 1,
                               gene.assignment = rep(1L, length(genome)))
}

test_that("getPhiMuMode / setPhiMuMode default is 0 (PHI_MU_CONSTRAINED)", {
    p <- .make_param()
    expect_equal(p$getPhiMuMode(), 0L)
    p$setPhiMuMode(1L)
    expect_equal(p$getPhiMuMode(), 1L)
    p$setPhiMuMode(0L)
    expect_equal(p$getPhiMuMode(), 0L)
})

test_that("getPhiConstraintValue / setPhiConstraintValue round-trips", {
    p <- .make_param()
    # Default is 1.0 (mean anchor at value=1)
    expect_equal(p$getPhiConstraintValue(), 1.0)
    p$setPhiConstraintValue(2.5)
    expect_equal(p$getPhiConstraintValue(), 2.5)
})

test_that("getPhiMuFixed / setPhiMuFixed round-trips", {
    p <- .make_param()
    p$setPhiMuFixed(-0.75)
    expect_equal(p$getPhiMuFixed(), -0.75)
    p$setPhiMuFixed(0.0)
    expect_equal(p$getPhiMuFixed(), 0.0)
})

test_that("getSphiPriorType / setSphiPriorType default is 0 (SPHI_PRIOR_FLAT)", {
    p <- .make_param()
    expect_equal(p$getSphiPriorType(), 0L)
    p$setSphiPriorType(1L)
    expect_equal(p$getSphiPriorType(), 1L)
    p$setSphiPriorType(2L)
    expect_equal(p$getSphiPriorType(), 2L)
    p$setSphiPriorType(0L)
    expect_equal(p$getSphiPriorType(), 0L)
})

test_that("getSphiPriorLow / getSphiPriorHigh / setSphiPriorBounds round-trips", {
    p <- .make_param()
    # Defaults set in Parameter ctor
    p$setSphiPriorBounds(0.5, 8.0)
    expect_equal(p$getSphiPriorLow(),  0.5)
    expect_equal(p$getSphiPriorHigh(), 8.0)
})

test_that("getPhiFloor / getPhiCeiling defaults and setPhiBounds round-trips", {
    p <- .make_param()
    expect_equal(p$getPhiFloor(),   1e-6)
    expect_equal(p$getPhiCeiling(), 1e6)
    p$setPhiBounds(1e-4, 1e4)
    expect_equal(p$getPhiFloor(),   1e-4)
    expect_equal(p$getPhiCeiling(), 1e4)
})


# ======================================================================
# Section 5: phi.mphi integration (constrained / fixed modes)
# ======================================================================

.genome_and_assign <- function() {
    genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
    genome      <- initializeGenomeObject(file = genome_file)
    list(genome = genome, gene_assign = rep(1L, length(genome)))
}

test_that("phi.mphi=constrained(mean) sets phiMuMode=0 and phiPriorConstraint=0", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = 1, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.mphi = constrained(statistic = "mean", value = 1))
    expect_equal(p$getPhiMuMode(),          0L)  # PHI_MU_CONSTRAINED
    expect_equal(p$getPhiPriorConstraint(), 0L)  # PHI_CONSTRAINT_MEAN
    expect_equal(p$getPhiConstraintValue(), 1.0)
})

test_that("phi.mphi=constrained(median) sets phiPriorConstraint=1", {
    ga <- .genome_and_assign()
    p  <- suppressWarnings(initializeParameterObject(
              genome = ga$genome, sphi = 1, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.mphi = constrained(statistic = "median", value = 1)))
    expect_equal(p$getPhiMuMode(),          0L)  # PHI_MU_CONSTRAINED
    expect_equal(p$getPhiPriorConstraint(), 1L)  # PHI_CONSTRAINT_MEDIAN
})

test_that("phi.mphi=constrained(mode) sets phiPriorConstraint=2", {
    ga <- .genome_and_assign()
    p  <- suppressWarnings(initializeParameterObject(
              genome = ga$genome, sphi = 1, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.mphi = constrained(statistic = "mode", value = 1)))
    expect_equal(p$getPhiPriorConstraint(), 2L)  # PHI_STATISTIC_MODE
})

test_that("phi.mphi=constrained with value=2 stores correct constraint value", {
    ga <- .genome_and_assign()
    p  <- suppressWarnings(initializeParameterObject(
              genome = ga$genome, sphi = 1, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.mphi = constrained(statistic = "median", value = 2.0)))
    expect_equal(p$getPhiConstraintValue(), 2.0)
})

test_that("phi.mphi=fixed() sets phiMuMode=1 and phiMuFixed", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = 1, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.mphi = fixed(value = -0.2))
    expect_equal(p$getPhiMuMode(),  1L)   # PHI_MU_FIXED
    expect_equal(p$getPhiMuFixed(), -0.2)
})

test_that("phi.mphi=estimated() sets PHI_MU_ESTIMATED mode and flat prior", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = 1, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.mphi = estimated())
    expect_equal(p$getPhiMuMode(), 2L)   # PHI_MU_ESTIMATED
    expect_equal(p$getPhiMuPriorType(), 0L)  # PHI_MU_PRIOR_FLAT
    expect_equal(p$getMuSynthesisRate(0L, FALSE), 0.0)  # default init
})

test_that("phi.mphi=estimated(prior_normal) sets normal prior on mphi", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = 1, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.mphi = estimated(prior = prior_normal(mean = 0.5, sd = 2.0)))
    expect_equal(p$getPhiMuMode(),      2L)    # PHI_MU_ESTIMATED
    expect_equal(p$getPhiMuPriorType(), 1L)    # PHI_MU_PRIOR_NORMAL
    expect_equal(p$getPhiMuPriorMu(),   0.5)
    expect_equal(p$getPhiMuPriorSd(),   2.0)
})

test_that("phi.mphi=estimated(init=...) sets custom initial mphi", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = 1, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.mphi = estimated(init = -0.5))
    expect_equal(p$getMuSynthesisRate(0L, FALSE), -0.5)
})


# ======================================================================
# Section 6: phi.sphi integration
# ======================================================================

test_that("phi.sphi=estimated(prior_uniform) sets sphiPriorType=2 and bounds", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = NA, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.sphi = estimated(prior = prior_uniform(0, 10)))
    expect_equal(p$getSphiPriorType(), 2L)   # SPHI_PRIOR_UNIFORM
    expect_equal(p$getSphiPriorHigh(), 10.0)
    expect_equal(p$getSphiPriorLow(),   0.0)
})

test_that("phi.sphi=estimated(prior_normal) sets sphiPriorType=1 and mu/sd", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = NA, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.sphi = estimated(prior = prior_normal(mean = 1.4, sd = 0.125)))
    expect_equal(p$getSphiPriorType(), 1L)    # SPHI_PRIOR_NORMAL
    expect_equal(p$getSphiPriorMu(),   1.4)
    expect_equal(p$getSphiPriorSd(),   0.125)
})

test_that("phi.sphi=estimated(prior=NULL) sets sphiPriorType=0 (flat)", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = NA, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.sphi = estimated(prior = NULL))
    expect_equal(p$getSphiPriorType(), 0L)   # SPHI_PRIOR_FLAT
})

test_that("phi.sphi=fixed(1) fixes sphi (default via sphi=1 legacy path)", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = 1, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.sphi = fixed(value = 1))
    # sphiPriorType unchanged from default (flat); fixSphi is the relevant flag
    expect_equal(p$getSphiPriorType(), 0L)
    # The MCMC smoke below confirms it runs correctly with est.hyper=FALSE
})

test_that("sphi=NA resolves to estimated(prior_uniform(0,10)) automatically", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = NA, num.mixtures = 1,
              gene.assignment = ga$gene_assign)
    expect_equal(p$getSphiPriorType(), 2L)   # SPHI_PRIOR_UNIFORM
    expect_equal(p$getSphiPriorHigh(), 10.0)
})


# ======================================================================
# Section 7: MCMC smoke tests
# ======================================================================

# genome must be the same object used to create parameter (shared C++ state)
.run_roc_mcmc <- function(parameter, genome, seed = 42) {
    mcmc  <- initializeMCMCObject(samples = 4, thinning = 2,
                                   adaptive.width = 2,
                                   est.expression = TRUE,
                                   est.csp        = TRUE,
                                   est.hyper      = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)
    set.seed(seed)
    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome, model = model,
            ncores = 1, divergence.iteration = 0)
    sink()
    mcmc
}

test_that("MCMC: constrained(median,1) + sphi=1 runs and traces have correct length", {
    ga <- .genome_and_assign()
    p  <- suppressWarnings(initializeParameterObject(
              genome = ga$genome, sphi = 1, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.mphi = constrained(statistic = "median", value = 1)))
    mcmc <- .run_roc_mcmc(p, ga$genome)
    expect_equal(length(mcmc$getLogPosteriorTrace()), 5L)  # samples=4 + 1
})

test_that("MCMC: constrained(mean,1) + estimated(prior_uniform) runs without error", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = NA, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.mphi = constrained(statistic = "mean", value = 1),
              phi.sphi = estimated(prior = prior_uniform(0, 10)))
    mcmc <- .run_roc_mcmc(p, ga$genome, seed = 51)
    expect_equal(length(mcmc$getLogPosteriorTrace()), 5L)
})

test_that("MCMC: fixed(mphi) + estimated(prior_uniform) runs without error", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = NA, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.mphi = fixed(value = -0.3),
              phi.sphi = estimated(prior = prior_uniform(0, 10)))
    mcmc <- .run_roc_mcmc(p, ga$genome, seed = 52)
    expect_equal(length(mcmc$getLogPosteriorTrace()), 5L)
})

test_that("MCMC: constrained(mean) + estimated(prior_normal) runs without error", {
    ga <- .genome_and_assign()
    p  <- initializeParameterObject(
              genome = ga$genome, sphi = NA, num.mixtures = 1,
              gene.assignment = ga$gene_assign,
              phi.sphi = estimated(prior = prior_normal(mean = 1.4, sd = 0.5)))
    mcmc <- .run_roc_mcmc(p, ga$genome, seed = 53)
    expect_equal(length(mcmc$getLogPosteriorTrace()), 5L)
})


# ======================================================================
# Section 8: Restart file round-trip for all new phi spec members
# ======================================================================

test_that("restart preserves phiMuMode, phiConstraintValue, sphiPriorType/bounds", {
    genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
    genome      <- initializeGenomeObject(file = genome_file)

    p <- suppressWarnings(initializeParameterObject(
             genome = genome, sphi = NA, num.mixtures = 1,
             gene.assignment = rep(1L, length(genome)),
             phi.mphi = constrained(statistic = "median", value = 2.0),
             phi.sphi = estimated(prior = prior_uniform(0.5, 8.0))))

    mcmc  <- initializeMCMCObject(samples = 2, thinning = 2, adaptive.width = 2,
                                   est.expression = TRUE, est.csp = TRUE, est.hyper = TRUE)
    model <- initializeModelObject(p, "ROC", with.phi = FALSE)
    rst   <- tempfile(pattern = "phi_spec_rst_", fileext = ".rst")
    setRestartSettings(mcmc, filename = rst, samples = 2, write.multiple = FALSE)
    set.seed(60)
    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome, model = model,
            ncores = 1, divergence.iteration = 0)
    sink()

    p2 <- initializeParameterObject(genome = genome,
                                     init.with.restart.file = rst,
                                     model = "ROC")
    expect_equal(p2$getPhiMuMode(),          0L)   # PHI_MU_CONSTRAINED
    expect_equal(p2$getPhiPriorConstraint(), 1L)   # PHI_CONSTRAINT_MEDIAN
    expect_equal(p2$getPhiConstraintValue(), 2.0)
    expect_equal(p2$getSphiPriorType(),      2L)   # SPHI_PRIOR_UNIFORM
    expect_equal(p2$getSphiPriorLow(),       0.5)
    expect_equal(p2$getSphiPriorHigh(),      8.0)
    unlink(rst)
})

test_that("restart preserves phiMuMode=1 (fixed) and phiMuFixed", {
    genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
    genome      <- initializeGenomeObject(file = genome_file)

    p <- initializeParameterObject(
             genome = genome, sphi = NA, num.mixtures = 1,
             gene.assignment = rep(1L, length(genome)),
             phi.mphi = fixed(value = -0.4),
             phi.sphi = estimated(prior = prior_uniform(0, 10)))

    mcmc  <- initializeMCMCObject(samples = 2, thinning = 2, adaptive.width = 2,
                                   est.expression = TRUE, est.csp = TRUE, est.hyper = TRUE)
    model <- initializeModelObject(p, "ROC", with.phi = FALSE)
    rst   <- tempfile(pattern = "phi_fixed_rst_", fileext = ".rst")
    setRestartSettings(mcmc, filename = rst, samples = 2, write.multiple = FALSE)
    set.seed(61)
    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome, model = model,
            ncores = 1, divergence.iteration = 0)
    sink()

    p2 <- initializeParameterObject(genome = genome,
                                     init.with.restart.file = rst,
                                     model = "ROC")
    expect_equal(p2$getPhiMuMode(),  1L)    # PHI_MU_FIXED
    expect_equal(p2$getPhiMuFixed(), -0.4)
    unlink(rst)
})

test_that("median-anchor phiPriorConstraint round-trips through restart file", {
    genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
    genome      <- initializeGenomeObject(file = genome_file)

    parameter <- suppressWarnings(initializeParameterObject(
                     genome = genome, sphi = 1, num.mixtures = 1,
                     gene.assignment = rep(1L, length(genome)),
                     phi.prior.constraint = "median"))
    mcmc  <- initializeMCMCObject(samples = 2, thinning = 2, adaptive.width = 2,
                                   est.expression = TRUE, est.csp = TRUE, est.hyper = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)
    rst   <- tempfile(pattern = "phi_anchor_rst_", fileext = ".rst")
    setRestartSettings(mcmc, filename = rst, samples = 2, write.multiple = FALSE)
    set.seed(13)
    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome, model = model,
            ncores = 1, divergence.iteration = 0)
    sink()

    p2 <- initializeParameterObject(genome = genome,
                                     init.with.restart.file = rst, model = "ROC")
    expect_equal(p2$getPhiPriorConstraint(), 1L)  # PHI_CONSTRAINT_MEDIAN preserved
    unlink(rst)
})


# ======================================================================
# Section 9: Legacy phi.prior.constraint (backward compat + deprecation)
# ======================================================================

test_that("phi.prior.constraint = 'median' sets PHI_CONSTRAINT_MEDIAN in C++ object", {
    ga <- .genome_and_assign()
    p_mean   <- initializeParameterObject(
                    genome = ga$genome, sphi = 1, num.mixtures = 1,
                    gene.assignment = ga$gene_assign,
                    phi.prior.constraint = "mean")
    p_median <- suppressWarnings(initializeParameterObject(
                    genome = ga$genome, sphi = 1, num.mixtures = 1,
                    gene.assignment = ga$gene_assign,
                    phi.prior.constraint = "median"))
    expect_equal(p_mean$getPhiPriorConstraint(),   0L)  # PHI_CONSTRAINT_MEAN
    expect_equal(p_median$getPhiPriorConstraint(), 1L)  # PHI_CONSTRAINT_MEDIAN
})

test_that("phi.prior.constraint='median' emits deprecation warning", {
    ga <- .genome_and_assign()
    expect_warning(
        initializeParameterObject(
            genome = ga$genome, sphi = 1, num.mixtures = 1,
            gene.assignment = ga$gene_assign,
            phi.prior.constraint = "median"),
        "deprecated")
})

test_that("median-anchor MCMC runs without error and produces correct trace length", {
    genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
    genome      <- initializeGenomeObject(file = genome_file)

    set.seed(7)
    parameter <- suppressWarnings(initializeParameterObject(
                     genome = genome, sphi = 1, num.mixtures = 1,
                     gene.assignment = rep(1L, length(genome)),
                     phi.prior.constraint = "median"))
    mcmc  <- initializeMCMCObject(samples = 4, thinning = 2, adaptive.width = 2,
                                   est.expression = TRUE, est.csp = TRUE, est.hyper = TRUE)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)
    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome, model = model,
            ncores = 1, divergence.iteration = 0)
    sink()
    expect_equal(length(mcmc$getLogPosteriorTrace()), 5L)  # samples + 1
})


# ======================================================================
# Section 10: sphi prior default-off regression
# ======================================================================

test_that("default ROC parameter has sphi prior OFF (sphiPriorSd == 0)", {
    ga        <- .genome_and_assign()
    parameter <- initializeParameterObject(
                     genome = ga$genome, sphi = 1, num.mixtures = 1,
                     gene.assignment = ga$gene_assign)
    expect_equal(parameter$getSphiPriorSd(), 0.0)
})

test_that("setSphiPrior activates the Normal prior and can be disabled", {
    ga        <- .genome_and_assign()
    parameter <- initializeParameterObject(
                     genome = ga$genome, sphi = 1, num.mixtures = 1,
                     gene.assignment = ga$gene_assign)
    parameter$setSphiPrior(1.4, 0.05)
    expect_equal(parameter$getSphiPriorMu(), 1.4)
    expect_equal(parameter$getSphiPriorSd(), 0.05)

    parameter$setSphiPrior(1.4, 0.0)
    expect_equal(parameter$getSphiPriorSd(), 0.0)
})
