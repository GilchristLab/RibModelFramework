library(testthat)
library(AnaCoDa)

# ============================================================================
# Tests for roc_sphi_est.stan and initializeStan() via cmdstanr.
#
# Skipped automatically if cmdstanr is not installed or cmdstan is not found.
# These tests exercise mathematical identities (no sampling required beyond
# a single optimize(iter=1) call to instantiate a fit object for log_prob()).
#
# Tests ported and adapted from scripts/test_reparam_inertness_stan.R and
# scripts/test_reparam_stan.R, updated for the new phi spec data block.
#
# Key identity: with noncentered=1 and deta_scale_anchor=1,
#   dEta_eff * phi = dEta * exp(ref + sphi*z)
# is INDEPENDENT of mphi (phi_mphi_fixed).  With phi_mphi_mode=1 there is
# also no prior on phi_mphi_fixed, so the full log-prob is zero-invariant.
# With scale_anchor=0, dEta*phi = dEta*exp(mphi + sphi*z) DEPENDS on mphi.
# That contrast is the discriminating test.
# ============================================================================

# ---- availability guard ------------------------------------------------------

has_cmdstan <- local({
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(FALSE)
  path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)
  !is.null(path) && nzchar(path)
})

# ---- compile and build base data once ----------------------------------------

stan_mod  <- NULL
stan_base <- NULL
genome_g  <- NULL   # genome object (shared)
scuo_g    <- NULL   # SCUO vector (shared)
K_g       <- NULL   # number of non-reference codons
G_g       <- NULL   # number of genes

if (has_cmdstan) {
  stan_file <- system.file("stan/roc_sphi_est.stan", package = "AnaCoDa")
  fasta     <- system.file("extdata", "genome.fasta", package = "AnaCoDa")

  genome_g <- initializeGenomeObject(file = fasta)
  scuo_g   <- calculateSCUO(genome_g)$SCUO

  stan_mod <- tryCatch(
    cmdstanr::cmdstan_model(stan_file,
                            cpp_options             = list(stan_threads = TRUE),
                            compile_model_methods   = TRUE,
                            force_recompile         = TRUE,
                            quiet                   = TRUE),
    error = function(e) NULL
  )

  if (!is.null(stan_mod)) {
    # stan_base pins the phi spec to the package defaults shared by BOTH backends:
    #   phi.mphi = constrained("mean", 1)          matches the MCMC native default
    #                                              (PHI_MU_CONSTRAINED, MEAN, 1.0).
    #   phi.sphi = estimated(prior_uniform(0, 10)) matches the MCMC default reached
    #                                              via initializeParameterObject(sphi = NA),
    #                                              which resolves to the same Uniform(0,10).
    # These are the apples-to-apples cross-backend defaults, not arbitrary choices.
    # (The legacy unbounded-flat sphi prior is still available via estimated(NULL).)
    stan_base <- suppressWarnings(
      initializeStan(genome_g, scuo = scuo_g,
                     phi.mphi = constrained("mean", 1),
                     phi.sphi = estimated(prior_uniform(0, 10)),
                     noncentered       = 0L,
                     deta_scale_anchor = 0L,
                     deta_anchor_ref   = 0.0,
                     deta_phi_center   = 0.0)
    )
    K_g <- stan_base$data$K
    G_g <- stan_base$data$G
  }
}

# ---- shared helpers ----------------------------------------------------------

# Run optimize(iter=1) to bind data to a fit object, then evaluate log_prob
# at a fixed unconstrained theta (jacobian=TRUE = standard HMC convention).
.lp_at <- function(data, theta) {
  fit <- stan_mod$optimize(
    data = data, iter = 1L, algorithm = "lbfgs",
    threads = 1L, refresh = 0, show_messages = FALSE
  )
  fit$log_prob(unconstrained_variables = theta, jacobian = TRUE)
}

# Build a fixed unconstrained theta.
# Parameter order: dM[K], dEta[K], latent_phi[G], sphi_unconstrained, mphi_param.
# mphi_param is ALWAYS in the parameters block (post-PR#50, fixes #47); in modes
# 0/1 it is pinned via std_normal() and contributes -0.5*mphi_param^2 to the
# log-prob, so two log-probs compared at the SAME mphi_param value cancel that term.
# sphi Uniform(0,10): unconstrained = logit(sphi/10) = log(sphi/(10-sphi)).
# dEta must be NON-ZERO so the likelihood depends on mphi when scale_anchor=0.
.make_theta <- function(K, G, seed = 42L, sphi = 1.0, sphi_high = 10.0, mphi = 0.0) {
  set.seed(seed)
  sphi_u <- log(sphi / (sphi_high - sphi))
  c(rnorm(K, 0, 0.3),   # dM
    rnorm(K, 0, 0.5),   # dEta (non-zero)
    rnorm(G, 0, 1.0),   # latent_phi
    sphi_u,             # sphi (unconstrained)
    mphi)               # mphi_param
}

# Build data with phi_mphi_mode=1 (fixed mphi), noncentered=nc, deta_scale_anchor=sa.
.d_fixed_mphi <- function(mphi_val, noncentered = 1L, scale_anchor = 1L) {
  suppressWarnings(
    initializeStan(genome_g, scuo = scuo_g,
                   phi.mphi = fixed(mphi_val),
                   phi.sphi = estimated(prior_uniform(0, 10)),
                   noncentered       = noncentered,
                   deta_scale_anchor = scale_anchor,
                   deta_anchor_ref   = 0.0,
                   deta_phi_center   = 0.0)
  )$data
}

# Build data with phi_mphi_mode=2 (ESTIMATED mphi, FLAT prior), noncentered=nc,
# deta_scale_anchor=sa.  A flat prior (estimated(NULL)) puts NO prior term on
# mphi_param, so the log-prob is a pure function of the likelihood + sphi prior.
# This is the configuration in which deta_scale_anchor actually collapses the
# dEta-phi ridge (the original purpose of issue #47).
.d_est_mphi <- function(noncentered = 1L, scale_anchor = 1L) {
  suppressWarnings(
    initializeStan(genome_g, scuo = scuo_g,
                   phi.mphi = estimated(NULL),                 # flat prior on mphi
                   phi.sphi = estimated(prior_uniform(0, 10)),
                   noncentered       = noncentered,
                   deta_scale_anchor = scale_anchor,
                   deta_anchor_ref   = 0.0,
                   deta_phi_center   = 0.0)
  )$data
}

# Replace the mphi_param element (last entry) of an unconstrained theta vector.
.set_mphi_param <- function(theta, val) {
  theta[length(theta)] <- val
  theta
}


# ============================================================================
# Section 1: compilation
# ============================================================================

context("roc_sphi_est.stan: compilation")

test_that("roc_sphi_est.stan compiles via cmdstanr", {
  skip_if(!has_cmdstan, "cmdstan not available")
  # Re-attempt compile here so the error message is visible in the test output
  stan_file <- system.file("stan/roc_sphi_est.stan", package = "AnaCoDa")
  expect_true(nzchar(stan_file), label = "inst/stan/roc_sphi_est.stan not found in installed package")
  mod_test <- cmdstanr::cmdstan_model(stan_file,
                                      cpp_options           = list(stan_threads = TRUE),
                                      compile_model_methods = TRUE,
                                      quiet                 = TRUE)
  expect_false(is.null(mod_test), label = "cmdstan_model() compilation failed")
})


# ============================================================================
# Section 2: log-prob sanity with initializeStan() data
# ============================================================================

context("roc_sphi_est.stan: log-prob sanity with initializeStan() data")

test_that("optimize(iter=1) at initializeStan() init gives finite log-prob", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  fit <- stan_mod$optimize(
    data = stan_base$data, init = list(stan_base$init),
    iter = 1L, algorithm = "lbfgs",
    threads = 1L, refresh = 0, show_messages = FALSE
  )
  expect_true(is.finite(fit$lp()))
})

test_that("phi spec data fields have correct types for default constrained('mean',1)", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  d <- stan_base$data
  expect_equal(d$phi_mphi_mode,      0L)   # constrained
  expect_equal(d$phi_mphi_statistic, 0L)   # mean
  expect_equal(d$phi_mphi_value,     1.0)
  expect_equal(d$sphi_low,           0.0)
  expect_equal(d$sphi_high,          10.0)
  expect_equal(d$sphi_prior_type,    0L)   # uniform (bounds only)
})

test_that("initializeStan() init list includes mphi_param (always in parameters block)", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  # Post-PR#50 (fixes #47), mphi_param is always declared in the Stan parameters
  # block (free when estimated, pinned via std_normal() otherwise), so the init
  # list must always supply it.
  expected <- c("dM", "dEta", "latent_phi", "sphi", "mphi_param")
  expect_true(all(expected %in% names(stan_base$init)))
})


# ============================================================================
# Section 3: MAP estimate quality
# ============================================================================

context("roc_sphi_est.stan: MAP estimate (optimize)")

test_that("optimize(iter=500) gives finite, valid MAP estimates", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  fit <- stan_mod$optimize(
    data = stan_base$data, init = list(stan_base$init),
    iter = 500L, algorithm = "lbfgs",
    threads = 1L, refresh = 0, show_messages = FALSE
  )
  draws    <- fit$draws(format = "df")
  sphi_map <- draws[["sphi"]]

  expect_true(is.finite(fit$lp()),    label = "lp() is finite")
  expect_true(is.finite(sphi_map),    label = "sphi MAP is finite")
  expect_gt(sphi_map, 0,              label = "sphi MAP > 0")
  expect_lt(sphi_map, 10,             label = "sphi MAP < upper bound")

  dM_cols <- paste0("dM[", 1:K_g, "]")
  dM_vals <- as.numeric(unlist(lapply(dM_cols, function(col) draws[[col]])))
  expect_true(all(is.finite(dM_vals)), label = "all dM MAP values are finite")
})


# ============================================================================
# Section 4: deta_anchor_ref inertness (scale_anchor = 0)
# ============================================================================

context("roc_sphi_est.stan: deta_anchor_ref inertness")

test_that("log-prob is invariant to deta_anchor_ref when deta_scale_anchor=0", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  theta <- .make_theta(K_g, G_g)

  d_r0 <- stan_base$data; d_r0$deta_anchor_ref <- 0.0
  d_r5 <- stan_base$data; d_r5$deta_anchor_ref <- 5.0

  lp0 <- .lp_at(d_r0, theta)
  lp5 <- .lp_at(d_r5, theta)

  expect_true(is.finite(lp0))
  expect_lt(abs(lp0 - lp5), 1e-9,
            label = "deta_anchor_ref does not affect log-prob when scale_anchor=0")
})


# ============================================================================
# Section 5: scale-anchor mphi-independence (noncentered=1, phi_mphi_mode=1)
# ============================================================================

context("roc_sphi_est.stan: scale-anchor mphi-independence")

test_that("scale_anchor=1, noncentered=1: log-prob is zero-invariant to phi_mphi_fixed", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  theta <- .make_theta(K_g, G_g)

  lp1 <- .lp_at(.d_fixed_mphi(-0.2, noncentered = 1L, scale_anchor = 1L), theta)
  lp2 <- .lp_at(.d_fixed_mphi( 0.3, noncentered = 1L, scale_anchor = 1L), theta)

  expect_true(is.finite(lp1) && is.finite(lp2))
  expect_lt(abs(lp1 - lp2), 1e-9,
            label = "log-prob does not change with phi_mphi_fixed when scale_anchor=1, noncentered=1")
})

test_that("scale_anchor=0, noncentered=1: log-prob IS sensitive to phi_mphi_fixed (contrast)", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  theta <- .make_theta(K_g, G_g)

  lp1 <- .lp_at(.d_fixed_mphi(-0.2, noncentered = 1L, scale_anchor = 0L), theta)
  lp2 <- .lp_at(.d_fixed_mphi( 0.3, noncentered = 1L, scale_anchor = 0L), theta)

  expect_true(is.finite(lp1) && is.finite(lp2))
  expect_gt(abs(lp1 - lp2), 0.01,
            label = "log-prob changes with phi_mphi_fixed when scale_anchor=0 (dEta*phi depends on mphi)")
})

# --- estimated() mphi: the configuration deta_scale_anchor was built for -------
#
# With phi_mphi_mode=2 and a FLAT prior, mphi = mphi_param is a free parameter
# with no prior term.  These tests vary the mphi_param ELEMENT of theta (not a
# data field), which is the genuine test of whether the scale-anchor decorrelates
# the sampled mphi from the rest of the posterior (issue #47).

test_that("estimated() flat prior, scale_anchor=1, noncentered=1: log-prob is zero-invariant to mphi_param value", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  d     <- .d_est_mphi(noncentered = 1L, scale_anchor = 1L)
  theta <- .make_theta(K_g, G_g)

  lp1 <- .lp_at(d, .set_mphi_param(theta, -0.2))
  lp2 <- .lp_at(d, .set_mphi_param(theta,  0.3))

  expect_true(is.finite(lp1) && is.finite(lp2))
  expect_lt(abs(lp1 - lp2), 1e-9,
            label = "free mphi_param does not affect log-prob when scale_anchor=1, noncentered=1, flat prior")
})

test_that("estimated() flat prior, scale_anchor=0, noncentered=1: log-prob IS sensitive to mphi_param value (contrast)", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  d     <- .d_est_mphi(noncentered = 1L, scale_anchor = 0L)
  theta <- .make_theta(K_g, G_g)

  lp1 <- .lp_at(d, .set_mphi_param(theta, -0.2))
  lp2 <- .lp_at(d, .set_mphi_param(theta,  0.3))

  expect_true(is.finite(lp1) && is.finite(lp2))
  expect_gt(abs(lp1 - lp2), 0.01,
            label = "free mphi_param changes log-prob when scale_anchor=0 (dEta*phi depends on mphi)")
})


# ============================================================================
# Section 6: computeMPhi formula in Stan == R formula (constrained mode)
# ============================================================================
#
# At sphi=1, constrained("mean", 1) gives mphi = log(1) - 0.5*1^2 = -0.5,
# which is identical to fixed(-0.5).  Same check for median=1 → mphi=0.
# Uses a theta pinned at sphi=1 (unconstrained for Uniform(0,10)).

context("roc_sphi_est.stan: computeMPhi formula consistency")

# theta pinned at sphi=1 (Uniform(0,10) unconstrained = logit(1/10) = log(1/9))
# Guarded: K_g/G_g are NULL when compilation failed; evaluated lazily in tests.
.theta_sphi1 <- if (!is.null(K_g) && !is.null(G_g)) {
  set.seed(42L)
  sphi_u <- log(1.0 / 9.0)   # logit(1/10) = log(1/9)
  # trailing 0.0 is mphi_param (always last in the parameters block post-PR#50)
  c(rnorm(K_g, 0, 0.3), rnorm(K_g, 0, 0.5), rnorm(G_g, 0, 1.0), sphi_u, 0.0)
} else {
  NULL
}

test_that("constrained('mean', 1) == fixed(-0.5) at sphi=1", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  d_c <- suppressWarnings(
    initializeStan(genome_g, scuo = scuo_g,
                   phi.mphi = constrained("mean", 1),
                   phi.sphi = estimated(prior_uniform(0, 10)),
                   noncentered = 0L)
  )$data

  d_f <- suppressWarnings(
    initializeStan(genome_g, scuo = scuo_g,
                   phi.mphi = fixed(-0.5),   # log(1) - 0.5*1^2 = -0.5
                   phi.sphi = estimated(prior_uniform(0, 10)),
                   noncentered = 0L)
  )$data

  lp_c <- .lp_at(d_c, .theta_sphi1)
  lp_f <- .lp_at(d_f, .theta_sphi1)

  expect_true(is.finite(lp_c) && is.finite(lp_f))
  expect_lt(abs(lp_c - lp_f), 1e-9,
            label = "constrained('mean',1) and fixed(-0.5) give identical log-prob at sphi=1")
})

test_that("constrained('median', 1) == fixed(0.0) at sphi=1", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  d_c <- suppressWarnings(
    initializeStan(genome_g, scuo = scuo_g,
                   phi.mphi = constrained("median", 1),
                   phi.sphi = estimated(prior_uniform(0, 10)),
                   noncentered = 0L)
  )$data

  d_f <- suppressWarnings(
    initializeStan(genome_g, scuo = scuo_g,
                   phi.mphi = fixed(0.0),   # log(1) = 0
                   phi.sphi = estimated(prior_uniform(0, 10)),
                   noncentered = 0L)
  )$data

  lp_c <- .lp_at(d_c, .theta_sphi1)
  lp_f <- .lp_at(d_f, .theta_sphi1)

  expect_true(is.finite(lp_c) && is.finite(lp_f))
  expect_lt(abs(lp_c - lp_f), 1e-9,
            label = "constrained('median',1) and fixed(0) give identical log-prob at sphi=1")
})

test_that("constrained('mode', 1) == fixed(1.0) at sphi=1", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  d_c <- suppressWarnings(
    initializeStan(genome_g, scuo = scuo_g,
                   phi.mphi = constrained("mode", 1),
                   phi.sphi = estimated(prior_uniform(0, 10)),
                   noncentered = 0L)
  )$data

  d_f <- suppressWarnings(
    initializeStan(genome_g, scuo = scuo_g,
                   phi.mphi = fixed(1.0),   # log(1) + 1^2 = 0 + 1 = 1
                   phi.sphi = estimated(prior_uniform(0, 10)),
                   noncentered = 0L)
  )$data

  lp_c <- .lp_at(d_c, .theta_sphi1)
  lp_f <- .lp_at(d_f, .theta_sphi1)

  expect_true(is.finite(lp_c) && is.finite(lp_f))
  expect_lt(abs(lp_c - lp_f), 1e-9,
            label = "constrained('mode',1) and fixed(1.0) give identical log-prob at sphi=1")
})


# ============================================================================
# Section 7: Normal sphi prior pulls MAP toward prior mean
# ============================================================================
#
# With a tight Normal(5, 0.1) prior, optimize() should return sphi near 5
# regardless of the data.  With Uniform(0,10), sphi is data-driven and for
# the small test genome (G=8) will be far from 5.

context("roc_sphi_est.stan: Normal sphi prior pulls MAP toward prior mean")

test_that("tight Normal(5, 0.1) prior anchors sphi MAP near 5", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  s_norm <- suppressWarnings(
    initializeStan(genome_g, scuo = scuo_g,
                   phi.mphi = constrained("mean", 1),
                   phi.sphi = estimated(prior_normal(5.0, 0.1)),
                   noncentered = 0L,
                   sphi.init   = 5.0)
  )

  fit <- stan_mod$optimize(
    data = s_norm$data, init = list(s_norm$init),
    iter = 500L, algorithm = "lbfgs",
    threads = 1L, refresh = 0, show_messages = FALSE
  )
  sphi_map <- fit$draws(format = "df")[["sphi"]]

  expect_true(is.finite(sphi_map))
  expect_gt(sphi_map, 4.0, label = "tight N(5,0.1) pulls sphi above 4")
  expect_lt(sphi_map, 6.0, label = "tight N(5,0.1) keeps sphi below 6")
})

test_that("Uniform(0,10) prior gives different (data-driven) sphi MAP", {
  skip_if(!has_cmdstan, "cmdstan not available")
  skip_if(is.null(stan_mod), "compilation failed")

  fit_u <- stan_mod$optimize(
    data = stan_base$data, init = list(stan_base$init),
    iter = 500L, algorithm = "lbfgs",
    threads = 1L, refresh = 0, show_messages = FALSE
  )
  sphi_unif <- fit_u$draws(format = "df")[["sphi"]]

  s_norm <- suppressWarnings(
    initializeStan(genome_g, scuo = scuo_g,
                   phi.mphi = constrained("mean", 1),
                   phi.sphi = estimated(prior_normal(5.0, 0.1)),
                   noncentered = 0L,
                   sphi.init   = 5.0)
  )
  fit_n <- stan_mod$optimize(
    data = s_norm$data, init = list(s_norm$init),
    iter = 500L, algorithm = "lbfgs",
    threads = 1L, refresh = 0, show_messages = FALSE
  )
  sphi_norm <- fit_n$draws(format = "df")[["sphi"]]

  # Tight Normal prior near 5 should produce substantially different sphi from data-driven estimate
  expect_gt(abs(sphi_norm - sphi_unif), 0.5,
            label = "Normal(5,0.1) and Uniform(0,10) priors give meaningfully different MAP sphi")
})
