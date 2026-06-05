# PANSE Stan phi prior: behavioral smoke tests.
#
# Runs 1 chain x 10 warmup / 10 sampling for all four phi prior combinations.
# Checks: no error, finite draws, correct parameter structure.
# Does NOT check convergence or scientific correctness.
library(testthat)

context("panse-phi-prior: smoke fits for all four combinations")

.smoke_config <- function(fx, phi_cfg = NULL) {
    cfg <- list(
        name   = "phi-prior-smoke",
        genome = list(input.dir = fx, pattern = "rfp_counts.csv"),
        fit    = list(
            model              = "sphi-est",
            parameterization   = "noncentered",
            with.phi           = TRUE,
            sphi.prior.sd      = 2.5,
            init.mode          = "fixed",
            init.partition.function = "auto",
            init.alpha.files   = list(file.path(fx, "truth_alpha.csv")),
            init.lambda.files  = list(file.path(fx, "truth_lambda.csv")),
            init.nserate.files = list(file.path(fx, "truth_nse.csv")),
            init.phi.file      = file.path(fx, "truth_phi.csv"),
            alpha.prior.lower  = 1e-3, alpha.prior.upper  = 100,
            lambda.prior.lower = 1e-3, lambda.prior.upper = 100,
            nserate.prior = list(type = "Log-Uniform",
                                 uniform.lower = 1e-7, uniform.upper = 0.1)
        ),
        stan = list(chains = 1L, parallel_chains = 1L,
                    warmup = 10L, sampling = 10L,
                    threads_per_chain = 1L, adapt_delta = 0.8,
                    max_treedepth = 6L, metric = "diag_e",
                    grainsize = 1L, seed = 42L)
    )
    if (!is.null(phi_cfg)) cfg$phi <- phi_cfg
    cfg
}

.data_spec <- function(pfx, col) {
    paste0("data(", file.path(pfx, "phi_prior_all.csv"), ", col=", col, ")")
}

test_that("smoke: population/population runs and produces finite draws", {
    skip_if_no_cmdstan()
    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    out <- tempfile("phi-prior-smoke-pp-")
    fit <- AnaCoDa::fit_panse_stan(.smoke_config(fx), out_dir = out,
                                   no_log_lik = TRUE, verbose = FALSE)
    draws <- fit$draws(format = "df")
    expect_true(all(is.finite(draws[["sphi"]])))
    expect_true(all(is.finite(draws[["log_phi[1]"]])))
})

test_that("smoke: data/population runs and produces finite draws", {
    skip_if_no_cmdstan()
    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    pfx <- testthat::test_path("fixtures", "panse-stan", "phi-prior")
    phi_cfg <- list(gene = list(prior.mu    = .data_spec(pfx, "log_mu"),
                                prior.sigma = "population"))
    out <- tempfile("phi-prior-smoke-dp-")
    fit <- AnaCoDa::fit_panse_stan(.smoke_config(fx, phi_cfg), out_dir = out,
                                   no_log_lik = TRUE, verbose = FALSE)
    draws <- fit$draws(format = "df")
    expect_true(all(is.finite(draws[["sphi"]])))
    expect_true(all(is.finite(draws[["log_phi[1]"]])))
})

test_that("smoke: population/data runs and produces finite draws", {
    skip_if_no_cmdstan()
    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    pfx <- testthat::test_path("fixtures", "panse-stan", "phi-prior")
    phi_cfg <- list(gene = list(prior.mu    = "population",
                                prior.sigma = .data_spec(pfx, "sigma")))
    out <- tempfile("phi-prior-smoke-pd-")
    fit <- AnaCoDa::fit_panse_stan(.smoke_config(fx, phi_cfg), out_dir = out,
                                   no_log_lik = TRUE, verbose = FALSE)
    draws <- fit$draws(format = "df")
    expect_true(all(is.finite(draws[["sphi"]])))
    expect_true(all(is.finite(draws[["log_phi[1]"]])))
})

test_that("smoke: data/data runs and produces finite draws", {
    skip_if_no_cmdstan()
    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    pfx <- testthat::test_path("fixtures", "panse-stan", "phi-prior")
    phi_cfg <- list(gene = list(prior.mu    = .data_spec(pfx, "log_mu"),
                                prior.sigma = .data_spec(pfx, "sigma")))
    out <- tempfile("phi-prior-smoke-dd-")
    fit <- AnaCoDa::fit_panse_stan(.smoke_config(fx, phi_cfg), out_dir = out,
                                   no_log_lik = TRUE, verbose = FALSE)
    draws <- fit$draws(format = "df")
    expect_true(all(is.finite(draws[["log_phi[1]"]])))
    # sphi may be absent if population block not needed for data/data;
    # just check the fit ran and phi draws are finite
})

test_that("smoke: population/population draws match current model (regression)", {
    skip_if_no_cmdstan()
    # Run both old-style config (model: sphi-est) and new-style phi block config.
    # With identical seeds and init, the draws must be bit-for-bit identical.
    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    out_old <- tempfile("phi-prior-smoke-old-")
    out_new <- tempfile("phi-prior-smoke-new-")

    cfg_old <- .smoke_config(fx)
    cfg_new <- .smoke_config(fx)
    cfg_new$phi <- list(gene = list(prior.mu = "population",
                                    prior.sigma = "population"))

    fit_old <- AnaCoDa::fit_panse_stan(cfg_old, out_dir = out_old,
                                       no_log_lik = TRUE, verbose = FALSE)
    fit_new <- AnaCoDa::fit_panse_stan(cfg_new, out_dir = out_new,
                                       no_log_lik = TRUE, verbose = FALSE)

    draws_old <- fit_old$draws(format = "df")
    draws_new <- fit_new$draws(format = "df")

    expect_equal(draws_old[["sphi"]], draws_new[["sphi"]],
                 tolerance = 1e-10,
                 info = "old and new population/population must produce identical sphi draws")
    expect_equal(draws_old[["log_phi[1]"]], draws_new[["log_phi[1]"]],
                 tolerance = 1e-10,
                 info = "old and new population/population must produce identical phi draws")
})
