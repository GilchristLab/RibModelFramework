# PANSE Stan phi prior: Stan-layer correctness.
#
# Verifies that the generalized Stan model (panse_sphi_est_noncentered.stan)
# accepts phi_prior_mu / phi_prior_sigma / phi_use_data in the data block and
# produces finite lp__ for both population/population and data/data configs.
# Also verifies that the two configs yield distinct stan_data (phi_use_data
# differs), which implies distinct Stan log_prob evaluations.
#
# All Stan tests skip if CmdStan is not available.
library(testthat)

context("panse-phi-prior: Stan model correctness")

.base_config <- function(fx) {
    list(
        name   = "phi-prior-logprob-test",
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
        )
    )
}

test_that("population/population stan_data has phi_use_data=0 and finite lp__", {
    skip_if_no_cmdstan()
    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    cfg <- .base_config(fx)
    sd  <- AnaCoDa::build_panse_stan_data(cfg, verbose = FALSE)

    expect_equal(sd$phi_use_data, 0L,
                 info = "population/population must set phi_use_data=0")

    mod <- cmdstanr::cmdstan_model(
        AnaCoDa:::.panse_stan_file("panse_sphi_est_noncentered.stan"),
        dir = tempdir(), quiet = TRUE
    )
    # Quick 1-chain mini-sample to verify the model accepts the new data fields
    suppressMessages(
        fit <- mod$sample(data = sd, chains = 1L,
                          iter_warmup = 5L, iter_sampling = 1L,
                          refresh = 0, show_messages = FALSE, seed = 42L)
    )
    lp <- as.numeric(fit$lp())
    expect_true(any(is.finite(lp)),
                info = "population/population must produce finite lp__")
})

test_that("data/data stan_data has phi_use_data=1 and finite lp__", {
    skip_if_no_cmdstan()
    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    pfx <- testthat::test_path("fixtures", "panse-stan", "phi-prior")
    cfg <- .base_config(fx)
    cfg$phi <- list(
        gene = list(
            prior.mu    = paste0("data(", file.path(pfx, "phi_prior_all.csv"),
                                 ", col=log_mu)"),
            prior.sigma = paste0("data(", file.path(pfx, "phi_prior_all.csv"),
                                 ", col=sigma)")
        )
    )
    sd  <- AnaCoDa::build_panse_stan_data(cfg, verbose = FALSE)

    expect_equal(sd$phi_use_data, 1L,
                 info = "data/data must set phi_use_data=1")

    mod <- cmdstanr::cmdstan_model(
        AnaCoDa:::.panse_stan_file("panse_sphi_est_noncentered.stan"),
        dir = tempdir(), quiet = TRUE
    )
    suppressMessages(
        fit <- mod$sample(data = sd, chains = 1L,
                          iter_warmup = 5L, iter_sampling = 1L,
                          refresh = 0, show_messages = FALSE, seed = 42L)
    )
    lp <- as.numeric(fit$lp())
    expect_true(any(is.finite(lp)),
                info = "data/data must produce finite lp__")
})

test_that("data/data stan_data differs from population/population (phi_use_data)", {
    # Pure R test: the two configs produce different phi_use_data and
    # phi_prior_mu values, which means the Stan model computes different log_prob.
    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    pfx <- testthat::test_path("fixtures", "panse-stan", "phi-prior")

    cfg_pop <- .base_config(fx)
    cfg_dat <- .base_config(fx)
    cfg_dat$phi <- list(
        gene = list(
            prior.mu    = paste0("data(", file.path(pfx, "phi_prior_all.csv"),
                                 ", col=log_mu)"),
            prior.sigma = paste0("data(", file.path(pfx, "phi_prior_all.csv"),
                                 ", col=sigma)")
        )
    )
    sd_pop <- AnaCoDa::build_panse_stan_data(cfg_pop, verbose = FALSE)
    sd_dat <- AnaCoDa::build_panse_stan_data(cfg_dat, verbose = FALSE)

    expect_equal(sd_pop$phi_use_data, 0L)
    expect_equal(sd_dat$phi_use_data, 1L)
    # phi_prior_mu differs: population uses rep(0,G), data uses file values
    expect_false(isTRUE(all.equal(sd_pop$phi_prior_mu, sd_dat$phi_prior_mu)),
                 info = "population and data phi_prior_mu must differ")
    # phi_prior_sigma differs: population uses rep(1,G), data uses file values
    expect_false(isTRUE(all.equal(sd_pop$phi_prior_sigma, sd_dat$phi_prior_sigma)),
                 info = "population and data phi_prior_sigma must differ")
})
