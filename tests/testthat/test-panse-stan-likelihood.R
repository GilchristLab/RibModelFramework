# PANSE Stan port -- likelihood agreement against an independent R reference.
# Deterministic (generated quantities at truth; no MCMC). Covers test-plan
# items 2 (likelihood agreement) and 3 (recovery, via the committed fixture).
library(testthat)

context("panse-stan: likelihood")

test_that("panse-stan: generated-quantities log_lik matches the R reference at truth", {
    skip_if_no_cmdstan()
    fx  <- load_panse_stan_fixture()
    mod_name <- "panse_sphi_est_noncentered_sharednse"
    sd  <- panse_stan_data_for(fx, mod_name)
    mod <- panse_stan_model(mod_name)
    dr  <- panse_stan_truth_draws(fx)

    gq <- mod$generate_quantities(fitted_params = dr, data = sd,
                                  threads_per_chain = 1, parallel_chains = 1)
    ll <- as.numeric(gq$draws("log_lik"))

    expect_equal(length(ll), sd$P)
    expect_true(all(is.finite(ll)))
    expect_equal(sum(ll), sum(fx$ref_log_lik), tolerance = 1e-4)
    expect_lt(max(abs(ll - fx$ref_log_lik)), 1e-5)
})
