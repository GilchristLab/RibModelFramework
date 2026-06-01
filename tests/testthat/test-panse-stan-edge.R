# PANSE Stan port -- edge cases: zeros, masked positions, single-codon gene.
# Deterministic via generated quantities at truth.
library(testthat)

context("panse-stan: edge cases")

test_that("panse-stan: all-zero counts in a gene give finite log_lik", {
    skip_if_no_cmdstan()
    fx  <- load_panse_stan_fixture()
    sd  <- panse_stan_data_for_fit(fx)
    mod <- panse_stan_model("panse_sphi_est_noncentered_sharednse")
    dr  <- panse_stan_truth_draws(fx)

    p0 <- sd$gene_offset[1]; p1 <- sd$gene_offset[2] - 1
    sd$y[p0:p1] <- 0L
    gq <- mod$generate_quantities(fitted_params = dr, data = sd,
                                  threads_per_chain = 1, parallel_chains = 1)
    ll <- as.numeric(gq$draws("log_lik"))
    expect_true(all(is.finite(ll)))
    expect_true(all(ll[p0:p1] <= 0))
})

test_that("panse-stan: masked positions contribute exactly 0 to log_lik", {
    skip_if_no_cmdstan()
    fx  <- load_panse_stan_fixture()
    sd  <- panse_stan_data_for_fit(fx)
    mod <- panse_stan_model("panse_sphi_est_noncentered_sharednse")
    dr  <- panse_stan_truth_draws(fx)

    p0 <- sd$gene_offset[1]; p1 <- sd$gene_offset[2] - 1
    sd$like_mask[p0:p1] <- 0L
    sd$all_unmasked <- 0L
    gq <- mod$generate_quantities(fitted_params = dr, data = sd,
                                  threads_per_chain = 1, parallel_chains = 1)
    ll <- as.numeric(gq$draws("log_lik"))
    expect_true(all(is.finite(ll)))
    expect_equal(ll[p0:p1], rep(0, p1 - p0 + 1))   # masked -> log_lik 0
    expect_true(all(ll[(p1 + 1):sd$P] <= 0))        # unmasked rest still contributes
})

test_that("panse-stan: a single-codon gene is handled (survival has no prior position)", {
    skip_if_no_cmdstan()
    fx  <- load_panse_stan_fixture()
    sd  <- panse_stan_data_for_fit(fx)
    mod <- panse_stan_model("panse_sphi_est_noncentered_sharednse")
    dr  <- panse_stan_truth_draws(fx)

    # Truncate gene 1 to a single position by masking all but its first.
    p0 <- sd$gene_offset[1]; p1 <- sd$gene_offset[2] - 1
    if (p1 > p0) { sd$like_mask[(p0 + 1):p1] <- 0L; sd$all_unmasked <- 0L }
    gq <- mod$generate_quantities(fitted_params = dr, data = sd,
                                  threads_per_chain = 1, parallel_chains = 1)
    ll <- as.numeric(gq$draws("log_lik"))
    expect_true(is.finite(ll[p0]))          # first position: survival = 1 (log 0)
    expect_lte(ll[p0], 0)
})
