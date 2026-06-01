# PANSE Stan port -- survival-term (log_psuccess) correctness + #24 regression.
library(testthat)

context("panse-stan: survival term")

test_that("panse-stan: log_psuccess stays <= 0 across the prior grid incl. the q>q* breach (#24)", {
    skip_if_no_cmdstan()
    mod <- panse_stan_functions()

    # Grid spans published codon range AND well past the breach q* ~ 0.78-0.99,
    # where the OLD 2nd-order Taylor went positive (psuccess > 1). The hybrid
    # must switch to the closed form and remain bounded <= 0 everywhere.
    grid <- expand.grid(
        alpha  = c(1e-2, 1, 1.5, 2, 5, 50),
        lambda = c(0.035, 0.1, 0.5),
        nse    = c(1e-7, 1e-5, 1e-3, 1e-2, 0.1, 0.3)
    )
    breach_seen <- FALSE
    for (i in seq_len(nrow(grid))) {
        a <- grid$alpha[i]; l <- grid$lambda[i]; n <- grid$nse[i]
        q <- a * n / l
        stan_v <- mod$functions$log_psuccess_hybrid(a, l, n)
        r_v    <- ref_log_psuccess(a, l, n)
        expect_lte(stan_v, 1e-9)                       # psuccess <= 1, physical
        expect_equal(stan_v, r_v, tolerance = 1e-6,
                     info = sprintf("alpha=%g lambda=%g nse=%g q=%g", a, l, n, q))
        if (q > 1) breach_seen <- TRUE
    }
    expect_true(breach_seen)  # grid actually exercises the breach region
})

test_that("panse-stan: log_upper_incomplete_gamma matches the R reference (incl. s<=0)", {
    skip_if_no_cmdstan()
    mod <- panse_stan_functions()
    cases <- list(c(1 - 1.5, 3.5), c(1 - 2.0, 10), c(0, 5), c(1 - 0.5, 0.2), c(1 - 50, 700))
    for (cs in cases) {
        s <- cs[1]; x <- cs[2]
        expect_equal(mod$functions$log_upper_incomplete_gamma(s, x),
                     .panse_log_uig(s, x), tolerance = 1e-6,
                     info = sprintf("s=%g x=%g", s, x))
    }
})
