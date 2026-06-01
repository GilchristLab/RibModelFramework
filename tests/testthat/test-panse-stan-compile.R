# PANSE Stan port -- every production model compiles and runs a few iterations.
library(testthat)

context("panse-stan: compile + smoke")

test_that("panse-stan: all production models compile and sample a few iters without error", {
    skip_if_no_cmdstan()
    fx <- load_panse_stan_fixture()

    for (nm in PANSE_STAN_MODELS) {
        sd  <- panse_stan_data_for(fx, nm)
        mod <- panse_stan_model(nm)
        fit <- tryCatch(
            mod$sample(data = sd, chains = 1, iter_warmup = 5, iter_sampling = 5,
                       threads_per_chain = 1, refresh = 0,
                       show_messages = FALSE, show_exceptions = FALSE, seed = 1),
            error = function(e) e
        )
        expect_false(inherits(fit, "error"),
                     info = sprintf("model %s failed to sample: %s", nm,
                                    if (inherits(fit, "error")) conditionMessage(fit) else ""))
        if (!inherits(fit, "error"))
            expect_gt(posterior::ndraws(fit$draws()), 0)
    }
})
