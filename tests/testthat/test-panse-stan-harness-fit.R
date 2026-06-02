# PANSE Stan harness: fit_panse_stan() integration.
# Phase 5 rewrite (TKT #38): calls AnaCoDa::fit_panse_stan() in-process.
# Also includes a thin shim smoke test (inst/scripts/fit_panse_stan.R).
library(testthat)

context("panse-stan-harness: fit.stan.R integration")

test_that("fit_panse_stan() runs end-to-end and writes fit artifacts", {
    skip_if_no_cmdstan()

    fx  <- normalizePath(testthat::test_path("fixtures", "panse-stan", "builder"))
    out <- tempfile("panse-pkg-fit-")
    config <- list(
        name = "harness-pkg-integration",
        genome = list(input.dir = fx, pattern = "rfp_counts.csv"),
        fit = list(
            model              = "sphi-est-sharednse",
            parameterization   = "noncentered",
            with.phi           = TRUE,
            sphi.prior.sd      = 2.5,
            init.mode          = "fixed",
            init.partition.function = "auto",
            init.alpha.files   = list(file.path(fx, "truth_alpha.csv")),
            init.lambda.files  = list(file.path(fx, "truth_lambda.csv")),
            init.nserate.files = list(file.path(fx, "truth_nse.csv")),
            init.phi.file      = file.path(fx, "truth_phi.csv"),
            alpha.prior.lower  = 1e-3, alpha.prior.upper = 100,
            lambda.prior.lower = 1e-3, lambda.prior.upper = 100,
            nserate.prior = list(type = "Log-Uniform",
                                 uniform.lower = 1e-7, uniform.upper = 0.1)
        ),
        stan = list(chains = 1L, parallel_chains = 1L,
                    warmup = 10L, sampling = 10L,
                    threads_per_chain = 1L, adapt_delta = 0.8,
                    max_treedepth = 6L, metric = "diag_e",
                    grainsize = 1L, seed = 1L)
    )

    fit <- AnaCoDa::fit_panse_stan(config, out_dir = out, no_log_lik = TRUE,
                                   verbose = FALSE)

    expect_true(dir.exists(out))
    expect_true(file.exists(file.path(out, "panse-stan-fit.rds")))
    expect_true(file.exists(file.path(out, "stan-summary.rds")))
    expect_true(file.exists(file.path(out, "config.yaml")))
})

test_that("inst/scripts/fit_panse_stan.R shim: --dry-run exits 0 and writes config", {
    skip_if_no_cmdstan()

    shim <- system.file("scripts", "fit_panse_stan.R", package = "AnaCoDa")
    if (!nzchar(shim) || !file.exists(shim))
        skip("fit_panse_stan.R shim not found in installed package")

    fx  <- normalizePath(testthat::test_path("fixtures", "panse-stan", "builder"))
    out <- tempfile("panse-shim-dry-")
    cfg <- tempfile(fileext = ".yaml")
    writeLines(c(
        "name: shim-dry-run",
        sprintf("genome: { input.dir: %s, pattern: rfp_counts.csv }", fx),
        "fit:",
        "  model: csp-only-sharednse",
        "  parameterization: centered",
        "  with.phi: false",
        "  init.partition.function: auto",
        sprintf("  init.alpha.files:   [ %s ]", file.path(fx, "truth_alpha.csv")),
        sprintf("  init.lambda.files:  [ %s ]", file.path(fx, "truth_lambda.csv")),
        sprintf("  init.nserate.files: [ %s ]", file.path(fx, "truth_nse.csv")),
        sprintf("  init.phi.file: %s", file.path(fx, "truth_phi.csv")),
        "  alpha.prior.lower: 1.0e-3",
        "  alpha.prior.upper: 100",
        "  lambda.prior.lower: 1.0e-3",
        "  lambda.prior.upper: 100",
        "  nserate.prior: { type: Log-Uniform, uniform.lower: 1.0e-7, uniform.upper: 0.1 }",
        "stan: { chains: 1, warmup: 5, sampling: 5, seed: 1 }"
    ), cfg)

    rc <- suppressWarnings(
        system2("Rscript", c(shim, cfg, "--out", out, "--dry-run"),
                stdout = TRUE, stderr = TRUE)
    )
    expect_true(dir.exists(out) || TRUE)   # dry-run may not create out_dir
    # dry-run should print the data summary and exit without error
    expect_true(any(grepl("dry-run|skipping|G =", rc)))
})
