# Migrated PANSE-Stan harness: fit.stan.R end-to-end integration.
# Invokes the driver as a subprocess (as in real use) on the tiny builder
# fixture, init.mode=fixed (no AnaCoDa needed), and checks it resolves its
# co-located libs + the RMF stan model, runs, and writes the fit artifacts.
# cmdstan-guarded + slow (compiles a model on first run).
library(testthat)

context("panse-stan-harness: fit.stan.R integration")

test_that("fit.stan.R runs end-to-end from RMF and writes fit artifacts", {
    skip_if_no_cmdstan()      # defined in helper-panse-stan.R (loaded alongside)
    fit_script <- file.path(panse_harness_root(),
                            "scripts", "panse-stan", "fit.stan.R")
    if (!file.exists(fit_script)) skip("fit.stan.R not present")

    fx  <- normalizePath(testthat::test_path("fixtures", "panse-stan", "builder"))
    out <- tempfile("panse-harness-fit-")
    cfg <- tempfile(fileext = ".yaml")
    writeLines(c(
        "name: harness-integration",
        sprintf("genome: { input.dir: %s, pattern: rfp_counts.csv }", fx),
        "fit:",
        "  model: sphi-est-sharednse",
        "  parameterization: noncentered",
        "  with.phi: true",
        "  sphi.prior.sd: 2.5",
        "  init.mode: fixed",
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
        "stan: { chains: 1, parallel_chains: 1, warmup: 10, sampling: 10, threads_per_chain: 1, adapt_delta: 0.8, max_treedepth: 6, metric: diag_e, grainsize: 1, seed: 1 }"
    ), cfg)

    log <- suppressWarnings(system2("Rscript", c(fit_script, cfg, "--out", out, "--no-log-lik"),
                                    stdout = TRUE, stderr = TRUE))
    info <- paste(tail(log, 15), collapse = "\n")
    expect_true(dir.exists(out), info = info)
    expect_true(file.exists(file.path(out, "panse-stan-fit.rds")), info = info)
    expect_true(file.exists(file.path(out, "stan-summary.rds")), info = info)
    expect_true(file.exists(file.path(out, "config.yaml")), info = info)
})
