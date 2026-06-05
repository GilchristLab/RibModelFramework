# PANSE Stan phi prior: config parsing layer.
#
# Verifies that build_panse_stan_data() correctly populates phi_prior_mu[G]
# and phi_prior_sigma[G] for all four source combinations:
#   population/population -- both from hierarchical hyperparameters
#   data/population       -- gene-specific mu from file, shared sphi
#   population/data       -- shared mphi, gene-specific sigma from file
#   data/data             -- both gene-specific from file
#
# These are pure R-layer tests; no Stan compilation needed.
library(testthat)

context("panse-phi-prior: build_panse_stan_data phi prior vectors")

.base_config <- function(fx) {
    list(
        name   = "phi-prior-build-test",
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

test_that("population/population: phi_prior_mu and sigma are G-length shared vectors", {
    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    cfg <- .base_config(fx)
    # population/population is the default when no phi block overrides
    sd  <- AnaCoDa::build_panse_stan_data(cfg, verbose = FALSE)

    expect_true(!is.null(sd$phi_prior_mu),
                info = "phi_prior_mu missing from stan_data")
    expect_true(!is.null(sd$phi_prior_sigma),
                info = "phi_prior_sigma missing from stan_data")
    expect_length(sd$phi_prior_mu,    sd$G)
    expect_length(sd$phi_prior_sigma, sd$G)

    # All entries equal: shared hyperparameters broadcast to G genes
    expect_equal(length(unique(sd$phi_prior_mu)),    1L,
                 info = "population mu should be identical across genes")
    expect_equal(length(unique(sd$phi_prior_sigma)), 1L,
                 info = "population sigma should be identical across genes")

    # Sigma entries are positive
    expect_true(all(sd$phi_prior_sigma > 0))
})

test_that("population/population: model compiles and samples with new data fields", {
    skip_if_no_cmdstan()

    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    cfg <- .base_config(fx)
    sd  <- AnaCoDa::build_panse_stan_data(cfg, verbose = FALSE)

    # phi_use_data=0 path must not break the model: compile and run a mini sample
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
                info = "population/population must produce finite lp__")
})

test_that("data/population: gene-specific mu from file; sigma shared", {
    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    pfx <- testthat::test_path("fixtures", "panse-stan", "phi-prior")
    cfg <- .base_config(fx)
    cfg$phi <- list(
        gene = list(
            prior.mu    = paste0("data(", file.path(pfx, "phi_prior_all.csv"),
                                 ", col=log_mu)"),
            prior.sigma = "population"
        )
    )
    sd <- AnaCoDa::build_panse_stan_data(cfg, verbose = FALSE)

    expect_length(sd$phi_prior_mu,    sd$G)
    expect_length(sd$phi_prior_sigma, sd$G)

    # mu values should vary across genes (gene-specific from file)
    expect_gt(length(unique(sd$phi_prior_mu)), 1L,
              label = "data mu should differ across genes")

    # sigma should be shared (population)
    expect_equal(length(unique(sd$phi_prior_sigma)), 1L,
                 label = "population sigma should be identical across genes")

    # Check specific values from fixture: g1=-0.347, g2=0.693, g3=-0.693
    gene_ids <- attr(sd, "gene_ids")
    expect_equal(sd$phi_prior_mu[gene_ids == "g1"], -0.347, tolerance = 1e-6)
    expect_equal(sd$phi_prior_mu[gene_ids == "g2"],  0.693, tolerance = 1e-6)
    expect_equal(sd$phi_prior_mu[gene_ids == "g3"], -0.693, tolerance = 1e-6)
})

test_that("population/data: shared mu; gene-specific sigma from file", {
    fx  <- testthat::test_path("fixtures", "panse-stan", "builder")
    pfx <- testthat::test_path("fixtures", "panse-stan", "phi-prior")
    cfg <- .base_config(fx)
    cfg$phi <- list(
        gene = list(
            prior.mu    = "population",
            prior.sigma = paste0("data(", file.path(pfx, "phi_prior_all.csv"),
                                 ", col=sigma)")
        )
    )
    sd <- AnaCoDa::build_panse_stan_data(cfg, verbose = FALSE)

    expect_length(sd$phi_prior_mu,    sd$G)
    expect_length(sd$phi_prior_sigma, sd$G)

    # mu shared
    expect_equal(length(unique(sd$phi_prior_mu)), 1L)

    # sigma varies (gene-specific from file)
    expect_gt(length(unique(sd$phi_prior_sigma)), 1L)

    # Check specific values: g1=0.8, g2=0.5, g3=1.2
    gene_ids <- attr(sd, "gene_ids")
    expect_equal(sd$phi_prior_sigma[gene_ids == "g1"], 0.8, tolerance = 1e-6)
    expect_equal(sd$phi_prior_sigma[gene_ids == "g2"], 0.5, tolerance = 1e-6)
    expect_equal(sd$phi_prior_sigma[gene_ids == "g3"], 1.2, tolerance = 1e-6)
    expect_true(all(sd$phi_prior_sigma > 0))
})

test_that("data/data: both gene-specific from file", {
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
    sd <- AnaCoDa::build_panse_stan_data(cfg, verbose = FALSE)

    expect_length(sd$phi_prior_mu,    sd$G)
    expect_length(sd$phi_prior_sigma, sd$G)

    gene_ids <- attr(sd, "gene_ids")
    expect_equal(sd$phi_prior_mu[gene_ids == "g1"],    -0.347, tolerance = 1e-6)
    expect_equal(sd$phi_prior_sigma[gene_ids == "g1"],  0.8,   tolerance = 1e-6)
    expect_equal(sd$phi_prior_mu[gene_ids == "g2"],     0.693, tolerance = 1e-6)
    expect_equal(sd$phi_prior_sigma[gene_ids == "g2"],  0.5,   tolerance = 1e-6)
    expect_gt(length(unique(sd$phi_prior_mu)),    1L)
    expect_gt(length(unique(sd$phi_prior_sigma)), 1L)
    expect_true(all(sd$phi_prior_sigma > 0))
})

test_that("gene order in phi_prior_mu/sigma matches gene_ids attribute", {
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
    sd <- AnaCoDa::build_panse_stan_data(cfg, verbose = FALSE)

    # phi file has genes in g1, g2, g3 order; stan_data gene order may differ.
    # The R layer must match by gene_id, not by row position.
    gene_ids <- attr(sd, "gene_ids")
    phi_tbl  <- read.csv(file.path(pfx, "phi_prior_all.csv"),
                         stringsAsFactors = FALSE)
    phi_map  <- setNames(phi_tbl$log_mu, phi_tbl$gene_id)
    expected <- unname(phi_map[gene_ids])
    expect_equal(sd$phi_prior_mu, expected, tolerance = 1e-6)
})
