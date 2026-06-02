# PANSE Stan harness: build_panse_stan_data() data-list contract.
# Phase 5 rewrite (TKT #38): calls AnaCoDa::build_panse_stan_data() directly.
library(testthat)

context("panse-stan-harness: build_panse_stan_data")

test_that("build_panse_stan_data produces the expected Stan data contract", {
    fx <- testthat::test_path("fixtures", "panse-stan", "builder")
    config <- list(
        genome = list(input.dir = fx, pattern = "rfp_counts.csv"),
        fit = list(
            init.alpha.files   = list(file.path(fx, "truth_alpha.csv")),
            init.lambda.files  = list(file.path(fx, "truth_lambda.csv")),
            init.nserate.files = list(file.path(fx, "truth_nse.csv")),
            init.phi.file      = file.path(fx, "truth_phi.csv"),
            alpha.prior.lower = 1e-3, alpha.prior.upper = 100,
            lambda.prior.lower = 1e-3, lambda.prior.upper = 100,
            nserate.prior = list(type = "Log-Uniform",
                                 uniform.lower = 1e-7, uniform.upper = 0.1),
            sphi.prior.sd = 2.5, init.partition.function = "auto"
        )
    )
    sd <- AnaCoDa::build_panse_stan_data(config, verbose = FALSE)

    # dimensions
    expect_equal(sd$G, 3L)
    expect_equal(sd$C, 4L)
    expect_equal(sd$P, 12L)

    # CSR layout: gene_offset is 1-indexed, monotone, length G+1, spans P
    expect_length(sd$gene_offset, 4L)
    expect_equal(sd$gene_offset, c(1L, 5L, 9L, 13L))
    expect_equal(sd$gene_offset[1], 1L)
    expect_equal(sd$gene_offset[sd$G + 1L], sd$P + 1L)

    # codon ids in range; counts preserved
    expect_true(all(sd$codon_at_pos >= 1L & sd$codon_at_pos <= sd$C))
    expect_length(sd$y, 12L)
    expect_equal(sum(sd$y), 36L)

    # codon ordering from alpha CSV
    expect_equal(names(attr(sd, "init_alpha")), c("TTT", "TTC", "TAT", "TAC"))
    expect_equal(unname(attr(sd, "init_alpha")[["TTT"]]), 1.5)

    # partition function positive; no Mixture column -> all positions in likelihood
    expect_gt(sd$U, 0)
    expect_true(all(sd$like_mask == 1L))
    expect_true(sd$all_unmasked == 1L)
})
