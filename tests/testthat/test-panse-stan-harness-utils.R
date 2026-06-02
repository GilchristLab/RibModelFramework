# PANSE Stan harness: shared init/codon-bias utility functions.
# Phase 5 rewrite (TKT #38): calls AnaCoDa::* directly; no source_panse_harness().
library(testthat)

context("panse-stan-harness: shared utils")

test_that("build_aa_codon_map: 19 groups with Ser split into S(TCN)/Z(AGY)", {
    m <- AnaCoDa::build_aa_codon_map()
    expect_equal(length(m), 19L)
    expect_false(any(c("M", "W", "X") %in% names(m)))
    expect_setequal(m[["S"]], c("TCT", "TCC", "TCA", "TCG"))
    expect_setequal(m[["Z"]], c("AGT", "AGC"))
    expect_setequal(m[["R"]], c("CGT","CGC","CGA","CGG","AGA","AGG"))
    expect_setequal(m[["L"]], c("TTA","TTG","CTT","CTC","CTA","CTG"))
    expect_setequal(m[["F"]], c("TTT", "TTC"))
    all_codons <- unlist(m, use.names = FALSE)
    expect_equal(anyDuplicated(all_codons), 0L)
    expect_false(any(c("TAA","TAG","TGA") %in% all_codons))
})

test_that("derive_null_from_genome: per-AA freqs sum to 1; uniform fallback on zero", {
    m <- AnaCoDa::build_aa_codon_map()
    codons <- unlist(m, use.names = FALSE)
    counts <- matrix(1L, nrow = 1, ncol = length(codons), dimnames = list("g1", codons))
    counts[1, "TTC"] <- 9L; counts[1, "TTT"] <- 1L
    counts[1, m[["P"]]] <- 0L
    null <- AnaCoDa::derive_null_from_genome(counts, m)
    for (aa in names(null)) expect_equal(sum(null[[aa]]), 1, tolerance = 1e-12)
    expect_equal(unname(null[["F"]][c("TTT","TTC")]), c(0.1, 0.9), tolerance = 1e-12)
    expect_equal(unname(null[["P"]]), rep(1/length(m[["P"]]), length(m[["P"]])))
})

test_that("scuo_to_log_phi: rank-order-preserving, length-preserving, deterministic", {
    set.seed(1); scuo <- runif(50)
    lp1 <- AnaCoDa::scuo_to_log_phi(scuo, sphi_seed = 1.0)
    lp2 <- AnaCoDa::scuo_to_log_phi(scuo, sphi_seed = 1.0)
    expect_length(lp1, 50)
    expect_true(all(is.finite(lp1)))
    expect_identical(lp1, lp2)
    expect_equal(order(lp1), order(scuo))
    lp_wide <- AnaCoDa::scuo_to_log_phi(scuo, sphi_seed = 2.0)
    expect_gt(sd(lp_wide), sd(lp1))
})

test_that("mle_to_init_list: groups indexed params into vectors in index order", {
    df <- data.frame(`dM[2]` = 20, `dM[1]` = 10, sphi = 1.5, junk = 99,
                     check.names = FALSE)
    out <- AnaCoDa::mle_to_init_list(df, param.names = c("dM", "sphi"))
    expect_setequal(names(out), c("dM", "sphi"))
    expect_equal(out$dM, c(10, 20))
    expect_equal(out$sphi, 1.5)
    expect_false("junk" %in% names(out))
})
