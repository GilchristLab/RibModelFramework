# Migrated PANSE-Stan harness: shared init/codon-bias utility functions.
# Pure R (no cmdstan/AnaCoDa). Validates the functions fit.stan.R links to.
library(testthat)

context("panse-stan-harness: shared utils")

test_that("build_aa_codon_map: 19 groups with Ser split into S(TCN)/Z(AGY)", {
    source_panse_harness("scripts/lib/codon_bias_metrics.R")
    m <- build_aa_codon_map()
    expect_equal(length(m), 19L)                                  # 18 std AAs - Ser + S + Z
    expect_false(any(c("M", "W", "X") %in% names(m)))             # no Met/Trp/stop groups
    expect_setequal(m[["S"]], c("TCT", "TCC", "TCA", "TCG"))      # TCN
    expect_setequal(m[["Z"]], c("AGT", "AGC"))                    # AGY
    expect_setequal(m[["R"]], c("CGT","CGC","CGA","CGG","AGA","AGG"))  # 6-fold, NOT split
    expect_setequal(m[["L"]], c("TTA","TTG","CTT","CTC","CTA","CTG")) # 6-fold, NOT split
    expect_setequal(m[["F"]], c("TTT", "TTC"))
    # every group's codons are unique and the union has no stop codons
    all_codons <- unlist(m, use.names = FALSE)
    expect_equal(anyDuplicated(all_codons), 0L)
    expect_false(any(c("TAA","TAG","TGA") %in% all_codons))
})

test_that("derive_null_from_genome: per-AA freqs sum to 1; uniform fallback on zero", {
    source_panse_harness("scripts/lib/codon_bias_metrics.R")
    m <- build_aa_codon_map()
    codons <- unlist(m, use.names = FALSE)
    # one gene, all counts = 1 except Phe heavily TTC-biased; Pro all zero
    counts <- matrix(1L, nrow = 1, ncol = length(codons), dimnames = list("g1", codons))
    counts[1, "TTC"] <- 9L; counts[1, "TTT"] <- 1L
    counts[1, m[["P"]]] <- 0L
    null <- derive_null_from_genome(counts, m)
    for (aa in names(null)) expect_equal(sum(null[[aa]]), 1, tolerance = 1e-12)
    expect_equal(unname(null[["F"]][c("TTT","TTC")]), c(0.1, 0.9), tolerance = 1e-12)
    expect_equal(unname(null[["P"]]), rep(1/length(m[["P"]]), length(m[["P"]])))  # uniform fallback
})

test_that("scuo_to_log_phi: rank-order-preserving, length-preserving, deterministic", {
    source_panse_harness("scripts/lib/init_helpers.R")
    set.seed(1); scuo <- runif(50)
    lp1 <- scuo_to_log_phi(scuo, sphi_seed = 1.0)
    lp2 <- scuo_to_log_phi(scuo, sphi_seed = 1.0)
    expect_length(lp1, 50)
    expect_true(all(is.finite(lp1)))
    expect_identical(lp1, lp2)                              # deterministic
    expect_equal(order(lp1), order(scuo))                  # monotone in SCUO rank
    # larger sphi_seed -> wider spread of log_phi
    lp_wide <- scuo_to_log_phi(scuo, sphi_seed = 2.0)
    expect_gt(sd(lp_wide), sd(lp1))
})

test_that("mle_to_init_list: groups indexed params into vectors in index order", {
    source_panse_harness("scripts/lib/init_helpers.R")
    df <- data.frame(`dM[2]` = 20, `dM[1]` = 10, sphi = 1.5, junk = 99,
                     check.names = FALSE)
    out <- mle_to_init_list(df, param.names = c("dM", "sphi"))
    expect_setequal(names(out), c("dM", "sphi"))
    expect_equal(out$dM, c(10, 20))      # ordered dM[1], dM[2]
    expect_equal(out$sphi, 1.5)
    expect_false("junk" %in% names(out)) # non-param dropped
})
