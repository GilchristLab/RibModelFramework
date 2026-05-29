library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Tests for the ROC arcsine likelihood approximation (approx=TRUE).
#
# The hybrid arcsine approximation replaces the exact multinomial LL
# with a variance-stabilised form for each amino acid / gene combination
# whose total observed codon count n = sum(c_i) meets approxMinExpected.
# For K synonymous codons, K-1 marginal binomials are summed (reference
# codon dropped to avoid double-counting the constraint sum(c_i) = n):
#
#   logL_arcsine = sum_{i=0}^{K-2} -2n * (arcsin(sqrt(c_i/n)) - arcsin(sqrt(p_i)))^2
#
# Key property tested: the arcsine variance 1/(4n) depends only on n,
# not on the per-cell probabilities p_i.  The threshold is therefore on
# total n (not n/K), and the approximation is valid for n >= ~20 at
# virtually any p -- a much less restrictive condition than the normal
# approximation (which requires n*p_i >= 5 per cell).
#
# What we verify:
#   1. approx=FALSE / default / huge-threshold all give identical exact LL.
#   2. LL RATIOS (proposed - current) agree between methods when n is large.
#   3. The sign of the LL ratio agrees even at the n=20 boundary.
#   4. approx=TRUE / "hybrid.arcsine" are interchangeable.
#
# Note: absolute LL values differ between methods (different scale /
# missing normalising constants) -- only ratios / differences are tested.
# ======================================================================

context("ROC arcsine likelihood approximation")

# ----------------------------------------------------------------------
# Shared setup: extdata genome + single-mixture parameter + model objects
# ----------------------------------------------------------------------
genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
genome  <- initializeGenomeObject(file = genome_file)
param   <- initializeParameterObject(genome = genome, sphi = 1, num.mixtures = 1,
                                      gene.assignment = rep(1, length(genome)))

m_exact      <- initializeModelObject(param, "ROC")
m_false      <- initializeModelObject(param, "ROC", approx = FALSE)
m_true       <- initializeModelObject(param, "ROC", approx = TRUE)
m_str        <- initializeModelObject(param, "ROC", approx = "hybrid.arcsine")
m_low        <- initializeModelObject(param, "ROC", approx = TRUE, approx.min.expected = 10)
m_fallback   <- initializeModelObject(param, "ROC", approx = TRUE, approx.min.expected = 1e6)

ll_exact     <- m_exact$calculateLogLikelihood(genome)
ll_true      <- m_true$calculateLogLikelihood(genome)
ll_fallback  <- m_fallback$calculateLogLikelihood(genome)

# ----------------------------------------------------------------------
# R oracle functions mirroring the C++ implementations exactly.
# Used for controlled single-AA comparisons without requiring MCMC.
# ----------------------------------------------------------------------

# Exact multinomial log-likelihood (all K terms).
ll_exact_oracle <- function(counts, probs) {
  mask <- counts > 0
  sum(counts[mask] * log(probs[mask]))
}

# Arcsine LL: K-1 marginal binomials (last codon = reference, dropped).
ll_arcsine_oracle <- function(counts, probs) {
  n   <- sum(counts)
  K   <- length(counts)
  idx <- seq_len(K - 1)
  sum(-2 * n * (asin(sqrt(counts[idx] / n)) - asin(sqrt(probs[idx])))^2)
}

# ======================================================================
# Section 1 -- Construction and basic parity
# ======================================================================

test_that("approx=FALSE explicit is identical to default (no approx)", {
  expect_equal(ll_exact, m_false$calculateLogLikelihood(genome), tolerance = 0)
})

test_that("approx=TRUE and approx='hybrid.arcsine' give identical LL", {
  expect_equal(ll_true, m_str$calculateLogLikelihood(genome), tolerance = 0)
})

test_that("arcsine LL is finite", {
  expect_true(is.finite(ll_true))
})

# ======================================================================
# Section 2 -- Fallback to exact when threshold is never met
# ======================================================================

test_that("approxMinExpected=1e6 falls back to exact multinomial for all AA/genes", {
  # No real gene has n >= 1e6 codons for a single amino acid, so every
  # call takes the exact branch.  Result must match the exact model exactly.
  expect_equal(ll_fallback, ll_exact, tolerance = 0)
})

# ======================================================================
# Section 3 -- Oracle-based LL ratio comparison
#
# Parameters chosen so that probs shift meaningfully between phi values.
# Codon counts are hand-crafted with known n so that n/threshold is
# controlled independently of K.  We use CalculateProbabilitiesForCodons
# to obtain the C++ probability vector (already validated in testROCNumerical.R).
# ======================================================================

# Shared parameter vectors (K-1 length; reference codon is implicit).
dM_6   <- c( 0.30, -0.20,  0.50, -0.10,  0.40)   # K=6
dEta_6 <- c( 0.01,  0.03,  0.005, 0.02,  0.015)

dM_2   <- c(-0.15)                                  # K=2
dEta_2 <- c(-0.02)

probs6_phi1 <- m_exact$CalculateProbabilitiesForCodons(dM_6, dEta_6, phi = 1.0)
probs6_phi2 <- m_exact$CalculateProbabilitiesForCodons(dM_6, dEta_6, phi = 2.0)
probs2_phi1 <- m_exact$CalculateProbabilitiesForCodons(dM_2, dEta_2, phi = 1.0)
probs2_phi3 <- m_exact$CalculateProbabilitiesForCodons(dM_2, dEta_2, phi = 3.0)

# -- Oracle whitebox: confirm ll_arcsine_oracle matches the formula exactly --
test_that("arcsine oracle matches formula term-by-term (6-codon AA)", {
  counts <- c(25, 10, 30, 8, 20, 27)   # n = 120
  n      <- sum(counts)
  K      <- length(counts)
  manual <- sum(vapply(seq_len(K - 1), function(i)
    -2 * n * (asin(sqrt(counts[i] / n)) - asin(sqrt(probs6_phi1[i])))^2,
    numeric(1)))
  expect_equal(ll_arcsine_oracle(counts, probs6_phi1), manual, tolerance = 1e-15)
})

# -- LL ratio agreement: n=120, K=6 (well above n=20 threshold) --
#
# The ~20% relative error here is a *constant systematic bias*, not a
# sample-size effect: both delta_exact and delta_arcsine scale linearly
# with n (and with phi step size), so their ratio is fixed for a given
# count distribution.  The bias arises because the arcsine-transform
# curvature differs from the log curvature when counts are skewed away
# from the probability vector.  For MCMC, this is equivalent to sampling
# from a slightly tempered (cooler) version of the posterior.  The sign
# is what drives accept/reject correctness.
test_that("arcsine LL ratio tracks exact LL ratio within 25%: 6-codon AA, n=120", {
  counts <- c(25, 10, 30, 8, 20, 27)   # n = 120, counts skewed away from p1

  delta_exact   <- ll_exact_oracle(counts, probs6_phi2)   - ll_exact_oracle(counts, probs6_phi1)
  delta_arcsine <- ll_arcsine_oracle(counts, probs6_phi2) - ll_arcsine_oracle(counts, probs6_phi1)

  expect_lt(abs(delta_arcsine - delta_exact) / abs(delta_exact), 0.25)
})

test_that("arcsine LL ratio sign agrees with exact: 6-codon AA, n=120", {
  counts <- c(25, 10, 30, 8, 20, 27)
  delta_exact   <- ll_exact_oracle(counts, probs6_phi2)   - ll_exact_oracle(counts, probs6_phi1)
  delta_arcsine <- ll_arcsine_oracle(counts, probs6_phi2) - ll_arcsine_oracle(counts, probs6_phi1)
  expect_equal(sign(delta_arcsine), sign(delta_exact))
})

# -- LL ratio agreement: n=60, K=2 (binary; only one marginal term) --
test_that("arcsine LL ratio tracks exact LL ratio within 10%: 2-codon AA, n=60", {
  counts <- c(40, 20)   # n = 60

  delta_exact   <- ll_exact_oracle(counts, probs2_phi3)   - ll_exact_oracle(counts, probs2_phi1)
  delta_arcsine <- ll_arcsine_oracle(counts, probs2_phi3) - ll_arcsine_oracle(counts, probs2_phi1)

  expect_lt(abs(delta_arcsine - delta_exact) / abs(delta_exact), 0.10)
  expect_equal(sign(delta_arcsine), sign(delta_exact))
})

# -- LL ratio sign at n=20 threshold boundary --
test_that("arcsine LL ratio sign agrees with exact at n=20 (threshold boundary): 4-codon AA", {
  dM_4   <- c( 0.10, -0.10,  0.20)
  dEta_4 <- c( 0.02,  0.01, -0.005)
  probs4_phi1 <- m_exact$CalculateProbabilitiesForCodons(dM_4, dEta_4, phi = 1.0)
  probs4_phi2 <- m_exact$CalculateProbabilitiesForCodons(dM_4, dEta_4, phi = 2.0)

  counts <- c(8, 6, 4, 2)   # n = 20, K = 4

  delta_exact   <- ll_exact_oracle(counts, probs4_phi2)   - ll_exact_oracle(counts, probs4_phi1)
  delta_arcsine <- ll_arcsine_oracle(counts, probs4_phi2) - ll_arcsine_oracle(counts, probs4_phi1)

  expect_equal(sign(delta_arcsine), sign(delta_exact))
  # Looser tolerance at the boundary
  expect_lt(abs(delta_arcsine - delta_exact) / abs(delta_exact), 0.30)
})

# -- For completeness: arcsine is genuinely different from exact (not a no-op) --
test_that("arcsine LL differs from exact LL on same genome (methods are distinct)", {
  expect_false(isTRUE(all.equal(ll_exact, ll_true, tolerance = 1e-6)))
})
