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

# ======================================================================
# Section 4 -- genomeToStanData() shape and consistency checks
# ======================================================================

context("genomeToStanData: Stan data list preparation")

stan_d <- genomeToStanData(genome)

test_that("genomeToStanData returns the required fields", {
  required <- c("G", "A", "K", "aa_start", "aa_end", "y_k", "N_ga",
                "approx_min_n",
                "dM_prior_mean", "dM_prior_sd",
                "dEta_prior_mean", "dEta_prior_sd",
                "sphi_prior_mean", "sphi_prior_sd",
                "noncentered", "anchor_phi", "mphi_prior_sd",
                "grainsize")
  expect_true(all(required %in% names(stan_d)))
})

test_that("A = 19 and K = 40 (standard genetic code, ROCParameter::groupList)", {
  expect_equal(stan_d$A, 19L)
  expect_equal(stan_d$K, 40L)
})

test_that("y_k has shape G x K and is integer", {
  expect_equal(dim(stan_d$y_k), c(stan_d$G, 40L))
  expect_true(is.integer(stan_d$y_k))
})

test_that("N_ga has shape G x A and is integer", {
  expect_equal(dim(stan_d$N_ga), c(stan_d$G, 19L))
  expect_true(is.integer(stan_d$N_ga))
})

test_that("aa_start[1] = 1 and aa_end[A] = K = 40", {
  expect_equal(stan_d$aa_start[1], 1L)
  expect_equal(stan_d$aa_end[stan_d$A], 40L)
})

test_that("aa ranges are non-empty and cover [1, K] without gaps", {
  # Each AA has at least one non-reference codon
  expect_true(all(stan_d$aa_end >= stan_d$aa_start))
  # Ranges tile [1, K] exactly (no overlaps, no gaps)
  expect_equal(sum(stan_d$aa_end - stan_d$aa_start + 1L), 40L)
})

test_that("non-ref codon counts are <= total AA counts (conservation)", {
  K   <- stan_d$K
  A   <- stan_d$A
  for (a in seq_len(A)) {
    ks       <- stan_d$aa_start[a]:stan_d$aa_end[a]
    nonref_g <- rowSums(stan_d$y_k[, ks, drop = FALSE])
    expect_true(all(nonref_g <= stan_d$N_ga[, a]),
                info = sprintf("AA index %d: non-ref counts exceed N_ga", a))
  }
})

test_that("prior parameter vectors are length K", {
  K <- stan_d$K
  expect_equal(length(stan_d$dM_prior_mean),   K)
  expect_equal(length(stan_d$dM_prior_sd),     K)
  expect_equal(length(stan_d$dEta_prior_mean), K)
  expect_equal(length(stan_d$dEta_prior_sd),   K)
})

test_that("approx_min_n defaults to 20L (integer)", {
  expect_identical(stan_d$approx_min_n, 20L)
})

test_that("scalar prior expansion works correctly", {
  d2 <- genomeToStanData(genome, dM_prior_sd = 0.5)
  expect_true(all(d2$dM_prior_sd == 0.5))
  expect_equal(length(d2$dM_prior_sd), 40L)
})

test_that("length-K prior vector is accepted as-is", {
  custom_mean <- seq(-1, 1, length.out = 40L)
  d3 <- genomeToStanData(genome, dM_prior_mean = custom_mean)
  expect_equal(d3$dM_prior_mean, custom_mean)
})

test_that("wrong-length prior vector raises an error", {
  expect_error(genomeToStanData(genome, dM_prior_sd = rep(1, 5)),
               regexp = "length")
})

# ======================================================================
# Section 5 -- dM prior from SCUO and genomeToStanInit()
# ======================================================================

context("genomeToStanData: SCUO-based dM prior")

# Pre-compute SCUO once (reused across tests)
scuo_vals <- calculateSCUO(genome)$SCUO

test_that("dM.prior='scuo' (default) gives non-zero dM_prior_mean", {
  # At least some dM priors should be non-zero when estimated from data
  expect_false(all(stan_d$dM_prior_mean == 0))
  expect_equal(length(stan_d$dM_prior_mean), 40L)
  expect_true(all(is.finite(stan_d$dM_prior_mean)))
})

test_that("dM.prior='flat' gives all-zero dM_prior_mean", {
  d_flat <- genomeToStanData(genome, dM.prior = "flat")
  expect_true(all(d_flat$dM_prior_mean == 0))
})

test_that("explicit dM_prior_mean overrides dM.prior", {
  custom <- rep(0.5, 40L)
  d_ov   <- genomeToStanData(genome, dM_prior_mean = custom, dM.prior = "scuo")
  expect_equal(d_ov$dM_prior_mean, custom)
})

test_that("pre-computed scuo arg gives same result as internal computation", {
  d_pre  <- genomeToStanData(genome, scuo = scuo_vals)
  d_auto <- genomeToStanData(genome)
  expect_equal(d_pre$dM_prior_mean, d_auto$dM_prior_mean)
})

test_that("SCUO-based dM prior has correct sign: rare codons get positive dM", {
  # dM_k = log(count_ref) - log(count_k).
  # Reference codons (e.g. GCT for Ala) tend to be most common in low-expression
  # genes, so most non-ref codons should have positive dM (= lower usage than ref).
  # Not all will be positive (mutation patterns vary by AA), but the median should be.
  expect_gt(median(stan_d$dM_prior_mean), 0)
})

context("genomeToStanInit: phi initialisation")

stan_init <- genomeToStanInit(genome, stan_d, scuo = scuo_vals)

test_that("genomeToStanInit returns the required parameter fields", {
  required <- c("dM", "dEta", "latent_phi", "sphi", "mphi_param")
  expect_true(all(required %in% names(stan_init)))
})

test_that("latent_phi has length G", {
  expect_equal(length(stan_init$latent_phi), stan_d$G)
})

test_that("dM and dEta init are zero vectors of length K", {
  expect_equal(stan_init$dM,   rep(0, stan_d$K))
  expect_equal(stan_init$dEta, rep(0, stan_d$K))
})

test_that("SCUO init is NOT all zeros (uses actual codon bias)", {
  expect_false(all(stan_init$latent_phi == 0))
})

test_that("SCUO init latent_phi is centered near 0 (mean close to 0)", {
  expect_lt(abs(mean(stan_init$latent_phi)), 0.5)
})

test_that("SCUO init sd ~ sphi.init (default 1.0)", {
  # With G=8 the sd won't be exactly 1, but should be in a reasonable range
  expect_gt(sd(stan_init$latent_phi), 0.3)
  expect_lt(sd(stan_init$latent_phi), 3.0)
})

test_that("phi.init='uniform' gives latent_phi = 0 everywhere", {
  init_u <- genomeToStanInit(genome, stan_d, phi.init = "uniform")
  expect_equal(init_u$latent_phi, rep(0.0, stan_d$G))
})

test_that("numeric phi.init vector is accepted and log-transformed", {
  phi_vec  <- rep(2.0, stan_d$G)   # all genes at phi=2
  init_phi <- genomeToStanInit(genome, stan_d, phi.init = phi_vec)
  # All genes same phi -> after centering, latent_phi = 0
  expect_equal(init_phi$latent_phi, rep(0.0, stan_d$G), tolerance = 1e-10)
})

test_that("pre-computed scuo gives same init as internal computation", {
  init_pre  <- genomeToStanInit(genome, stan_d, scuo = scuo_vals)
  init_auto <- genomeToStanInit(genome, stan_d)
  expect_equal(init_pre$latent_phi, init_auto$latent_phi)
})

test_that("sphi.init controls the spread of latent_phi", {
  init_narrow <- genomeToStanInit(genome, stan_d, scuo = scuo_vals, sphi.init = 0.5)
  init_wide   <- genomeToStanInit(genome, stan_d, scuo = scuo_vals, sphi.init = 2.0)
  expect_lt(sd(init_narrow$latent_phi), sd(init_wide$latent_phi))
})


# ======================================================================
# Section 6: adviToWarmStart() -- init + inv_metric from ADVI fit
# ======================================================================
# Build a synthetic ADVI-like draws data frame to test adviToWarmStart()
# without running actual Stan. We need a data.frame with columns matching
# what mod$variational()$draws(format="df") would return.

local({
  G_s <- stan_d$G   # 8
  K_s <- stan_d$K   # 40
  N   <- 50L        # synthetic draw count

  set.seed(42L)
  # check.names=FALSE preserves "[" in column names (data.frame() would
  # otherwise mangle them to e.g. "dM.1." via make.names())
  synth_draws <- data.frame(
    matrix(rnorm(N * K_s, mean = 0.5, sd = 0.3), nrow = N,
           dimnames = list(NULL, paste0("dM[", 1:K_s, "]"))),
    matrix(rnorm(N * K_s, mean = 0.0, sd = 1.5), nrow = N,
           dimnames = list(NULL, paste0("dEta[", 1:K_s, "]"))),
    matrix(rnorm(N * G_s, mean = 0.0, sd = 0.8), nrow = N,
           dimnames = list(NULL, paste0("latent_phi[", 1:G_s, "]"))),
    sphi       = abs(rnorm(N, mean = 1.3, sd = 0.15)),
    mphi_param = rnorm(N, mean = -0.5, sd = 0.1),
    check.names = FALSE
  )

  # Mock fit object with $draws() method
  mock_fit <- list(draws = function(format = "df") synth_draws)

  ws <<- adviToWarmStart(mock_fit, stan_d)
})

test_that("adviToWarmStart returns init and inv_metric", {
  expect_named(ws, c("init", "inv_metric"))
})

test_that("init has all required parameter fields", {
  expect_named(ws$init, c("dM", "dEta", "latent_phi", "sphi", "mphi_param"),
               ignore.order = TRUE)
})

test_that("init field lengths match data dimensions", {
  expect_length(ws$init$dM,         stan_d$K)
  expect_length(ws$init$dEta,       stan_d$K)
  expect_length(ws$init$latent_phi, stan_d$G)
  expect_length(ws$init$sphi,       1L)
  expect_length(ws$init$mphi_param, 1L)
})

test_that("inv_metric length is 2K + G + 2", {
  expect_length(ws$inv_metric, 2L * stan_d$K + stan_d$G + 2L)
})

test_that("all inv_metric values are strictly positive", {
  expect_true(all(ws$inv_metric > 0))
})

test_that("all init values are finite", {
  expect_true(all(is.finite(ws$init$dM)))
  expect_true(all(is.finite(ws$init$dEta)))
  expect_true(all(is.finite(ws$init$latent_phi)))
  expect_true(is.finite(ws$init$sphi))
  expect_true(is.finite(ws$init$mphi_param))
})

test_that("init$sphi is positive", {
  expect_gt(ws$init$sphi, 0)
})

test_that("inv_metric dEta block (K values) is wider than dM block", {
  K_s    <- stan_d$K
  dM_var   <- ws$inv_metric[seq_len(K_s)]
  dEta_var <- ws$inv_metric[K_s + seq_len(K_s)]
  # dEta draws had sd=1.5, dM draws had sd=0.3 -> dEta vars should be larger
  expect_gt(mean(dEta_var), mean(dM_var))
})

test_that("sphi inv_metric entry is var(log(sphi)), not var(sphi)", {
  # The sphi entry is at position 2K + G + 1
  K_s <- stan_d$K; G_s <- stan_d$G
  sphi_metric_idx <- 2L * K_s + G_s + 1L
  # Should be small (log-scale variance), not on the raw scale
  expect_lt(ws$inv_metric[sphi_metric_idx], 1.0)
})
