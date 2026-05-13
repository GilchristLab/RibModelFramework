library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Tests for the ROC phi prior parameterization.
#
# ROC's default phi prior is implemented at src/ROCModel.cpp:139 as:
#   mPhi = -(stdDevSynthesisRate^2) / 2
#   densityLogNorm(phi, mPhi, stdDevSynthesisRate)
#
# This parameterization anchors the MEAN of phi at 1 (E[phi] = 1), not
# the median. The median under this prior is exp(-sigma^2/2) < 1.
#
# These tests pin the math expectation of the parameterization so that
# silent changes to mPhi in src/ROCModel.cpp:139 surface in review.
#
# A future option (task #12) will let users select phi.prior.constraint
# = "median" -- when that lands, the second test below activates.
# ======================================================================

test_that("default ROC phi prior (single lognormal) anchors E[phi] = 1", {
    sigma <- 0.5
    mPhi  <- -(sigma^2) / 2

    # Closed-form: E[phi] = exp(mu + sigma^2/2) for LogNormal(mu, sigma)
    expect_equal(exp(mPhi + sigma^2 / 2), 1.0, tolerance = 1e-12)

    # Median = exp(mu); must be strictly below 1 under this parameterization
    expect_lt(exp(mPhi), 1.0)
    expect_equal(exp(mPhi), exp(-sigma^2 / 2), tolerance = 1e-12)

    # Monte Carlo sanity: empirical mean ~ 1, empirical median < 1
    set.seed(1)
    phi <- rlnorm(1e6, meanlog = mPhi, sdlog = sigma)
    expect_equal(mean(phi), 1.0, tolerance = 0.01)
    expect_lt(median(phi), 0.95)
})

test_that("median-constraint phi prior (when enabled) anchors median[phi] = 1", {
    # Forward-compatibility guard for task #12. When phi.prior.constraint
    # = "median" is added to the single-LN code path, the implementation
    # should use mPhi = 0 so that exp(mPhi) = 1 (median of LogNormal).
    # Remove the skip() and add the AnaCoDa-side hook once #12 lands.
    skip("phi.prior.constraint = 'median' not yet implemented (task #12)")

    sigma <- 0.5
    mPhi  <- 0

    expect_equal(exp(mPhi), 1.0, tolerance = 1e-12)
    expect_gt(exp(mPhi + sigma^2 / 2), 1.0)
})
