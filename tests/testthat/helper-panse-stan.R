# ============================================================================
# helper-panse-stan.R -- shared helpers for the PANSE *Stan port* test suite
# (test-panse-stan-*.R).  DISTINCT from the native PANSE MCMC tests; do not mix.
#
# Provides: skip guard, stan/ locator, compile+cache, fixture loader, and an
# independent R reference for the per-codon survival term (copied verbatim from
# the analysis-side ground-truth test, so the Stan model is checked against an
# implementation it does not share).
# ============================================================================

# ---- skip guard: Stan tests need cmdstanr + a working cmdstan toolchain -----
skip_if_no_cmdstan <- function() {
    testthat::skip_on_cran()
    testthat::skip_if_not_installed("cmdstanr")
    testthat::skip_if_not_installed("posterior")
    v <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                  error = function(e) NA_character_)
    if (is.na(v) || is.null(v)) testthat::skip("cmdstan toolchain not installed")
}

# ---- locate the stan/ directory by walking up from the test dir -------------
panse_stan_dir <- function() {
    d <- normalizePath(testthat::test_path("."), mustWork = FALSE)
    for (i in 1:6) {
        cand <- file.path(d, "stan")
        if (file.exists(file.path(cand, "panse_csp_only.stan")))
            return(normalizePath(cand))
        d <- dirname(d)
    }
    testthat::skip("cannot locate stan/ directory (panse_csp_only.stan not found)")
}

# ---- compile + cache models once per test-run -------------------------------
.panse_stan_models <- new.env(parent = emptyenv())
panse_stan_model <- function(name) {
    if (!is.null(.panse_stan_models[[name]])) return(.panse_stan_models[[name]])
    stan_file <- file.path(panse_stan_dir(), paste0(name, ".stan"))
    m <- cmdstanr::cmdstan_model(stan_file,
                                 cpp_options = list(stan_threads = TRUE),
                                 quiet = TRUE)
    .panse_stan_models[[name]] <- m
    m
}

PANSE_STAN_MODELS <- c(
    "panse_csp_only", "panse_csp_only_sharednse",
    "panse_basic", "panse_basic_sharednse",
    "panse_sphi_est_centered_sharednse",
    "panse_sphi_est_noncentered_sharednse"
)

# ---- fixture ----------------------------------------------------------------
load_panse_stan_fixture <- function() {
    readRDS(testthat::test_path("fixtures", "panse-stan", "likelihood_fixture.rds"))
}

# Stan data list ready for the sphi-est models: add sphi_prior_sd, drop the
# extra (undeclared) phi field that the data-builder attaches.
panse_stan_data_for_fit <- function(fx, sphi_prior_sd = 2.5) {
    sd <- fx$stan_data
    sd$sphi_prior_sd <- sphi_prior_sd
    sd$phi <- NULL
    sd
}

# Build a 1-draw draws_array of the truth parameters for the noncentered
# shared-NSE model (for generate_quantities).
panse_stan_truth_draws <- function(fx) {
    tr <- fx$truth; G <- fx$meta$G; C <- fx$meta$C
    sphi  <- tr$sphi
    z_phi <- (log(tr$phi) + 0.5 * sphi^2) / sphi
    nm <- c(paste0("log_alpha[", 1:C, "]"),
            paste0("log_lambdaPrime[", 1:C, "]"),
            "log_NSERate_shared", "sphi",
            paste0("z_phi[", 1:G, "]"))
    vals <- c(log(tr$alpha), log(tr$lambda), log(tr$nse_shared), sphi, z_phi)
    posterior::as_draws_array(matrix(vals, nrow = 1, dimnames = list(NULL, nm)))
}

# ---- independent R reference for the survival term --------------------------
# (verbatim from panse/.../scripts/sim/test_log_psuccess.R)
.panse_uig_cf <- function(s, x, depth = 10000L) {
    rv <- depth / x
    for (i in depth:1) {
        if (i %% 2 == 0) rv <- (i %/% 2) / (x + rv)
        else             rv <- ((i %/% 2) + 1 - s) / (1 + rv)
    }
    x + rv
}
.panse_log_uig  <- function(s, x) s * log(x) - x - log(.panse_uig_cf(s, x))
.panse_closed   <- function(a, l, n) { v <- 1/n; lv <- l*v; a*log(lv) + lv + .panse_log_uig(1 - a, lv) }
.panse_taylor2  <- function(a, l, n) { v <- 1/n; alv <- a/(l*v); -alv + alv/(l*v) + 0.5*alv*alv }
ref_log_psuccess <- function(a, l, n) {
    if (a * n / l < 0.005) .panse_taylor2(a, l, n) else .panse_closed(a, l, n)
}
