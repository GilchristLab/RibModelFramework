# ============================================================================
# helper-panse-stan.R -- shared helpers for the PANSE *Stan port* test suite
# (test-panse-stan-*.R).  DISTINCT from the native PANSE MCMC tests; do not mix.
#
# Provides: skip guard, stan/ locator, compile+cache, fixture loader, and an
# independent R reference for the per-codon survival term (copied verbatim from
# the analysis-side ground-truth test, so the Stan model is checked against an
# implementation it does not share).
#
# RUNNING THESE TESTS MANUALLY
# ----------------------------
# The skip guard calls testthat::skip_on_cran(), which causes all Stan tests
# to be skipped when the suite is run via test_check("AnaCoDa") (the standard
# CI harness).  To run them locally set NOT_CRAN=true:
#
#   NOT_CRAN=true Rscript -e '
#     library(AnaCoDa)
#     library(testthat)
#     Sys.setenv(NOT_CRAN = "true")
#     test_dir("tests/testthat", filter = "panse-stan")
#   '
#
# CmdStan and the R packages cmdstanr + posterior must be installed.
# Stan models are compiled on first run and cached by cmdstanr.
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

# ---- locate the stan/ directory (inst/stan preferred, top-level fallback) ---
panse_stan_dir <- function() {
    # 1. Installed package path (inst/stan/ after R CMD INSTALL).
    pkg_dir <- system.file("stan", package = "AnaCoDa")
    if (nzchar(pkg_dir) && file.exists(file.path(pkg_dir, "panse_csp_only.stan")))
        return(normalizePath(pkg_dir))
    # 2. In-source fallback: walk up to the repo-root stan/ (worktree / devtools).
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

# Stan data list tailored to each model family's data block:
#   csp_only*  -> phi is FIXED input data
#   basic*     -> phi is sampled; sphi is FIXED input data
#   sphi_est*  -> phi & sphi are parameters; needs sphi_prior_sd
panse_stan_data_for <- function(fx, model, sphi_prior_sd = 2.5) {
    sd <- fx$stan_data
    sd$phi <- NULL
    if (grepl("csp_only", model)) {
        sd$phi <- as.numeric(fx$truth$phi)
    } else if (grepl("basic", model)) {
        sd$sphi <- fx$truth$sphi
    } else {                                  # sphi_est_*
        sd$sphi_prior_sd <- sphi_prior_sd
    }
    sd
}

# A separate model compiled WITH standalone functions, for expose_functions().
# (expose_functions fails on a sampling exe loaded from cache.)
.panse_stan_fn_models <- new.env(parent = emptyenv())
panse_stan_functions <- function(name = "panse_sphi_est_noncentered_sharednse") {
    if (!is.null(.panse_stan_fn_models[[name]])) return(.panse_stan_fn_models[[name]])
    # Compile from a temp copy so there is no cached sampling exe to load as
    # "pre-compiled" (which blocks expose_functions); forces a fresh build that
    # retains the C++ needed to export the functions block.
    src <- file.path(panse_stan_dir(), paste0(name, ".stan"))
    tmp <- file.path(tempdir(), paste0(name, "__fns.stan"))
    file.copy(src, tmp, overwrite = TRUE)
    m <- cmdstanr::cmdstan_model(tmp, compile_standalone = TRUE, quiet = TRUE)
    m$expose_functions(verbose = FALSE)
    .panse_stan_fn_models[[name]] <- m
    m
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
