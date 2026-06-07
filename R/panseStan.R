# panseStan.R -- resolver for installed PANSE Stan model files.
#
# The six panse_*.stan models live in inst/stan/ and are installed into the
# package as read-only files.  cmdstanr compiles them to a WRITABLE cache dir
# (never next to the read-only source).
#
# Path resolution order (internal .panse_stan_file):
#   1. system.file("stan", basename, package="AnaCoDa")  [installed package]
#   2. file.path(dirname(dirname(..)),  "stan", basename) [in-source/worktree]
# Stops with an informative error if neither resolves.


# Valid PANSE Stan model basenames.
.PANSE_STAN_MODELS <- c(
    "panse_csp_only.stan",
    "panse_csp_only_sharednse.stan",
    "panse_basic.stan",
    "panse_basic_sharednse.stan",
    "panse_sphi_est_noncentered.stan",
    "panse_sphi_est_noncentered_pointval.stan",
    "panse_sphi_est_noncentered_singular.stan",   # dwell.dist=singular (Poisson + point-value)
    "panse_sphi_est_noncentered_invgamma.stan",   # dwell.dist=invgamma (Sichel surrogate + exact surv)
    "panse_sphi_est_sumzero.stan",
    "panse_sphi_est_sumzero_geomean.stan",
    "panse_sphi_est_centered.stan",
    "panse_sphi_est_centered_sharednse.stan",
    "panse_sphi_est_noncentered_sharednse.stan",
    "panse_sphi_est_noncentered_aanse.stan"
)


# Internal resolver (used by fit_panse_stan and tests).
.panse_stan_file <- function(basename) {
    if (!basename %in% .PANSE_STAN_MODELS)
        stop(".panse_stan_file: unknown model '", basename, "'; valid: ",
             paste(.PANSE_STAN_MODELS, collapse = ", "))

    # 1. Installed package path (post R CMD INSTALL).
    pkg_path <- system.file("stan", basename, package = "AnaCoDa")
    if (nzchar(pkg_path) && file.exists(pkg_path)) return(pkg_path)

    # 2. In-source fallback (running from a worktree or devtools::load_all()).
    src_path <- system.file("../stan", basename, package = "AnaCoDa")
    if (nzchar(src_path) && file.exists(src_path)) return(src_path)

    stop(".panse_stan_file: could not resolve '", basename, "' from the ",
         "installed package or in-source stan/. Run R CMD INSTALL or ",
         "devtools::load_all() from the package root.")
}


#' Return the path to an installed PANSE Stan model file
#'
#' Resolves the absolute path to one of the six \code{panse_*.stan} models
#' shipped with \pkg{AnaCoDa}.  Useful for A-RMF analysis tools that need to
#' inspect or manually compile a model.  For fitting, prefer
#' \code{\link{fit_panse_stan}} which handles compilation automatically.
#'
#' @param basename Character; one of the six valid PANSE Stan model basenames:
#'   \code{"panse_csp_only.stan"}, \code{"panse_csp_only_sharednse.stan"},
#'   \code{"panse_basic.stan"}, \code{"panse_basic_sharednse.stan"},
#'   \code{"panse_sphi_est_centered_sharednse.stan"},
#'   \code{"panse_sphi_est_noncentered_sharednse.stan"}.
#' @return Absolute path string.
#' @export
panse_stan_model_path <- function(basename) {
    .panse_stan_file(basename)
}


#' Compile a PANSE Stan model with cmdstanr
#'
#' Looks up the installed \code{panse_*.stan} source via
#' \code{\link{panse_stan_model_path}}, then compiles it to a writable cache
#' directory so the read-only \code{inst/stan/} source is never modified.
#'
#' @param basename PANSE Stan model basename (see \code{\link{panse_stan_model_path}}).
#' @param build_dir Writable cache directory for the compiled binary.  Defaults
#'   to \code{tools::R_user_dir("AnaCoDa", "cache")}.
#' @param threads Logical; compile with \code{stan_threads = TRUE} (default).
#' @param ... Additional arguments forwarded to \code{cmdstanr::cmdstan_model()}.
#' @return A \code{CmdStanModel} object.
#' @export
panse_compile_model <- function(basename,
                                build_dir = tools::R_user_dir("AnaCoDa", "cache"),
                                threads   = TRUE,
                                ...) {
    if (!requireNamespace("cmdstanr", quietly = TRUE))
        stop("panse_compile_model requires cmdstanr. ",
             "Install with: install.packages('cmdstanr', repos='https://mc-stan.org/r-packages/')")

    stan_file <- .panse_stan_file(basename)
    if (!dir.exists(build_dir)) dir.create(build_dir, recursive = TRUE)

    cpp_opts <- if (threads) list(stan_threads = TRUE) else list()
    cmdstanr::cmdstan_model(stan_file, dir = build_dir, cpp_options = cpp_opts, ...)
}
