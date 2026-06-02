# panseStanFinalize.R -- shared post-processing for PANSE Stan fits.
#
# Called by fit_panse_stan() after sampling and by recover_postproc.R when
# the wrapper died after sampling completed.
#
# Depends on posterior (via fit$summary()) -- pre-existing cmdstanr dependency.

#' Diagnose, save, and summarize a PANSE Stan fit
#'
#' Saves the fit RDS, computes and saves the posterior summary RDS, and prints
#' convergence diagnostics.  Designed to be called immediately after
#' \code{CmdStanModel$sample()} completes.
#'
#' @param fit A \code{CmdStanMCMC} object from \code{cmdstanr}.
#' @param config Parsed YAML config list.
#' @param stan_data Stan data list passed to \code{$sample()}.
#' @param out_dir Output directory for RDS files.
#' @param wall_sec Sampling wall time in seconds (\code{NA_real_} when
#'   recovering from existing CSV output).
#' @param git_sha Short git SHA string (default \code{"unknown"}).
#' @param diagnose_fn Function called for diagnostics; receives \code{fit}.
#'   Default calls \code{fit$diagnostic_summary()}.
#' @return Invisibly, the summary data frame.
#' @export
panse_stan_finalize <- function(fit, config, stan_data, out_dir,
                                wall_sec = NA_real_, git_sha = "unknown",
                                diagnose_fn = function(fit) print(fit$diagnostic_summary())) {

    cat("\n--- Diagnostics ---\n")
    diagnose_fn(fit)

    cat("\n--- Saving fit ---\n")
    saveRDS(list(
        fit       = fit,
        config    = config,
        stan_data = stan_data,
        wall_sec  = wall_sec,
        git_sha   = git_sha
    ), file.path(out_dir, "panse-stan-fit.rds"))
    cat("[ok] saved panse-stan-fit.rds\n")

    # Pre-compute summary (CSP parameters + phi when sampled).
    # Use stan_variables (parameter stems) not variables (indexed expansion).
    stan_vars    <- fit$metadata()$stan_variables
    summary_vars <- c("alpha", "lambdaPrime")
    summary_vars <- c(summary_vars,
                      if ("NSERate_shared" %in% stan_vars) "NSERate_shared" else "NSERate")
    if ("log_phi" %in% stan_vars) summary_vars <- c(summary_vars, "log_phi")
    if ("sphi"    %in% stan_vars) summary_vars <- c(summary_vars, "sphi")
    cat("[summary] variables:", paste(summary_vars, collapse = ", "), "\n")

    sm <- fit$summary(variables = summary_vars)
    saveRDS(sm, file.path(out_dir, "stan-summary.rds"))
    cat("[ok] saved stan-summary.rds\n")

    cat("\nR-hat   range:", range(sm$rhat,     na.rm = TRUE), "\n")
    cat("ESS bulk range:", range(sm$ess_bulk, na.rm = TRUE), "\n")
    cat("ESS tail range:", range(sm$ess_tail, na.rm = TRUE), "\n")
    cat("# params w/ R-hat > 1.01:", sum(sm$rhat > 1.01, na.rm = TRUE), "/",
        nrow(sm), "\n")

    cat("\n[ok] Output ->", normalizePath(out_dir), "\n")
    invisible(sm)
}
