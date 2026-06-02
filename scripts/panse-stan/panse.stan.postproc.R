# panse.stan.postproc.R -- shared post-processing for PANSE Stan fits
#
# Provides panse_stan_finalize(), called by fit.stan.R after sampling and by
# recover_postproc.R when the wrapper died after sampling completed.
#
# Depends on posterior (via fit$summary()) -- pre-existing cmdstanr dependency.

# Diagnose, save fit RDS, compute and save summary RDS, print convergence stats.
#
# fit         -- CmdStanMCMC object (from mod$sample() or as_cmdstan_fit())
# config      -- parsed YAML config list
# stan_data   -- list passed to mod$sample()
# out_dir     -- output directory for RDS files
# wall_sec    -- sampling wall time in seconds (NA_real_ when recovering from CSV)
# git_sha     -- short git SHA string
# diagnose_fn -- function(fit) called for diagnostics; caller supplies the
#                appropriate method:
#                  fit.stan.R:        function(fit) fit$cmdstan_diagnose()
#                  recover_postproc:  function(fit) print(fit$diagnostic_summary())
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

    # Pre-compute summary (CSP plus phi when sampled).
    # Use stan_variables (parameter STEMS) not variables (indexed expansion).
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
