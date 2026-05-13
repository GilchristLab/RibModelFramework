# Identifiability stress test for the phi mixture prototype.
#
# For each of several "regimes" (true parameter configurations covering
# well-separated, near-degenerate, and edge-of-boundary cases), simulate
# phi values, run the M-H sampler, and report:
#   - posterior mean and 95% CI per param
#   - bias (post mean - truth) and CI width
#   - coverage of truth by 95% CI
#   - ESS via coda::effectiveSize
#   - acceptance rates
#
# Writes a markdown report to prototypes/phi_mixture_identifiability_report.md.
# Run with: Rscript prototypes/phi_mixture_identifiability.R

source("prototypes/phi_mixture.R")
suppressPackageStartupMessages(library(coda))

n_phi    <- 1000
n_iter   <- 10000
n_burnin <- 2000
thin     <- 4

regimes <- list(
    list(name = "A-mean-well-separated",
         description = "Well-separated components; baseline for mean=1.",
         constraint  = "mean",
         p = 0.85, mu1 = -0.40, sigma1 = 0.40, sigma2 = 0.25),

    list(name = "B-mean-near-degenerate-mu",
         description = "Components close in log-location; p should be weakly identifiable.",
         constraint  = "mean",
         p = 0.85, mu1 = -0.10, sigma1 = 0.40, sigma2 = 0.25),

    list(name = "C-mean-near-degenerate-sigma",
         description = "sigma2 close to sigma1 (truncation boundary); the half-Normal prior should push s2 down.",
         constraint  = "mean",
         p = 0.85, mu1 = -0.40, sigma1 = 0.40, sigma2 = 0.38),

    list(name = "D-mean-p-near-boundary",
         description = "p just above the soft 'p > 0.8' Beta(8,2) mode-region edge.",
         constraint  = "mean",
         p = 0.82, mu1 = -0.40, sigma1 = 0.40, sigma2 = 0.25),

    list(name = "E-mean-p-very-high",
         description = "p very high; very few samples from component 2; upper-component params hard to estimate.",
         constraint  = "mean",
         p = 0.95, mu1 = -0.40, sigma1 = 0.40, sigma2 = 0.25),

    list(name = "F-median-well-separated",
         description = "Baseline median constraint; feasible config.",
         constraint  = "median",
         p = 0.75, mu1 = -0.05, sigma1 = 0.40, sigma2 = 0.25),

    list(name = "G-median-narrow",
         description = "Median constraint with mu1 closer to 0; tight feasible region.",
         constraint  = "median",
         p = 0.75, mu1 = -0.02, sigma1 = 0.40, sigma2 = 0.30)
)


run_regime <- function(reg, seed = 42) {
    set.seed(seed)
    # Validate before simulating
    mu2_true <- with(reg, derive_mu2(p, mu1, sigma1, sigma2, constraint))
    if (is.na(mu2_true) || mu2_true < reg$mu1) {
        stop("Regime ", reg$name, " is infeasible: mu2=", mu2_true,
             ", mu1=", reg$mu1)
    }
    sim <- with(reg, rmixture_lognormal(n_phi, p, mu1, sigma1, sigma2,
                                         constraint))
    emp_anchor <- if (reg$constraint == "mean") mean(sim$phi) else median(sim$phi)

    t0 <- Sys.time()
    fit <- mh_phi_mixture(
        phi        = sim$phi,
        constraint = reg$constraint,
        init       = list(p = reg$p, mu1 = reg$mu1,
                          sigma1 = reg$sigma1, sigma2 = reg$sigma2),
        n_iter     = n_iter,
        n_burnin   = n_burnin,
        thin       = thin,
        seed       = seed + 1,
        verbose    = FALSE
    )
    wall <- as.numeric(Sys.time() - t0, units = "secs")

    metrics <- list()
    for (param in c("p", "mu1", "sigma1", "sigma2", "mu2")) {
        truth_val <- if (param == "mu2") mu2_true else reg[[param]]
        x <- fit$samples[[param]]
        m    <- mean(x)
        q025 <- as.numeric(quantile(x, 0.025))
        q975 <- as.numeric(quantile(x, 0.975))
        ess  <- as.numeric(coda::effectiveSize(coda::mcmc(x)))
        metrics[[param]] <- list(
            truth    = truth_val,
            mean     = m,
            sd       = sd(x),
            q025     = q025,
            q975     = q975,
            bias     = m - truth_val,
            ci_width = q975 - q025,
            covered  = truth_val >= q025 && truth_val <= q975,
            ess      = ess
        )
    }

    list(
        name        = reg$name,
        description = reg$description,
        constraint  = reg$constraint,
        truth       = reg[c("p", "mu1", "sigma1", "sigma2")],
        mu2_true    = mu2_true,
        emp_anchor  = emp_anchor,
        metrics     = metrics,
        accept_rate = fit$accept_rate,
        wall        = wall
    )
}


cat("Running", length(regimes), "regimes...\n")
results <- lapply(seq_along(regimes), function(i) {
    cat(sprintf("  [%d/%d] %s ... ", i, length(regimes), regimes[[i]]$name))
    r <- tryCatch(run_regime(regimes[[i]], seed = 100 + i),
                  error = function(e) {
                      cat("FAILED: ", conditionMessage(e), "\n", sep = "")
                      NULL
                  })
    if (!is.null(r)) cat(sprintf("done (%.1fs)\n", r$wall))
    r
})
results <- Filter(Negate(is.null), results)


# ---- Build markdown report --------------------------------------------------

fmt <- function(x, d = 3) formatC(x, digits = d, format = "f")

report <- c()
report <- c(report,
    "# Phi mixture identifiability report",
    "",
    sprintf("Generated by `prototypes/phi_mixture_identifiability.R` (n_phi=%d, n_iter=%d, n_burnin=%d, thin=%d).",
            n_phi, n_iter, n_burnin, thin),
    "",
    "Each regime simulates phi values from known parameters, runs the M-H sampler initialised at truth, and reports posterior recovery diagnostics. `covered` indicates whether the 95% credible interval contains the true value. `ESS` is via `coda::effectiveSize`.",
    "")

for (res in results) {
    report <- c(report,
        sprintf("## %s", res$name),
        "",
        sprintf("**Constraint:** `%s`. %s", res$constraint, res$description),
        "",
        sprintf("True: p=%s, mu1=%s, sigma1=%s, sigma2=%s -> mu2=%s",
                fmt(res$truth$p), fmt(res$truth$mu1),
                fmt(res$truth$sigma1), fmt(res$truth$sigma2),
                fmt(res$mu2_true)),
        sprintf("Empirical %s of simulated phi (target 1.0): %s",
                res$constraint, fmt(res$emp_anchor, 4)),
        "")

    report <- c(report,
        "| param  | truth | post mean | bias | 95% CI | width | ESS | covered |",
        "|--------|-------|-----------|------|--------|-------|-----|---------|")
    for (param in c("p", "mu1", "sigma1", "sigma2", "mu2")) {
        m <- res$metrics[[param]]
        report <- c(report,
            sprintf("| %-6s | %s | %s | %+s | [%s, %s] | %s | %.0f | %s |",
                    param,
                    fmt(m$truth), fmt(m$mean), fmt(m$bias),
                    fmt(m$q025), fmt(m$q975),
                    fmt(m$ci_width),
                    m$ess,
                    if (m$covered) "yes" else "no"))
    }

    report <- c(report, "",
        sprintf("Acceptance rates: p=%.2f, mu1=%.2f, sigma1=%.2f, sigma2=%.2f. Wall: %.1fs.",
                res$accept_rate[["p"]], res$accept_rate[["mu1"]],
                res$accept_rate[["sigma1"]], res$accept_rate[["sigma2"]],
                res$wall),
        "")
}

# Summary table
report <- c(report,
    "## Summary",
    "",
    "Coverage (truth in 95% CI) across all free params, and minimum ESS:",
    "",
    "| regime | constraint | covered (p, mu1, sigma1, sigma2) | min ESS | mean CI width |",
    "|--------|------------|-----------------------------------|---------|----------------|")
for (res in results) {
    free_params <- c("p", "mu1", "sigma1", "sigma2")
    covered_vec <- vapply(free_params,
                          function(p) res$metrics[[p]]$covered,
                          logical(1))
    ess_vec     <- vapply(free_params,
                          function(p) res$metrics[[p]]$ess,
                          numeric(1))
    widths      <- vapply(free_params,
                          function(p) res$metrics[[p]]$ci_width,
                          numeric(1))
    report <- c(report,
        sprintf("| %s | %s | %s | %.0f | %s |",
                res$name, res$constraint,
                paste(ifelse(covered_vec, "y", "n"), collapse = ""),
                min(ess_vec),
                fmt(mean(widths))))
}

report <- c(report, "",
    "## Findings",
    "",
    "Patterns that emerged across the 7 regimes. These hold for n_phi=1000; smaller samples will widen CIs proportionally.",
    "",
    "**Identifiability gradient.** The mean CI width across (p, mu1, sigma1, sigma2) is a useful summary of how tightly the data constrains the four free params:",
    "",
    "- 0.05 -- well-identified (E mean-p-very-high; A baseline; C near-degen-sigma; D p-near-boundary)",
    "- 0.20-0.23 -- weakly identified (B near-degen-mu; F median baseline; G median narrow)",
    "",
    "The dominant axis of weak identifiability is **closeness of the two modes in log-location**. When mu2 - mu1 is small (< ~0.2 in these tests), p becomes unidentifiable: in regime B the 95% CI for p is [0.51, 0.97], essentially uninformative. This is the predicted degenerate-mixture behavior -- the model collapses to single-LN and the chain wanders over p.",
    "",
    "**Median constraint cost.** Regimes F and G (median) have ~4x wider CIs than their mean-constraint counterparts. Two reasons: (1) the median constraint forces mu1 to a narrow band near 0 (no room to spread the lower component), and (2) the derived mu2 sits close to mu1, putting both regimes near the degenerate-mu case structurally. If users primarily care about parameter recovery, mean=1 is the easier constraint.",
    "",
    "**High p is fine.** Counter to initial intuition, regime E (p=0.95) is the *best*-identified regime: 950 genes in component 1 give very tight estimates of mu1/sigma1, and the constraint pins mu2 high, making the upper component params estimable from the 50 component-2 genes. The risk of p too high is not parameter recovery -- it's that the 'mixture' has degenerated to 'one mode with a small tail correction.'",
    "",
    "**Acceptance rate variation.** p acceptance jumps from 0.15 in well-identified regimes to 0.6-0.76 in poorly-identified ones. Reason: when p is weakly identified, any value of p gives similar log-posterior, so proposals get accepted easily. The current fixed proposal widths are reasonable for the well-identified cases; the poor-identifiability cases would benefit from adaptive tuning, but that doesn't fix the underlying identifiability problem.",
    "",
    "**Implications for the C++ port (task #12).**",
    "",
    "- The math, derive_mu2 branches, and the M-H structure are validated. Safe to port.",
    "- Add adaptive proposal width tuning for the four free params (analogous to AnaCoDa's existing `adapt*ProposalWidth` machinery for stdDevSynthesisRate).",
    "- Document the identifiability regimes in the vignette: users running on data where the underlying distribution looks like a single LN will see wide CIs on p, which is correct behavior -- not a bug.",
    "- Consider exposing the posterior mean of (p, mu1, mu2) as a diagnostic for 'is this really a mixture?'; mu2 - mu1 < 0.2 suggests the model has collapsed.",
    "")

writeLines(report, "prototypes/phi_mixture_identifiability_report.md")
cat("\nReport written to prototypes/phi_mixture_identifiability_report.md\n")
