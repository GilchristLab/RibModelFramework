# Validation: simulate phi from a known mixture, run the M-H sampler, check
# parameter recovery. This is a single regime; multi-regime stress testing
# lives in task #14 (phi_mixture_identifiability.R).
#
# Run with: Rscript prototypes/phi_mixture_sampler_check.R

source("prototypes/phi_mixture.R")

ok <- function(msg, cond) {
    cat(if (isTRUE(cond)) "  ok    " else "  FAIL  ", msg, "\n", sep = "")
    if (!isTRUE(cond)) stop("check failed: ", msg)
}


# ---------------- Test 1: mean constraint, well-separated components --------

cat("=== test 1: mean=1, well-separated (p=0.85, mu1=-0.4, s1=0.4, s2=0.25) ===\n")

set.seed(101)
truth <- list(p = 0.85, mu1 = -0.4, sigma1 = 0.4, sigma2 = 0.25)
mu2_true <- with(truth, derive_mu2(p, mu1, sigma1, sigma2, "mean"))
cat(sprintf("  true: p=%.2f mu1=%.2f s1=%.2f s2=%.2f -> mu2=%.3f\n",
            truth$p, truth$mu1, truth$sigma1, truth$sigma2, mu2_true))

n_phi <- 1000
sim <- with(truth, rmixture_lognormal(n_phi, p, mu1, sigma1, sigma2, "mean"))
cat(sprintf("  simulated %d phi values; empirical mean=%.4f (target 1)\n",
            n_phi, mean(sim$phi)))

t0 <- Sys.time()
fit <- mh_phi_mixture(
    phi = sim$phi,
    constraint = "mean",
    init = list(p = 0.8, mu1 = -0.3, sigma1 = 0.5, sigma2 = 0.3),
    n_iter = 8000, n_burnin = 2000, thin = 4,
    seed = 102,
    verbose = TRUE
)
cat(sprintf("  sampler wall time: %.1fs\n", as.numeric(Sys.time() - t0)))

.q <- function(samples, prob) {
    sapply(samples, function(x) unname(quantile(x, prob)))
}
post_mean <- sapply(fit$samples[, c("p", "mu1", "sigma1", "sigma2", "mu2")],
                    mean)
post_q025 <- .q(fit$samples[, c("p", "mu1", "sigma1", "sigma2", "mu2")], 0.025)
post_q975 <- .q(fit$samples[, c("p", "mu1", "sigma1", "sigma2", "mu2")], 0.975)

cat("\n  posterior summary:\n")
cat(sprintf("  %-6s  %8s   %8s  %8s   %8s\n",
            "param", "true", "post mean", "2.5%", "97.5%"))
for (param in c("p", "mu1", "sigma1", "sigma2")) {
    cat(sprintf("  %-6s  %8.4f   %8.4f  %8.4f   %8.4f\n",
                param, truth[[param]], post_mean[[param]],
                post_q025[[param]], post_q975[[param]]))
}
cat(sprintf("  %-6s  %8.4f   %8.4f  %8.4f   %8.4f\n",
            "mu2 (derived)", mu2_true, post_mean[["mu2"]],
            post_q025[["mu2"]], post_q975[["mu2"]]))

# Coverage check: do 95% CIs contain the true values?
for (param in c("p", "mu1", "sigma1", "sigma2")) {
    covered <- truth[[param]] >= post_q025[[param]] &&
               truth[[param]] <= post_q975[[param]]
    ok(sprintf("95%% CI covers true %s = %.3f", param, truth[[param]]), covered)
}

# Acceptance rates in a reasonable range (0.15 - 0.6 per param is fine for
# random walks)
for (param in c("p", "mu1", "sigma1", "sigma2")) {
    rate <- fit$accept_rate[[param]]
    ok(sprintf("acceptance for %s in (0.10, 0.70): %.3f", param, rate),
       rate > 0.10 && rate < 0.70)
}


# ---------------- Test 2: median constraint, feasible config ----------------

cat("\n=== test 2: median=1 (p=0.75, mu1=-0.05, s1=0.4, s2=0.25) ===\n")

set.seed(201)
truth2 <- list(p = 0.75, mu1 = -0.05, sigma1 = 0.4, sigma2 = 0.25)
mu2_true2 <- with(truth2, derive_mu2(p, mu1, sigma1, sigma2, "median"))
cat(sprintf("  true: p=%.2f mu1=%.2f s1=%.2f s2=%.2f -> mu2=%.3f\n",
            truth2$p, truth2$mu1, truth2$sigma1, truth2$sigma2, mu2_true2))

sim2 <- with(truth2, rmixture_lognormal(n_phi, p, mu1, sigma1, sigma2, "median"))
cat(sprintf("  simulated %d phi values; empirical median=%.4f (target 1)\n",
            n_phi, median(sim2$phi)))

fit2 <- mh_phi_mixture(
    phi = sim2$phi,
    constraint = "median",
    init = list(p = 0.75, mu1 = -0.05, sigma1 = 0.4, sigma2 = 0.25),
    n_iter = 8000, n_burnin = 2000, thin = 4,
    seed = 202,
    verbose = TRUE
)

post_mean2 <- sapply(fit2$samples[, c("p", "mu1", "sigma1", "sigma2")], mean)
post_q025_2 <- .q(fit2$samples[, c("p", "mu1", "sigma1", "sigma2")], 0.025)
post_q975_2 <- .q(fit2$samples[, c("p", "mu1", "sigma1", "sigma2")], 0.975)

cat("\n  posterior summary:\n")
cat(sprintf("  %-6s  %8s   %8s  %8s   %8s\n",
            "param", "true", "post mean", "2.5%", "97.5%"))
for (param in c("p", "mu1", "sigma1", "sigma2")) {
    cat(sprintf("  %-6s  %8.4f   %8.4f  %8.4f   %8.4f\n",
                param, truth2[[param]], post_mean2[[param]],
                post_q025_2[[param]], post_q975_2[[param]]))
}

for (param in c("p", "mu1", "sigma1", "sigma2")) {
    covered <- truth2[[param]] >= post_q025_2[[param]] &&
               truth2[[param]] <= post_q975_2[[param]]
    ok(sprintf("95%% CI covers true %s = %.3f", param, truth2[[param]]), covered)
}


cat("\nAll sampler recovery checks passed.\n")
