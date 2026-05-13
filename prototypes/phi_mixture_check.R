# Sanity checks for prototypes/phi_mixture.R
# Run with: Rscript prototypes/phi_mixture_check.R

source("prototypes/phi_mixture.R")

ok <- function(msg, cond) {
    cat(if (isTRUE(cond)) "  ok    " else "  FAIL  ", msg, "\n", sep = "")
    if (!isTRUE(cond)) stop("check failed: ", msg)
}
near <- function(a, b, tol = 1e-10) isTRUE(all.equal(a, b, tolerance = tol))


cat("=== derive_mu2 closed-form correctness ===\n")

# mean constraint: verify p*exp(mu1+s1^2/2) + (1-p)*exp(mu2+s2^2/2) == 1
for (cfg in list(
    list(p = 0.9, mu1 = -0.5, sigma1 = 0.4, sigma2 = 0.3),
    list(p = 0.5, mu1 = -1.0, sigma1 = 0.6, sigma2 = 0.2),
    list(p = 0.85, mu1 = -0.2, sigma1 = 0.5, sigma2 = 0.3)
)) {
    mu2 <- with(cfg, derive_mu2(p, mu1, sigma1, sigma2, "mean"))
    em  <- with(cfg, p * exp(mu1 + sigma1^2 / 2) +
                     (1 - p) * exp(mu2 + sigma2^2 / 2))
    ok(sprintf("mean constraint holds: p=%.2f mu1=%.1f s1=%.1f s2=%.1f -> E[phi]=%.15f",
               cfg$p, cfg$mu1, cfg$sigma1, cfg$sigma2, em),
       near(em, 1))
}

# median constraint: verify p*pnorm(-mu1/s1) + (1-p)*pnorm(-mu2/s2) == 0.5
# AND the ordering constraint mu2 >= mu1 holds.
#
# Feasibility under median=1 + ordering is narrower than under mean=1 alone:
# mu2 >= mu1 forces mu1 <= some upper bound that depends on (p, sigma1, sigma2).
# At p = 0.9 with these sigmas the bound is essentially mu1 <= 0; large positive
# mu1 produces a "feasible" mu2 that's below mu1 (label-switched).
for (cfg in list(
    list(p = 0.85, mu1 = -0.01, sigma1 = 0.4, sigma2 = 0.3),
    list(p = 0.5,  mu1 = -0.3,  sigma1 = 0.5, sigma2 = 0.2),
    list(p = 0.75, mu1 = -0.05, sigma1 = 0.5, sigma2 = 0.3)
)) {
    mu2 <- with(cfg, derive_mu2(p, mu1, sigma1, sigma2, "median"))
    ok(sprintf("feasibility for p=%.2f mu1=%.2f -> mu2=%s (ordering holds)",
               cfg$p, cfg$mu1, format(mu2, digits = 4)),
       !is.na(mu2) && mu2 >= cfg$mu1)
    lhs <- with(cfg, p * pnorm(-mu1 / sigma1) +
                     (1 - p) * pnorm(-mu2 / sigma2))
    ok(sprintf("median constraint holds: p=%.2f mu1=%.2f -> CDF(1)=%.15f",
               cfg$p, cfg$mu1, lhs),
       near(lhs, 0.5))
}


cat("\n=== label-switching guard (mu2 < mu1) ===\n")

# Two distinct cases produce mu2 < mu1:
#   (a) Constraint is satisfiable, but the derived mu2 happens to be smaller
#       than mu1 -- typical for "median" with high p and small positive mu1.
#   (b) For "mean" constraint: pick mu1 high enough that the lower component
#       alone overshoots E[phi]=1, forcing mu2 down.

# (a) median case: feasible derive_mu2 but ordering fails
mu2 <- derive_mu2(p = 0.9, mu1 = 0.05, sigma1 = 0.4, sigma2 = 0.3, "median")
ok(sprintf("median p=0.9 mu1=0.05 -> mu2=%.3f (feasible, but mu2 < mu1)", mu2),
   !is.na(mu2) && mu2 < 0.05)
ld <- dmixture_lognormal(c(0.5, 1, 2), p = 0.9, mu1 = 0.05,
                          sigma1 = 0.4, sigma2 = 0.3, "median")
ok("dmixture returns -Inf when ordering fails (median, label-switched)",
   all(is.infinite(ld) & ld < 0))

# (b) mean case: pick a config where derive_mu2 is feasible (m2 > 0) but
# mu2 ends up below mu1. p=0.3 mu1=0.5 gives m1 ~ 1.72; (1 - 0.3*1.72)/0.7
# > 0 so feasible, but mu2 ~ -0.42 < 0.5.
mu2 <- derive_mu2(p = 0.3, mu1 = 0.5, sigma1 = 0.3, sigma2 = 0.3, "mean")
ok(sprintf("mean p=0.3 mu1=0.5 -> mu2=%.3f (< mu1)", mu2),
   !is.na(mu2) && mu2 < 0.5)
ld <- dmixture_lognormal(c(0.5, 1, 2), p = 0.3, mu1 = 0.5,
                          sigma1 = 0.3, sigma2 = 0.3, "mean")
ok("dmixture returns -Inf when ordering fails (mean, label-switched)",
   all(is.infinite(ld) & ld < 0))


cat("\n=== infeasibility (constraint cannot be satisfied) ===\n")

# Mean constraint: if p*exp(mu1+s1^2/2) >= 1, no real mu2 works.
mu2 <- derive_mu2(p = 0.5, mu1 = 1.5, sigma1 = 1.0, sigma2 = 0.3, "mean")
ok("mean infeasibility returns NA from derive_mu2", is.na(mu2))
ld <- dmixture_lognormal(c(1), p = 0.5, mu1 = 1.5,
                          sigma1 = 1.0, sigma2 = 0.3, "mean")
ok("dmixture returns -Inf when mean-infeasible", is.infinite(ld) && ld < 0)

# Median constraint: if lower component already has > 0.5/p of its mass below
# phi=1, no mu2 >= mu1 can satisfy median(phi) = 1.
mu2 <- derive_mu2(p = 0.9, mu1 = -0.5, sigma1 = 0.4, sigma2 = 0.3, "median")
ok("median infeasibility returns NA from derive_mu2 (low mu1 + high p)",
   is.na(mu2))
ld <- dmixture_lognormal(c(1), p = 0.9, mu1 = -0.5,
                          sigma1 = 0.4, sigma2 = 0.3, "median")
ok("dmixture returns -Inf when median-infeasible",
   is.infinite(ld) && ld < 0)


cat("\n=== simulator -> empirical constraint ===\n")

set.seed(42)
# mean constraint: any mu1 works (no ordering tension)
sim_mean <- rmixture_lognormal(50000, p = 0.9, mu1 = -0.5,
                                sigma1 = 0.4, sigma2 = 0.3, "mean")
ok(sprintf("empirical mean ~ 1 (constraint=mean): mean=%.4f", mean(sim_mean$phi)),
   abs(mean(sim_mean$phi) - 1) < 0.02)

# median constraint: use config that satisfies both feasibility and ordering
sim_med <- rmixture_lognormal(50000, p = 0.85, mu1 = -0.01,
                               sigma1 = 0.4, sigma2 = 0.3, "median")
ok(sprintf("empirical median ~ 1 (constraint=median): median=%.4f",
           median(sim_med$phi)),
   abs(median(sim_med$phi) - 1) < 0.02)


cat("\n=== density integrates to 1 ===\n")

for (cfg in list(
    list(constraint = "mean",   p = 0.9,  mu1 = -0.5,  sigma1 = 0.4, sigma2 = 0.3),
    list(constraint = "median", p = 0.85, mu1 = -0.01, sigma1 = 0.4, sigma2 = 0.3)
)) {
    f <- function(phi) with(cfg, dmixture_lognormal(phi, p, mu1, sigma1, sigma2,
                                                     constraint, log = FALSE))
    I <- integrate(f, lower = 1e-6, upper = 50, rel.tol = 1e-8)$value
    ok(sprintf("integral over phi (constraint=%s) = %.6f", cfg$constraint, I),
       abs(I - 1) < 1e-4)
}


cat("\n=== nested case: p = 1 reduces to single LN at component 1 ===\n")

# At p exactly = 1, derive_mu2 returns NA (handled by guard at the top), so
# the density returns -Inf. p = 0.999 should be very close to single LN.
cfg <- list(p = 0.999, mu1 = -0.08, sigma1 = 0.4, sigma2 = 0.3)
phi_grid <- c(0.5, 1.0, 1.5, 2.0, 3.0)
mix_dens <- with(cfg, dmixture_lognormal(phi_grid, p, mu1, sigma1, sigma2,
                                          "mean", log = FALSE))
single_dens <- dlnorm(phi_grid, meanlog = cfg$mu1, sdlog = cfg$sigma1)
rel_err <- max(abs(mix_dens - single_dens) / single_dens)
ok(sprintf("p=0.999 mixture matches single LN: max relative error = %.5f",
           rel_err),
   rel_err < 0.05)


cat("\n=== priors: log_prior is finite at sensible defaults ===\n")

lp <- log_prior(p = 0.85, mu1 = -0.1, sigma1 = 0.5, sigma2 = 0.3)
ok(sprintf("log_prior(p=0.85, mu1=-0.1, s1=0.5, s2=0.3) = %.4f", lp),
   is.finite(lp))

lp <- log_prior(p = 0.85, mu1 = -0.1, sigma1 = 0.3, sigma2 = 0.5) # s2 > s1
ok("log_prior = -Inf when sigma2 > sigma1 (truncation wall)", is.infinite(lp))

# User-settable hyperparams: override just p prior
custom <- merge_hyperparams(list(p = list(alpha = 4, beta = 4)))
lp_default <- log_prior(0.5, 0, 0.4, 0.3)
lp_custom  <- log_prior(0.5, 0, 0.4, 0.3, hyperparams = custom)
ok("custom hyperparams change log_prior at p=0.5", lp_default != lp_custom)


cat("\nAll checks passed.\n")
