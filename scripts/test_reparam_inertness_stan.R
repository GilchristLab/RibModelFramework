# ============================================================================
# test_reparam_inertness_stan.R -- legacy mean(phi)=1 inertness/regression guard
#                                   for stan/roc_sphi_est.stan.
#
# NOTE (2026-06-02): This script uses the OLD genomeToStanData() API with
# anchor_phi=0 / mphi_param fields that have been replaced by the phi.mphi /
# phi.sphi spec in roc_sphi_est.stan.  It CANNOT run against the updated model.
# Update to initializeStan() (default constrained("mean",1)) before re-running.
#
# Previous ROC code used the mean(phi)=1 constraint (anchor_phi=0).  This test
# pins that legacy gauge so the dEta scale-anchor / phi-centering additions cannot
# silently change behavior when OFF.  Under the legacy config
#   anchor_phi=0 (mean), noncentered=0 (centered), deta_scale_anchor=0, deta_phi_center=0
# it asserts:
#
#   (1) INERTNESS -- log_prob at a fixed parameter point is INVARIANT to the
#       reparam-only data fields:
#         (a) deta_anchor_ref  (only used when deta_scale_anchor=1)
#         (b) mphi_param        (a phantom under anchor_phi=0: no prior, mphi=-sphi^2/2)
#       If either leaks into the log-density when "off", the reparam is not inert.
#
#   (2) REGRESSION -- log_prob at that fixed point equals a frozen GOLDEN value,
#       so any future change to the legacy log-density is flagged.  (GOLDEN is the
#       value measured on this branch; update deliberately if the model legitimately
#       changes, never to silence a surprise.)
#
# NOTE: the scale-anchor *correctness* test (test_reparam_stan.R) necessarily uses
# the median anchor (the scale-anchor is only defined under a free mphi).  This
# file is the complementary LEGACY-PARITY guard, in the mean=1 gauge.
#
# Input:  inst/extdata/genome.fasta (G=8).  Usage: Rscript scripts/test_reparam_inertness_stan.R [--threads N]
# Exit:   0 = pass, 1 = fail.   Reuses the cached roc_sphi_est exe (fast after first compile).
# AI Model: Claude Opus 4.8 (phi-deta session)
# ============================================================================

suppressMessages({ library(AnaCoDa); library(cmdstanr) })
source("R/stanDataHelpers.R")

args    <- commandArgs(trailingOnly = TRUE)
ti      <- which(args == "--threads")
threads <- if (length(ti)) as.integer(args[ti + 1L]) else 4L

# Frozen reference log_prob for the legacy config + fixed theta below.
# NA on first authoring -> script prints the measured value and SKIPS the
# regression assertion; bake the printed number here and re-run to lock it.
GOLDEN     <- -5330.449352
GOLDEN_TOL <- 1e-6
INERT_TOL  <- 1e-9
fails <- 0L
check <- function(name, ok) { cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", name)); if (!ok) fails <<- fails + 1L }

cat("============================================================\n")
cat("roc_sphi_est.stan  legacy mean=1 inertness / regression guard\n")
cat("cmdstan", as.character(cmdstan_version()), " threads", threads, "\n")
cat("============================================================\n")

# ---- data: LEGACY gauge (anchor_phi=0 mean, noncentered=0 centered) --------
fasta  <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
genome <- initializeGenomeObject(file = fasta)
scuo   <- calculateSCUO(genome)$SCUO
d <- genomeToStanData(genome, scuo = scuo,
                      noncentered = 0L,        # centered (legacy)
                      anchor_phi  = 0L,        # mean(phi)=1 (legacy constraint)
                      dM.prior    = "scuo",
                      sphi_prior_mean = 1.0, sphi_prior_sd = 2.0,
                      mphi_prior_sd   = 0.5, grainsize = 1L)
d$deta_scale_anchor <- 0L
d$deta_phi_center   <- 0.0
K <- d$K; G <- d$G
cat(sprintf("STAGE 1 data: G=%d A=%d K=%d  (mean gauge, centered)\n", G, d$A, K))

# ---- compile (reuses cached exe if source unchanged) -----------------------
te <- system.time(
  mod <- cmdstan_model("stan/roc_sphi_est.stan",
                       cpp_options = list(stan_threads = TRUE),
                       compile_model_methods = TRUE, force_recompile = TRUE,
                       quiet = TRUE)
)["elapsed"]
cat(sprintf("STAGE 2 compile: %.1fs\n", te))

# ---- fixed parameter point (order: dM[K], dEta[K], latent_phi[G], log(sphi), mphi_param)
set.seed(20260601)
dM_v <- rnorm(K, 0, 0.3); dEta_v <- rnorm(K, 0, 0.5); lat_v <- rnorm(G, 0, 1.0)
mk <- function(mphi) c(dM_v, dEta_v, lat_v, log(1.0), mphi)
theta <- mk(0.0)
stopifnot(length(theta) == 2L * K + G + 2L)

lp <- function(data, theta) {
  fit <- mod$optimize(data = data, iter = 1L, algorithm = "lbfgs",
                      threads = threads, refresh = 0, show_messages = FALSE)
  fit$log_prob(unconstrained_variables = theta, jacobian = TRUE)
}

# ---- (1a) inert to deta_anchor_ref when scale-anchor OFF -------------------
cat("STAGE 3: inertness checks (legacy config, flags off)\n")
d_r0 <- d; d_r0$deta_anchor_ref <- 0.0
d_r5 <- d; d_r5$deta_anchor_ref <- 5.0
lp_r0 <- lp(d_r0, theta); lp_r5 <- lp(d_r5, theta)
cat(sprintf("  log_prob(ref=0)=%.6f  log_prob(ref=5)=%.6f  |diff|=%.2e\n",
            lp_r0, lp_r5, abs(lp_r0 - lp_r5)))
check("inert to deta_anchor_ref (scale_anchor=0)", is.finite(lp_r0) && abs(lp_r0 - lp_r5) < INERT_TOL)

# ---- (1b) mphi_param enters ONLY via its std_normal phantom prior ----------
# Under anchor_phi=0 the model gives mphi_param a normal(0,1) PHANTOM prior
# (model line ~165: "never enters mphi or likelihood").  So changing mphi_param
# must move log_prob by EXACTLY the phantom-prior difference -- proving it does
# not leak into the data likelihood in the mean gauge.
lp_m0 <- lp_r0                       # theta has mphi_param=0
lp_m3 <- lp(d_r0, mk(3.0))           # same everything, mphi_param=3
phantom_expected <- (3.0^2 - 0.0^2) / (2 * 1.0^2)   # std_normal(0,1) -> 4.5
cat(sprintf("  log_prob(mphi=0)=%.6f  (mphi=3)=%.6f  diff=%.6f  expected(phantom)=%.4f\n",
            lp_m0, lp_m3, lp_m0 - lp_m3, phantom_expected))
check("mphi_param enters ONLY via std_normal phantom prior, not the likelihood",
      abs((lp_m0 - lp_m3) - phantom_expected) < 1e-6)
check("log_prob finite", is.finite(lp_r0))

# ---- (2) regression: frozen golden -----------------------------------------
cat("STAGE 4: regression vs frozen golden\n")
if (is.na(GOLDEN)) {
  cat(sprintf("  GOLDEN is NA -- measured log_prob = %.6f  (bake this into GOLDEN)\n", lp_r0))
} else {
  cat(sprintf("  measured=%.6f  golden=%.6f  |diff|=%.2e\n", lp_r0, GOLDEN, abs(lp_r0 - GOLDEN)))
  check("legacy log_prob matches frozen golden", abs(lp_r0 - GOLDEN) < GOLDEN_TOL)
}

cat(sprintf("\nfailures: %d\nOverall: %s\n", fails, if (fails == 0L) "PASS" else "FAIL"))
quit(status = if (fails == 0L && !is.na(GOLDEN)) 0L else if (is.na(GOLDEN)) 0L else 1L)
