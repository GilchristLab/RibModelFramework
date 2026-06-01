# ============================================================================
# test_reparam_stan.R -- correctness test for the dEta scale-anchor reparam
#                        in stan/roc_sphi_est.stan, via cmdstanr $log_prob().
#
# WHAT IT PROVES (no sampling; pure log-density identities):
#
#   The scale-anchor sets dEta_eff = dEta * exp(-(mphi - ref)) and, under
#   noncentered log_phi, phi = exp(mphi + sphi*z).  Then
#       dEta_eff * phi = dEta * exp(ref + sphi*z)
#   is INDEPENDENT of the global level mphi.  So with anchor_phi=1
#   (mphi = mphi_param, free), noncentered=1, deta_phi_center=0, the ONLY
#   mphi-dependent term left in the log-density is the prior
#   mphi_param ~ Normal(0, mphi_prior_sd).  Holding dM, dEta, z, sphi fixed:
#
#       log_prob(m1) - log_prob(m2) == (m2^2 - m1^2) / (2 * mphi_prior_sd^2)    (A)
#
#   With the scale-anchor OFF (and dEta != 0), the dEta*phi product DOES depend
#   on mphi, so the same difference deviates from (A) -- the contrast is the
#   test's discriminating power.
#
#   Two flaws this test deliberately avoids (would make it vacuous):
#     - dEta must be NON-ZERO (else the likelihood is mphi-invariant either way).
#     - noncentered must be 1 (under centered, phi does not carry +mphi, so the
#       scale-anchor does NOT cancel mphi and (A) would not hold).
#     - deta_phi_center must be 0 (else the dM_at0 prior reintroduces mphi).
#
# Stages: 1 data, 2 compile, 3 build theta, 4 anchor=1 invariance (assert == A),
#         5 anchor=0 non-invariance (assert != A), 6 summary + exit code.
#
# Input:  inst/extdata/genome.fasta (G=8) -- small is fine, this is an identity.
# Usage:  Rscript scripts/test_reparam_stan.R [--threads N]
# Exit:   0 = all assertions pass, 1 = failure.
#
# AI Model: Claude Opus 4.8 (programmer+statistician agent pair, phi-deta session)
# ============================================================================

suppressMessages({ library(AnaCoDa); library(cmdstanr) })
source("R/stanDataHelpers.R")

args    <- commandArgs(trailingOnly = TRUE)
ti      <- which(args == "--threads")
threads <- if (length(ti)) as.integer(args[ti + 1L]) else 4L
TOL_INV <- 1e-6     # tight: scale-anchor=1 must match (A) to this
TOL_DEV <- 1e-2     # scale-anchor=0 must deviate from (A) by at least this
sep <- function() cat(strrep("-", 64), "\n")
fails <- 0L
check <- function(name, ok) {
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", name))
  if (!ok) fails <<- fails + 1L
}

cat("============================================================\n")
cat("roc_sphi_est.stan  dEta scale-anchor reparam correctness test\n")
cat("cmdstan", as.character(cmdstan_version()), " threads", threads, "\n")
cat("============================================================\n")

# ---- Stage 1: data ---------------------------------------------------------
sep(); cat("STAGE 1: data assembly\n"); sep()
fasta  <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
genome <- initializeGenomeObject(file = fasta)
scuo   <- calculateSCUO(genome)$SCUO
MPHI_SD <- 1.0
d <- genomeToStanData(genome, scuo = scuo,
                      noncentered = 1L,        # REQUIRED for the invariance
                      anchor_phi  = 1L,        # REQUIRED: mphi is a free param
                      dM.prior    = "scuo",
                      sphi_prior_mean = 1.0, sphi_prior_sd = 2.0,
                      mphi_prior_sd   = MPHI_SD, grainsize = 1L)
d$deta_anchor_ref <- 0.0
d$deta_phi_center <- 0.0                       # REQUIRED: c=0 isolates the test
K <- d$K; G <- d$G
cat(sprintf("  G=%d  A=%d  K=%d  mphi_prior_sd=%.3f\n", G, d$A, K, MPHI_SD))

# ---- Stage 2: compile ------------------------------------------------------
sep(); cat("STAGE 2: compile\n"); sep()
te <- system.time(
  mod <- cmdstan_model("stan/roc_sphi_est.stan",
                       cpp_options = list(stan_threads = TRUE), quiet = TRUE)
)["elapsed"]
cat(sprintf("  compiled in %.1fs\n", te))

# ---- Stage 3: parameter vectors (NON-ZERO dEta and z) ----------------------
# Unconstrained order: dM[K], dEta[K], latent_phi[G], log(sphi), mphi_param.
sep(); cat("STAGE 3: build unconstrained theta (dEta != 0, z != 0)\n"); sep()
set.seed(20260601)
dM_v   <- rnorm(K, 0, 0.3)
dEta_v <- rnorm(K, 0, 0.5)            # NON-ZERO -- makes the likelihood mphi-sensitive when anchor off
z_v    <- rnorm(G, 0, 1.0)           # latent_phi = z (noncentered); NON-ZERO
sphi_v <- 1.0                         # log(sphi)=0
m1 <- -0.2; m2 <- 0.3
mk_theta <- function(m) c(dM_v, dEta_v, z_v, log(sphi_v), m)
th1 <- mk_theta(m1); th2 <- mk_theta(m2)
stopifnot(length(th1) == 2L * K + G + 2L)
expected_A <- (m2^2 - m1^2) / (2 * MPHI_SD^2)   # = log_prob(m1) - log_prob(m2)
cat(sprintf("  theta length %d (=2K+G+2)  m1=%.2f m2=%.2f\n", length(th1), m1, m2))
cat(sprintf("  expected (A): log_prob(m1)-log_prob(m2) = %.6f\n", expected_A))

# helper: log_prob at theta for a given deta_scale_anchor flag
lp_at <- function(scale_anchor, theta) {
  dd <- d; dd$deta_scale_anchor <- as.integer(scale_anchor)
  fit <- mod$optimize(data = dd, iter = 1L, algorithm = "lbfgs",
                      threads = threads, refresh = 0, show_messages = FALSE)
  fit$init_model_methods(verbose = FALSE)      # enable $log_prob (bridgestan)
  fit$log_prob(unconstrained_variables = theta, jacobian = TRUE)
}

# ---- Stage 4: scale-anchor = 1  -> must satisfy (A) ------------------------
sep(); cat("STAGE 4: deta_scale_anchor=1  (expect mphi-INVARIANT likelihood)\n"); sep()
d1m1 <- lp_at(1L, th1); d1m2 <- lp_at(1L, th2)
delta1 <- d1m1 - d1m2
cat(sprintf("  log_prob(m1)=%.6f  log_prob(m2)=%.6f\n", d1m1, d1m2))
cat(sprintf("  delta=%.6f   expected(A)=%.6f   |err|=%.2e\n",
            delta1, expected_A, abs(delta1 - expected_A)))
check("anchor=1: log_prob diff equals prior-only formula (A)",
      is.finite(delta1) && abs(delta1 - expected_A) < TOL_INV)

# ---- Stage 5: scale-anchor = 0  -> must DEVIATE from (A) -------------------
sep(); cat("STAGE 5: deta_scale_anchor=0  (expect mphi-DEPENDENT likelihood)\n"); sep()
d0m1 <- lp_at(0L, th1); d0m2 <- lp_at(0L, th2)
delta0 <- d0m1 - d0m2
cat(sprintf("  log_prob(m1)=%.6f  log_prob(m2)=%.6f\n", d0m1, d0m2))
cat(sprintf("  delta=%.6f   expected(A)=%.6f   |dev|=%.4e\n",
            delta0, expected_A, abs(delta0 - expected_A)))
check("anchor=0: log_prob diff DEVIATES from (A) (likelihood depends on mphi)",
      is.finite(delta0) && abs(delta0 - expected_A) > TOL_DEV)
check("contrast: anchor=1 and anchor=0 deltas differ (reparam changes geometry)",
      abs(delta1 - delta0) > TOL_DEV)

# ---- Stage 6: summary ------------------------------------------------------
sep(); cat("SUMMARY\n"); sep()
cat(sprintf("  anchor=1 delta %.6f (expect %.6f)\n", delta1, expected_A))
cat(sprintf("  anchor=0 delta %.6f (deviates by %.4e)\n", delta0, abs(delta0 - expected_A)))
cat(sprintf("  failures: %d\n", fails))
cat(sprintf("Overall: %s\n", if (fails == 0L) "PASS" else "FAIL"))
quit(status = if (fails == 0L) 0L else 1L)
