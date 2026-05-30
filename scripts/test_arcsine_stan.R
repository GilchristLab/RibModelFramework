# ============================================================================
# test_arcsine_stan.R -- End-to-end test of roc_arcsine.stan via cmdstanr.
#
# Stages:
#   1. Compile  roc_arcsine.stan with STAN_THREADS=true
#   2. Check    log_prob at default init is finite (gradient sanity)
#   3. optimize()    -- MAP estimate (primary warm-start use case)
#   4. variational() -- ADVI mean-field  (ELBO warm-start use case)
#
# Input:  /tmp/simROC_2000.fasta  (G=2000 simulated genes from benchmark_arcsine.R)
#         Falls back to inst/extdata/genome.fasta (G=8) if the sim file is absent.
#
# Output: Printed summary at each stage; exits 0 on success.
#
# Usage:
#   Rscript scripts/test_arcsine_stan.R
#   Rscript scripts/test_arcsine_stan.R --threads 8
# ============================================================================

suppressMessages({
  library(AnaCoDa)
  library(cmdstanr)
})
source("R/stanDataHelpers.R")

# --------------------------------------------------------------------------
# CLI argument: --threads N  (default 4)
# --------------------------------------------------------------------------
args     <- commandArgs(trailingOnly = TRUE)
thr_idx  <- which(args == "--threads")
threads  <- if (length(thr_idx)) as.integer(args[thr_idx + 1]) else 4L

sep <- function() cat(strrep("-", 60), "\n")

cat("============================================================\n")
cat("roc_arcsine.stan end-to-end test\n")
cat("cmdstan:", cmdstan_version(), "  threads:", threads, "\n")
cat("============================================================\n\n")

# --------------------------------------------------------------------------
# 1. Data preparation
# --------------------------------------------------------------------------
sep()
cat("STAGE 1: Data preparation\n")
sep()

sim_fasta   <- "/tmp/simROC_2000.fasta"
extdata     <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
genome_file <- if (file.exists(sim_fasta)) sim_fasta else extdata
cat("Genome:", genome_file, "\n")

genome <- initializeGenomeObject(file = genome_file)
G      <- genome$getGenomeSize(FALSE)
cat("Genes:", G, "\n")

# Compute SCUO once; reused for both dM prior and phi init
t_scuo <- system.time(scuo <- calculateSCUO(genome)$SCUO)["elapsed"]
cat(sprintf("SCUO computed in %.1f s  range=[%.3f, %.3f]\n",
            t_scuo, min(scuo), max(scuo)))

# genomeToStanData: uses SCUO-estimated dM_prior_mean (dM.prior="scuo" default)
d <- genomeToStanData(
  genome,
  scuo         = scuo,   # reuse pre-computed SCUO
  approx_min_n = 20L,
  noncentered  = 0L,
  anchor_phi   = 0L
)

cat("A =", d$A, "  K =", d$K, "  approx_min_n =", d$approx_min_n, "\n")
pct_arcsine <- mean(d$N_ga >= d$approx_min_n) * 100
cat(sprintf("Gene/AA pairs using arcsine branch: %.1f%%\n", pct_arcsine))
cat(sprintf("dM_prior_mean range: [%.3f, %.3f]  (from low-SCUO genes)\n",
            min(d$dM_prior_mean), max(d$dM_prior_mean)))
cat("\n")

# --------------------------------------------------------------------------
# 2. Compile
# --------------------------------------------------------------------------
sep()
cat("STAGE 2: Compile roc_arcsine.stan\n")
sep()

branch   <- system("git branch --show-current", intern = TRUE)
sha      <- system("git rev-parse --short HEAD", intern = TRUE)
build_dir <- file.path("stan", "build", paste0(branch, "-", sha, "-th"))
dir.create(build_dir, showWarnings = FALSE, recursive = TRUE)
exe_path <- file.path(build_dir, "roc_arcsine")

t_compile <- system.time({
  mod <- cmdstan_model(
    "stan/roc_arcsine.stan",
    cpp_options = list(stan_threads = TRUE),
    exe_file    = exe_path,
    quiet       = FALSE
  )
})["elapsed"]
cat(sprintf("Compiled in %.1f s  ->  %s\n\n", t_compile, exe_path))

# --------------------------------------------------------------------------
# 3. Log-prob sanity check at default init
# --------------------------------------------------------------------------
sep()
cat("STAGE 3: log_prob sanity check\n")
sep()

# Default Stan init: parameters drawn from Uniform(-2, 2) on unconstrained space.
# Just run optimize with 1 iteration to verify the log_prob is finite.
t_lp <- system.time({
  fit_lp <- tryCatch(
    mod$optimize(
      data        = d,
      threads     = threads,
      iter        = 1L,
      algorithm   = "lbfgs",
      refresh     = 0,
      show_messages = FALSE
    ),
    error = function(e) e
  )
})["elapsed"]

if (inherits(fit_lp, "error")) {
  cat("FAIL: log_prob check errored:", conditionMessage(fit_lp), "\n")
  quit(status = 1)
}
lp_init <- fit_lp$lp()
cat(sprintf("lp at iter 1: %.2f  (%.1f s)\n", lp_init, t_lp))
if (!is.finite(lp_init)) {
  cat("FAIL: lp is not finite at default init\n")
  quit(status = 1)
}
cat("PASS: lp is finite\n\n")

# --------------------------------------------------------------------------
# 4. optimize() -- MAP estimate
# --------------------------------------------------------------------------
sep()
cat("STAGE 4: optimize() -- MAP\n")
sep()

# SCUO-based init: latent_phi proportional to log(SCUO), centered and scaled.
# dM=0, dEta=0 gives uniform codon probs at the start; phi varies by gene from
# codon bias -- already in the right ballpark without any MCMC.
scuo_init <- genomeToStanInit(genome, d, scuo = scuo)

t_opt <- system.time({
  fit_opt <- tryCatch(
    mod$optimize(
      data        = d,
      init        = list(scuo_init),
      threads     = threads,
      iter        = 5000L,
      algorithm   = "lbfgs",
      refresh     = 200,
      tol_rel_grad = 1e-8,
      show_messages = TRUE
    ),
    error = function(e) e
  )
})["elapsed"]

if (inherits(fit_opt, "error")) {
  cat("FAIL:", conditionMessage(fit_opt), "\n")
  quit(status = 1)
}

lp_opt <- fit_opt$lp()
cat(sprintf("\nMAP lp*: %.2f  (%.1f s)\n", lp_opt, t_opt))

draws_opt <- fit_opt$draws(format = "df")
sphi_opt  <- draws_opt$sphi
dM_opt    <- as.numeric(draws_opt[, paste0("dM[", 1:d$K, "]")])
dEta_opt  <- as.numeric(draws_opt[, paste0("dEta[", 1:d$K, "]")])

cat(sprintf("sphi (MAP):  %.4f\n", sphi_opt))
cat(sprintf("dM   range:  [%.3f, %.3f]  median %.3f\n",
            min(dM_opt), max(dM_opt), median(dM_opt)))
cat(sprintf("dEta range:  [%.3f, %.3f]  median %.3f\n",
            min(dEta_opt), max(dEta_opt), median(dEta_opt)))

phi_opt   <- as.numeric(draws_opt[, paste0("phi[", 1:d$G, "]")])
cat(sprintf("phi  range:  [%.3f, %.3f]  median %.3f  geomean %.3f\n",
            min(phi_opt), max(phi_opt), median(phi_opt), exp(mean(log(phi_opt)))))

# Sanity checks
ok_opt <- TRUE
if (!is.finite(lp_opt))            { cat("FAIL: lp* not finite\n");        ok_opt <- FALSE }
if (sphi_opt <= 0)                  { cat("FAIL: sphi <= 0\n");             ok_opt <- FALSE }
if (any(!is.finite(dM_opt)))        { cat("FAIL: NaN/Inf in dM\n");        ok_opt <- FALSE }
if (any(!is.finite(dEta_opt)))      { cat("FAIL: NaN/Inf in dEta\n");      ok_opt <- FALSE }
if (any(!is.finite(phi_opt)))       { cat("FAIL: NaN/Inf in phi\n");       ok_opt <- FALSE }
if (any(phi_opt <= 0))              { cat("FAIL: phi <= 0\n");              ok_opt <- FALSE }
if (ok_opt) cat("PASS: optimize() returned sane MAP estimates\n")
cat("\n")

# --------------------------------------------------------------------------
# 5. variational() -- ADVI mean-field
# --------------------------------------------------------------------------
sep()
cat("STAGE 5: variational() -- ADVI mean-field\n")
sep()

t_vi <- system.time({
  fit_vi <- tryCatch(
    mod$variational(
      data           = d,
      init           = list(scuo_init),
      threads        = threads,
      algorithm      = "meanfield",
      iter           = 20000L,
      grad_samples   = 1L,
      elbo_samples   = 100L,
      tol_rel_obj    = 0.01,
      output_samples = 1000L,
      refresh        = 2000,
      show_messages  = TRUE
    ),
    error = function(e) e
  )
})["elapsed"]

if (inherits(fit_vi, "error")) {
  cat("FAIL:", conditionMessage(fit_vi), "\n")
  quit(status = 1)
}

cat(sprintf("\nADVI wall time: %.1f s\n", t_vi))

draws_vi  <- fit_vi$draws(format = "df")
sphi_vi   <- draws_vi$sphi
dM_vi     <- as.numeric(unlist(draws_vi[, paste0("dM[", 1:d$K, "]")]))
dEta_vi   <- as.numeric(unlist(draws_vi[, paste0("dEta[", 1:d$K, "]")]))
phi_vi    <- as.numeric(unlist(draws_vi[, paste0("phi[", 1:d$G, "]")]))

cat(sprintf("sphi (ADVI):  mean=%.4f  sd=%.4f  95%%CI=[%.3f,%.3f]\n",
            mean(sphi_vi), sd(sphi_vi),
            quantile(sphi_vi, 0.025), quantile(sphi_vi, 0.975)))
cat(sprintf("dM   (ADVI):  median=%.3f  sd=%.3f\n",
            median(dM_vi), sd(dM_vi)))
cat(sprintf("dEta (ADVI):  median=%.3f  sd=%.3f\n",
            median(dEta_vi), sd(dEta_vi)))
cat(sprintf("phi  geomean (ADVI): %.3f  (should be ~1 for anchor_phi=0)\n",
            exp(mean(log(phi_vi[phi_vi > 0])))))

# Compare MAP vs ADVI means
cat(sprintf("\nsphi: MAP=%.4f  ADVI_mean=%.4f\n", sphi_opt, mean(sphi_vi)))

ok_vi <- TRUE
if (any(!is.finite(sphi_vi)))  { cat("FAIL: NaN/Inf in sphi ADVI draws\n"); ok_vi <- FALSE }
if (any(sphi_vi <= 0))          { cat("FAIL: sphi <= 0 in ADVI draws\n");    ok_vi <- FALSE }
if (any(!is.finite(dM_vi)))    { cat("FAIL: NaN/Inf in dM ADVI draws\n");   ok_vi <- FALSE }
if (ok_vi) cat("PASS: variational() returned sane draws\n")

# --------------------------------------------------------------------------
# 6. sample() -- HMC gradient check (1 chain, short run)
# --------------------------------------------------------------------------
sep()
cat("STAGE 6: sample() -- HMC gradient check\n")
sep()
cat("1 chain, warmup=300, sampling=200, threads_per_chain=", threads, "\n")

t_hmc <- system.time({
  fit_hmc <- tryCatch(
    mod$sample(
      data             = d,
      init             = list(scuo_init),
      chains           = 1L,
      iter_warmup      = 300L,
      iter_sampling    = 200L,
      threads_per_chain = threads,
      refresh          = 100,
      show_messages    = TRUE
    ),
    error = function(e) e
  )
})["elapsed"]

if (inherits(fit_hmc, "error")) {
  cat("FAIL:", conditionMessage(fit_hmc), "\n")
  quit(status = 1)
}

cat(sprintf("\nHMC wall time: %.1f s\n", t_hmc))

diag       <- fit_hmc$diagnostic_summary(quiet = TRUE)
n_div      <- sum(diag$num_divergent)
n_max_td   <- sum(diag$num_max_treedepth)
ebfmi      <- diag$ebfmi

cat(sprintf("Divergences:      %d  (want 0 or very few)\n", n_div))
cat(sprintf("Max treedepth:    %d\n", n_max_td))
cat(sprintf("E-BFMI:           %.3f  (want > 0.3)\n", ebfmi))

draws_hmc  <- fit_hmc$draws(format = "df")
sphi_hmc   <- draws_hmc$sphi
phi_hmc    <- as.numeric(unlist(draws_hmc[, paste0("phi_out[", 1:min(d$G, 10), "]")]))

sphi_summ  <- fit_hmc$summary("sphi")
cat(sprintf("sphi:  mean=%.4f  sd=%.4f  rhat=%.3f  ess_bulk=%.0f\n",
            sphi_summ$mean, sphi_summ$sd,
            sphi_summ$rhat, sphi_summ$ess_bulk))
cat(sprintf("phi geomean (HMC sample): %.3f\n",
            exp(mean(log(phi_hmc[phi_hmc > 0])))))

ok_hmc <- TRUE
if (n_div > 10)              { cat("WARN: many divergences\n") }
if (ebfmi < 0.3)             { cat("WARN: low E-BFMI (", round(ebfmi,3), ")\n") }
if (sphi_summ$rhat > 1.1)    { cat("WARN: high R-hat for sphi\n") }
if (any(!is.finite(sphi_hmc))){ cat("FAIL: NaN/Inf in sphi HMC draws\n"); ok_hmc <- FALSE }
if (any(sphi_hmc <= 0))       { cat("FAIL: sphi <= 0 in HMC draws\n");    ok_hmc <- FALSE }
if (ok_hmc) cat("PASS: sample() completed with finite draws\n")
cat("\n")

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
sep()
cat("SUMMARY\n")
sep()
cat(sprintf("  Genome:       %s  (G=%d)\n", basename(genome_file), G))
cat(sprintf("  Arcsine %%:    %.1f%%\n", pct_arcsine))
cat(sprintf("  Compile:      %.1f s\n", t_compile))
cat(sprintf("  optimize():   %.1f s   lp*=%.1f   sphi=%.4f\n", t_opt, lp_opt, sphi_opt))
cat(sprintf("  variational():%.1f s   sphi mean=%.4f  sd=%.4f\n",
            t_vi, mean(sphi_vi), sd(sphi_vi)))
cat(sprintf("  sample():     %.1f s   sphi mean=%.4f  divs=%d  E-BFMI=%.3f\n",
            t_hmc, mean(sphi_hmc), n_div, ebfmi))
cat(sprintf("  Threads:      %d\n", threads))

all_ok <- ok_opt && ok_vi && ok_hmc
cat(sprintf("\nOverall: %s\n", if (all_ok) "PASS" else "FAIL"))
if (!all_ok) quit(status = 1)
