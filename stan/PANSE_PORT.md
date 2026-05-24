# PANSE Stan port

Stan/HMC implementations of the PANSE (PAusing time + NonSense Error)
ribosome-footprint model from AnaCoDa's C++ MCMC.  Companion to
`src/PANSEModel.cpp` (the math reference) and
`PANSE.data.analyses/s.cerevisiae/adapter.dev/` (the R driver +
validation infrastructure).

Branch: `feat/hmc-stan-panse`.  Dedicated worktree at
`~/Repositories/RibModelFramework-hmc-stan-panse/` to avoid HEAD slips
shared with the ROC Stan port.

## Why HMC?

Per-AA block-Metropolis in the C++ implementation handles each codon's
(alpha, lambda', NSE) update conditional on phi.  But these parameters
couple to per-gene phi through the partition function U and the
multiplicative ridge `mu = alpha * phi * sigma / (U * lambda')`.
Block-MH can't traverse the resulting cross-block correlations
efficiently; the ROC->Stan port (2026-05-21+) saw 100-400x ESS-per-iter
improvements from HMC + reduce_sum.  PANSE has analogous structure and
shows similar gains in this port.

## Likelihood

For each codon `c` of gene `g` at position `p`:

    mu[g, p] = alpha[c] * phi[g] * sigma[g, p] / (U * lambda'[c])
    y[g, p]  ~ NegBinomial2(mu[g, p], r = alpha[c])

where:

- `sigma[g, p]` = cumulative survive-probability through positions
  1..p-1 (off-by-one: position 1 has `sigma = 1`, position 2 multiplies
  by `prob_successful[c_1]`, etc.).
- `prob_successful[c]` = `exp(log_psuccess[c])` where
  `log_psuccess[c] = -a/(lv) + a/(lv)^2 + 0.5*(a/(lv))^2`,
  the 2nd-order Taylor expansion of the upper incomplete gamma
  integrand around the wait-time mean.  Matches
  `PANSEModel.cpp::elongationUntilIndexApproximation2ProbabilityLog`.
- `U = Z(k) / Y`, the partition function computed once at data-prep
  time from the init phi/alpha/lambda'/NSE and baked into the data
  list.  See `lib/build_panse_stan_data.R` for the sigma-aware Z
  formula (CRITICAL: must use realized per-position sigma, not assume
  sigma=1 -- see [[panse-builder-sigma-aware-Z]] memory).

## mean(phi) = 1 anchor

PANSE has a multiplicative degeneracy: `mu` is invariant under
`phi -> k*phi, lambda' -> k*lambda'`.  AnaCoDa breaks this by enforcing
`mean(phi) = 1` strictly (subtract mean(log_phi) every iteration).
The Stan port emulates this with a soft anchor in the model block:

    target += -0.5 * square((mean(log_phi) + 0.5 * sphi^2) / 0.01);

The 0.01 SD is much tighter than the per-gene-prior natural SE on
`mean(log_phi)` (= `sphi/sqrt(G)`), so it effectively pins the
ensemble mean.

The anchor IS necessary -- without it, the chain drifts along the
multiplicative ridge (R-hat 1.09-1.14 on lambda, lambda log-bias
several percent).  With it, lambda recovers to within 0.001 log-bias.
But the soft anchor creates a side-effect: a posterior coupling
`cor(sphi, mean_log_phi) = -0.91` that diagonal mass-matrix HMC can't
traverse efficiently.  See "sphi mixing" below.

## Model variants

| File | model:parameterization | What's sampled | When to use |
|---|---|---|---|
| `panse_csp_only.stan` | `csp-only` | alpha, lambda', NSE per-codon (phi as data) | v0 sanity check; data-builder + likelihood validation without the multiplicative ridge |
| `panse_csp_only_sharednse.stan` | `csp-only-sharednse` | alpha, lambda' per-codon; single shared NSE scalar | v0 with shared NSE assumption |
| `panse_basic.stan` | `basic` | + log_phi per-gene (sphi fixed) | v1 baseline; diagnostic for the multiplicative-ridge issue |
| `panse_basic_sharednse.stan` | `basic-sharednse` | same + shared NSE | v1 baseline (shared NSE) |
| `panse_sphi_est_centered.stan` | `sphi-est:centered` | + sphi | v2 when data per gene is so informative sphi posterior is essentially a delta |
| `panse_sphi_est_centered_sharednse.stan` | `sphi-est-sharednse:centered` | same + shared NSE | rare; centered with sharednse |
| `panse_sphi_est_noncentered.stan` | `sphi-est:noncentered` | (same as v2 centered, reparam'd) | per-codon NSE production |
| `panse_sphi_est_noncentered_sharednse.stan` | `sphi-est-sharednse:noncentered` | (same as v2nc, sharednse) | **PRODUCTION DEFAULT** for shared-NSE fits |
| `panse_sphi_est_sumzero_sharednse.stan` | `sphi-est-sharednse:sumzero` | (same as v2nc, sumzero) | when dense_e infeasible (param count > few thousand) AND sphi bias verified harmless at your G |

All variants:

- accept `nse_log_uniform` flag (1 = Log-Uniform prior on NSE, 0 = Natural-Uniform; AnaCoDa default is Log-Uniform)
- accept `emit_log_lik` flag (1 = compute per-position log_lik in
  generated quantities for WAIC/LOO; 0 = skip, saves ~24 GB on full
  Weinberg fits and ~30% wall time)
- hard-bound `log_alpha` and `log_lambdaPrime` (HMC stability fix
  2026-05-24; without bounds, warmup hits exp() overflow with -nan
  log-likelihood)
- use `reduce_sum` over the gene index with `grainsize` configurable
  in YAML

## Validation results

### Compact sim (100 genes x 500 positions x 22 codons; truth NSE shared = 1e-5; truth sphi = 0.5)

| Fit | lambda log-bias | lambda R-hat | lambda ESS bulk | sphi R-hat | sphi ESS | sphi covers truth | Sampling efficiency (lambda) |
|---|---|---|---|---|---|---|---|
| csp-only-sharednse                                | -0.001 | 0.999-1.009 | 2441-4128 | N/A | N/A | N/A | 8.6 ESS/s |
| basic-sharednse (sphi=0.5 fixed)                  | +0.046 | 1.035-1.050 | 71-115 | (fixed) | (fixed) | (fixed) | 0.30 ESS/s |
| sphi-est-sharednse:centered + soft anchor         | +0.002 | 1.033-1.045 | 69-153 | 1.109 | 28 | IN | 0.69 ESS/s |
| sphi-est-sharednse:noncentered + soft anchor      | +0.0001 | 1.012-1.030 | 136-301 | 1.093 | 42 | IN | 0.95 ESS/s |
| sphi-est-sharednse:noncentered + anchor + **dense_e** | +0.001 | **0.999-1.004** | **1427-1936** | **1.006** | **748** | IN | **6.7 ESS/s** |
| sphi-est-sharednse:sumzero (no anchor)            | -0.008 | 1.008-1.017 | 172-226 | 1.035 | 69 | **OUT** (post 0.562) | 0.38 ESS/s |
| sphi-est:noncentered (per-codon NSE) + anchor     | -0.008 | 1.030-1.058 | 50-124 | 1.209 | 14 | IN | 0.26 ESS/s |

Conclusions:

1. **`sphi-est-sharednse:noncentered` + `metric: dense_e` is the
   production winner.**  7x lambda sampling efficiency, 16x sphi
   sampling efficiency over diag_e at the same wall time; all R-hat
   < 1.01; all 90% CIs cover truth.
2. **Per-codon NSE works** but pays ~3x ESS penalty on lambda/phi
   because 22 weakly-identified NSE params add posterior dimensions.
   Use sharednse by default; switch to per-codon only if per-codon
   variation is the question being asked.
3. **`sumzero` is algebraically cleaner** (no anchor needed) but has
   a small sphi bias at G=100 (post 0.562, truth 0.500; suspect
   interaction with implicit (G-1)/G variance reduction).  Probably
   fine at G > 500 -- needs verification.

### sphi mixing diagnostic

When dense_e is unavailable (memory budget at very large G), the
sphi mixing pathology with diag_e + soft anchor is diagnosable by:

```r
draws <- fit$draws(variables = c("sphi", "lp__"), format = "draws_df")
log_phi_mean <- rowMeans(as.matrix(
    fit$draws(variables = "log_phi", format = "draws_matrix")))
draws$mean_log_phi <- log_phi_mean
print(round(cor(draws[, c("sphi", "mean_log_phi", "lp__")]), 3))
```

If `cor(sphi, mean_log_phi)` magnitude > 0.8 with no divergences,
the soft anchor's SD is tighter than the data's natural SE and the
ridge is built in by hand.  Fix: dense_e, sum_to_zero_vector, or
loosen anchor.  See
`PANSE.data.analyses/notes/stan.hmc.mixing.insights.md` for the full
story.

### Full Weinberg 0-ramp (G=1781, C=61, P=746K)

v0 (csp-only) fit DONE 2026-05-24 with sigma-aware-Z builder.
Clean (R-hat 0.999-1.007, 0 divergences).  Production v2 fit on
real data pending.

## Performance + scaling

- Per-position compute is dominated by the `neg_binomial_2_log_lpmf`
  call inside `partial_sum`.  Vectorized over position within gene
  when `all_unmasked == 1` (the typical case); falls back to
  per-position loop when any `like_mask[p] == 0`.
- `reduce_sum` parallelizes over genes.  `grainsize: 1` is the safe
  default; can tune up to reduce overhead on fits with many short
  genes.
- Compile time ~30-60s.  Cached at
  `stan/build/feat-hmc-stan-panse-<git-sha>/`; re-compiles only when
  the SHA changes.

## Build + run

```bash
cd ~/Repositories/PANSE.data.analyses/s.cerevisiae/adapter.dev
# Auto-compile if SHA changed; auto-load if cached:
Rscript scripts/fit.stan.R runs/<your-config>.yaml \
    --out Output/<your-run>-stan \
    --no-log-lik
```

The driver:
1. Selects the Stan model based on `fit.model` + `fit.parameterization`
2. Builds the CSR data tensor (gene_offset, codon_at_pos, y, like_mask)
   + computes sigma-aware U
3. Inits parameters from CSV files specified in YAML
4. Auto-compiles if no cached binary matches current SHA
5. Samples with 4 chains x 2 threads (configurable)
6. Saves `panse-stan-fit.rds` (full fit + wall_sec) and
   `stan-summary.rds` (parameter summary)

For the cleanest cross-comparison, always pass `--no-log-lik` unless
you actually need WAIC/LOO (saves storage + 30% wall).

## Known issues

1. **Warmup `phi[1] is -nan` info-messages** -- common during step-size
   adaptation; Stan rejects the leapfrog step and recovers.  Watch
   for these turning into actual divergences (different message).

2. **WaitingTime CSV from RMF has unclear unit relationship to NSE**
   -- to recover NSE point estimates from an RMF fit, prefer loading
   `R_objects/parameter.Rda` via AnaCoDa rather than inverting
   WaitingTime in the CSV.

3. **`sum_to_zero_vector` sphi bias at small G** -- noted above;
   needs verification at G > 500 before using in production.

## Related

- `~/Repositories/PANSE.data.analyses/notes/stan.hmc.mixing.insights.md`
  -- the diagnostic + fix story for the soft-anchor / dense_e /
  sum_to_zero trio.  Applies to ROC port too.
- `PANSE.data.analyses/s.cerevisiae/adapter.dev/CLAUDE.md`
  -- production recipe + driver invocation.
- `PANSE.data.analyses/s.cerevisiae/adapter.dev/runs/template-production.yaml`
  -- starting template for new species/datasets.
