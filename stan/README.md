# Stan models for ROC

> **Branch status: RETIRED 2026-05-26.**
> All content merged into `upstream-develop` (GilchristLab/RibModelFramework).
> This branch (`feat/hmc-stan-reduce-sum`) is preserved read-only for history.
> Active development continues on `upstream-develop`.

HMC implementations of the ROC codon-usage model, intended as an alternative
sampling backend to the in-house AnaCoDa native MCMC (Rcpp).  Same model;
different sampler.

**Mixture-LN design rationale:** see `MIXTURE_LN_DESIGN.md` in this
directory for cross-method (ROC / PANSE / FONSE) Stan mixture-LN
prior structure, failure-mode catalog, SCUO init mechanics, anchor
convention discussion, and recommended starting architecture.

## When to use Stan vs native

Stan's per-iteration mixing is 100-400x better than native+ar-and-cov on the
sim S288c benchmark (G=1000).  The phi-dEta posterior correlation that
bottlenecks native (dEta mixes 2.4x slower than dM, rare-codon ESS = 10-30)
is handled natively by HMC gradient-based exploration.

Use Stan when:
- You want accurate posteriors for slow-mixing codons.
- You're estimating sphi or the phi-mixture hyperparameters and need
  reliable ESS on those scalars.

Use native when:
- You need AnaCoDa-format outputs (Parameter / Trace objects, .Rdata, .rst)
  for downstream tools that assume that state.
- The phi-mixture C++ branch is required for a chunked / restart workflow.

The two backends are alternatives, not replacements.

## File naming

**Active (unified) files:**

| File                                      | Sphi    | Notes |
|-------------------------------------------|---------|-------|
| `roc_basic.stan`                          | fixed   | centered only (no funnel when sphi fixed) |
| `roc_sphi_est.stan`                       | sampled | unified: centered + NC via `noncentered` data flag |
| `roc_mixture_sphi_geomean_soft.stan`      | mixture | unified: centered + NC via `noncentered` data flag; geomean(phi)=1 |

**Legacy files (retained for reference, not called by fit.stan.R):**

| File                                            | Superseded by |
|-------------------------------------------------|---------------|
| `roc_sphi_est_centered.stan`                    | `roc_sphi_est.stan` with `noncentered=0` |
| `roc_sphi_est_noncentered.stan`                 | `roc_sphi_est.stan` with `noncentered=1` |
| `roc_mixture_sphi_geomean_soft_noncentered.stan`| `roc_mixture_sphi_geomean_soft.stan` with `noncentered=1` |
| `roc_mixture_sphi_centered.stan`                | geomean_soft variant preferred |
| `roc_mixture_sphi_noncentered.stan`             | geomean_soft variant preferred |
| `roc_mixture_sphi_ordered.stan`                 | geomean_soft variant preferred |
| `roc_mixture_sphi_centered_smooth.stan`         | geomean_soft variant preferred |
| `roc_mixture_sphi_geomean.stan`                 | geomean_soft variant preferred |

### Unified binary: noncentered data flag

Both `roc_sphi_est.stan` and `roc_mixture_sphi_geomean_soft.stan` cover
centered and non-centered parameterizations in a single compiled binary via:

```stan
data { int<lower=0,upper=1> noncentered; ... }
parameters { vector[G] latent_phi; ... }
transformed parameters {
    vector[G] log_phi = noncentered ? (mphi + sphi * latent_phi) : latent_phi;
}
```

When `noncentered=0`: `latent_phi` IS `log_phi` (centered).
When `noncentered=1`: `latent_phi` is z-score; `log_phi` reconstructed in
transformed parameters.

No recompile needed to switch parameterizations -- just change the YAML.

### Centered vs non-centered: which one to pick

**Use centered (default).** It is the empirically-validated choice at
production G.

At G=1000 (sim S288c) centered wins on all parameters:

| At G=1000 sphi=1 (sim S288c) | dM median ESS | dEta median ESS | sphi ESS_bulk |
|---|---:|---:|---:|
| centered                | 3097          | 365             | 91            |
| non-centered            | 7             | 7               | 7             |

At G=2000 (real S288c mixture-LN cross-validation, 2026-05-25/26):
- Centered: E-BFMI < 0.3 on all 4 chains (sigma1 funnel).
- Non-centered: E-BFMI < 0.3 on 1/4 chains (improved), but chain 1
  mode-flipped entirely (R-hat 1.44-1.76, sep CI includes negative).

Conclusion: Stan HMC cannot reliably cross-validate the ROC mixture-LN
model regardless of parameterization.  The posterior is genuinely
multi-modal (label-switching + phi-dEta coupling); NC fixes the sigma1
funnel but introduces mode-flip sensitivity.  Native MCMC is authoritative
for the mixture model.

Non-centered is appropriate only in the data-sparse regime (very small G,
few codons per gene, or known small sphi posterior).

## Threading (reduce_sum)

All models have the per-gene loop wrapped in `reduce_sum` for within-chain
parallelism via Stan's TBB backend.  **All binaries are always compiled with
`STAN_THREADS=true`** -- `threads_per_chain` is a runtime parameter, not a
compile-time choice.

```r
mod <- cmdstan_model("roc_sphi_est.stan",
                     cpp_options = list(stan_threads = TRUE),
                     exe_file = "stan/build/feat-hmc-stan-reduce-sum-XXXXXXXX-th/roc_sphi_est")

fit <- mod$sample(data = data.list,
                  threads_per_chain = 4,
                  chains = 4, ...)
```

Expected speedup: ~0.7-0.85x linear (TBB scheduling overhead ~1-2% of compute
at G=1000).  Threads=4 typically gives ~3x, threads=8 typically gives ~5-6x.
Diminishing returns past 8 threads.

Choose between more chains and more threads:
- More chains -> more total ESS, same wall.  Plateau at ~4 chains for
  convergence diagnostics; further chains are decorative.
- More threads -> less wall per chain.  Only way to actually wait less.

For G in thousands at fixed core budget: 4 chains x N threads is usually
the sweet spot.

`grainsize` controls reduce_sum partition size.  Default 1 lets TBB choose
automatically (recommended).  Explicit `G / (2 * threads_per_chain)` is
sometimes faster for small G.

## Branch-tagged binaries

The wrapper compiles each model to:

```
stan/build/<branch>-<short_sha>-th/<model_name>
```

so multiple branch / commit versions can coexist.  Build paths are
git-ignored.  The `-th` suffix is always present (STAN_THREADS always on).

## Hyperprior conventions

Mirrors C++ defaults (`src/Parameter.cpp:phiMixtureHyper_*`):

```yaml
fit:
  mixture:
    p.alpha: 8           # Beta(8, 2) -> E[p] = 0.8
    p.beta:  2
    mu1.prior.mean: 0    # N(0, 10) on mu1
    mu1.prior.sd:   10
    sigma1.prior.scale: 1   # half-normal(0, 1) on sigma1
    sigma2.prior.scale: 1   # half-normal(0, 1) on sigma2
```

mu2 is DERIVED from the geomean(phi) = 1 constraint (PHI_CONSTRAINT_GEOMEAN),
not free.  Label-switching guard `mu2 >= mu1` is in the model block.

## Data layout

Per-AA ragged: for AA `a`, non-ref codons sit at indices
`aa_start[a]..aa_end[a]` (1-indexed inclusive) within the dM and dEta
vectors.  See `roc_basic.stan` header for the full shape contract.
