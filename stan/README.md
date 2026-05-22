# Stan models for ROC

HMC implementations of the ROC codon-usage model, intended as an alternative
sampling backend to the in-house AnaCoDa native MCMC (Rcpp).  Same model;
different sampler.

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

| File                                | Sphi    | log_phi parameterization |
|-------------------------------------|---------|--------------------------|
| `roc_basic.stan`                    | fixed   | centered (no choice; sphi fixed -> no funnel) |
| `roc_sphi_est_centered.stan`        | sampled | log_phi as direct parameter |
| `roc_sphi_est_noncentered.stan`     | sampled | log_phi = -sphi^2/2 + sphi * z_phi |
| `roc_mixture_sphi_centered.stan`    | mixture | log_phi as direct parameter |
| `roc_mixture_sphi_noncentered.stan` | mixture | log_phi = mu1 + sigma1 * z_phi (component-1 anchored) |

### Centered vs non-centered: which one to pick

**Use centered (default) for the typical workload.**

At G in the thousands with normal codon counts per gene, the data strongly
informs each log_phi[g].  Centered mixes well in that data-informative
regime.  Non-centered creates a (sphi, z_phi) coupling the chain has to
navigate and -- empirically -- hurts dM/dEta ESS dramatically:

| At G=1000 sphi=1 (sim S288c) | dM median ESS | dEta median ESS | sphi ESS_bulk |
|---|---:|---:|---:|
| centered                | 3097          | 365             | 91            |
| non-centered            | 7             | 7               | 7             |

Non-centered earns its keep only in the **data-sparse** regime where the
prior dominates and Neal's funnel kicks in: very small G, genes with few
codons, or known small sphi posterior.

The wrapper (`fit.stan.R`) defaults to centered.  Set
`fit.parameterization: noncentered` in the YAML to override.

## Threading (reduce_sum)

All five models have the per-gene loop wrapped in `reduce_sum` for
within-chain parallelism via Stan's TBB backend.  Compile with
`STAN_THREADS=true` and pass `threads_per_chain > 1` to enable:

```r
mod <- cmdstan_model("roc_basic.stan",
                     cpp_options = list(stan_threads = TRUE),
                     exe_file = "stan/build/feat-hmc-stan-reduce-sum-XXXXXXXX-th/roc_basic")

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
stan/build/<branch>-<short_sha>[-dirty][-th]/<model_name>
```

so multiple branch / commit versions can coexist.  Build paths are
git-ignored.

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

mu2 is DERIVED from the mean(phi) = 1 constraint (PHI_CONSTRAINT_MEAN),
not free.  Label-switching guard `mu2 >= mu1` is in the model block.

## Data layout

Per-AA ragged: for AA `a`, non-ref codons sit at indices
`aa_start[a]..aa_end[a]` (1-indexed inclusive) within the dM and dEta
vectors.  See `roc_basic.stan` header for the full shape contract.
