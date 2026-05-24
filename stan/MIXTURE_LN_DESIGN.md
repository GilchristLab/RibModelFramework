# Mixture-LN Stan port: design rationale

Method-level guidance for the Stan HMC port of any AnaCoDa codon-usage
model that uses a mixture-of-lognormals prior on per-gene phi.
Currently applies to ROC mixture variants in this directory; the same
geometric principles apply to FONSE and PANSE mixture extensions when
they happen.

Empirical observations behind these recommendations come from the ROC
mixture-LN sim work in 2026-05; for organism-specific numbers (yeast
WAIC, Loki sphi values, simulation tables) see the analysis repos.

## TL;DR for new mixture ports

Mixture-LN identifiability under Stan/HMC is fragile in ways that the
native AnaCoDa Metropolis-Hastings sampler is not.  Three structural
fixes are load-bearing:

1. **Asymmetric informative priors on the component locations**:
   mu1 ~ N(-1, 0.5) and mu2 ~ N(+2, 0.5) under mean(phi)=1 anchor.
   Loose symmetric priors (e.g. N(0, 10)) let the posterior collapse
   to a single broad LN.
2. **sigma2 << sigma1 prior asymmetry**: e.g. Half-N(0, 0.3) for
   sigma2 vs Half-N(0, 1.0) for sigma1.  Biologically motivated --
   highly-expressed genes (component 2) cluster tightly under strong
   translational selection; bulk genes (component 1) have wide
   expression variance.  Also serves as identifiability help: prevents
   the sigma2-inflation collapse mode.
3. **Avoid derived-parameter constraints under HMC**: deriving mu2
   from (mu1, p, sigma1, sigma2) via an equality (mean=1 or
   geomean=1) creates a curved codimension-1 manifold that NUTS
   leapfrog steps routinely leave.  Mass-matrix adaptation then locks
   in the wrong basin during warmup.  Use soft Gaussian constraints
   with both mu1 and mu2 free instead.

## Failure modes (with mechanism)

### Derived-parameter manifold pathology (`centered_smooth` style)

A Stan model that derives mu2 via the mean(phi)=1 constraint:

    mu2 = log( max(1 - p*exp(mu1 + sigma1^2/2), eps) / (1 - p) )
          - sigma2^2 / 2

makes the feasible (p, mu1, sigma1, sigma2) region a curved
codimension-1 manifold in 4D hyperparameter space.  Even initialized
AT TRUE hyperparameters on a bimodal sphi sim, chains drift to a
single-broad-LN basin during warmup and never return.  Trace plots
show post-warmup iter 1 already off truth: mu1 ~ -2.5 vs truth -1.15;
sigma2 ~ 2.0 vs truth 0.4.

The mechanism: NUTS leapfrog steps perturb all parameters
simultaneously; trajectories routinely cross the manifold; safe_numer
floors + soft penalties keep the lp finite but the chain doesn't
recover the original geometry.  Mass-matrix adaptation averages over
the off-manifold excursions and locks in a metric tuned to the wrong
basin.

**Rule:** if a Stan mixture model needs to anchor a moment of the
mixture distribution (mean / geomean / median = 1), do NOT derive
hyperparameters from each other.  Use a soft Gaussian penalty on the
deviation: `log(observed_anchor) ~ normal(0, anchor_sd)`.

### Denominator singularity (`geomean` style)

A linear constraint geomean(phi)=1 with `mu2 = -p*mu1/(1-p)` has a
singularity at p -> 1.  HMC pushes the unconstrained logit-p into
large values where (1-p) ~ 1e-5; exp(mu2) overflows in a large
fraction of post-warmup draws.  Safe-floor on the LOG argument
doesn't fix this -- it bounds log(.), not the denominator.

**Rule:** any constraint with (1-p) or (1-q) in a denominator is
unsafe under HMC, even with floored intermediates.  Prefer soft
Gaussian on a LINEAR combination:

    log_anchor = p*mu1 + (1-p)*mu2;
    log_anchor ~ normal(0, anchor_sd);

with both mu1 and mu2 as free parameters.  No singularity, no
overflow risk.

### Soft-constraint variant works but is underconstrained

A soft-constraint variant (`geomean_soft`-style: free mu1, mu2;
Gaussian soft prior on the linear combination) is numerically stable
and has no boundary issues.  But under loose hyperparameter priors,
sigma2 can still inflate to absorb both components into one broad
LN, and mu2 drifts toward mu1.

This is why the asymmetric mu priors AND tightened sigma2 prior (the
first two rules in the TL;DR) are needed.  Soft-constraint geometry
is necessary but not sufficient.

## Single-LN insights (apply to single-LN models too)

### SCUO init for log_phi is load-bearing

Single-LN MAP for log_phi without SCUO-based init lands in a LOW
basin (sphi ~ 0.06 in empirical work).  With SCUO-rank-based init
for log_phi, MAP lands in a HIGH basin (sphi ~ 1+ for typical
organisms).  Model-comparison metrics (WAIC etc.) strongly prefer
the HIGH basin in real-data fits.

The init mechanism: compute per-gene SCUO via AnaCoDa::calculateSCUO,
rank-normalize, map to standard-normal quantiles, scale by
`sphi_seed`.  Use as init values for `log_phi` for ALL chains (with
optional per-chain jitter).

**Rule:** default init_mode for any phi-based Stan model should be
SCUO-based.  Pure MAP from random init is a bistability trap.

### Bistability is a model property, not a data artifact

Single-LN ROC has a dual-basin posterior structure that appears
identically in sim and real data.  The LOW basin corresponds to
"selection signal is weak, all genes near phi=1."  The HIGH basin
corresponds to real biological expression heterogeneity.  Chain
initialization can land in either basin; transitions between them
are slow.

Any phi-based model (ROC / PANSE / FONSE single-LN OR mixture)
likely exhibits similar bistability.  Anticipate this when
validating new Stan ports; validate convergence to the biologically-
plausible basin, not just convergence to *some* basin.

### Identifiability rule: G * len_aa >= 500k AND G >= 250

For mixture-LN sphi identifiability under standard sim conditions:
total amino acid count G * mean(len_aa) >= 500,000 AND G >= 250 (so
the low-component sample is large enough).  Below this regime the
posterior is under-determined regardless of sampler.  Hyperparameter
posteriors are also bounded by the count of genes in the rarer
mixture component (typically 1-p ~ 0.05-0.10 of G, so 12-50 genes
at the lower bound).

For PANSE the relevant unit is total RFP-position counts rather than
AA counts; needs empirical recalibration when PANSE mixture comes
online.

## Anchor convention: mean / median / geomean = 1

These are not equivalent under a mixture.  Three conventions in use:

| Convention | Math | Used by |
|---|---|---|
| `median(phi)=1` | mPhi = -sphi^2/2 (single-LN); generalizes via per-component median for mixture | Native AnaCoDa (mature ROC) |
| `mean(phi)=1` | E[phi] = 1; mPhi = 0 under single-LN; sum constraint under mixture | Stan ROC port (current) |
| `geomean(phi)=1` | exp(E[log phi]) = 1; linear constraint p*mu1 + (1-p)*mu2 = 0 under mixture | Stan geomean variants |

These differ by translation in log space.  At sphi=1 the native vs
Stan single-LN scales differ by ~1.66x.  Mixture introduces
additional cross-convention translation.

**Rule for new ports:**
- Pick ONE convention; document it in the YAML / output metadata.
- When comparing Stan results to native results, translate
  hyperparameters explicitly before comparing.
- Per-gene phi values are sensitive to convention -- a multiplicative
  shift, not a re-ranking, so SCUO-correlation etc. is convention-
  invariant but absolute phi values are not.

## Recommended starting point for new mixture ports

Use a `geomean_soft`-style architecture as the template:

- Free hyperparameters: `p`, `mu1`, `mu2`, `sigma1`, `sigma2`
- Hyperpriors per the TL;DR (asymmetric mu, sigma2 < sigma1)
- Soft anchor: `log(anchor_phi) ~ normal(0, ~0.05)` where
  `anchor_phi` is your chosen mean/geomean
- log_phi init via SCUO ranking
- Threading via `reduce_sum` over per-gene log-likelihood
- Stage-1 SCUO init only (no Stage-2 MAP pre-fit needed; MAP without
  SCUO bias finds the LOW basin and corrupts hyperparameter init)

`roc_mixture_sphi_geomean_soft.stan` in this directory is the
canonical reference implementation.  `roc_mixture_sphi_ordered.stan`
is a label-switching-safe alternative using Stan's `ordered[2]`.

## Dense mass matrix (`metric: dense_e`) is high-leverage

PANSE Stan port (2026-05-24, sphi-est sim sweep, compact 100-gene
simulation, results in `PANSE.data.analyses/.../adapter.dev/Output/
compare_all_sim_fits.rds`) found:

| Variant | sphi ESS | sphi R-hat | log_phi ESS-min |
|---|---|---|---|
| v2-centered (diag_e) | 404 | 1.011 | 31 |
| v2-noncentered (diag_e) | 25 | 1.154 | 23 |
| v2-noncentered-anchor (diag_e) | 42 | 1.093 | 143 |
| **v2-noncentered-anchor + dense_e** | **748** | **1.006** | **1263** |

`dense_e` gave 17.8x sphi ESS over diag_e + noncentered + anchor, with
R-hat closing from 1.093 to 1.006.  log_phi ESS-min jumped 8.8x.
Non-centered ALONE without dense_e actually HURT sphi mixing on this
data (25 vs 404 for centered).  The dense mass matrix learns the
sphi <-> z_phi negative cross-correlation (-0.91 observed) that
diag_e cannot represent.

### When to use dense_e

- Posterior has known cross-block correlation (typical for mixture
  hyperparameters; sphi <-> log_phi in non-centered; mu1 <-> mu2 in
  geomean_soft).
- R-hat for a scalar hyperparameter fails to close even with longer
  warmup -- often indicates a correlation the diag metric can't
  navigate.
- Affordable cost: warmup roughly 2x for full mass-matrix estimation;
  per-iter cost is also higher.  Worth it when you'd otherwise need
  to double sample budget.

### When to avoid dense_e

- Per-parameter scale very heterogeneous AND no strong cross-
  correlation: diag_e is cheaper and equivalent.
- Very large parameter dim (G in tens of thousands): the metric grows
  as O(D^2) memory and warmup cost.  At G=5000+ the dense mass
  matrix is ~25M entries; dense_e becomes prohibitive.
- Quick-iteration prototyping: pay the dense_e warmup tax only when
  you're going for final results.

### Wiring

Cmdstanr exposes `metric = "dense_e"` via `mod$sample()`.  PANSE and
ROC adapter.dev wrappers expose a `fit.metric` YAML field + a
`--metric` CLI flag, both gating to "diag_e" / "dense_e" / "unit_e".

## Centered vs non-centered for log_phi (Neal's funnel)

For data-dense regimes (G in the thousands with ~1000+ AA positions
per gene), the **centered** parameterization for log_phi (sampling
log_phi directly) is the default and works well.

For data-sparse regimes (small G, low per-gene count, or just a
sparse rare-component), **non-centered** reparametrization of log_phi
removes Neal's funnel in the (log_phi, sigma) plane:

    log_phi[g] = mu1 + sigma1 * z_phi[g]      // transformed
    z_phi[g] ~ p * N(0, 1) + (1-p) * N((mu2-mu1)/sigma1, sigma2/sigma1)

Anchoring on component 1 fully removes the funnel for the bulk
component; the high component still sees sigma2 through the scale
ratio so its funnel is partially mitigated but not removed.

`roc_mixture_sphi_geomean_soft_noncentered.stan` is the non-centered
analog of `geomean_soft`.

### When to switch from centered to non-centered

Symptoms that suggest non-centered:
- R-hat for sigma1 fails to close even with longer warmup / more chains
- Chains stick in tight-sigma1 regions and explore wider-sigma1 only
  rarely
- Divergence rate spikes correlate with sigma1 excursions toward zero
- Per-gene log_phi posteriors near the bulk-component mean show wider
  posterior than the prior suggests they should

For mixture models on MAG-scale data (1500-3000 genes) with a small
rare component (e.g., 1-p = 0.05 -> ~75-150 high-expression genes),
the rare-component-side sigma is the most fragile dimension.  Non-
centering provides geometric improvement for the bulk but the rare
component remains a difficulty.

### Per-component non-centering (future work)

A more thorough fix would non-center BOTH components separately,
e.g. via mixture-model-aware reparametrization with explicit
component labels.  Stan's discrete-parameter limitation makes this
non-trivial (each gene's component assignment is marginalized in the
log_mix expression).  See e.g. the Stan mixture model documentation
for label-switching reparametrizations.

For now, the single-component non-centering above is the simplest
useful intervention.

## See also

- `stan/README.md` for the full inventory of Stan models in this directory.
- `~/Repositories/ROC.data.analyses.adaptive-mcmc-dev/s.cerevisiae/adapter.dev/scripts/fit.stan.R` -- canonical wrapper script (YAML-driven; SCUO init; scuo_low_csv dM prior; per-parameterization init translation).
- `~/Repositories/PANSE.data.analyses/notes/mixture-ln-port-insights-from-roc.md` -- organism-specific PANSE-side observations.
