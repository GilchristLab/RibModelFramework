# CHANGES IN AnaCoDa 0.1.10 (unreleased; in-progress)

## NEW FEATURES

- **FONSE elongation cost `a2` promoted to a full parameter.** The positional
  slope in the FONSE nonsense-error cost function `beta(k) = a1 + a2 * k` was
  previously hard-coded as the literal `4.0`. It is now a first-class
  `FONSEParameter` value (`a2`) mirroring the initiation cost `a1`: log-normal
  random-walk proposal, adaptive proposal width, trace, restart-file I/O
  (`elongation_cost` / `std_elongation_cost`), and R/Rcpp exposure. New R
  argument `init.elongation.cost` (default 4); new trace selector
  `what = "ElongationCost"` for `plot()` and `as.mcmc()`; new methods
  `parameter$estimateElongationCost()` / `parameter$fixedElongationCost()`.

- **BEHAVIOR CHANGE: FONSE `a1` and `a2` are now fixed at 4 by default.**
  Previously `a1` was estimated by default, which is weakly identified and could
  run away (observed reaching ~1e4 vs the biological ~4). Both coefficients are
  now held fixed at their initial value of 4 unless estimation is opted into via
  `parameter$estimateInitiationCost()` / `parameter$estimateElongationCost()`.
  This changes default FONSE results: `a1` no longer drifts, and the FONSE codon
  probabilities with default settings are bit-identical to the prior model with
  `a1 = a2 = 4`.

- **Pluggable CSP adaptive proposal-width schemes.** `Parameter` now
  holds a `std::unique_ptr<CSPAdaptationStrategy>` that controls how
  codon-specific parameter (CSP) proposal widths are tuned during MCMC.
  The default scheme (`native`) preserves the historical in-house
  0.8/1.2 multiplicative logic and is bit-identical to prior fits.
  A second scheme (`andrieu_thoms`) implements Andrieu and Thoms 2008
  (*Statistics and Computing* 18:343-373), Algorithm 4: a continuous
  Robbins-Monro update on `log(std_csp)` with a diminishing step
  schedule.  Selectable from R via `AdaptiveScheme.AndrieuThoms(...)`
  + `parameter$setCSPAdaptationScheme(...)`, or from a YAML config
  block in the v.3 Lokiarchaeota pipeline.  See
  `docs/csp-adaptation-api.md` for the full design.

- `schemes.available()`, `is.AdaptiveScheme()`, `print.AdaptiveScheme`,
  `format.AdaptiveScheme`, `adaptive.scheme.diagnostics()` exported.

- Rcpp methods on Parameter: `setCSPAdaptationScheme(name, params)` +
  `getCSPAdaptationSchemeName()`.

- **Stan ROC single-LN reparameterization flags for dEta mixing (closes #6).**
  `stan/roc_sphi_est.stan` gains two opt-in data flags, both OFF by default and
  byte-identical to the prior model when off:
  - `deta_scale_anchor` — samples `dEta` on a scale-anchored coordinate
    (`dEta = dS * exp(-(mphi - ref))`) so the global phi level `mphi` cancels from
    the `dEta * phi` product, breaking the multiplicative dEta-phi scale ridge that
    dominates dEta mixing.  Requires `anchor_phi = median` (a free, sampled `mphi`).
  - `deta_phi_center` — centers phi in the likelihood, applying the dM prior to the
    intercept-at-0 (unit-Jacobian shear).

  Paired with the non-centered `log_phi` parameterization (`noncentered = 1`), the
  scale-anchor makes the data likelihood exactly invariant to `mphi`, removing the
  Neal's-funnel that centering otherwise exposes.  On a G=2000 simulated genome this
  lifted dEta ESS/sec ~1.6x over the mean-anchored/centered baseline with the funnel
  gone and no regression on dM/log_phi.  Correctness/regression guards added:
  `scripts/test_reparam_stan.R` (likelihood invariance to `mphi`, verified to ~1e-13)
  and `scripts/test_reparam_inertness_stan.R` (legacy `mean(phi)=1` parity + frozen
  golden).  See the `phi-parameterization` vignette for guidance.

- **Composable phi prior specification API for the ROC single-LN phi prior
  (closes #28).** `initializeParameterObject()` accepts two new arguments,
  `phi.mphi` and `phi.sphi`, built from a set of composable mode and prior
  constructors:

  *Mode constructors* (control how mPhi or sphi is treated):
  - `constrained(statistic, value)` — derived so that the named distributional
    property of phi equals `value`.  `statistic` is one of `"mean"`,
    `"median"`, `"mode"`, `"variance"`, `"sd"`.  Default:
    `constrained("mean", 1)` which restores `mPhi = log(value) - sphi^2/2`,
    identical to the previous default behavior.
  - `fixed(value)` — holds mPhi (or sphi) at a constant.
  - `estimated(prior)` — mPhi/sphi is sampled during MCMC with an optional
    prior distribution.  Now implemented for both `phi.sphi` and `phi.mphi`.
    `phi.mphi = estimated(prior)` adds a Metropolis–Hastings update on mPhi
    each iteration, fixing issue #47 (deta_scale_anchor was a no-op when
    mphi was deterministic).

  *Prior distribution constructors* (used inside `estimated()`):
  - `prior_uniform(low, high)` — default for sphi; `prior_uniform(0, 10)` is
    the new system default, replacing the previous implicit improper flat prior.
  - `prior_normal(mean, sd)` — Normal prior on sphi.
  - `prior_student_t(df, mean, sd)` — Student-t prior (forward-compatible; not
    yet wired into the MCMC acceptance step in this release).
  - `prior_exponential(rate)` — Exponential prior (forward-compatible).
  - `NULL` — improper flat; retains backward compatibility with callers that
    passed no prior via `parameter$setSphiPrior()`.

  *Backward compatibility:* All previously valid calling patterns continue to
  work unchanged.  `sphi = 1` (numeric) maps to `phi.sphi = fixed(1)`;
  `sphi = NA` maps to `phi.sphi = estimated(prior = prior_uniform(0, 10))`;
  `phi.prior.constraint = "median"` maps to
  `phi.mphi = constrained("median", 1)` with a deprecation warning.

  *C++ internals:* `Parameter` gains `computeMPhi(sphi, statistic, value)` (a
  static helper with closed-form mPhi formulas for all five statistics),
  `phiMuMode` / `phiMuFixed` / `phiConstraintValue` members, and
  `sphiPriorType` / `sphiPriorLow` / `sphiPriorHigh` members.  All new
  members are persisted in `.rst` restart files (forward-compatible: old builds
  silently skip unknown keys).

- **Stan backend: `initializeStan()` with phi.mphi/phi.sphi spec.**
  `genomeToStanData()` + `genomeToStanInit()` are replaced by a single
  `initializeStan(genome, phi.mphi = ..., phi.sphi = ...)` that accepts
  the same `constrained()` / `fixed()` / `estimated()` / `prior_*()` spec
  objects as `initializeParameterObject()`.  Returns a two-slot list:
  `$data` (pass to `mod$sample(data = ...)`) and `$init` (pass to
  `mod$sample(init = list(...))`).  The matching Stan model
  `stan/roc_sphi_est.stan` is updated with:
  - Data-dependent bounds `real<lower=sphi_low, upper=sphi_high> sphi`
    replacing `real<lower=0> sphi`; a `Uniform(0, 10)` prior sets
    `sphi_low=0, sphi_high=10` with no prior density statement, matching
    the MCMC uniform convention exactly.
  - `phi_mphi_mode` / `phi_mphi_statistic` / `phi_mphi_value` / `phi_mphi_fixed`
    data fields replace the old `anchor_phi` / `mphi_param` / `mphi_prior_sd`
    flags.  All five `computeMPhi` formulas (mean, median, mode, variance, sd)
    are in the `transformed parameters` block.
  - `mphi_param` removed from the `parameters` block; mphi is always computed
    deterministically from the phi spec.
  `adviToWarmStart()` is updated to handle models with or without `mphi_param`.

## BUG FIXES

- **`Parameter` C++ ctor default for `sphiPriorType` changed from `SPHI_PRIOR_FLAT` (0)
  to `SPHI_PRIOR_UNIFORM` (0, 10) (closes #52, belt-and-suspenders).** The R
  initialization path (`initializeParameterObject()` → `.applyPhiSpec()`) always
  overwrites this default, so the change is only observable when constructing a
  `Parameter` object from C++ without going through the R path. Bounds `[0, 10]`
  are already set in the ctor.

- **MCMC RNG determinism (behavior change).** `runMCMC()` is now
  bit-reproducible at `ncores=1` under `set.seed()` for ROC, PA, and
  PANSE.  Two sources of non-determinism were removed from
  `MCMCAlgorithm.cpp`: (1) the codon-specific parameter update scan
  order was shuffled using a wall-clock seed (`std::chrono`), so every
  run used a different scan order; (2) per-gene Metropolis acceptance
  thresholds and mixture-assignment draws came from a
  `std::default_random_engine` that was default-seeded and shared across
  OpenMP threads (data race).  Both now draw from R's RNG via
  `Parameter::randUnif` / `Parameter::randExp`, which is controlled by
  `set.seed()`.  **Downstream impact:** the RNG stream shifts relative
  to all prior versions; any hardcoded numeric expectations or restart
  files calibrated against older code should be re-baselined.
  `ncores > 1` remains nondeterministic (OpenMP FP reduction order).

- **ROC memory leak.** `ROCModel::CalculateProbabilitiesForCodons` heap-
  allocated a `double[]` buffer that was copied into the return vector
  but never freed.  Fixed by adding `delete[] codonProb` before the
  return statement.

- `Parameter::operator=` now copies `lastIteration` and the five
  `restartFile*` build-info fields.  These were silently omitted by
  the explicit assignment operator; the implicit copy ctor handled
  them until the new `unique_ptr<CSPAdaptationStrategy>` member made
  the implicit copy ctor ill-formed.  Surfaced as a segfault in the
  evaluate pipeline (`ROCModel::getParameter()` returns by value;
  the copy's uninitialized `lastIteration=0` caused `traceLength=1`
  and out-of-bounds reads in `getEstimatedMixtureAssignmentForGene`).
  The new explicit copy ctor delegates to the default ctor before
  running `operator=` to safely initialize all primitives.

- **`sphi=<numeric>` regression: numeric init no longer freezes sphi (regression
  introduced in PR #37).** Since PR #37, passing a numeric `sphi` value to
  `initializeParameterObject()` silently resolved to `phi.sphi = fixed(value =
  mean(sphi))`, which called `fixSphi()` and set `fix_stdDevSynthesis = TRUE`.
  This made `proposeStdDevSynthesisRate()` a no-op, freezing sphi even when
  `est.hyper = TRUE`.  All package vignettes (`anacoda.Rmd`, `pa-panse.Rmd`,
  `fonse.Rmd`) pass a numeric sphi and were silently affected.
  Restored canonical pre-PR#37 behavior: a numeric `sphi=` is the per-mixture
  *initial value*; sphi is estimated by default with a Uniform(0, 10) prior.
  To freeze sphi, use `phi.sphi = fixed(value)` (modern) or `est.hyper = FALSE`
  (legacy loop toggle); both still work.

## REFACTORING

- `Parameter` member fields `adaptiveStepPrev` / `adaptiveStepCurr`
  renamed to `adaptiveSamplePrev` / `adaptiveSampleCurr` (they hold
  thinned-sample indices, not raw MCMC steps).  Parameter
  `lastIteration` on `adaptCodonSpecificParameterProposalWidth`
  renamed to `lastSample` (likewise).  Local `samples` inside that
  function renamed to `samplesSinceLastAdapt` to avoid shadowing the
  global `samples` concept.  Pure rename: bit-identical output.

# CHANGES IN AnaCoDa 0.1.2

## BUG FIXES
- fixed a bug were the scaling of observed phi values was used inconsitently, causing problems with estimates of Aphi and Sepsilon

## NEW FEATURES
- Added SCUO calculation and improved getCSPEstimates to include reference codons

# CHANGES IN AnaCoDa 0.1.1

## BUG FIXES
- fixed problem with getCSPEstimates where log scaling was falsely enabled

- fixed problem where the grouplist was not stored by writeParameterObject

## NEW FEATURES
- Added functions to calculate the Codon Adaptation Index, Effective Number of Codons and selection coefficients.

- Allow to set initial phi values based on observed phi values stored in genome object.


