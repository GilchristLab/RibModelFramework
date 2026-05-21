# CHANGES IN AnaCoDa 0.1.10 (unreleased; in-progress)

## NEW FEATURES

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

## BUG FIXES

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


