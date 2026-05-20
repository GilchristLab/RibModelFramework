# CSP adaptive-scheme test & doc coverage

## Tests

### R-level constructor / validation (tests/testthat/test-adaptive-scheme.R, 12 tests)

- `AdaptiveScheme.Native()` returns correct S3 class + empty params
- `AdaptiveScheme.AndrieuThoms()` defaults: target=0.234, alpha=0.7, c=1.0, t0=10
- Range rejection for each AT param:
  - `target`: rejects 0.0, 1.0, -0.1, 1.5, NA, NaN, Inf, length-2, character
  - `alpha`: rejects 0.5 (closed bound), 1.01, 0, NA, length-2
  - `c`: rejects 0, -0.5, NA
  - `t0`: rejects -1, NA, length-2
- Boundary cases that DO pass: alpha=1.0, t0=0, target=0.001 alpha=0.501
- `is.AdaptiveScheme()` returns TRUE for both constructed objects, FALSE for `list()`, NULL, numeric, character
- `print()` produces output matching expected text; `format()` returns char vector of expected shape
- `schemes.available()` returns `c("native", "andrieu_thoms")`
- `.scheme.name.to.constructor()` maps snake_case names to PascalCase constructors; errors on unknown names, empty string, NA, length-2

### Rcpp seam (test-adaptive-scheme.R, 11 tests)

- Default Parameter has `native` scheme bound
- `setCSPAdaptationScheme("andrieu_thoms", list(...))` switches the in-memory scheme; `getCSPAdaptationSchemeName()` reflects the change
- `setCSPAdaptationScheme("native", list())` switches back
- Unknown scheme name -> R error containing "unknown"
- Out-of-range params surface as R errors with the bound name in the message:
  - target=2.0 -> "target"
  - alpha=0.5 -> "alpha"
  - c=0.0 -> "c"
  - t0=-1 -> "t0"
- Extra param for native -> R error
- Missing required param for andrieu_thoms -> R error citing the missing key by name
- Extra param for andrieu_thoms -> R error
- Scheme name survives Parameter copy / operator=

**Total: 74 tests, all passing as of `at-integration-v1`.**

### Integration / end-to-end

- The 02v.4-AT.adapter sweep IS the end-to-end integration test: 21
  single-genome fits at 10000 stored samples each, all under
  `andrieu_thoms`.  As of this writing the sweep is running with all
  completed fits at rc=0.  Acceptance-rate patterns at end of round 2
  match the native sweep on the same genomes (same low-acceptance AAs
  on the same data), suggesting A-T and native converge to similar
  end states on this codebase.

## Tests we deliberately did NOT add

- **C++ Testing.cpp unit tests for `makeCSPAdapter`.**  The R-level
  Rcpp-seam tests exercise the same code paths via Rcpp::stop()
  translation of std::invalid_argument; adding a duplicate C++-only
  test would be redundant.

- **A "10-iter MCMC" smoke test in testthat.**  Setting up a fixture
  (genome, parameter, model, mcmc objects) for a runnable 10-iter
  test is a non-trivial amount of brittle code, and the 02v.4 sweep
  exercises the full path at scale.

- **AndrieuThomsCSPAdapter::update() unit test calling it directly
  with a synthetic CSPAdaptContext.**  Requires constructing a Trace
  object and CovarianceMatrix from scratch; the R-side
  setCSPAdaptationScheme + a real MCMC indirectly tests the same
  update logic and produces meaningful chains.

## Documentation

| File                                          | Audience          | Status |
|-----------------------------------------------|-------------------|--------|
| `docs/csp-adaptation-api.md`                  | API designers     | Done   |
| `docs/csp-adaptation-howto.md`                | End users         | Done   |
| `docs/csp-adaptation-test-coverage.md`        | Maintainers       | (this file) |
| `man/AdaptiveScheme.Native.Rd`                | R help (`?`)      | Done (roxygen-generated) |
| `man/AdaptiveScheme.AndrieuThoms.Rd`          | R help (`?`)      | Done |
| `man/schemes.available.Rd`                    | R help (`?`)      | Done |
| `man/is.AdaptiveScheme.Rd`                    | R help (`?`)      | Done |
| `man/adaptive.scheme.diagnostics.Rd`          | R help (`?`)      | Done |
| `man/print.AdaptiveScheme.Rd`                 | R help (`?`)      | Done |
| `man/format.AdaptiveScheme.Rd`                | R help (`?`)      | Done |
| `NEWS.md`                                     | Users / releases  | Done (0.1.10 unreleased) |
| Lokiarchaeota `02v.3/SCHEMA.md`               | YAML authors      | Done (fit.csp.adaptation block) |

## Documentation we deliberately did NOT add

- **Vignette section in `vignettes/anacoda.Rmd`.**  Existing vignettes
  use small embedded examples that run during R CMD INSTALL; adding
  a vignette section for adaptive-scheme selection would require a
  runnable MCMC example and would slow down the vignette build.
  The user-facing how-to (`docs/csp-adaptation-howto.md`) covers the
  same content without that build-time cost.

- **R Rcpp method roxygen for `setCSPAdaptationScheme` /
  `getCSPAdaptationSchemeName`.**  Rcpp module methods aren't picked
  up by roxygen2 automatically; they would need hand-written .Rd
  files.  The methods are documented in the design doc + how-to doc;
  R users mostly access them indirectly via `AdaptiveScheme.*`
  constructors which DO have roxygen help pages.

## Project memory (Lokiarchaeota)

- `project_pluggable_csp_adapter.md` -- design decisions
- `project_at_integration_2026_05_20.md` -- what landed + bugs found
- `feedback_mcmc_clean_slate.md` -- user feedback rule (clean-slate, not patch)
- `feedback_evidence_based_tweaks.md` -- user wants A/B tests for tweaks
- `reference_adaptive_ratio_semantics.md` -- the orthogonal [0,1] axis
