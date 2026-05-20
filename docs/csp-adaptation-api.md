# CSP adaptive proposal-width schemes: R API design

Status: design (2026-05-19). Implementation tracked in tasks #5-#9
of the Lokiarchaeota project session.

## Scope

This document defines the **R-level user contract** for pluggable
adaptive proposal-width schemes used during the codon-specific
parameter (CSP) update step of an MCMC fit. The C++ side serves the
R API; it is described only at the Rcpp seam (where it crosses into
R) and at the validation seam (where it must protect itself against
non-R callers).

In scope: the CSP proposal-width tuning step (the existing
`adaptCodonSpecificParameterProposalWidth` path in `Parameter.cpp`).
Out of scope (for this iteration): alternative proposal *mechanisms*
(HMC, MALA), and adapter abstractions for the hyperparameter and
synthesis-rate update steps.

## Goals

1. Multiple proposal-width-adaptation schemes selectable per fit, by
   name from R or from YAML, with documented per-scheme params.
2. The existing in-house scheme (renamed `native`) remains the default
   and is bit-identical when no scheme is explicitly specified.
3. Validation happens at both the R layer (clean error messages for
   interactive R users) and the C++ layer (so a standalone non-R
   caller still fails loudly with a clear `std::invalid_argument`).
4. Adding a new scheme requires touching one R file, one C++ class,
   one factory entry, and one set of tests -- no surgery on `Parameter`
   or other schemes.

Non-goals: per-AA scheme selection (single scheme per fit), within-run
scheme changes, R-defined schemes (callbacks into R from the MCMC
inner loop would tank performance).

## Background terms

- **CSP** = codon-specific parameters: per-codon mutation bias (`dM`)
  and selection efficacy (`dEta`).
- **Adapt fire** = one invocation of the proposal-width-tuning step.
  In AnaCoDa this happens every `adaptive.width` raw iterations
  (default `round(100 / thinning)`, so every 10 raw iters at
  `thinning=10`).
- **`adaptive.ratio`** = continuous fraction in `[0, 1]` controlling
  *what fraction of a round's iterations* the chosen scheme is allowed
  to adapt before being frozen for the rest of the round. Orthogonal
  to scheme choice. See `analysis/.../SCHEMA.md` for full semantics.
- **Scheme** = the algorithm that decides how proposal widths change
  on an adapt fire. Picked once per fit.

## Initial scheme set

Two schemes in v1:

| Name              | Source           | Params (defaults)                                |
|-------------------|------------------|--------------------------------------------------|
| `native`          | in-house (Cedric, Mike); shipped with AnaCoDa | none |
| `andrieu_thoms`   | Andrieu & Thoms 2008, *Statistics and Computing*, Algorithm 4 | `target=0.234, alpha=0.7, c=1.0, t0=10` |

### native (algorithmic summary, for reference only)

On each adapt fire, per AA:
- compute acceptance rate over the past `adaptive.width` raw iterations
- if acceptance < 0.225: multiply per-codon `std_csp` by 0.8, AND
  blend the per-AA proposal covariance toward a sample covariance
  computed from that adapt window's trace rows
- if acceptance > 0.325: multiply per-codon `std_csp` by 1.2
- otherwise: no change

Has been used in all RibModelFramework fits to date; reaches
stationarity in practice on Lokiarchaeota fits.

### andrieu_thoms (algorithmic summary)

Continuous Robbins-Monro update on `log(std_csp)`. Per AA:
```
log(std_csp[k])_{t+1} = log(std_csp[k])_t + gamma_t * (accept_t - target)
gamma_t                = c / (t + t0)^alpha
```
- `t` = adapt-fire count for that AA (per-AA counter, increments by 1
  each fire; starts at 0)
- `target` = single target acceptance rate (Gelman et al's
  high-dimensional optimum is 0.234 for d > ~5; default 0.234)
- `alpha in (0.5, 1]` -- diminishing-adaptation condition (Andrieu &
  Thoms 2008, Theorem 2)
- `c > 0` -- initial step magnitude
- `t0 >= 0` -- offset to avoid huge early steps

Proposal covariance structure is **not** updated (stays at its
initialization value). Only the scalar `std_csp[k]` (per codon)
moves.

## R API

### Scheme constructors

One constructor per scheme. Pattern: `AdaptiveScheme.<PascalName>(...)`.
Each constructor:
- validates its arguments (layer 2)
- returns an S3 object of class
  `c("AdaptiveScheme.<PascalName>", "AdaptiveScheme")`

```r
AdaptiveScheme.Native <- function() {
    structure(
        list(
            scheme = "native",
            params = list()
        ),
        class = c("AdaptiveScheme.Native", "AdaptiveScheme")
    )
}

AdaptiveScheme.AndrieuThoms <- function(target = 0.234,
                                        alpha  = 0.7,
                                        c      = 1.0,
                                        t0     = 10) {
    stopifnot(
        is.numeric(target), length(target) == 1L, is.finite(target),
        target > 0, target < 1,
        is.numeric(alpha),  length(alpha)  == 1L, is.finite(alpha),
        alpha > 0.5, alpha <= 1.0,
        is.numeric(c),      length(c)      == 1L, is.finite(c),
        c > 0,
        is.numeric(t0),     length(t0)     == 1L, is.finite(t0),
        t0 >= 0
    )
    structure(
        list(
            scheme = "andrieu_thoms",
            params = list(target = target, alpha = alpha, c = c, t0 = t0)
        ),
        class = c("AdaptiveScheme.AndrieuThoms", "AdaptiveScheme")
    )
}
```

Validation predicates (layer 2):

| Param   | Constraint            | Reason                                          |
|---------|-----------------------|-------------------------------------------------|
| target  | `0 < target < 1`      | acceptance rate                                 |
| alpha   | `0.5 < alpha <= 1.0`  | diminishing-adaptation theorem (A&T Theorem 2)  |
| c       | `c > 0`               | positive step size                              |
| t0      | `t0 >= 0`             | offset is non-negative                          |

All numeric params also require `length == 1`, `is.finite`, `is.numeric`
checks. Same predicates are re-implemented in C++ (layer 4).

### Discovery and inspection

```r
schemes.available()
# c("native", "andrieu_thoms")

is.AdaptiveScheme(x)
# TRUE / FALSE; tests inherits(x, "AdaptiveScheme")

print(at)
# AdaptiveScheme: andrieu_thoms
#   target = 0.234, alpha = 0.7, c = 1, t0 = 10

format(at)
# returns the character string used by print
```

`schemes.available()` returns the canonical names (lowercase
snake_case) suitable for YAML. The corresponding R constructor for
name `"foo_bar"` is `AdaptiveScheme.FooBar()` (snake -> Pascal).
A helper `.scheme.name.to.constructor()` is provided for the YAML
loader to do this mapping.

### Use in fits

```r
fitSingleRoundOfROC(
    ...,
    adaptive.scheme = NULL,   # NULL -> AdaptiveScheme.Native(), bit-identical
    ...
)
```

`adaptive.scheme` is a fit-level argument, not a round-level one. It
appears once in the function signature; all rounds in a fit share it.

Inside `fitSingleRoundOfROC`:
```r
if (is.null(adaptive.scheme)) adaptive.scheme <- AdaptiveScheme.Native()
stopifnot(is.AdaptiveScheme(adaptive.scheme))
parameter$setCSPAdaptationScheme(
    adaptive.scheme$scheme,
    adaptive.scheme$params
)
```

`adaptive.ratio` continues to live per round in `settings.tbl` (v.2)
or `rounds:` (v.3) and is unchanged by this work.

### Post-fit diagnostics

```r
adaptive.scheme.diagnostics(parameter, mcmc)
# returns a list:
# $scheme.name      : character(1)        -- the scheme used for this fit
# $params           : list                -- params at fit time
# $std_csp.trace    : array               -- (codon x adapt-call x round)
# $acceptance.trace : array               -- (AA    x adapt-call x round)
# $scheme.specific  : list                -- scheme-specific extras
```

For `native`: `$scheme.specific$which.branch` is a per-AA-per-call
character series in {`"shrink"`, `"none"`, `"grow"`}.
For `andrieu_thoms`: `$scheme.specific$gamma.trace` is a per-AA-per-call
double series of step sizes.

The `$std_csp.trace` and `$acceptance.trace` come from existing
`Trace` machinery in AnaCoDa; the `$scheme.specific` slot is new and
populated by each strategy via a `getDiagnostics()` virtual method
on the C++ side.

## YAML schema (v.3)

New optional block under `fit`:

```yaml
fit:
  ## ... existing fit fields ...
  csp.adaptation:                    # optional; missing -> native scheme
    scheme: andrieu_thoms            # one of schemes.available()
    target: 0.234                    # scheme-specific param
    alpha:  0.7
    c:      1.0
    t0:     10
```

Loader behavior (`lib/config.R`, layer 1):
- If the `csp.adaptation` block is missing or null: no R object is
  constructed; `fitSingleRoundOfROC` is called with
  `adaptive.scheme = NULL`, which the runner translates to
  `AdaptiveScheme.Native()`.
- If present, `scheme` is required and must be in
  `schemes.available()`. All other keys are scheme-specific params.
- Unknown keys (typos) inside `csp.adaptation` are rejected by the
  same unknown-key validator already used for the top-level config.
- The loader does NOT do per-param range checks. It does check that
  the keys present are the keys expected for the named scheme; it
  then constructs the R `AdaptiveScheme.*` object via the
  scheme-name-to-constructor mapping, which performs the range
  validation (layer 2). Failed validation raises an R error pointing
  to the offending YAML field.

Backward compatibility: existing YAMLs that don't set
`csp.adaptation` get the `native` scheme with the existing in-house
0.8/1.2 logic, bit-identical to today.

## C++ R-side surface (Rcpp seam)

Single new method on `Parameter`:

```cpp
void setCSPAdaptationScheme(const std::string& name,
                            const Rcpp::List& params);
```

Implementation:
- Convert `params` (named numeric list) to `std::map<std::string, double>`
- Call `makeCSPAdapter(name, params_map)` (layer 3, layer 4)
- Replace this Parameter's `cspAdapter` with the returned strategy

Errors are raised via `Rcpp::stop()` (which is caught and re-thrown
as an R error) so any validation failure surfaces cleanly to the R
user. The factory function `makeCSPAdapter` is also callable
directly from C++ standalone code -- it uses `std::invalid_argument`
for non-R callers, and `setCSPAdaptationScheme` translates that to
`Rcpp::stop` at the Rcpp seam.

## Tests (R)

In `tests/testthat/test-adaptive-scheme.R`:
- Constructor accepts valid params
- Constructor rejects invalid params (one assertion per bound, both
  edges)
- Length-not-1 inputs rejected
- Non-finite (NA, NaN, Inf) inputs rejected
- `schemes.available()` returns at least `c("native", "andrieu_thoms")`
- `is.AdaptiveScheme` recognises both scheme classes
- `print` does not error

In `tests/testthat/test-adaptive-scheme-csp.R` (slower; tagged `@slow`):
- Bit-identical regression: a short B21-or-similar fit with
  `AdaptiveScheme.Native()` produces identical `std_csp` trace to the
  same fit without the `adaptive.scheme` arg (legacy code path)
- `AdaptiveScheme.AndrieuThoms()` runs to completion on the same
  fit with no errors; `std_csp` moves; `gamma` decreases monotonically
- `adaptive.scheme.diagnostics()` returns the documented shape for
  both schemes

C++ tests in `src/Testing.cpp`:
- `makeCSPAdapter("native", {})` returns a non-null strategy
- `makeCSPAdapter("andrieu_thoms", {target: 0.234, ...})` returns
  a non-null strategy
- `makeCSPAdapter("does_not_exist", {})` throws
- `makeCSPAdapter("andrieu_thoms", {target: 2.0, ...})` throws
- Each scheme's `update()` is callable and modifies `std_csp` only
  via the documented mechanism

## Help pages (roxygen)

Required:
- `?AdaptiveScheme.Native` -- algorithm + history + no params
- `?AdaptiveScheme.AndrieuThoms` -- algorithm + math (above) + param
  meanings and bounds + paper citation
- `?schemes.available`
- `?adaptive.scheme.diagnostics` -- output shape, per-scheme
  `$scheme.specific` documentation

## File layout (proposed)

R side:
- `R/adaptiveScheme.R` -- constructors, `schemes.available`,
  `is.AdaptiveScheme`, print/format S3 methods, name-to-constructor
  mapping, `adaptive.scheme.diagnostics`

C++ side (under `src/`):
- `CSPAdaptationStrategy.h` -- abstract base
- `NativeCSPAdapter.h` / `.cpp` -- existing logic extracted
- `AndrieuThomsCSPAdapter.h` / `.cpp` -- new
- `CSPAdaptationFactory.h` / `.cpp` -- `makeCSPAdapter` factory
- `Parameter.cpp` -- delegates `adaptCodonSpecificParameterProposalWidth`
  to the strategy; adds `setCSPAdaptationScheme` Rcpp method
- `Parameter.h` -- holds `std::unique_ptr<CSPAdaptationStrategy>`,
  default-constructed to `NativeCSPAdapter`

Tests:
- `tests/testthat/test-adaptive-scheme.R` -- R-side unit tests
- `tests/testthat/test-adaptive-scheme-csp.R` -- end-to-end regression
- `src/Testing.cpp` -- C++ unit additions for factory and validators

## Open questions

None blocking implementation. Possible follow-ups:
- Should `adaptive.scheme.diagnostics()` also return the timing
  / wall-clock cost of each adapt fire? Useful for benchmarking
  schemes against each other.
- Per-AA scheme selection (`scheme: andrieu_thoms` globally with
  `scheme.overrides: {ARG: native, MET: native}` -- speculative,
  not required for v1).
- Adapter state serialization in restart files: not in v1; each
  round starts with a fresh adapter state. Worth revisiting if the
  cross-round scheme reset is observably worse than continuity.

## See also

- Lokiarchaeota project memory: `pluggable-csp-adapter`
- `analysis/02v.3_fit.Sarina.More.Complete.Genomes/SCHEMA.md` --
  v.3 YAML schema; will gain a `fit.csp.adaptation` section when
  task #8 lands.
- Andrieu, C., & Thoms, J. (2008). A tutorial on adaptive MCMC.
  *Statistics and Computing*, 18(4), 343-373.
  [https://doi.org/10.1007/s11222-008-9110-y](https://doi.org/10.1007/s11222-008-9110-y)
