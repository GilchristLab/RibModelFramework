# How to choose a CSP adaptive proposal-width scheme

This is a short user-facing guide.  See `csp-adaptation-api.md` for the
full API design and `NEWS.md` for the release-level changelog.

## When to read this

You're running ROC (or FONSE, etc.) MCMC fits with AnaCoDa and want to
either (a) get the historical default behavior (= do nothing -- it's the
default) or (b) try the new Andrieu-Thoms scheme as an alternative
proposal-width tuner.

## Which scheme should I pick?

Two ship in this build:

| Scheme | When to use it |
|--------|----------------|
| **`native`** (default) | The in-house scheme that produced all historical fits.  Reaches stationarity in practice.  No tuning required; takes no parameters. |
| **`andrieu_thoms`** | An alternative from Andrieu & Thoms 2008 (*Statistics and Computing* 18:343-373) Algorithm 4.  Continuous Robbins-Monro update on `log(std_csp)` with a diminishing step schedule.  Mathematically principled (satisfies the diminishing-adaptation theorem for `alpha in (0.5, 1.0]`).  Try this if you're investigating mixing behavior on a chain that's slow to converge under native -- particularly if proposal widths are drifting unstably or acceptance rates aren't settling. |

In short: if you don't have a specific reason to switch, stay on
`native`.  If you do, try `andrieu_thoms` with the literature defaults
first (see below).

## How to select a scheme

### From R, interactively

```r
library(AnaCoDa)

# Use the default native scheme (no setup needed)
parameter <- initializeParameterObject(genome = ..., sphi = 1, ...)
# parameter$getCSPAdaptationSchemeName()  # -> "native"

# Switch to Andrieu-Thoms with literature defaults
at <- AdaptiveScheme.AndrieuThoms()        # target=0.234, alpha=0.7, c=1.0, t0=10
parameter$setCSPAdaptationScheme(at$scheme, at$params)

# Verify
parameter$getCSPAdaptationSchemeName()      # -> "andrieu_thoms"
```

### From R, with custom A-T parameters

The `target` defaults to 0.234 (Gelman et al's optimum-d optimum).  If
you have a reason to target a different acceptance rate -- e.g. a
low-dimensional model where a higher target is appropriate -- pass it:

```r
at <- AdaptiveScheme.AndrieuThoms(target = 0.30, alpha = 0.7)
parameter$setCSPAdaptationScheme(at$scheme, at$params)
```

Allowed ranges:

| Param   | Range                | Default | Meaning |
|---------|----------------------|---------|---------|
| target  | `(0, 1)`             | 0.234   | Target acceptance rate. |
| alpha   | `(0.5, 1.0]`         | 0.7     | Step-size decay exponent.  Smaller = adapt longer; larger = freeze sooner. |
| c       | `> 0`                | 1.0     | Initial step magnitude in the SD update. |
| t0      | `>= 0`               | 10      | Schedule offset (`gamma_t = c / (t + t0)^alpha`); larger values shrink early-iteration step magnitude. |

Out-of-range values raise an R error at the constructor; the C++ side
also revalidates so standalone (non-R) callers fail loudly too.

### From a YAML config (Lokiarchaeota v.3 pipeline)

Add an optional `csp.adaptation` block under `fit:` in your run YAML:

```yaml
fit:
  sphi: 1
  dM.prior.method: gcBias
  dM.prior.sd: 0.8
  # ... other fit fields ...
  csp.adaptation:
    scheme: andrieu_thoms
    target: 0.234
    alpha:  0.7
    c:      1.0
    t0:     10
```

Omitting the block is equivalent to `scheme: native` and gives you the
historical default.

## How adaptation runs across rounds

The scheme is fit-level (one choice per fit, applied to all rounds).
Within each round, whether the scheme actually *runs* is gated by that
round's `adaptive.ratio` field:

- `adaptive.ratio: 1.0` -- scheme adapts the entire round
- `adaptive.ratio: 0.5` -- scheme adapts the first 50% of iterations
  in that round, then freezes for the rest
- `adaptive.ratio: 0.0` -- scheme is frozen for the whole round
  (proposal widths are whatever they were at the end of the previous
  round)

So the typical "adapt early, sample late" pattern in a 2-round YAML is:

```yaml
rounds:
  - { samples: 5000, adaptive.ratio: 1.0 }   # scheme adapts here
  - { samples: 5000, adaptive.ratio: 0.0 }   # scheme frozen, sample
```

Both rounds use the same scheme; only the per-round "ratio" controls
when it acts.

## How to inspect what's in effect

```r
parameter$getCSPAdaptationSchemeName()
# -> "native"  or  "andrieu_thoms"

adaptive.scheme.diagnostics(parameter)
# list with $scheme.name and (in future versions) traces / params
```

## What does NOT change

- Posterior summaries, trace files, restart files: format unchanged.
- The acceptance-rate diagnostic the MCMC prints every adapt fire:
  format unchanged across schemes.
- Random-number streams: native is bit-identical to pre-2026-05-20
  builds when you don't set a scheme.  Andrieu-Thoms uses the same
  RNG draws -- only the acceptance-criterion-conditioned scaling of
  `std_csp` differs.
