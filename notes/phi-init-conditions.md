# Phi initial conditions for ROC (and the other models)

Dev note on how AnaCoDa sets the starting synthesis-rate (phi) vector, what
the available "init modes" actually do, and why a downstream pipeline may
prefer a deterministic SCUO seed over the built-in C++ default.

Written 2026-06-12. Companion to `vignettes/phi-parameterization.Rmd` (which
covers the phi *scale gauge*: mean(phi)=1 <=> mphi=-sphi^2/2, etc.). This note
is specifically about *initialization*, not the posterior parameterization.

## Where init happens

R: `initializeParameterObject()` -> `initializeROCParameterObject()`
(`R/parameterObject.R`). For ROC the branch is:

```r
if (init.by.random)                                   # only if explicitly asked
  parameter$initializeSynthesisRateByRandom(sphi[1])
else if (is.null(expressionValues) && init.w.obs.phi == FALSE)   # <-- THE DEFAULT
  parameter$initializeSynthesisRateByGenome(genome, mean(sphi))
else if (init.w.obs.phi == TRUE && is.null(expressionValues))
  ... initializeSynthesisRateByList(observed.phi)     # from observed phi set
else if (!is.null(expressionValues))
  parameter$initializeSynthesisRateByList(expressionValues)
```

So with no `initial.expression.values`, no observed phi, and
`init.by.random = FALSE` (all defaults), ROC takes the **`ByGenome`** path.
The same structure holds for FONSE/PA/PANSE except their default branch is
`initializeSynthesisRateByRandom` (no SCUO) -- ROC is the one that seeds from
SCUO by default.

## What each C++ initializer does (`src/Parameter.cpp`)

1. **`InitializeSynthesisRate(Genome&, sd_phi)`  ("ByGenome", the ROC default)**
   - For every gene i: compute `SCUO_i = calculateSCUO(gene_i)` and draw a
     random `expr_i = randLogNorm(-sd_phi^2/2, sd_phi)`.
   - Sort genes by SCUO; independently sort the random `expr` draws.
   - Assign the k-th smallest `expr` draw to the gene with the k-th smallest
     SCUO (rank-matching), and set `std_phi = 0.1`.
   - Net effect: **phi is SCUO-rank-ordered**, with magnitudes given by the
     *order statistics of G random LN(-sd^2/2, sd) draws*. Stochastic -- the
     exact values depend on the RNG state. `sd_phi = mean(sphi)`. Mean-anchored
     (`-sd^2/2`).

2. **`InitializeSynthesisRate(sd_phi)`  ("ByRandom")**
   - phi_i = `randLogNorm(-sd_phi^2/2, sd_phi)`, NO SCUO ordering. This is the
     only genuinely *uninformative* (codon-usage-agnostic) seed. Reached only
     via `init.by.random = TRUE`.

3. **`InitializeSynthesisRate(std::vector<double> expression)`  ("ByList")**
   - phi taken directly from a supplied vector (observed phi, or any
     externally computed seed such as the R SCUO/ENC' helpers below).

## The deterministic R seed (`scuo_to_log_phi`, R/initHelpers.R)

A pipeline that wants a *reproducible* SCUO seed can compute it in R and pass
it through `ByList`:

```r
r <- rank(scuo, ties.method = "average")
q <- r / (G + 1)
log_phi <- qnorm(q, mean = -0.5 * sphi_seed^2, sd = sphi_seed)   # theoretical quantiles
```

This is the *same idea* as ByGenome (SCUO rank order, LN(-s^2/2, s) spread,
mean-anchored) but replaces the noisy Monte-Carlo order statistics with the
smooth inverse-CDF quantiles `qnorm(rank/(G+1), ...)`. Differences vs ByGenome:

| aspect            | C++ ByGenome (`default`)          | R `scuo_to_log_phi` (`scuo`)        |
|-------------------|-----------------------------------|-------------------------------------|
| SCUO rank order   | yes                               | yes                                 |
| magnitude source  | sorted RANDOM LN draws (stochastic, RNG-dependent) | qnorm theoretical quantiles (deterministic) |
| spread parameter  | `sd_phi = mean(sphi)`             | `sphi_seed` (caller's choice)       |
| anchor            | mean, `-sd^2/2`                   | mean, `-sphi_seed^2/2`              |
| reproducible      | no (depends on RNG state)         | yes                                 |
| cross-backend     | C++ only                          | identical to the Stan ROC seed      |

An ENC'-based variant (`encp`) reuses the *same* `scuo_to_log_phi` mapping but
feeds it negated ENC' (Novembre 2002) ranks instead of SCUO ranks.

## Practical guidance / recommendation

- The ROC built-in `default` is **not** an agnostic baseline -- it is a
  *stochastic SCUO-rank* seed. Do not describe `default` runs as
  "uninformative"; only `init.by.random = TRUE` (ByRandom) is that.
- For reproducible, cross-backend-consistent fits prefer the deterministic
  `scuo` (or `encp`) seed via `initial.expression.values` / `ByList`. It is
  RNG-independent, lets you set the seed spread explicitly (`sphi_seed`), and
  matches the Stan ROC initialization exactly.
- Keep `default` (C++ ByGenome) only as a *reference arm* when you explicitly
  want to compare against AnaCoDa's historical init, e.g. an init-condition
  study.
- A true uninformative control requires exposing `init.by.random` (ByRandom);
  it is not reachable through the SCUO/ENC' seed options.

## Pointers

- `R/parameterObject.R` -- `initializeROCParameterObject()` dispatch.
- `src/Parameter.cpp` -- `InitializeSynthesisRate(...)` overloads (ByGenome /
  ByRandom / ByList) and `calculateSCUO()`.
- `R/initHelpers.R` -- `scuo_to_log_phi()` (deterministic SCUO seed).
- `vignettes/phi-parameterization.Rmd` -- the phi *scale gauge* (anchoring).
