# Autonomous task: fix `sphi=<numeric>` silent-freeze regression

Branch: `fix/sphi-numeric-estimated-regression` (off origin/main @ 1ee783a)
Worktree: `~/Repositories/RibModelFramework-sphi-numeric-fix`
Author: Claude (autonomous run, 2026-06-06), delegated by M. Gilchrist who is away.

NOTE: This file is a scratch/decision log. Delete it (or move to a PR comment)
before the final PR if it shouldn't live in the repo. It is committed so the
reasoning survives context compaction.

---

## Working protocol (per user)

- **Commit early and often.** Each logical step = its own commit. Do not batch.
- **Slow and careful > fast and sloppy.** Verify each step before moving on.
- Stop at an **open PR** against GilchristLab/main. **Do NOT merge / force-push / delete.**
- Log every judgment call here as it happens.
- If a genuine ambiguity can't be resolved safely, STOP and document it here; do not guess.
- Run slow builds/tests with `run_in_background: true`; the session re-invokes on completion.

## The bug (confirmed empirically on main, 2026-06-06)

`initializeParameterObject(genome, sphi = <numeric>)` resolves (since PR #37,
commit bfa682c) to `phi.sphi <- fixed(value = mean(sphi))`, which calls
`parameter$fixSphi()` -> `fix_stdDevSynthesis = true`. That makes
`proposeStdDevSynthesisRate()` propose sphi UNCHANGED, so sphi is frozen even
when `est.hyper = TRUE`. The resolution runs before model dispatch, so it hits
ROC, FONSE, PA, PANSE alike. All package vignettes pass numeric sphi, so they
all now silently freeze sphi.

Empirical proof (ROC, est.hyper=TRUE, 200 samples, test genome):
```
sphi=1  (numeric)    n_unique=2    max=1.0    <- FROZEN
sphi=NA (estimated)  n_unique=160  max=2.1    <- estimated/moving
```

Canonical pre-PR#37 behavior: `sphi` = initial value only; estimation toggled
by `est.hyper`; `fixSphi()` was never called from the R path; `is.na(sphi)`
never existed. So numeric sphi was historically ESTIMATED (seeded at the value).

## Locked design (no further user input needed)

1. `sphi=<numeric>` -> `phi.sphi <- estimated(prior_uniform(0,10))`, initial
   value = the numeric (per mixture). Does NOT call fixSphi(). sphi is estimated.
2. `sphi=NA` -> unchanged: estimated(prior_uniform(0,10)), init 1.0.
3. Fixing sphi remains available via `phi.sphi=fixed(value)` (modern) or
   `est.hyper=FALSE` (legacy loop toggle). Both must still work.
4. Prior = Uniform(0,10) for the numeric case (matches NA branch + the
   apples-to-apples MCMC/Stan parity preference; NOT the old unbounded-flat).
5. Fix roxygen for `sphi`: it is the initial value, estimated by default; point
   to phi.sphi=fixed() / est.hyper=FALSE for fixing.
6. Keep the change minimal: only the resolution block in
   `R/parameterObject.R` (the `else { phi.sphi <- fixed(...) }` branch). The
   `sphi` init plumbing downstream is untouched.

## Implementation plan (commit after each)

- [x] C1: edit `R/parameterObject.R` resolution block (numeric -> estimated + init).
- [x] C2: update roxygen for `sphi` (+ regenerate man/ if doc build is set up).
- [x] C3: regression tests (new file `tests/testthat/testSphiNumericEstimated.R`):
      (a) sphi=1 + est.hyper=TRUE -> sphi trace moves (n_unique large);
      (b) phi.sphi=fixed(1) -> sphi frozen (n_unique==1 among samples);
      (c) sphi=NA -> moves (baseline);
      (d) smoke: all four models initialize without error via initializeParameterObject.
- [x] C4: empirical multi-model verification script + record results in decision log.
      ROC/FONSE: confirmed moving (n_unique~85, range>0.8, prior_type=2).
      PA/PANSE: init smoke OK (mini_rfp.csv, prior_type=2). MCMC not run (too small).
- [x] C5: NEWS.md entry (regression fix; restores estimation for numeric sphi).
- [x] C6: build + full relevant test run; all green.
      Tested: testSphiNumericEstimated.R, testPhiPrior.R, testPhiMixturePrior.R,
      testROCLogLikelihood.R, testROCNumerical.R, testPAModelBasic.R, testMCMCPANSE.R.
- [ ] C7: also Issue #52 — SEPARATE branch/PR: flip C++ ctor
      `sphiPriorType = SPHI_PRIOR_FLAT` -> `SPHI_PRIOR_UNIFORM` (both ctors,
      src/Parameter.cpp ~lines 103, 185; bounds already [0,10]); update any
      testPhiPrior.R assertion of the flat ctor default; NEWS line. Do this only
      after the sphi PR is green.

## Build notes

- R-only change for the sphi fix: prefer `devtools::load_all()` or
  `R CMD INSTALL --no-docs` to avoid full C++ recompile when iterating on R/tests.
- The worktree starts with stale build artifacts inherited from the checkout;
  a clean `R CMD INSTALL` may be needed once before tests load the new R code.
- #52 touches C++ -> requires a full recompile + install.

## Open risks / things to watch

- PA/PANSE empirical check needs RFP-count genome data, not FASTA. Find package
  test data (inst/extdata, tests/testthat data); if none usable, document and
  limit empirical proof to ROC+FONSE, note PA/PANSE covered by init smoke test.
- Confirm no existing test asserts numeric->fixed (earlier grep:
  testPhiPrior.R:429 tests `phi.sphi=fixed(1)` explicitly — that's the EXPLICIT
  fixed() path and should still pass; it does NOT assert numeric->fixed).
- Make sure `mean(sphi)` semantics: for multi-mixture numeric sphi, init should
  be the per-mixture vector, not collapsed to the mean. Verify the init plumbing
  passes the full vector (the old fixed() used mean(); estimated init should keep
  the vector).

## Decision log (append as work proceeds)

- 2026-06-06: worktree + branch created off origin/main 1ee783a. Spec locked.
- 2026-06-07: C1-C4 complete. Empirical multi-model verification (100 samples, set.seed(42)):
  ROC  sphi=1 numeric: n_unique=86, range=0.807, prior_type=2 (UNIFORM) -- MOVING
  ROC  sphi=NA       : n_unique=86, range=0.906, prior_type=2 (UNIFORM) -- MOVING
  FONSE sphi=1 numeric: n_unique=84, range=2.401, prior_type=2 (UNIFORM) -- MOVING
  FONSE sphi=NA       : n_unique=83, range=2.627, prior_type=2 (UNIFORM) -- MOVING
  PA    sphi=1.5 init : OK, prior_type=2 (UNIFORM) -- smoke only, no MCMC
  PANSE sphi=1.5 init : OK, prior_type=2 (UNIFORM) -- smoke only, no MCMC
  (PA/PANSE RFP-count genomes require specific fixture; mini_rfp.csv used for init smoke)
- 2026-06-07: testSphiNumericEstimated.R trace[1]=0 sentinel found; excluded in tests.
- 2026-06-07: testPhiPrior.R getSphiPriorType default assertion updated: now 2 (UNIFORM)
  not 0 (FLAT), because numeric sphi= now resolves through estimated(uniform) not fixed().
