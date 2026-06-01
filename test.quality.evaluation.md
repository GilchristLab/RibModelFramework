# Test-Suite Quality Evaluation -- AnaCoDa / RibModelFramework

> **Status:** This is the *baseline* assessment captured before the work on
> branch `feat/...` (the `eval-testcases` line). Several gaps below have since
> been addressed in this branch: the ROC/PA MCMC RNG non-determinism is fixed
> (P1/P4), the ROC `CalculateProbabilitiesForCodons` leak is fixed, the ROC
> log-space likelihood now has tests, and PA has active tests + a modern RFP
> fixture (P2, partial). Remaining open: C++ model-likelihood unit tests beyond
> ROC, PANSE coverage, and the deferred C++ `testPAParameter` revival. Treat the
> findings below as the "before" picture.

- **Scope:** test cases on the main line (`origin/main`, GilchristLab, commit
  `57f4914`). This is local `main` (`3591a24`) plus the merged Stan/arcsine work.
- **Evaluated:** 23 R test files (~4,500 lines) under `tests/testthat/`, the
  3,756-line C++ unit-test driver `src/Testing.cpp`, and the Stan-layer code
  under test `R/stanDataHelpers.R`.
- **Method:** static review of test design, coverage, and assertion strength
  (suite was not executed in this pass).
- **Date:** 2026-05-30  Host: watauga  Worktree: eval-testcases

---

## Headline verdict

**Overall: B- / C+, but strongly bimodal.**

The suite splits cleanly into two populations:

- **Recent R-level tests (2026 work) -- genuinely good (A-/B+).**
  Oracle-based numerical checks, regression tests tied to specific commits/bug
  numbers, documented fixtures, exhaustive boundary validation. This is
  publication-grade test engineering.
- **Legacy C++ unit layer + the main MCMC integration tests -- weak (D/C).**
  Almost entirely data-structure plumbing (getters/setters/IO round-trips);
  the scientific core (model likelihoods, proposals, MCMC steps) is untested at
  the C++ level; several files are disabled, smoke-only, or non-functional.

The trend line is clearly upward -- the newest files show exactly the right
patterns -- but the oldest and most safety-critical paths remain the least
tested.

---

## Per-file scorecard

| File | Lines | State | Grade |
|---|---|---|---|
| `src/Testing.cpp` (C++ layer) | 3756 | accessor/IO plumbing, model math absent | **D+** |
| testROCNumerical.R | 191 | active, R-oracle vs C++ @1e-12 | **B+** |
| testMCMCROC.R | 321 | active, but regression-snapshot @1% slack | **C** |
| testSimulateGenome.R | 255 | active, structural/IO only | **C+** |
| testMCMCROCSnapshot.R | 44 | repaired but no stored snapshot (no-op) | **D** |
| testFONSE.R | 253 | active, R-oracle + MCMC structural | **B** |
| testPAModel.R | 203 | entirely commented out | **F** |
| testGenome.R | 103 | C++ runner commented; thin R remains | **C** |
| testGene.R | 67 | active (C++ runner + R behavioral) | **B-** |
| testSequenceSummary.R | 156 | active, exhaustive table enumeration | **A-** |
| testFilterRFP.R | 191 | active, per-flag + integration | **A-** |
| testutility.R | 8 | C++ smoke passthrough only | **D** |
| testConvergence.R | 332 | active, all diagnostic APIs | **B+** |
| testAssertChainsMoving.R | 124 | active, true/false-positive paths | **B** |
| test-adaptive-scheme.R | 417 | active, exhaustive boundary validation | **A-** |
| testParameterPersistence.R | 259 | active, bug #388 + concat regressions | **A-** |
| testRestartFileLoad.R | 148 | active, segfault-path regressions | **B** |
| testRestartFilenameStyle.R | 125 | active, both naming conventions | **B+** |
| testParameter.R | 121 | hybrid (C++ black box + real R) | **B** |
| testCovarianceMatrix.R | 22 | C++ smoke; R-side all TODO-commented | **D** |
| testMCMCAlgorithm.R | 46 | getter/setter round-trips only | **C** |
| testPhiPrior.R | 93 | active, sphi-default-off regression | **B+** |
| testPhiMixturePrior.R | 524 | active, math+API+MCMC+recovery | **A-** |
| testArcsineApprox.R | 457 | active, arcsine oracle + Stan-data | **B+** |

---

## Test architecture

Two layers:

1. **C++ unit tests** (`src/Testing.cpp`) exposed via the `Test_mod` Rcpp module
   and invoked by thin R wrappers: `testUtility`, `testSequenceSummary`,
   `testGene`, `testGenome`, `testParameter`, `testCovarianceMatrix`,
   `testMCMCAlgorithm`.
   - The failure mechanism is **correct**: each function latches a `globalError`
     flag on any internal check and `return globalError`, so
     `expect_equal(testX(), 0)` genuinely catches regressions in tested code
     (verified by tracing `testParameter`, `Testing.cpp:2285-3195`).
   - BUT it is a **black box** to R: a failure reports only `0 vs nonzero`, never
     which sub-assertion broke; and numeric comparisons use `==`/`!=` on doubles
     with **no floating-point tolerance** (platform-fragile).

2. **Pure-R integration / numerical tests** -- the larger, newer files.

---

## Strengths (patterns worth propagating)

1. **Independent R oracles vs C++.** `testROCNumerical.R:37-50`,
   `testFONSE.R` (calculateCodonProbabilityVector mirror), and
   `testArcsineApprox.R:65-70` re-implement the model math in R and compare to
   C++ at `1e-12`, including both algebraic branches, limit cases, and
   numerical-overflow regimes. This is the gold standard in the suite.
2. **Named regression tests tied to real bugs/commits.** bug #388 first-element
   drop (`testParameterPersistence.R:210-258`, asserts `tr[1]` and overlap-skip
   `tr[6]/tr[7]`); sphi-prior-default-OFF (`testPhiPrior.R:75`,
   `getSphiPriorSd()==0.0`); restart-file segfault path
   (`testRestartFileLoad.R`); totalRFPCount preservation (`Testing.cpp:1924-2001`).
3. **Documented, intent-revealing fixtures.** `testFilterRFP.R:7-20` states which
   gene passes/fails which filter; each filter tested in isolation, then combined.
4. **Exhaustive boundary validation.** `test-adaptive-scheme.R` rejects every
   out-of-range scheme parameter at both the R and C++ seams (146 assertions).
5. **Complete lookup-table coverage.** `testSequenceSummary.R` enumerates
   `AAToCodon()` for all 22 AA codes in both split-serine modes -- a high-value
   table that would silently break on a codon-order change.
6. **Cholesky math check** (`Testing.cpp:3284-3315`) independently re-implements
   the decomposition -- the one genuine numerical test in the C++ layer.

---

## Weaknesses (priority order)

### P1 -- The scientific core is untested at the C++ level
`ROCModel`, `FONSEModel`, `PAModel`, `PANSEModel` are never touched by
`Testing.cpp`. No likelihood-ratio, proposal-distribution, Metropolis-acceptance,
or adaptation-step unit test exists. The ROC/FONSE codon-probability math is
covered only via R oracles of the **non-log** path; the log-space variant used
in the actual MCMC likelihood (`calculateLogCodonProbabilityVector`,
`calculateLogLikelihoodRatioPerAA`) has no unit test at all.

### P2 -- PA / PANSE essentially have no tests
`testPAModel.R` is 100% commented out (was `skip()`-wrapped even before).
PANSE appears in **no** test file -- not the prob_successful accumulation, the
`lgamma_rfp_alpha` table, nor the three NSE prior options. `testPAParameter`/
`testPATrace` are commented out of the C++ module (`Testing.cpp:3367,3404`),
so Trace initialization has zero coverage for any model.

### P3 -- Disabled / no-op tests that read as "covered"
- `testGenome()` C++ runner is commented out in its wrapper
  (`testGenome.R:6-8`) -- all FASTA/RFP IO and serialization round-trips are
  unreachable from `R CMD check`.
- `testCovarianceMatrix.R:10-22` -- every R-side assertion is in a `#TODO` block.
- `testMCMCROCSnapshot.R` -- syntax repaired, but **no `_snaps/` reference
  exists**, so it cannot fail on a regression (silently records or errors
  "snapshot not found" on first CI run).
- `readObservedPhiValues` block commented out inside `testGenome()`
  (`Testing.cpp:2157-2270`).

### P4 -- Loose integration tolerances + non-determinism
`testMCMCROC.R` pins the final log-posterior to hard-coded round numbers
(`-945000/-946000`) at `tolerance=0.01` -- 1% of ~9.5e5 is **+/-9,450 units** of
slack, chosen to absorb the known thread-unsafe-RNG-in-OpenMP run-to-run variance
(`testMCMCROC.R:54-58`). A real numerical regression of a few thousand units
would pass silently. Root cause (no fixed-seed determinism under OpenMP) is
flagged but never fixed; three trailing MCMC runs (`:153-220`) have zero
assertions.

### P5 -- C++ black-box passthrough + no FP tolerance
`testParameter.R:7`, `testCovarianceMatrix.R:7`, `testMCMCAlgorithm.R:7` expose
only a pass/fail integer. `Testing.cpp` compares doubles with `==`/`!=`.
`testParameter` exercises only the `allUnique` mutation/selection state; the
other three states carry `TODO` markers and are never run
(`Testing.cpp:2567+`).

### P6 -- Untested newer Stan branches
`genomeToStanInit(noncentered=1)` and `genomeToStanData(anchor_phi=1)` code paths
are only checked for field *presence*, never invoked
(`stanDataHelpers.R:221-226`). `.dMPriorFromSCUO` Laplace pseudocount and
`scuo.low.frac` non-default are untested.

---

## Coverage map

| Area | Coverage |
|---|---|
| Codon-probability math (ROC, FONSE, non-log) | Good -- R oracles @1e-12 |
| Log-space likelihood / per-AA LL ratio | **None** |
| ROC MCMC integration | Weak -- loose snapshot, nondeterministic |
| FONSE MCMC | Structural only (no param recovery) |
| PA model | **None** (commented out) |
| PANSE model | **None** (nowhere) |
| simulateGenome distributional correctness | **None** (structural only) |
| Sequence/codon tables, Genome IO (R level) | Good |
| Genome IO (C++ level) | Disabled |
| Convergence diagnostics (Geweke/Gelman/as.mcmc) | Good API coverage; no numeric thresholds |
| Adaptive proposal scheme config | Excellent (config); behavior untested |
| Parameter / restart / MCMC persistence | Good -- real regression guards |
| Phi prior + phi-mixture prior | Excellent |
| Arcsine / Stan-data helpers | Good; two branches untested |
| Covariance matrix (R level) | None (TODO) |
| withPhi measurement-error model | Smoke only |

---

## Meta-finding: CLAUDE.md is stale

`CLAUDE.md` "Test Coverage Gaps" claims testGene.R and testPAModel.R are
"entirely commented out" and "no test files exist for FONSE or PANSE." On the
current main line: testGene.R and testFONSE.R are **active** with real
assertions; only testPAModel.R is still commented; PANSE is genuinely untested.
The doc predates the 2026 test-development work and should be refreshed.

---

## Recommendations (high-leverage, in order)

1. **Add C++ likelihood unit tests** for at least ROC and FONSE: feed known
   parameters + a tiny gene, assert `calculateLogLikelihoodRatioPerAA` against a
   hand-computed value. Closes the single largest gap (P1).
2. **Make ROC MCMC deterministic** (fixed seed, single-thread test build) so
   `testMCMCROC.R` can tighten from 1% to a real regression bound (P4).
3. **Re-enable `testGenome()`** in its wrapper and add FP tolerance to C++
   numeric comparisons (P3, P5).
4. **Stand up a minimal PA/PANSE test** (short fixed-seed chain + recovery on
   simulated RFP data); resurrect `testPAModel.R` (P2).
5. **Commit a reference snapshot** for `testMCMCROCSnapshot.R` or delete the file
   (P3).
6. **Cover `noncentered=1` / `anchor_phi=1`** Stan-init branches (P6).
7. **Refresh CLAUDE.md** coverage section to match reality.
