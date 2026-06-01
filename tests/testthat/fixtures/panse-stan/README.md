# PANSE Stan test fixtures

These fixtures are for the **PANSE Stan port** test suite (`test-panse-stan-*.R`,
which compile and check the `stan/panse_*.stan` models via cmdstanr). They are
distinct from the native **PANSE MCMC** tests (developed separately); do not
mix the two.

## `likelihood_fixture.rds`

A small, deterministic fixture for the PANSE Stan likelihood/regression tests
(no MCMC; the model's generated-quantities `log_lik` is checked against an
independent R reference at the truth parameters).

Contents (a single list):
- `stan_data`  -- Stan data list (G=15, C=22, P=900) built by
  `build_panse_stan_data()` from a tiny `nb-2o-approx` simulation
  (shared NSE = 1e-5, sphi = 1.0, seed = 424242), init = truth.
- `truth`      -- truth parameters: `alpha[C]`, `lambda[C]`, `nse_shared`,
  `phi[G]`, `sphi`, `U`, `codon_order`.
- `ref_log_lik`-- length-P reference `log_lik` at truth, computed independently
  in R (NB2 mean `phi*alpha/(U*lambda')*survival`, size `alpha`; survival via
  the exact forward-Lentz upper-incomplete-gamma hybrid). Total = -3059.30.
- `meta`       -- G/C/P, seed, note.

The fixture embeds **no absolute paths** and is platform-independent.

## Provenance / regeneration

Generated from a tiny simulation produced by the PANSE simulator. The
generator currently lives in the analysis repo (A-RMF
`panse/s.cerevisiae/scripts/sim/simulate_panse_dataset.R` +
`scripts/lib/build_panse_stan_data.R`); it will move into this repo with the
harness consolidation (GilchristLab/RibModelFramework#26), after which a
repo-local `make_fixtures.R` will regenerate this file. Until then the
committed `.rds` is the source of truth for the tests.

The reference `log_lik` deliberately re-implements the likelihood independently
(it is NOT produced by the Stan model), so the tests cross-check the Stan
implementation against it.
