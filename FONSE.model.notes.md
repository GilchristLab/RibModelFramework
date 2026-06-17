# FONSE model: math notes for the dEta + Stan changes (feat/fonse-native)

Context for reviewers of this branch. The full derivation lives in the FONSE
paper repo (`fonse_paper/semppr_math/eta.analyses.tex`, §"Defining η" and
§"Approximations"); this is the short version needed to review the code.

## The model

FONSE (First-Order approximation of NonSense Error) treats codon usage as the
evolutionary outcome of selection against nonsense (premature-termination)
errors during translation. The per-codon log-probability is an exact
categorical logit over the synonymous codons of each amino acid:

    η_i = −dM_i − dEta_i·φ − dOmega_i·φ·β(pos),     β(pos) = a1 + a2·(pos−1)
    P(codon i | pos) = softmax(η)_i,                reference codon: η = 0

- **dM_i**      mutation (GC/nucleotide) bias; position-independent.
- **dEta_i**    position-INDEPENDENT elongation selection — identical to ROC's
                selection coefficient (ribosome-sequestration / elongation cost).
                **New in this branch.**
- **dOmega_i**  position-DEPENDENT nonsense-error selection (the FONSE term);
                ω_i = b/c_i, the nonsense-to-elongation odds.
- **φ**         gene-specific protein synthesis rate (latent, estimated).
- reference codon = last alphabetical synonym, fixed at 0; softmax normalises
  within each amino acid.

So **FONSE = ROC (dM + dEta) + the nonsense term (dOmega)**, and it reduces
exactly to ROC as dOmega → 0 (nonsense rate b → 0).

## Where the exact model lives (why this is "first-order")

The exact protein-production cost (Gilchrist & Wagner) carries a survival
product:

    η = Σ_i β_i (b/c_i) ∏_{k>i}(1 + b/c_k) + β_{n+1},   p_k = c_k/(c_k+b)

- **FONSE** = the first-order Taylor expansion of this around b = 0 (the product
  → 1). **This branch implements first-order only** (native, Stan, and the R
  logistic prototype). The 2nd-order ("SONSE") and full forms are future work.
- **PANSE** already implements the exact `c/(c+b)` survival (incomplete-gamma)
  — but for ribosome-footprint data, not genomic codon usage
  (`inst/stan/panse_*.stan`).

## What this branch changes

1. **dEta as a third codon-specific (CSP) category.** `FONSEParameter::dEtaCSP =
   2` — a distinct index, because the inherited `Parameter::dEta = 1` *aliases*
   `dOmega` in FONSE and cannot be reused. dEta uses `numSelectionCategories`;
   the per-AA covariance grows to `(numMut + 2·numSel)·numCodons`. dEta defaults
   **fixed at 0** → byte-identical to prior FONSE; `estimateDEta()` / `--est_eta`
   opts in (mirrors the a1/a2 fixed-by-default decision).

2. **Segfault fix.** `Trace::getCodonSpecificCategory` returned the unset
   `nse = −1` sentinel as `unsigned` for FONSE's paramType-2 retrieval →
   out-of-bounds → segfault on any Eta posterior/variance/quantile call. Now
   falls back to the selection category when `nse < 0` (PANSE's nse path
   unchanged).

3. **Stan backend** (`inst/stan/fonse_sphi_est.stan` + `R/fonseStan.R`). Because
   β(pos) sits inside the exponent, the per-AA normaliser is position-dependent,
   so — unlike ROC/PANSE — the likelihood **cannot collapse to per-gene codon
   counts**. The data is therefore **per-position** (one row per sense codon),
   grouped by gene for `reduce_sum`. Reuses ROC's phi.mphi / phi.sphi spec and
   the centered/non-centered switch verbatim.

## Suggested reviewer focus

- **Covariance-block indexing** in `FONSEParameter::proposeCodonSpecificParameter`
  — the dEta proposal offset `(numMut+numSel)·numCodons + j` and the iid draw
  count `numCodons·(numMut + 2·numSel)` must stay consistent with the covariance
  matrix size, or proposals corrupt.
- The **paramType-2 category fallback** in `Trace.cpp` and that it leaves PANSE's
  NSERate path untouched.
- **Restart I/O backward-compat**: old `.rst` files lacking the
  `>currentEtaParameter:` block must still load (dEta defaults to 0).
- The Stan **per-position data layout** and `reduce_sum` grouping by gene
  (`build_fonse_stan_data` / `fonse_genome_to_obs`).
