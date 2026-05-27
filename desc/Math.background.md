# Mathematical background for RibModelFramework

The models implemented here (ROC, FONSE, PA, PANSE) rest on a longer
analytical lineage of derivations around the elongation cost function
`eta(c)` and its moments.  The foundational math notes live in a
separate repo:

    https://github.com/mikegilchrist/semppr_math
    (clone: git@github.com:mikegilchrist/semppr_math.git)

That repo (renamed 2026-05-21 from `eta_analysis`, original SEMPPR-era
location `/home/semppr/Projects/SEMPPR/Analysis/eta.analyses.tex`)
contains a single ~1300-line LaTeX document, `eta.analyses.tex`, with
~13 years of derivations:

- Distribution `f(eta)` from amino-acid sequence + codon elongation
  table; mean `E(eta)`, variance `Var(eta)` (Russ Zaretzki's corrections
  via covariance).
- First- and second-order Taylor expansion of `eta` around `b = 0`
  -- the math that motivates the **FONSE** approximation.
- Cross-position cov terms -- foundation for **PANSE**'s per-position
  formulation.
- Fixation probability of `eta` under Gamma / Normal distributions
  (the population-genetic foundation of **ROC**).
- `Delta eta` for codon substitutions (Laura Salter).
- Multinomial-parameter analysis and an alternative `eta`
  approximation (Christopher Oballe, REU 2013).
- Missense error model.
- FONSE model notes (Jeremy Rogers, 2016).
- Additive formulation, SELAC cross-references.

If you are reading the C++ source for a model and want to know where a
particular equation comes from, `semppr_math/eta.analyses.tex` is
likely the place.  The document is not actively maintained (content
frozen since Feb 2017) but its derivations are the source-of-truth
that the paper PDFs ultimately reference.

## Also referenced as a submodule

The `fonse_paper` repo
([`mikegilchrist/fonse_paper`](https://github.com/mikegilchrist/fonse_paper))
includes `semppr_math` as a git submodule at `semppr_math/`.  The
broader nonsense-error-rates paper repo
([`mikegilchrist/Nonsense_error_rates`](https://github.com/mikegilchrist/Nonsense_error_rates))
does not currently reference it as a submodule, though it could.
