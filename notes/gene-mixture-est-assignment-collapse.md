# Gene-mixture (est.mix) collapse: analysis + upstream issue draft

> Promoted from Analyses-RibModelFramework `roc/notes/` @ 0476b94 (2026-06-19).
> Filed upstream as acope3/RibModelFramework#431.

**Date:** 2026-06-17
**Context:** ROC native MCMC, `num.mixtures>1` + `mixture.definition="allUnique"`,
estimated gene assignments (`est.mix=TRUE`). Observed: the two components collapse
onto identical dM/dEta. Investigated whether this is a bug.

## Finding (code-verified)

The collapse is **inherent finite-mixture non-identifiability**, not a bug. The
machinery is statistically correct:

- Gene reassignment is a correct Gibbs draw `p(z_g=k) prop w_k * L_k(gene)`
  (`src/MCMCAlgorithm.cpp:186-254`).
- The per-gene assignment likelihood uses each component's OWN dM/dEta
  (`src/ROCModel.cpp:197-218`; allUnique => `categories[k].delM=k, delEta=k`,
  `src/Parameter.cpp:1170-1172`).
- The CSP MH update partitions genes by their current assignment correctly
  (`src/ROCModel.cpp:288-301`).

With identical exchangeable priors on both components and a codon-usage
likelihood that is only weakly separating once phi is conditioned on, the
symmetric state (both components equal, genes split ~50/50) is a stable fixed
point. Fixed assignments work only because the partition is held externally
(`getCSPbyLogit` seeds components from gene subsets, `R/parameterObject.R:1419`).

## Real, fixable wart: Dirichlet mixing weights have no pseudocounts

`src/MCMCAlgorithm.cpp:401-411` sets the Dirichlet parameters to the raw
assignment counts with no base concentration; `randDirichlet`
(`src/Parameter.cpp:2991-3017`) draws `y_k ~ Gamma(count_k, 1)`. If a component
empties, `Gamma(0,1)=0` deterministically -> its weight is exactly 0 and the
state is **absorbing** (a dead component can never revive). This is independent
of the non-identifiability but makes degeneration fast and irreversible.

### Proposed upstream issue (acope3/RibModelFramework)

**Title:** Gene-mixture Dirichlet weight update lacks concentration pseudocounts -> empty component is an absorbing state

**Body:**
In the mixture-proportion update (`MCMCAlgorithm.cpp:401-411`), the Dirichlet
parameters are the raw per-category assignment counts with no base concentration
added. Combined with `randDirichlet`'s `Gamma(shape=count,1)` draw
(`Parameter.cpp:2991-3017`), a category that transiently empties gets shape 0 ->
weight exactly 0 -> it can never recover, even if some genes would prefer it.
This makes `est.mix=TRUE` runs prone to irreversible component death.

Suggested fix: add a small concentration `alpha` (e.g. 1, or 1/K) to each
category count before the Dirichlet draw:
`dirichletParameters[k] = alpha + count_k`. This keeps empty components
recoverable (a symmetric Dirichlet(alpha) prior on the weights). Optionally
expose `alpha` as a setting.

Note: this does NOT fix the deeper non-identifiability of `allUnique` gene
mixtures under identical per-component priors (that needs component-distinct
priors, an ordering constraint, or fixed/anchored assignments) -- it only removes
the absorbing-empty-component failure mode.

## Make estimated-assignment mixtures usable (beyond the wart)

1. Component-distinct priors -- the per-component prior already accepts a
   per-category matrix (`R/parameterObject.R:522,538`); pass distinct means.
2. Ordering/identifiability constraint on the gene-mixture categories.
3. Dirichlet pseudocounts (above).
4. Semi-fixed / anchored assignments (pin markers, let the rest move).
