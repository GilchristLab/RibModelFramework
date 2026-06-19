# Handoff: unified NSEROC model family (GLM-EM + Stan) -- merge into main

**Date:** 2026-06-19  **Status:** plan for a fresh session (this worktree); design
settled in conversation, not yet built. **Worktree:** this one
(`RibModelFramework-nseroc-glmm`, branch `feat/nseroc-glmm-unified`, off RMF main
@ 394ee47, which already carries the validated ROC backend trio).

## Goal
ONE unified codon-model family in which the named models are special cases of a
single parameterization. Land it in main as the **GLM-EM (R) fitter + Stan models**
now; the **AnaCoDa C++ model unification is a deferred follow-up** (the production
native models stay separate for now).

  ROC      = NSEROC with the NSE term OFF
  FONSE/SONSE/NSE = NSEROC with the ROC term OFF ("ROC=0"), at first/second/full order
  NSEROC (both terms on, full order) is the UNIVERSAL PARENT.

## The unified model (one equation, two orthogonal switches)
```
eta_c(g,pos) = log P(c | g,pos)/P(ref) = -dM_c - phi_g * ( I_eta * dEta_c  +  I_om * dOmega_c * W_k(order) )
```
- `dM_c`   mutation (always present; CSP intercept).
- `dEta_c` ROC selection (position-INDEPENDENT). Present iff `I_eta` (est.dEta).
- `dOmega_c` NSE selection (position-DEPENDENT via the known weight `W_k`).
  Present iff `I_om` (est.dOmega).
- `W_k(order)` known per-position weight (elongation waiting times u + scalar
  nonsense rate b); set by `order in {first, second, full}` -- the validated
  `.compute_W` from fonse-orders-glm.R / fit.fonse.glm.em.R.
- NATIVE AnaCoDa sign (positive dEta/dM/dOmega = cost); report dM=-a, dEta=-b,
  dOmega=-g from the raw GLM coefs (a,b,g), so all models agree with native + ROC.

| model     | I_eta (ROC) | I_om (NSE) | order  |
|-----------|-------------|------------|--------|
| ROC       | on          | OFF        | --     |
| FONSEROC  | on          | on         | first  |
| SONSEROC  | on          | on         | second |
| NSEROC    | on          | on         | full   |
| FONSE     | OFF         | on         | first  |
| SONSE     | OFF         | on         | second |
| NSE       | OFF         | on         | full   |

Mechanically it is ONE per-AA multinomial GLM with 1-3 covariate columns:
intercept = dM always; covariate `phi` -> dEta_c when I_eta; covariate `phi*W_k`
-> dOmega_c when I_om. The phi-step (NB/Poisson-Gamma Fisher-scoring, median-gauged
Gamma(sphi) prior) runs on the composite `(I_eta*dEta_c + I_om*dOmega_c*W_k)`. For
SONSE/full add the 1-D plug-in EM step for `b` (recompute W_k); first-order skips it.

## What already exists vs the gap
- `fit.roc.glm.em.R` (RMF prototypes/, validated): the `I_eta=1, I_om=0` corner.
- `fit.fonse.glm.em.R` v0.2.0 (in the OTHER session's worktree
  `Analyses-RibModelFramework-fonse-glmm-gamma`, branch feat/fonse-glmm-gamma, 7
  ahead, UNMERGED): already has `--dOmega T/F` (off => ROC reduction) and
  `--order FONSEROC|SONSEROC|NSEROC`, i.e. the WHOLE `I_eta=1` row.
- **THE GAP:** no `est.dEta` toggle -> the pure-NSE cases (FONSE/SONSE/NSE with
  ROC=0, `I_eta=0`) cannot be expressed yet. Add that second switch.

## Plan (this session)
1. **Port + consolidate (R GLM-EM):** copy the validated `fit.fonse.glm.em.R` +
   `fonse-orders-glm.R` `.compute_W` here as the unified `prototypes/fit.nseroc.glm.em.R`,
   ADD the `est.dEta` toggle (drop the `phi` covariate when off), and expose a
   `--model {ROC,FONSEROC,SONSEROC,NSEROC,FONSE,SONSE,NSE}` preset over
   `{est.deta, est.domega, order}`. PORT ONLY after the fonse-glmm-gamma session is
   done + retired (Mike retires it once its tests/validations pass) -- do NOT edit
   that worktree.
2. **Stan family:** generalize `stan/roc_basic.stan` to `stan/nseroc_basic.stan`
   carrying the `dOmega * W_k` covariate (W_k passed as data, since it is known);
   ROC = `dOmega` block off / zeroed. Drive it from the clean RMF Stan wrapper
   pattern (`prototypes/fit.roc.stan.R`). C++ AnaCoDa unification deferred.
3. **Checks:** extend the per-backend recovery + cross-backend agreement checks to
   the family (GLM-EM vs Stan), at the canonical fixed sphi=1.4.

## Validation ladder (do in order)
1. `I_om=0` (ROC) MUST reproduce the validated RMF ROC trio anchor on S288c
   (phi ~0.97/0.996, dEta/dM |r| ~0.99). Same anchor as fit.roc.glm.em.R.
2. Add dOmega (FONSEROC first-order): match native FONSE + the plug-in
   `fonse-orders-glm.R` estimates.
3. Pure-NSE (`I_eta=0`): sanity that dOmega alone is recoverable; and that ROC vs
   NSE are separable when BOTH are on -- this leans entirely on the position
   leverage in W_k (the dEta/dOmega correlation / `b_c*u_c` product is the known
   identifiability risk; TEST it, do not assume).
4. Cross-backend GLM-EM vs Stan on the NSE family (mirror roc_cross_backend_check.R).

## Conventions to carry
- Canonical FIXED **sphi=1.4** for benchmarks (sphi is a gauge; see the ROC notes +
  the sphi-prior-default memory). simulate AND fit at 1.4.
- Gamma phi prior = the analytic-phi-step enabler (closed-form Laplace transform;
  lognormal has none -> forces MCMC). Median-gauged.
- Native extraction gauge gotcha: getCSPEstimates(relative.to.optimal.codon=FALSE)
  to compare to the parameter-object truth.

## Pointers
- ROC anchor + theory: this repo `prototypes/fit.roc.glm.em.R`, `prototypes/fit.roc.stan.R`,
  `notes/sphi-Ne-metrics-gamma-theory.md` (Sec 9b = cross-backend result).
- W_k / orders source (port once retired): fonse-glmm-gamma worktree
  `fonse/fit.fonse.glm.em.R` v0.2.0, `fonse/fonse-orders-glm.R`,
  `fonse/notes/fonse-glmm-gamma-handoff.md` + `terminology.md`.
- Bigger payoff (from the fonse handoff): seed fonse-stan with the fast joint-phi
  fit + Fisher SEs; reconstruct PANSE footprint-derived estimates from sequence-only
  F/SONSE.
