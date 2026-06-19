# ROC: what sphi/dEta mean (Ne, gBGC), the metric bridge, and Gamma analytics

> Promoted from Analyses-RibModelFramework `roc/notes/` @ 0476b94 (2026-06-19).
> Companion prototype: `prototypes/fit.roc.glm.em.R`.

**Date:** 2026-06-17
**Status:** theory synthesis from a working conversation (Mike Gilchrist + Claude).
Method-level; not Loki-specific. Companion to the Loki convergence work
(analysis/03x.../convergence.summary.md) and [[dEta-as-Ne-signal]].

The model throughout: for an AA, gene g uses synonymous codon c with
p_c(g) = softmax(-(ΔM_c + ΔEta_c·φ_g)) -- mutation term ΔM (gene-independent)
plus selection term ΔEta·φ (gene-specific via expression φ).

---

## 1. sphi is identified only through dEta·φ

The codon data depend on φ ONLY through the product dEta·φ. So they inform the
distribution of dEta·φ, not φ itself. sphi = SD(log φ) is identifiable only via
the selection channel: roughly sd(dEta·φ) ~ |dEta|·sphi. Shrink |dEta| -> the
data constrain |dEta|·sphi small, and the dEta prior + near-flat φ likelihood
pull recovered sphi DOWN. As dEta -> 0, φ is fully unidentified and recovered
sphi -> 0 regardless of the true value.

Implication: recovered sphi is a measure of the STRENGTH of selection on the
expression gradient (Ne·s magnitude), not a clean expression-variance estimate.
Low recovered sphi can mean weak codon selection (small Ne), not narrow
expression. (Testable by simulation: true sphi fixed at 1.5, dEta dialed down ->
recovered sphi falls. Queued.)

## 2. sphi vs Ne is NON-MONOTONIC (a hump)

A codon shows a usable across-gene gradient only when γ_g = dEta·φ_g ~ O(1)
straddles the φ range (some genes drift-dominated |γ|<1, some selection-
dominated |γ|>1). Both Ne extremes kill the gradient:
- Ne·s -> 0: all γ << 1, mutation equilibrium, no gradient -> recovered sphi -> 0.
- Ne·s -> inf: all γ >> 1, every gene saturated at the optimum, no gradient ->
  recovered sphi -> 0.
So recovered-sphi is a humped, saturating function of Ne, NOT a monotone ruler:
two species at opposite Ne extremes can both show low sphi. (Under expression-
DEPENDENT selection the saturation erodes the gradient top-down; full uniformity
needs Ne·s·φ_min >> 1.)

## 3. The dEta scale bundles Ne, generation time, and expression scale

s is a per-generation fitness effect. Per-generation cost = (per-event cost) x
(translations per generation) = dEta' · φ_rate · T (T = generation time). With
φ gauged to mean/median = 1 (relative phi = φ_raw/C, C = mean or median absolute
expression):

  dEta = Ne · dEta' · C   and   C = <φ_rate> · T   =>   dEta = Ne · dEta' · <φ_rate> · T.

So dEta is a ruler for **Ne · T · <φ_rate> · dEta'**, four things the codon data
cannot separate. Only dEta' (intrinsic per-event cost) is plausibly species-
constant. Cross-species dEta_i/dEta_j = (Ne·T·<φ_rate>)_i / (...)_j -- NOT Ne
alone. A slow deep-sea archaeon's long T can offset its small Ne. sphi is clean
of all this (= SD(log φ_raw), gauge- and T-invariant); only dEta carries the
composite. This sharpens [[dEta-as-Ne-signal]]: the unidentified scale on dEta
is Ne·T·<φ_rate>, and "dEta ratio = Ne ratio" needs T and absolute expression
scale species-constant (needs external absolute-expression data to pin).

## 3b. tRNA gene copy number as an external anchor (and null test) for dEta

We might further constrain dEta values -- and hence the scale that ties sphi to
dEta (Sec 1, 3) -- by regressing dEta against **tRNA gene copy number**, the
standard proxy for cognate tRNA abundance. The degenerate-AA structure tells us
where this can and cannot work:

- **2-codon AAs are a NULL, not a probe.** They are traditionally read by a
  SINGLE shared anticodon (one codon cognate, the other via wobble). If dEta is
  driven SOLELY by cognate/wobble tRNA reading efficiency, then increasing that
  one tRNA's gene copy number scales the translation rate of BOTH codons together
  -- so the WITHIN-AA dEta *difference* (the only thing ROC estimates, the
  reference codon being pinned) should NOT shift with copy number. So 2-codon AAs
  give no leverage to predict the dEta contrast from tRNA counts under the simplest
  model. Usefully, this also makes them a FALSIFIABLE control: if a 2-codon dEta
  contrast DOES move with the shared tRNA's copy number across genes/species, then
  dEta is NOT solely cognate-wobble driven (other forces -- e.g. wobble-pair
  modification, mRNA secondary structure, demand-vs-supply effects -- are in play).
- **4- and 6-codon AAs are where the signal lives.** These read with MULTIPLE
  tRNA species / anticodons whose gene copy numbers can differ within the AA, so
  the per-codon dEta contrasts CAN be regressed on (differences in) cognate tRNA
  copy number. That variation is the external handle for anchoring the dEta
  magnitude.

This is a route to the dEta-scale prior the GH analysis says we need (v0.6
below / issue #14): instead of (or alongside) an arbitrary ridge SD on dEta, the
tRNA-copy-number regression on the 4/6-codon AAs could supply an INFORMATIVE,
mechanistic prior on the dEta scale, with the 2-codon AAs as the consistency
check. Connects to [[dEta-prior-from-trna]] and [[dEta-as-Ne-signal]]. Caveat:
needs per-genome tRNA-gene annotations (tRNAscan-SE), which the metagenomic Loki
bins may not cleanly provide.

## 4. gBGC: phi-independence is the discriminator

gBGC is selection-like and Ne-dependent (strength ~ Ne·b·r) but depends on
recombination rate r, not expression φ.
- Uniform recombination -> gBGC is a constant GC push == dM (absorbed in mutation).
- ROC's dEta·φ term is phi-DEPENDENT; gBGC is phi-INDEPENDENT. So ROC
  AUTOMATICALLY separates translational selection (phi-dependent) from
  mutation+gBGC (phi-flat). This is exactly what the intergenic-dM prior
  exploits: non-coding GC = mutation+gBGC, no translational selection.
- ROC's dM is genome-UNIFORM, so VARIABLE-recombination gBGC cannot be absorbed
  by dM; it shows as per-region residual and only contaminates dEta if local r
  correlates with φ. To separate it: per-region/local intergenic-dM, or a
  recombination-rate covariate (yeast: Mancera 2008 map, or population rho).
- Caveat for the Ne program: gBGC strength ~ Ne·r·b scales with the SAME Ne we
  want from dEta -- controlling for it sharpens the signal but r and Ne aren't
  independent knobs.

## 5. The two Ne extremes are disambiguated by ENCp level

Both extremes flatten the across-gene VARIATION (sphi -> 0), but differ in the
ABSOLUTE bias LEVEL:
- High Ne: ENCp uniformly very low (extreme bias beyond mutation), low Nc.
- Low Ne: ENCp near the mutational null (bias ~= dM only), high Nc unless dM
  extreme. gBGC ~ 0 here too (drift-limited).
So sphi (variation) says "you're at an extreme"; ENCp (level, vs the dM null)
says WHICH. For drift-dominated low-Ne MAGs, low sphi AND ENCp near the
mutational null should agree.

## 6. Metrics are scalar projections of the same softmax (logit linearity)

For a 2-codon AA the model is logistic, and in log-odds it is LINEAR in φ:

  logit p_1(g) = log(p_1/p_2) = ΔM + ΔEta·φ_g     (intercept ΔM, slope ΔEta).

Every metric is a nonlinear reduction of that line, differing by REFERENCE and
Renyi ORDER:

| metric | reference | isolates | order |
|--------|-----------|----------|-------|
| ENC    | uniform   | total bias | 2 (Simpson, F=sum p^2) |
| SCUO   | uniform   | total bias | 1 (Shannon; H_max-H = KL(p\|\|uniform)) |
| ENCp   | mutational null softmax(-ΔM) | SELECTION (ΔEta·φ) | 2 |
| CAI    | high-φ optimum | adaptation ~ φ | geom-mean / order-1-ish |

Exact for the 2-codon AA:
- ENC: F = p1^2+p2^2 = cosh(x)/(1+cosh x) => **Nc = 1 + sech(x)**, x = ΔM+ΔEta·φ
  (sech 0 =1 -> Nc=2 even; sech inf =0 -> Nc=1 one codon).
- ENCp: logit(observed) - logit(mutational null at φ=0) = (ΔM+ΔEta·φ) - ΔM =
  **ΔEta·φ**. ENCp's background subtraction == removing the ΔM intercept; its
  residual ~ sech of ΔEta·φ. (ENCp's correction is the metric-world analog of
  ROC's dM term: background = dM, residual = dEta·φ.)
- SCUO: 1 - H_b(σ(x))/log2 (binary entropy of the same x).

So the model and metrics are the same equation read three ways: ENC/SCUO see the
FULL logit, ENCp sees logit MINUS intercept, the slope itself is dEta.

Mean ENCp across genes ~ dEta·C (selection at the reference gene = the Ne·T·
<φ_rate>·dEta' composite); across-gene SPREAD of ENCp/CAI ~ dEta·(spread of φ)
= the gradient ROC formalizes as sphi. So sphi == across-gene variance of the
selection metric; both Ne extremes flatten it (Sec 2/5).

Operational link already in the pipeline: φ is seeded by scuo_to_log_phi(
-encp_raw) ([[encp-init-pattern]]) -- the metric ranks genes by selective bias
(monotone in φ), the model then factors that one scalar into ΔM, ΔEta, φ.

## 7. CAI critique (falls out of Sec 6)

CAI takes high-φ genes' OBSERVED usage as "optimal." In the line picture that
reference is p1(φ_high) = σ(ΔM + ΔEta·φ_high): a point at FINITE φ on an
UNSATURATED sigmoid that still carries the ΔM intercept. So CAI's optimum is
wrong twice:
1. Mutation-contaminated: includes ΔM, can't tell "optimal" from "mutation-
   favored." True selective optimum is the φ->inf asymptote (pure ΔEta ranking,
   ΔM-free). When ΔM and ΔEta favor DIFFERENT codons, CAI's reference codon !=
   ROC's ΔEta-optimal codon -- a concrete testable divergence.
2. Unsaturated: codon freq still RISING with φ at the top (observed in yeast for
   some AAs) means ΔEta·φ_max ~ O(1), not >>1 -- even the highest-expression
   genes sit below the ceiling, so the reference is a moving, sub-asymptotic
   target. (Also: this "still climbing" is direct evidence yeast is in the
   informative INTERMEDIATE Ne regime, not saturation.)
ROC escapes both: estimates slope (dEta) and intercept (dM) separately and reads
the asymptote; CAI has only finite-φ frequencies and no mutation reference.

## 8. Gamma analytics for a per-gene phi estimate

From the line, the bare per-AA MLE is closed form (invert the logit):
  **φ_hat_a = (log(n1/n2) - ΔM_a) / ΔEta_a.**
Joint over a gene's 2-codon AAs: concave 1-D log-lik; Gaussian/Laplace form is an
inverse-variance-weighted average of the per-AA inversions, weights = Fisher
info I_a = n_a p_a(1-p_a) ΔEta_a^2. Fast, analytical, given dM/dEta -- the basis
of a coordinate-ascent/EM ROC (φ-step = this; CSP-step = weighted logistic
regression of codon log-odds on φ; cheaper than AnaCoDa per-gene MH).

Breakdown: a fixed codon (n2=0) -> logit = inf -> φ_hat -> inf (logistic
"perfect separation"), exactly the high-φ genes. A prior on φ is needed.

WHY Gamma helps -- and the sharpening: a Gamma prior placed DIRECTLY on φ gives
closed-form prior-predictive codon frequencies/moments via its Laplace transform
E[e^{-sφ}] = (1+s/b)^{-a} (lognormal has NO closed-form Laplace transform, which
is why it forces MCMC), and it regularizes the separation blow-up; but the full
posterior of φ is NOT closed form (the Poisson e^{-lambda} = e^{-kappa e^{-dEta φ}}
double-exponential resists). The EXACTLY conjugate object is the Gamma prior on
the SELECTION MULTIPLIER psi = e^{-dEta·φ} (the multiplicative codon-preference
fold-change), in the rare-codon/Poisson limit (softmax normalizer ~ 1):
  n2 ~ Poisson(kappa·psi), psi ~ Gamma(a,b)  =>  psi | n2 ~ Gamma(a+n2, b+kappa)
  (Poisson-Gamma conjugacy), marginal n2 ~ NegBinom, all moments closed form.
psi = e^{-dEta·φ}  =>  φ = -log(psi)/dEta, so a Gamma posterior on psi is a
LOG-GAMMA posterior on φ (positive support, right tail -- a lognormal cousin).
So the user's "Gamma on φ" is close; the analytically-blessed version is Gamma
on e^{-dEta·φ}.

## 8b. Median-gauged Gamma reparameterization (the sphi analog)

To use phi ~ Gamma as a median-gauged sibling of the lognormal, reparameterize so
ONE parameter sets the gauge (median phi = 1) and the OTHER controls the spread,
with the spread = sphi = SD(log phi) -- the SAME quantity as the lognormal sphi.
Two facts make it clean:
- Var(log phi) = trigamma(shape k), shape-ONLY (independent of rate/scale).
- median = m_k / rate, with m_k = qgamma(0.5, k, rate=1) = median of Gamma(k,1).
So gauge (rate) and spread (shape) are ORTHOGONAL, mirroring (mphi, sphi) for the
lognormal. The reparameterization:

  shape k    = inv_trigamma(sphi^2)          # spread; sphi = SD(log phi)
  rate  beta = qgamma(0.5, k, rate = 1)      # gauge;  median(phi) = 1

| sphi | shape k | rate beta | median | SD(logphi) | CV    | mean(phi) |
|------|---------|-----------|--------|------------|-------|-----------|
| 0.25 | 16.495  | 16.163    | 1.000  | 0.250      | 0.246 | 1.021     |
| 0.50 | 4.479   | 4.151     | 1.000  | 0.500      | 0.473 | 1.079     |
| 1.00 | 1.426   | 1.110     | 1.000  | 1.000      | 0.837 | 1.285     |
| 1.50 | 0.811   | 0.512     | 1.000  | 1.500      | 1.110 | 1.584     |
| 2.00 | 0.566   | 0.286     | 1.000  | 2.000      | 1.329 | 1.983     |

So sphi is a UNIVERSAL spread parameter: lognormal sets sigma=sphi; Gamma sets
trigamma(k)=sphi^2; both median-gauged at 1, so the two models compare apples to
apples. Note: for sphi >~ 1 the shape k < 1.5 (and < 1 for sphi > ~1.4) -> the
median-gauged Gamma is monotone-decreasing (mode at 0). This is also the
median-gauge regularizer used by the GLM-EM fitter (fit.roc.glm.em.R).

CAVEAT (feeds Sec 9): matching BOTH median AND sphi still leaves the LINEAR-scale
spread very different at high sphi -- CV at matched sphi:
  sphi 0.25: logN 0.254 vs Gamma 0.246  (interchangeable)
  sphi 1.00: logN 1.31  vs Gamma 0.84
  sphi 1.50: logN 2.91  vs Gamma 1.11   (lognormal far heavier-tailed; mean 3.08 vs 1.58)
Near-identical at low sphi, sharply divergent at high sphi.

## 9. Gamma vs lognormal behavior (and the normal limit)

Both are positive, right-skewed; both -> Normal as the coefficient of variation
CV = SD/mean -> 0 ("E(X) >> Var(X)" = concentrated relative to mean). At matched
small CV their skewness differs by a constant:
- Gamma(shape a): CV = 1/sqrt(a), skew = 2/sqrt(a) = **2·CV**, exponential right tail.
- Lognormal(σ): CV = sqrt(e^{σ^2}-1) ~ σ, skew ~ **3·CV** (=3σ+σ^3), heavy
  sub-exponential right tail.
So at the same CV, **lognormal is ~1.5x more skewed and heavier-tailed** ->
approaches Normal SLOWER, and models the extreme high-expression tail better
(expression is multiplicatively generated -> lognormal is the natural CLT).
Near 0 they also differ: Gamma(a<1) piles mass at 0 (mode at 0); lognormal
density -> 0 at the origin always.
Practical consequence: the two are nearly interchangeable (both ~ Normal) ONLY
in the LOW-sphi / low-CV regime -- i.e. exactly the low-sphi Loki MAGs, where the
Gamma analytics (Sec 8) are an excellent approximation. For high-sphi yeast
(σ ~ 1.5, large CV) they diverge sharply and Gamma is a poor stand-in for the
heavy lognormal tail. So: use Gamma for cheap analytics where CV is small;
keep lognormal where the tail matters.

## GLM-EM proof of concept (2026-06-18) -- VALIDATED on S288c

`Analyses-RibModelFramework/roc/scripts/fit.roc.glm.em.R`. ROC fit as a
multinomial GLMM by coordinate ascent: CSP step = per-AA multinomial logistic
GLM (VGAM::vglm, codon counts ~ phi, reference codon matched to AnaCoDa from the
published CSV); phi step = per-gene bounded 1-D maximize of the multinomial
loglik + the median-gauged Gamma(sphi=1) prior; NO separate gauge step (the prior
pins the scale). Deterministic, monotone log-post, converges in ~12 iterations
(seconds, vs MCMC hours).

S288c (1000 genes, all 19 AAs / 59 codons) vs the native dual-LN ROC fit +
published dM/dEta:
- **phi: Pearson 0.974, Spearman 0.995** (was 0.62 with only the 9 two-fold AAs
  -> the multinomial CSP over all AAs is what closes it).
- **dEta: |r| = 0.999** (40 non-ref codons); **dM: |r| = 0.992**. (Sign is flipped
  = my log-odds(c/ref) convention vs AnaCoDa's energy-form -(dM+dEta*phi).)
So the GLM-EM reproduces the native MCMC fit to within noise, in seconds, with no
MH / funnel / mode-flip stochasticity (the PC1/encp init anchors the mode). This
is the "GLM route" paying off concretely (Sec 8). Caveats: MAP not full posterior
(Laplace/bootstrap for SEs); prior sphi held fixed at 1.0 (recovered SD(log phi)=0.89
unprompted); the closed-form NB/Poisson-Gamma phi-step (warm start) not yet wired
(optimize() stands in); mild prior shrinkage in the low-phi tail.

### v0.3 (2026-06-18): analytic NB phi-step + joint sphi
- phi-step is now analytic damped Fisher-scoring (the NB / Poisson-Gamma update),
  vectorised, returning the per-gene Fisher info. FASTER than optimize: default
  (fixed sphi=1) converges in 6 iterations (vs 12), phi Pearson 0.979 / Spearman
  0.996, dEta/dM -0.998. The analytic step is correct + cheaper, and yields the
  Hessian for Laplace SEs.
- Joint sphi (--est-sphi) empirically CONFIRMS sphi is only weakly identified
  from codon data: a plain MAP-EM update RUNS AWAY (phi-modes spread to the
  bounds once the prior loosens to k<1 -- positive feedback). Capping sphi at the
  k=1 boundary (sqrt(trigamma(1))~1.28) stabilises it, but it then PEGS at the cap
  (1.28) vs the native 0.916 -- MAP-EM over-estimates. Tellingly, Spearman stays
  0.994 throughout: the phi RANKING (selection ordering) is well-identified;
  only the linear-scale Pearson degrades as sphi over-spreads. So the SCALE/sphi
  is the weakly-identified piece (exactly Sec 1-2). A proper joint sphi needs the
  NB MARGINAL likelihood (integrate phi out via Gamma-Poisson, estimate the
  dispersion = sphi) or an informative hyperprior; the point-estimate route
  cannot. Default remains fixed-sphi (best phi recovery).

### v0.4-v0.5 (2026-06-18): Laplace SEs + the NB-marginal sphi attempt
- v0.4: Laplace/GLM SEs. Per-gene SE(log phi) = 1/sqrt(Fisher info) from the
  phi-step; per-codon dM/dEta SEs from the vglm vcov. On S288c: SE(log phi)
  median 0.26, SE(dEta) 0.018, SE(dM) 0.026 (matches the published fit SDs).
  Exported to <out>_phi.csv / _csp.csv.
- v0.5: NB / Laplace MARGINAL sphi (--est-sphi). Integrate phi out gene-by-gene
  under Gamma(sphi) via Laplace and profile over sphi. RESULT: the full Laplace
  marginal is UNRELIABLE here -- it decreases monotonically (prefers sphi->0).
  Decomposing it explains why: the data-fit term (multinomial LL) gains only
  +44.7 over sphi 0.4->1.8 and PLATEAUS by sphi~1.0 (the dEta*phi scale
  degeneracy -- beyond ~1.0 extra phi-spread is absorbed by dEta), while the
  Gamma prior-density term dGammaLL falls ~ -950 and the Laplace volume term
  fails to offset it (a known Laplace failure when data are informative
  per-gene). So the marginal is dominated by the prior-density drop, not the
  data. The RELIABLE signal is the data-fit plateau: --est-sphi takes sphi-hat =
  smallest sphi achieving 95% of the data-fit gain = **1.00 on S288c (native
  ~0.92)**, with full recovery preserved (phi r=0.974, dEta/dM -0.999/-0.992).
  CAVEATS: (1) data-fit-plateau is a heuristic, not a likelihood max; a proper
  marginal needs per-gene Gauss-Hermite QUADRATURE (not Laplace) -- done in v0.6,
  which CORRECTS the interpretation below. (2) the final fit must warm-start from
  the clean PC1 init, not the spread high-sphi grid point (which falls into a bad
  basin).
  NOTE: v0.5 guessed the monotone marginal was a Laplace APPROXIMATION failure.
  v0.6 (Gauss-Hermite) shows that guess was WRONG -- see below.

### v0.6 (2026-06-18): Gauss-Hermite marginal -- the monotone marginal is REAL, not a Laplace artifact
- Implemented adaptive Gauss-Hermite quadrature (20 nodes, Golub-Welsch in base R)
  for the per-gene phi marginal: integrate in th=log(phi), centre each gene at the
  integrand mode (th.hat + 1/info.th, the Jacobian-corrected phi-MAP), scale by
  1/sqrt(info.th), reweight nodes by exp(x^2). Holds the converged CSP alpha/beta
  fixed (profile marginal). Function gh_logmarginal() in fit.roc.glm.em.R.
- RESULT on S288c: logML_GH tracks logML_Laplace to within ~+1.0..+1.3 log-units
  across the WHOLE sphi grid (0.4..1.8) -- the small, stable, positive offset is
  exactly the right-skew correction expected when the per-gene log-phi posteriors
  are near-Gaussian. BOTH marginals decrease monotonically and peak at the sphi=0.4
  boundary. So the Laplace was NOT artifactual; the monotone "prefers small sphi"
  is a REAL property of this (unpenalized-CSP) marginal.
- WHY (confirmed by decomposition): the Gamma prior-NORMALIZATION term
  k*log(bg)-lgamma(k), summed over G=1000 genes, swings by THOUSANDS of log-units
  across the grid (+6420 at sphi=0.4 down to -1000 at sphi=1.8), while the data-fit
  moves only +44.7 and plateaus by ~1.0. Because the CSP step is UNPENALIZED (no
  prior on dEta), the dEta*phi ridge lets the model rescale phi to fit any sphi
  equally well, so nothing in the data counteracts the prior-volume Occam gradient
  -> the marginal slides to small sphi. This is the dEta*phi scale degeneracy seen
  directly through a proper marginal.
- IMPLICATION: an informative sphi hyperprior N(1.4,0.5) is TOO WEAK to fix this
  (its ~4/unit^2 curvature cannot offset a ~750/unit marginal-likelihood slope).
  Native MCMC identifies sphi~0.92 because it has a PRIOR on dEta/dM that BREAKS
  the ridge -- but that is NOT a fix to emulate: native's sphi is set by ITS
  arbitrary CSP prior SD, so it has the SAME identifiability problem (its 0.92 is
  not a true value). An L2/ridge prior on dEta would just MOVE the arbitrary gauge
  choice from "fix sphi" to "fix tau_E" (sphi-hat is a monotone function of tau_E).
  RESOLUTION: FIX sphi (~1.4) explicitly -- the honest form of the gauge choice --
  for production GLM-EM fits; the phi RANKING and dEta/dM recovery are unaffected
  by the sphi choice (r=0.97 / |0.99| throughout). The ridge prior is NOTED but
  not pursued (see Open threads); only an EXTERNALLY-informed dEta prior (Sec 3b)
  would be non-circular.
- So the corrected count is: MAP-EM, the data-fit plateau, the Laplace marginal AND
  the proper GH marginal all agree sphi is weakly/not identified from codon data
  under an unpenalized CSP -- and GH pins down the MECHANISM (prior-volume vs a flat
  dEta*phi ridge), pointing to a dEta prior as the structural fix.

### v0.7 (2026-06-19): code-review robustness hardening (no method change)
Review-driven fixes; bit-for-bit identical on S288c (phi r=0.974, dEta/dM
-0.999/-0.992). Row-stable log-sum-exp in multinom_ll/loglik_total (were overflow-
prone, unlike the per-gene version); info.th now recomputed at the FINAL phi-step th
(was a pre/post-step mix when the last Newton step is non-trivial); per-AA vglm
wrapped in tryCatch + a column-count assertion (a sparse/aliased AA is skipped, not
silently recycled into mis-assigned dM/dEta) -- matters for sparse genomes (Loki
bins), not dense yeast; sphi.plateau falls back to the data-fit argmax instead of NA
when mll.v is non-increasing; GH parabola vertex clamped to its bracket; CLI hardened
(--help/--dry-run/--verbose, input-file checks, validated --sphi/--est-sphi parsing,
getarg trailing-flag guard + --flag=value).

## Handling sphi: fixed vs soft prior, and centering (esp. in Stan)

Practical consequence of the identifiability finding for HOW to put sphi into a
fit. `N(1.4, 0.01)` is just the sigma_p -> 0 limit of fixing, approached -- the
question is what it costs to keep sphi as a parameter while the dEta*phi
degeneracy is unbroken.

**The number is nearly the same either way.** With the degeneracy unbroken the
sphi marginal likelihood is a steep, near-linear DOWNSLOPE (~ -120 log-units per
unit sphi near 1.4 on S288c; no interior peak -- the v0.6 GH result). So a soft
prior's mode lands at approximately

  sphi_post  ~=  1.4 + slope * sigma_p^2  =  1.4 - 120 * sigma_p^2

  sigma_p=0.01 -> ~1.39 (SD ~0.01)   sigma_p=0.05 -> ~1.10   sigma_p=0.5 -> collapses.

So a soft prior loose enough to matter is dragged toward the small-sphi artifact,
and one tight enough not to be is just fixing with an extra parameter. No sweet
spot -- the steep slope makes sigma_p hypersensitive.

**In the GLM-EM / native MCMC: just FIX it.** Fixing rigidly pins the gauge ->
the dEta scale is exactly determined, the worst-conditioned eigen-direction of the
posterior (the sphi<->dEta ridge, the source of native's slow phi-dEta mixing) is
REMOVED, there is no sigma_p to tune, and the gauge choice is honest/explicit
(a reported "1.39 +/- 0.01" is just the prior echoed back). For the cross-species
goal, fixing the SAME sphi for every genome is the cleanest shared gauge, so dEta
ratios carry the Ne*T*<phi_rate> signal ([[dEta-as-Ne-signal]]).

**In Stan/HMC the advantage of fixing is BIGGER, and it is about geometry.** The
sphi<->dEta scale degeneracy is a curved (Neal's-funnel) direction, and HMC is
acutely sensitive to funnels: one global step size + diagonal mass matrix cannot
handle the stiff narrow neck and wide mouth at once -> DIVERGENCES, low E-BFMI,
sphi and the dEta scale wandering with high R-hat. This is the pathology the
earlier Stan ROC runs hit.
- **Pin sphi -> pass it as `data`.** Removes the funnel-generating dimension
  entirely; cleanest geometry, zero divergences from that source. Strictly better
  than `N(1.4,0.01)`, which KEEPS the dimension and a stiff direction (curvature
  1/0.01^2 = 1e4 sitting next to O(1) directions), forcing the global step size
  down and slowing every other parameter -- for no informational gain. A tiny soft
  prior is the WORST of the three options in HMC.
- **The real lever is the PARAMETERIZATION, not the prior width.** If sphi is to
  be a parameter at all, NON-CENTER phi: `phi_raw ~ normal(0,1);
  log_phi = mphi + sphi*phi_raw`. This decouples each phi from sphi and tames the
  funnel -- exactly what flipped [[mphi1-implicit-s288c-trial-results]] from
  E-BFMI 0.20 (centered) to 0.60 (non-centered). With the centered form even a
  tight prior can diverge; with NC a MODERATE prior is samplable.

**Rule of thumb.** While the degeneracy is unbroken: fix sphi (~1.4) -- as `data`
in Stan, as a constant in the GLM-EM/native. Once a dEta-scale prior (the #14
follow-up) or genuinely informative phi data BREAKS the degeneracy, the sphi
likelihood gains a real interior peak; THEN switch to a moderate soft prior
`N(1.4, 0.1-0.2)` on a non-centered phi (ideally hierarchical across genomes, to
TEST rather than assert sphi conservation). See [[sphi-prior-default]].

## Open threads / next
- **dEta/dM ridge prior in the CSP step (issue #14) -- NOTED, NOT PURSUED.** An
  L2 prior on dEta with SD tau_E would break the dEta*phi degeneracy and make sphi
  "identifiable", BUT the SD is arbitrary and sphi-hat is just a monotone function
  of tau_E -- so this only RELOCATES the gauge choice from "fix sphi" to "fix
  tau_E", it does not discover a true sphi. Matching native's sphi~0.92 is NOT a
  target: native's value is itself set by ITS (arbitrary) CSP prior SD and suffers
  the same identifiability problem. So fixing sphi explicitly is the honest form of
  the same gauge choice; the arbitrary ridge is not worth building. The ONLY
  non-circular version is an EXTERNALLY-INFORMED dEta prior (Sec 3b: tRNA gene
  copy number on the 4/6-codon AAs, 2-codon AAs as the null) -- that injects
  independent information rather than an arbitrary scale, and is the only route
  worth exploring later.
- s288c plot: logit(codon freq) vs φ for real 2-codon AAs -- line + the
  unsaturated high-φ climb (Sec 6/7). Running.
- Simulation (Sec 1/2): true sphi=1.5 (median gauge), realistic AT-biased dM,
  dEta dialed weak->strong, recover sphi vs dEta. HELD pending the one
  load-bearing choice: WIDE sphi prior (so data, not prior, set recovered sphi).
- gBGC/recombination (Sec 4): GC3-residual vs recombination on s288c (Mancera
  2008 / population rho), once a map is ingested.
- EM/coordinate-ascent ROC via the Gamma φ-step (Sec 8) -- a cheaper fitter.
