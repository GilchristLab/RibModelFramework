/* ============================================================================
 * NativeCSPAdapter.cpp -- implementation of the in-house adaptive scheme.
 *
 * Originally code moved verbatim from
 * Parameter::adaptCodonSpecificParameterProposalWidth inner loop (2016-era
 * scheme; cov-blend fired only inside the low-acceptance branch).  The
 * asymmetric "blend only when low, rescale only when high" structure was
 * found empirically (2026-05-20, sibling session on phi-mixture S288c
 * chains) to leave the cov matrix shape stale for AAs that drift into the
 * high-acceptance band -- ESS for Selection on 4-/6-codon AAs (Arg, Leu,
 * Ala, Gly, Val on A-ending codons specifically) was 5-10x lower than for
 * 2-codon AAs.  See feat/symmetric-csp-cov-adaptation 541f59e for the
 * minimal-patch evidence; this adapter incorporates the symmetry as the
 * permanent NativeCSPAdapter behaviour.
 *
 * Current structure:
 *   1. Gate on off-target acceptance (low OR high) as before.
 *   2. Blend cov toward sample cov in BOTH branches -- shape updates fire
 *      whenever the chain leaves the target band, regardless of direction.
 *   3. Apply asymmetric scale (0.8 if low, 1.2 if high) -- preserves the
 *      damping behaviour Alex relied on to keep high-acceptance AAs from
 *      running away.
 *   4. Re-decompose cov for the proposal generator.
 * ============================================================================ */

#include "include/base/NativeCSPAdapter.h"

void NativeCSPAdapter::update(const CSPAdaptContext& ctx) {
    // Off-target gate.  When acceptance is inside [low, high] the chain is
    // mixing well at its current proposal; do nothing (no scale, no shape
    // change, no Cholesky).
    const bool below = ctx.acceptanceLevel < acceptanceTargetLow;
    const bool above = ctx.acceptanceLevel > acceptanceTargetHigh;
    if (!below && !above) return;

    // 1) Shape update: blend the current cov toward the sample cov over the
    //    most recent adaptation window.  Fires for BOTH low- and high-
    //    acceptance off-target conditions (symmetric in direction).
    CovarianceMatrix covcurr(ctx.covarianceMatrix.getNumVariates());
    covcurr.calculateSampleCovariance(
        *ctx.traces.getCodonSpecificParameterTrace(),
        ctx.aa, ctx.samplesSinceLastAdapt, ctx.lastSample);
    CovarianceMatrix covprev = ctx.covarianceMatrix;
    covprev = (covprev * 0.6);
    covcurr = (covcurr * 0.4);
    ctx.covarianceMatrix = covprev + covcurr;

    // 2) Asymmetric scale update: shrink on low acceptance, grow on high.
    //    Direction-asymmetric on purpose (Alex's observation: the only-when-
    //    high arm needs the 1.2 multiplier to damp runaway acceptance).
    const double adjustFactor = below ? adjustFactorLow : adjustFactorHigh;
    for (unsigned k = ctx.aaStart; k < ctx.aaEnd; k++) {
        ctx.std_csp[k] *= adjustFactor;
    }
    ctx.covarianceMatrix *= adjustFactor;

    // 3) Re-decompose for the proposal generator.
    ctx.covarianceMatrix.choleskyDecomposition();
}
