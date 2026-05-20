/* ============================================================================
 * NativeCSPAdapter.cpp -- implementation of the in-house adaptive scheme.
 *
 * Code moved verbatim from Parameter::adaptCodonSpecificParameterProposalWidth
 * inner loop.  No behavior change vs pre-refactor; a bit-identical regression
 * test against HEAD guards this.
 * ============================================================================ */

#include "include/base/NativeCSPAdapter.h"

void NativeCSPAdapter::update(const CSPAdaptContext& ctx) {
    // Match the pre-refactor branch structure exactly: the off-target gate
    // wraps everything, including the final Cholesky.  When on-target,
    // the cov matrix is left untouched and NOT re-decomposed.
    if (!(ctx.acceptanceLevel < acceptanceTargetLow
          || ctx.acceptanceLevel > acceptanceTargetHigh)) {
        return;
    }

    double adjustFactor = 1.0;

    if (ctx.acceptanceLevel < acceptanceTargetLow) {
        adjustFactor = adjustFactorLow;

        // Update cov matrix toward the sample cov over the adapt window.
        //
        // Historical context (preserved from Parameter.cpp comments at HEAD):
        // In Cedric's original code, this fired only when acceptance was low.
        // Mike updated it to fire every off-target step.  Alex noted that
        // this drove acceptance into 0.4-0.5 for some AAs.  The current
        // branch (low-only) was restored in commit ec63bb21a1e9 (2016).
        CovarianceMatrix covcurr(ctx.covarianceMatrix.getNumVariates());
        covcurr.calculateSampleCovariance(
            *ctx.traces.getCodonSpecificParameterTrace(),
            ctx.aa, ctx.samples, ctx.lastIteration);
        CovarianceMatrix covprev = ctx.covarianceMatrix;
        covprev = (covprev * 0.6);
        covcurr = (covcurr * 0.4);
        ctx.covarianceMatrix = covprev + covcurr;
    } else /* acceptanceLevel > acceptanceTargetHigh */ {
        adjustFactor = adjustFactorHigh;
    }

    if (adjustFactor != 1.0) {
        for (unsigned k = ctx.aaStart; k < ctx.aaEnd; k++) {
            ctx.std_csp[k] *= adjustFactor;
        }
        ctx.covarianceMatrix *= adjustFactor;
    }

    // Re-decompose cov for the proposal generator (matches pre-refactor:
    // only fires when we entered the off-target branch).
    ctx.covarianceMatrix.choleskyDecomposition();
}
