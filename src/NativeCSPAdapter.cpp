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
#include <stdexcept>
#include <sstream>

NativeCSPAdapter::NativeCSPAdapter(double a, double w,
                                   double tD2, double tD4to6, double tD7plus,
                                   double bw)
    : aggressiveness(a),
      adjustFactorLow(1.0 - a),
      adjustFactorHigh(1.0 + a),
      prevWeight(w),
      arTargetD2(tD2),
      arTargetD4to6(tD4to6),
      arTargetD7plus(tD7plus),
      arBandHalfWidth(bw)
{
    if (!(a > 0.0 && a < 1.0)) {
        std::ostringstream o;
        o << "NativeCSPAdapter: aggressiveness must be in (0, 1); got " << a;
        throw std::invalid_argument(o.str());
    }
    if (!(w > 0.0 && w < 1.0)) {
        std::ostringstream o;
        o << "NativeCSPAdapter: prev.weight must be in (0, 1); got " << w;
        throw std::invalid_argument(o.str());
    }
    if (!(tD2 > 0.0 && tD2 < 1.0)) {
        std::ostringstream o;
        o << "NativeCSPAdapter: target.ar.d2 must be in (0, 1); got " << tD2;
        throw std::invalid_argument(o.str());
    }
    if (!(tD4to6 > 0.0 && tD4to6 < 1.0)) {
        std::ostringstream o;
        o << "NativeCSPAdapter: target.ar.d4to6 must be in (0, 1); got " << tD4to6;
        throw std::invalid_argument(o.str());
    }
    if (!(tD7plus > 0.0 && tD7plus < 1.0)) {
        std::ostringstream o;
        o << "NativeCSPAdapter: target.ar.d7plus must be in (0, 1); got " << tD7plus;
        throw std::invalid_argument(o.str());
    }
    if (!(bw > 0.0 && bw < 0.5)) {
        std::ostringstream o;
        o << "NativeCSPAdapter: ar.band.half.width must be in (0, 0.5); got " << bw;
        throw std::invalid_argument(o.str());
    }
}

void NativeCSPAdapter::update(const CSPAdaptContext& ctx) {
    // Per-AA target band keyed off the joint proposal dimension.  d here
    // is the number of variates the cov matrix proposes jointly
    // (= numCodons * (numMutCat + numSelCat); for single-mixture ROC,
    // 2*(n_codons - 1) per AA).  See targetBandFor doc in the header.
    const unsigned d =
        static_cast<unsigned>(ctx.covarianceMatrix.getNumVariates());
    std::pair<double, double> band = targetBandFor(d);
    const double acceptanceTargetLow  = band.first;
    const double acceptanceTargetHigh = band.second;

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
    covprev = (covprev * prevWeight);
    covcurr = (covcurr * (1.0 - prevWeight));
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
