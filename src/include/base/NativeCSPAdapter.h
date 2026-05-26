#ifndef NATIVE_CSP_ADAPTER_H
#define NATIVE_CSP_ADAPTER_H

/* ============================================================================
 * NativeCSPAdapter -- the in-house CSP adaptive proposal-width scheme.
 *
 * Wraps the hand-tuned logic that has shipped with RibModelFramework since
 * the early days (originally Cedric's code, later edits by Mike with
 * commentary from Alex).  Default adapter; chosen when no scheme is
 * explicitly requested from R or YAML.
 *
 * Algorithm (per AA per fire):
 *
 *   d = covarianceMatrix.getNumVariates()
 *   (lo, hi) = targetBandFor(d)              # dim-dependent target band
 *   if acceptanceLevel < lo:
 *       std_csp[k] *= (1 - aggressiveness)    for k in [aaStart, aaEnd)
 *       covarianceMatrix *= (1 - aggressiveness)
 *       blend toward sample cov (prevWeight * prev + (1 - prevWeight) * sample_cov)
 *       choleskyDecomposition()
 *   else if acceptanceLevel > hi:
 *       std_csp[k] *= (1 + aggressiveness)
 *       covarianceMatrix *= (1 + aggressiveness)
 *       blend toward sample cov  (symmetric shape update, 2026-05-20)
 *       choleskyDecomposition()
 *   else:
 *       no change (and no Cholesky)
 *
 * The target band is dimensionality-dependent (Roberts-Gelman-Gilks 1997
 * / Roberts & Rosenthal 2001) instead of a flat [0.225, 0.325].  The
 * per-AA joint proposal dimension is
 *   d = numCodons * (numMutationCategories + numSelectionCategories)
 * which for single-mixture ROC is 2*(n_codons - 1) per AA:
 *
 *   AA class             codons  d   optimal AR   default band (bw=0.05)
 *   2-codon (C,D,E,...)  2      2   0.35         [0.30, 0.40]
 *   3-codon (Ile)        3      4   0.27*        [0.22, 0.32]
 *   4-codon (A,G,...)    4      6   0.27         [0.22, 0.32]
 *   6-codon (L,R)        6      10  0.234        [0.19, 0.28]
 *
 * *Ile target merged from 0.28 to 0.27 (same as 4-codon): the d=4/d=6
 * gap is only 0.01, below any practical half-width, so keeping separate
 * targets creates unavoidable band overlap when bw < 0.005.
 *
 * Background: the legacy flat [0.225, 0.325] band was correct for d~4-5
 * but mis-targeted d=2 AAs (their optimum is 0.35, well above the high
 * edge), so R1 adaptation drove their proposal too wide.  Observed
 * empirically 2026-05-21 on the 02v.5-chunked sweep where 2-codon AAs
 * sat at AR 0.50-0.57 in the post-fix-cov frozen R2.
 *
 * The scale factors are user-tunable via `aggressiveness` a in (0, 1):
 *   adjustFactorLow  = 1 - a
 *   adjustFactorHigh = 1 + a
 * Recommended values: 0.1 (gentle), 0.2 (default, == legacy 0.8/1.2),
 * 0.3 (aggressive).
 *
 * The band half-width is user-tunable via `band.half.width` bw in (0, 0.5):
 *   band = [target - bw, target + bw]
 * Default 0.05 preserves the original +-0.05 calibration.  Tighter values
 * (e.g. 0.015) enforce stricter AR targeting; looser values tolerate more
 * deviation before a cov update fires.
 *
 * The cov-blend mixture is user-tunable via `prev.weight` (alias prevWeight)
 * w in (0, 1):
 *   covarianceMatrix <- w * covarianceMatrix + (1 - w) * sample_cov
 * Default 0.6 preserves the legacy 0.6 / 0.4 blend.  Smaller w => faster
 * cov-shape adaptation (more weight on the most recent sample cov) at the
 * cost of higher shape variance; larger w => smoother shape estimate.
 * Independent of `aggressiveness`: the two control distinct facets of the
 * adapter (scale vs shape).
 * ============================================================================ */

#include "CSPAdaptationStrategy.h"
#include <utility>

class NativeCSPAdapter : public CSPAdaptationStrategy {
public:
    // aggressiveness in (0, 1); default 0.2 preserves legacy 0.8/1.2 behavior.
    // prevWeight in (0, 1); default 0.6 preserves legacy 0.6/0.4 cov blend.
    // bandHalfWidth in (0, 0.5); default 0.05 preserves original +-0.05 bands.
    explicit NativeCSPAdapter(double aggressiveness  = 0.2,
                              double prevWeight      = 0.6,
                              double bandHalfWidth   = 0.05);
    ~NativeCSPAdapter() override = default;

    void update(const CSPAdaptContext& ctx) override;
    std::string name() const override { return "native"; }

    double getAggressiveness()  const { return aggressiveness; }
    double getPrevWeight()      const { return prevWeight; }
    double getBandHalfWidth()   const { return bandHalfWidth; }

    std::unique_ptr<CSPAdaptationStrategy> clone() const override {
        return std::unique_ptr<CSPAdaptationStrategy>(
            new NativeCSPAdapter(aggressiveness, prevWeight, bandHalfWidth));
    }

private:
    // Dimension-dependent target acceptance band.  d is the per-AA joint
    // proposal dimension (numCodons * (numMutCat + numSelCat); for
    // single-mixture ROC that is 2*(n_codons - 1)).  Target centers from
    // Roberts-Gelman-Gilks 1997 / Roberts & Rosenthal 2001; half-width
    // is user-configurable via bandHalfWidth (default 0.05).
    // Ile (d=4) merged with 4-codon target (d=6) at 0.27: gap is only 0.01.
    std::pair<double, double> targetBandFor(unsigned d) const {
        double target;
        if      (d <= 1) target = 0.44;   // d=1   target ~0.44
        else if (d == 2) target = 0.35;   // d=2   target ~0.35
        else if (d == 3) target = 0.32;   // d=3   target ~0.32 (rare in ROC)
        else if (d <= 6) target = 0.27;   // d=4,5,6 merged target ~0.27
        else             target = 0.234;  // d>=7  target ~0.234
        return std::make_pair(
            std::max(0.0, target - bandHalfWidth),
            std::min(1.0, target + bandHalfWidth));
    }

    double aggressiveness;       // a in (0, 1)
    double adjustFactorLow;      // 1 - a
    double adjustFactorHigh;     // 1 + a
    double prevWeight;           // w in (0, 1); legacy 0.6
    double bandHalfWidth;        // bw in (0, 0.5); default 0.05
};

#endif // NATIVE_CSP_ADAPTER_H
