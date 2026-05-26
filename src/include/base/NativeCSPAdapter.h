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
 * User-tunable parameters (YAML keys / R AdaptiveScheme.Native() args):
 *
 *   aggressiveness    (0, 1)     scale factor (1-a low, 1+a high); default 0.2
 *   prev.weight       (0, 1)     cov-blend weight on prior cov; default 0.6
 *   target.ar.d2      (0, 1)     optimal AR for d=2 (2-codon AAs); default 0.35
 *   target.ar.d4to6   (0, 1)     optimal AR for d=4..6 (Ile + 4-codon); default 0.27
 *   target.ar.d7plus  (0, 1)     optimal AR for d>=7 (6-codon L,R); default 0.234
 *   ar.band.half.width (0, 0.5)  half-width of the AR target band; default 0.05
 *
 * band = [target - ar.band.half.width, target + ar.band.half.width]
 * Tighter values (e.g. 0.015) enforce stricter AR targeting so the adapter
 * fires more often; looser values tolerate more deviation before firing.
 *
 * Default values are canonical class constants (kDefault*) used by both the
 * constructor and CSPAdaptationFactory to guarantee a single source of truth.
 *
 * The cov-blend mixture is user-tunable via `prev.weight` w in (0, 1):
 *   covarianceMatrix <- w * covarianceMatrix + (1 - w) * sample_cov
 * Default 0.6 preserves the legacy 0.6 / 0.4 blend.  Smaller w => faster
 * cov-shape adaptation; larger w => smoother shape estimate.
 * ============================================================================ */

#include "CSPAdaptationStrategy.h"
#include <utility>

class NativeCSPAdapter : public CSPAdaptationStrategy {
public:
    // Canonical defaults -- single source of truth for constructor and factory.
    static constexpr double kDefaultAggressiveness   = 0.2;
    static constexpr double kDefaultPrevWeight       = 0.6;
    static constexpr double kDefaultARTargetD2       = 0.35;
    static constexpr double kDefaultARTargetD4to6    = 0.27;
    static constexpr double kDefaultARTargetD7plus   = 0.234;
    static constexpr double kDefaultARBandHalfWidth  = 0.05;

    explicit NativeCSPAdapter(
        double aggressiveness   = kDefaultAggressiveness,
        double prevWeight       = kDefaultPrevWeight,
        double arTargetD2       = kDefaultARTargetD2,
        double arTargetD4to6    = kDefaultARTargetD4to6,
        double arTargetD7plus   = kDefaultARTargetD7plus,
        double arBandHalfWidth  = kDefaultARBandHalfWidth);
    ~NativeCSPAdapter() override = default;

    void update(const CSPAdaptContext& ctx) override;
    std::string name() const override { return "native"; }

    double getAggressiveness()    const { return aggressiveness; }
    double getPrevWeight()        const { return prevWeight; }
    double getARTargetD2()        const { return arTargetD2; }
    double getARTargetD4to6()     const { return arTargetD4to6; }
    double getARTargetD7plus()    const { return arTargetD7plus; }
    double getARBandHalfWidth()   const { return arBandHalfWidth; }

    std::unique_ptr<CSPAdaptationStrategy> clone() const override {
        return std::unique_ptr<CSPAdaptationStrategy>(
            new NativeCSPAdapter(aggressiveness, prevWeight,
                                 arTargetD2, arTargetD4to6, arTargetD7plus,
                                 arBandHalfWidth));
    }

private:
    // Returns [target - arBandHalfWidth, target + arBandHalfWidth] for
    // dimension d.  Three user-tunable groups cover all ROC/PANSE AA classes.
    // d<=1 and d==3 are hardcoded (never fire in standard single-mixture ROC).
    std::pair<double, double> targetBandFor(unsigned d) const {
        double target;
        if      (d <= 1) target = 0.44;         // hardcoded; unused in ROC
        else if (d == 2) target = arTargetD2;
        else if (d == 3) target = 0.32;         // hardcoded; unused in ROC
        else if (d <= 6) target = arTargetD4to6;
        else             target = arTargetD7plus;
        return std::make_pair(
            std::max(0.0, target - arBandHalfWidth),
            std::min(1.0, target + arBandHalfWidth));
    }

    double aggressiveness;    // a in (0,1)
    double adjustFactorLow;   // 1 - a
    double adjustFactorHigh;  // 1 + a
    double prevWeight;        // w in (0,1)
    double arTargetD2;        // optimal AR for d=2
    double arTargetD4to6;     // optimal AR for d in [4,6]
    double arTargetD7plus;    // optimal AR for d>=7
    double arBandHalfWidth;   // bw in (0,0.5)
};

#endif // NATIVE_CSP_ADAPTER_H
