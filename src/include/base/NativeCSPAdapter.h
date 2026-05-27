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
 *   aggressiveness     (0, 1)    scale factor (1-a low, 1+a high); default 0.2
 *   prev.weight        (0, 1)    cov-blend weight on prior cov; default 0.6
 *   ar.band.half.width (0, 0.5)  half-width of the AR target band; default 0.05
 *
 * The per-dimension target AR values are theory-driven (Roberts-Gelman-Gilks
 * 1997 / Roberts & Rosenthal 2001) and not user-tunable -- use targetARFor(d).
 * band = [targetARFor(d) - ar.band.half.width, targetARFor(d) + ar.band.half.width]
 * Tighter values (e.g. 0.015) make the adapter fire more often; looser values
 * tolerate more AR deviation before a cov update fires.
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
    static constexpr double kDefaultARTarget2codon   = 0.35;   // RGG 1997 d=2
    static constexpr double kDefaultARTarget4codon   = 0.27;   // d=4,6 merged
    static constexpr double kDefaultARTarget6codon   = 0.234;  // RGG 1997 d>=7
    static constexpr double kDefaultARBandHalfWidth  = 0.05;

    explicit NativeCSPAdapter(
        double aggressiveness   = kDefaultAggressiveness,
        double prevWeight       = kDefaultPrevWeight,
        double arTarget2codon   = kDefaultARTarget2codon,
        double arTarget4codon   = kDefaultARTarget4codon,
        double arTarget6codon   = kDefaultARTarget6codon,
        double arBandHalfWidth  = kDefaultARBandHalfWidth);
    ~NativeCSPAdapter() override = default;

    void update(const CSPAdaptContext& ctx) override;
    std::string name() const override { return "native"; }

    double getAggressiveness()   const { return aggressiveness; }
    double getPrevWeight()       const { return prevWeight; }
    double getARTarget2codon()   const { return arTarget2codon; }
    double getARTarget4codon()   const { return arTarget4codon; }
    double getARTarget6codon()   const { return arTarget6codon; }
    double getARBandHalfWidth()  const { return arBandHalfWidth; }

    std::unique_ptr<CSPAdaptationStrategy> clone() const override {
        return std::unique_ptr<CSPAdaptationStrategy>(
            new NativeCSPAdapter(aggressiveness, prevWeight,
                                 arTarget2codon, arTarget4codon, arTarget6codon,
                                 arBandHalfWidth));
    }

private:
    // Returns the optimal AR for the codon class implied by proposal dimension d.
    // Targets are user-configurable (YAML: ar.target.2codon etc.); defaults from
    // Roberts-Gelman-Gilks 1997 / Roberts & Rosenthal 2001.
    // Ile (3-codon, d=4) merged with 4-codon group: gap of 0.01 is below any
    // practical half-width.  d<=1 and d==3 are hardcoded (unused in ROC/PANSE).
    double targetARFor(unsigned d) const {
        if      (d <= 1) return 0.44;          // unused in standard ROC
        else if (d == 2) return arTarget2codon;
        else if (d == 3) return 0.32;          // unused in standard ROC
        else if (d <= 6) return arTarget4codon; // Ile (d=4) + 4-codon (d=6)
        else             return arTarget6codon; // 6-codon Leu/Arg (d=10)
    }

    // Returns [targetARFor(d) - arBandHalfWidth, targetARFor(d) + arBandHalfWidth].
    std::pair<double, double> targetBandFor(unsigned d) const {
        double t = targetARFor(d);
        return std::make_pair(std::max(0.0, t - arBandHalfWidth),
                              std::min(1.0, t + arBandHalfWidth));
    }

    double aggressiveness;    // a in (0,1)
    double adjustFactorLow;   // 1 - a
    double adjustFactorHigh;  // 1 + a
    double prevWeight;        // w in (0,1)
    double arTarget2codon;    // optimal AR for 2-codon AAs
    double arTarget4codon;    // optimal AR for Ile + 4-codon AAs
    double arTarget6codon;    // optimal AR for 6-codon AAs (Leu, Arg)
    double arBandHalfWidth;   // half-width of [target-bw, target+bw]
};

#endif // NATIVE_CSP_ADAPTER_H
