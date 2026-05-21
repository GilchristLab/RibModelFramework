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
 *       std_csp[k] *= 0.8    for k in [aaStart, aaEnd)
 *       covarianceMatrix *= 0.8
 *       blend toward sample cov (0.6 * prev + 0.4 * sample_cov)
 *       choleskyDecomposition()
 *   else if acceptanceLevel > hi:
 *       std_csp[k] *= 1.2
 *       covarianceMatrix *= 1.2
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
 *   AA class           codons  d   optimal AR   band
 *   2-codon (C,D,E,...) 2      2   0.35         [0.30, 0.40]
 *   3-codon (Ile)       3      4   0.28         [0.23, 0.33]
 *   4-codon (A,G,...)   4      6   0.27         [0.22, 0.32]
 *   6-codon (L,R)       6      10  0.234        [0.19, 0.28]
 *
 * Background: the legacy flat [0.225, 0.325] band was correct for d~4-5
 * but mis-targeted d=2 AAs (their optimum is 0.35, well above the high
 * edge), so R1 adaptation drove their proposal too wide.  Observed
 * empirically 2026-05-21 on the 02v.5-chunked sweep where 2-codon AAs
 * sat at AR 0.50-0.57 in the post-fix-cov frozen R2.
 *
 * The scale factors are user-tunable via a single `aggressiveness`
 * scalar a in (0, 1):
 *   adjustFactorLow  = 1 - a
 *   adjustFactorHigh = 1 + a
 * Recommended values: 0.1 (gentle), 0.2 (default, == legacy 0.8/1.2),
 * 0.3 (aggressive).  Larger a converges faster but with more thrash;
 * smaller a is steadier but slower.  Target band is NOT user-tunable
 * (theory-driven from optimal-AR scaling).
 * ============================================================================ */

#include "CSPAdaptationStrategy.h"
#include <utility>

class NativeCSPAdapter : public CSPAdaptationStrategy {
public:
    // aggressiveness in (0, 1); default 0.2 preserves legacy 0.8/1.2 behavior.
    explicit NativeCSPAdapter(double aggressiveness = 0.2);
    ~NativeCSPAdapter() override = default;

    void update(const CSPAdaptContext& ctx) override;
    std::string name() const override { return "native"; }

    double getAggressiveness() const { return aggressiveness; }

    std::unique_ptr<CSPAdaptationStrategy> clone() const override {
        return std::unique_ptr<CSPAdaptationStrategy>(
            new NativeCSPAdapter(aggressiveness));
    }

private:
    // Dimension-dependent target acceptance band.  d is the per-AA joint
    // proposal dimension (numCodons * (numMutCat + numSelCat); for
    // single-mixture ROC that is 2*(n_codons - 1)).  Calibration:
    // Roberts-Gelman-Gilks 1997 / Roberts & Rosenthal 2001 optimal AR
    // for d-dim Gaussian, +/-0.05 around the optimum to define the band.
    static std::pair<double, double> targetBandFor(unsigned d) {
        if (d <= 1) return std::make_pair(0.39, 0.49);   // d=1   target ~0.44
        if (d == 2) return std::make_pair(0.30, 0.40);   // d=2   target ~0.35
        if (d == 3) return std::make_pair(0.27, 0.37);   // d=3   target ~0.32
        if (d == 4) return std::make_pair(0.23, 0.33);   // d=4   target ~0.28
        if (d <= 6) return std::make_pair(0.22, 0.32);   // d=5,6 target ~0.27
        return        std::make_pair(0.19, 0.28);        // d>=7  target ~0.234
    }

    double aggressiveness;       // a in (0, 1)
    double adjustFactorLow;      // 1 - a
    double adjustFactorHigh;     // 1 + a
};

#endif // NATIVE_CSP_ADAPTER_H
