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
 *   if acceptanceLevel < 0.225:
 *       std_csp[k] *= 0.8    for k in [aaStart, aaEnd)
 *       covarianceMatrix *= 0.8
 *       blend toward sample cov (0.6 * prev + 0.4 * sample_cov)
 *       choleskyDecomposition()
 *   else if acceptanceLevel > 0.325:
 *       std_csp[k] *= 1.2
 *       covarianceMatrix *= 1.2
 *       choleskyDecomposition()
 *   else:
 *       no change (and no Cholesky)
 *
 * Constants are not user-tunable -- this is the "native" scheme as it
 * has historically existed.  Tunable schemes live in their own classes
 * (e.g. AndrieuThomsCSPAdapter).
 * ============================================================================ */

#include "CSPAdaptationStrategy.h"

class NativeCSPAdapter : public CSPAdaptationStrategy {
public:
    NativeCSPAdapter() = default;
    ~NativeCSPAdapter() override = default;

    void update(const CSPAdaptContext& ctx) override;
    std::string name() const override { return "native"; }

    std::unique_ptr<CSPAdaptationStrategy> clone() const override {
        // NativeCSPAdapter has no per-instance mutable state (only static
        // constexpr thresholds); a fresh instance is equivalent to a deep
        // copy of any other instance.
        return std::unique_ptr<CSPAdaptationStrategy>(new NativeCSPAdapter());
    }

private:
    // Hard-coded historical constants; see class doc.
    static constexpr double acceptanceTargetLow  = 0.225;
    static constexpr double acceptanceTargetHigh = 0.325;
    static constexpr double adjustFactorLow      = 0.8;
    static constexpr double adjustFactorHigh     = 1.2;
};

#endif // NATIVE_CSP_ADAPTER_H
