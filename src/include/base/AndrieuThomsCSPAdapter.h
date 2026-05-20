#ifndef ANDRIEU_THOMS_CSP_ADAPTER_H
#define ANDRIEU_THOMS_CSP_ADAPTER_H

/* ============================================================================
 * AndrieuThomsCSPAdapter -- continuous Robbins-Monro adaptation on log(SD).
 *
 * Algorithm 4 from Andrieu, C. and Thoms, J. (2008), "A tutorial on adaptive
 * MCMC", Statistics and Computing 18:343-373.  Per AA, on each adapt fire:
 *
 *     log(std_csp[k])_{t+1} = log(std_csp[k])_t + gamma_t * (accept_t - target)
 *     gamma_t                = c / (t + t0)^alpha
 *
 * where t is the per-AA adapt-fire count (0, 1, 2, ...).  The diminishing
 * step schedule satisfies the conditions of Andrieu-Thoms Theorem 2 for
 * alpha in (0.5, 1.0].
 *
 * Differences from the in-house "native" scheme:
 *   - continuous update (multiplicative on SD via exp(delta))
 *     vs discrete x0.8 / x1.2 multiplicative factor
 *   - single target acceptance (0.234, Gelman optimal-d default)
 *     vs threshold band [0.225, 0.325]
 *   - magnitude shrinks via gamma_t -> 0
 *     vs constant 20% per fire
 *   - covariance matrix structure NOT touched (only scalar SD adapts)
 *     vs cov-matrix reblend toward sample cov when acceptance is low
 *
 * Tunable parameters with documented bounds (enforced by ctor):
 *   target in (0, 1)
 *   alpha  in (0.5, 1.0]      -- diminishing-adaptation theorem
 *   c      > 0                 -- initial step size
 *   t0     >= 0                -- offset, avoids huge early steps
 *
 * See docs/csp-adaptation-api.md for the R-facing API contract.
 * ============================================================================ */

#include "CSPAdaptationStrategy.h"
#include <vector>

class AndrieuThomsCSPAdapter : public CSPAdaptationStrategy {
public:
    AndrieuThomsCSPAdapter(double target, double alpha, double c, double t0);
    ~AndrieuThomsCSPAdapter() override = default;

    void update(const CSPAdaptContext& ctx) override;
    std::string name() const override { return "andrieu_thoms"; }
    std::unique_ptr<CSPAdaptationStrategy> clone() const override;

    // Accessors for tests / diagnostics
    double getTarget() const { return target; }
    double getAlpha()  const { return alpha;  }
    double getC()      const { return c;      }
    double getT0()     const { return t0;     }

private:
    double target;
    double alpha;
    double c;
    double t0;
    // Per-AA adapt-fire counter; lazily sized in update() on first use
    // per AA index.  Cloned alongside hyperparameters via clone().
    std::vector<unsigned> aaFireCount;
};

#endif // ANDRIEU_THOMS_CSP_ADAPTER_H
