#ifndef VIHOLA_2012_CSP_ADAPTER_H
#define VIHOLA_2012_CSP_ADAPTER_H

/* ============================================================================
 * Vihola2012CSPAdapter -- Robust Adaptive Metropolis (RAM) per Vihola 2012.
 *
 * Vihola, M. (2012), "Robust adaptive Metropolis algorithm with coerced
 * acceptance rate", Statistics and Computing 22:997-1008, Algorithm 1.
 *
 * Per-AA, per-MH-step Cholesky rank-1 update of the proposal covariance to
 * drive the per-step acceptance rate toward a fixed target alpha_star.
 * Unlike NativeCSPAdapter (windowed multiplicative scale) and
 * AndrieuThomsCSPAdapter (scalar SD only), VAM adapts the proposal SHAPE.
 *
 * Algorithm (per AA, per MH step t):
 *   eta_t   = min(1, d * (t+1)^(-gamma))             # diminishing step
 *   sigma_t = sign(alpha_t - alpha_star) * eta_t * |alpha_t - alpha_star|
 *   v_t     = proposal direction applied at step t   # = L^T Z_t in this code
 *   C_{t+1} = C_t + sigma_t / |Z_t|^2 * v_t v_t^T    # rank-1 cov update
 *   L_{t+1} = chol(C_{t+1})                           # refreshed Cholesky
 *
 * The adapter is called once per adapt-fire (not per MH step) and consumes
 * the per-step buffers populated by Parameter::pushStepZ /
 * Model::recordCSPStepAlpha during the fire window.  The update is done in
 * cov-matrix space using rank-1 outer products, then the existing
 * CovarianceMatrix::choleskyDecomposition() refreshes the Cholesky factor.
 *
 * Within-window batching note: Vihola's per-step algorithm uses S_i (the
 * Cholesky factor at step i) when applying step i's update.  Here we keep
 * the L used by the proposal path fixed to L_0 (start of window) across
 * the W steps in the window, then apply all W rank-1 updates in
 * cov-matrix space, then refresh L for the next window.  Equivalent to
 * Vihola's algorithm with adaptive_width=1; an approximation otherwise.
 * For typical adaptive_width=20 the divergence is small in practice.
 *
 * Tunable parameters (enforced by ctor):
 *   target in (0, 1)             -- coerced acceptance rate (alpha_star)
 *   gamma  in (0.5, 1.0]         -- step-size decay; diminishing-adaptation
 *
 * See RAM_DESIGN.md in s.cerevisiae/adapter.dev/ for the design rationale
 * and the plumbing trade-offs (per-step Z + alpha capture in Parameter).
 * ============================================================================ */

#include "CSPAdaptationStrategy.h"
#include <vector>

class Vihola2012CSPAdapter : public CSPAdaptationStrategy {
public:
    Vihola2012CSPAdapter(double target, double gamma);
    ~Vihola2012CSPAdapter() override = default;

    void update(const CSPAdaptContext& ctx) override;
    std::string name() const override { return "vihola_2012"; }
    std::unique_ptr<CSPAdaptationStrategy> clone() const override;

    // Accessors for tests / diagnostics
    double getTarget() const { return target; }
    double getGamma()  const { return gamma;  }

private:
    double target;   // alpha_star, coerced AR target
    double gamma;    // step-size decay exponent
    // Per-AA absolute MH step counter; advances by ctx.samplesSinceLastAdapt
    // * thinning (= adaptationWidth raw iterations) each fire.  Lazily
    // sized in update() on first use per AA index.
    std::vector<unsigned long> aaStepCount;
};

#endif // VIHOLA_2012_CSP_ADAPTER_H
