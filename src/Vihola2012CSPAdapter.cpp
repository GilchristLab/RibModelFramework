/* ============================================================================
 * Vihola2012CSPAdapter.cpp -- implementation of Vihola 2012 RAM Algorithm 1.
 *
 * See header for algorithm reference and the within-window batching note.
 * Range checks here mirror those in the R-side AdaptiveScheme.Vihola2012()
 * constructor; layer 4 of the dual validation (R layer 2 + C++ layer 4)
 * per docs/csp-adaptation-api.md.
 * ============================================================================ */

#include "include/base/Vihola2012CSPAdapter.h"
#include "include/CovarianceMatrix.h"
#include <cmath>
#include <sstream>
#include <stdexcept>

Vihola2012CSPAdapter::Vihola2012CSPAdapter(double target, double gamma)
    : target(target), gamma(gamma)
{
    auto bad = [](const std::string& msg) {
        throw std::invalid_argument("Vihola2012CSPAdapter: " + msg);
    };

    if (!std::isfinite(target) || !(target > 0.0) || !(target < 1.0)) {
        std::ostringstream o;
        o << "target must be finite and in (0, 1); got " << target;
        bad(o.str());
    }
    if (!std::isfinite(gamma) || !(gamma > 0.5) || !(gamma <= 1.0)) {
        std::ostringstream o;
        o << "gamma must be finite and in (0.5, 1.0]; got " << gamma;
        bad(o.str());
    }
}


void Vihola2012CSPAdapter::update(const CSPAdaptContext& ctx)
{
    if (aaStepCount.size() <= ctx.aaIndex) {
        aaStepCount.resize(ctx.aaIndex + 1, 0ul);
    }

    const std::vector<std::vector<double>>& zBuf = ctx.stepZBuffer;
    const std::vector<double>&              aBuf = ctx.stepAlphaBuffer;

    // No-op if no steps recorded this window (e.g. AA not in subset).
    if (zBuf.empty() || aBuf.empty()) {
        return;
    }
    // Defensive: Z and alpha buffers must be parallel.  If mismatched,
    // truncate to the shorter -- something went wrong upstream but a
    // partial update is safer than UB.
    const std::size_t nsteps =
        (zBuf.size() < aBuf.size()) ? zBuf.size() : aBuf.size();

    CovarianceMatrix& C = ctx.covarianceMatrix;
    const int dInt = C.getNumVariates();
    if (dInt <= 0) return;
    const std::size_t d = static_cast<std::size_t>(dInt);

    // Read the Cholesky factor L (= what proposeCodonSpecificParameter uses
    // via transformIidNumbersIntoCovaryingNumbers).  Stored row-major n*n;
    // the proposal computes covarying[i] = sum_k L[k*n + i] * Z[k], i.e.
    // the actual proposal direction is v = L^T Z under this convention.
    // We hold L_0 fixed across the window (batch approximation) and update
    // the cov matrix in place; choleskyDecomposition() refreshes L at end.
    std::vector<double> L = *C.getCholeskyMatrix();
    std::vector<double>& cov = *C.getCovMatrix();
    if (L.size() != d * d || cov.size() != d * d) {
        // Storage shape unexpected -- skip rather than corrupt state.
        return;
    }

    // Loop over per-step Z + alpha pairs, accumulating rank-1 cov updates.
    // sigma_i = sign(alpha_i - target) * eta_i * |alpha_i - target|
    //         = eta_i * (alpha_i - target)
    // eta_i   = min(1, d * (t_global + 1)^{-gamma}) with t_global the
    //           absolute MH step count for this AA (advances per step).
    for (std::size_t i = 0; i < nsteps; ++i) {
        const std::vector<double>& z = zBuf[i];
        if (z.size() != d) continue;       // dim mismatch -- skip
        const double alpha_i = aBuf[i];

        // Per-step diminishing step size; +1 to avoid 0^(-gamma) at start.
        const double t_plus_one =
            static_cast<double>(aaStepCount[ctx.aaIndex] + i + 1);
        double eta_i = static_cast<double>(d) * std::pow(t_plus_one, -gamma);
        if (eta_i > 1.0) eta_i = 1.0;
        const double sigma_i = eta_i * (alpha_i - target);

        // ||Z||^2 -- denominator of the rank-1 update direction.
        double z2 = 0.0;
        for (std::size_t k = 0; k < d; ++k) z2 += z[k] * z[k];
        if (!(z2 > 0.0)) continue;          // degenerate Z; skip

        // v_i = L^T Z_i (this code's proposal direction, applied at step i)
        //     = sum_k L[k*d + j] * Z[k]   for j = 0..d-1
        std::vector<double> v(d, 0.0);
        for (std::size_t j = 0; j < d; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < d; ++k) {
                s += L[k * d + j] * z[k];
            }
            v[j] = s;
        }

        // Rank-1 cov update: cov += (sigma_i / ||Z||^2) * v v^T.
        const double scale = sigma_i / z2;
        for (std::size_t r = 0; r < d; ++r) {
            const double vr = v[r];
            for (std::size_t c2 = 0; c2 < d; ++c2) {
                cov[r * d + c2] += scale * vr * v[c2];
            }
        }
    }

    // Advance per-AA step counter by the number of consumed steps so the
    // next fire's eta picks up where this one left off.
    aaStepCount[ctx.aaIndex] += static_cast<unsigned long>(nsteps);

    // Refresh L from the updated cov.  If the rank-1 downdates drove cov
    // toward indefiniteness, choleskyDecomposition will produce NaNs on
    // the diagonal sqrt(negative); we tolerate that here (the next
    // proposal will be NaN-poisoned and a downstream isnan guard will
    // reject it).  Production hardening: clip negative diag of cov to a
    // small epsilon before refresh -- deferred to a Phase-2 robustness pass.
    C.choleskyDecomposition();

    // std_csp is the per-codon proposal SD scalar used by the legacy
    // proposal path BEFORE Cholesky factoring; native/A-T mutate it but
    // VAM owns shape via the cov matrix and treats std_csp as a constant
    // multiplier.  Leave std_csp untouched.
    (void) ctx.std_csp;
    (void) ctx.aaStart;
    (void) ctx.aaEnd;
}


std::unique_ptr<CSPAdaptationStrategy> Vihola2012CSPAdapter::clone() const
{
    auto p = std::unique_ptr<Vihola2012CSPAdapter>(
        new Vihola2012CSPAdapter(target, gamma));
    p->aaStepCount = aaStepCount;
    return p;
}
