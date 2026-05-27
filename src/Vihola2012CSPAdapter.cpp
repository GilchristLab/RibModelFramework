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

    std::vector<double>& cov = *C.getCovMatrix();
    if (cov.size() != d * d) return;       // unexpected shape; skip

    // Loop over per-step Z + alpha pairs, applying the rank-1 cov update
    // sequentially with L refreshed AFTER each step.  This mirrors Vihola's
    // per-MH-step algorithm: step i's update sees the cov shrunk/grown by
    // steps 0..i-1, so cumulative downdates cannot drive cov non-PSD beyond
    // the per-step PSD bound (|sigma_i| <= 1).
    //
    // The COST of the in-loop choleskyDecomposition is O(d^3) per step.
    // For d=10 and W=200 steps/window, that's ~200K flops per AA per fire,
    // negligible vs the MH inner loop.
    //
    // The ACTUAL chain used L_0 (start-of-window L) for all W proposals,
    // so this replay produces a slightly different cov than Vihola's true
    // per-step algorithm would.  But the proposal stationarity is set by
    // the actual proposals (L_0); the adapter's job is to evolve L for
    // future windows, and per-step replay does this in a PSD-preserving way.
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

        // v_i = L_i^T Z_i  using the CURRENT (evolving) Cholesky factor.
        // L_i is the choleskyMatrix as updated by the previous step's
        // choleskyDecomposition() call (or the initial L_0 for i=0).
        const std::vector<double>& Lcurr = *C.getCholeskyMatrix();
        if (Lcurr.size() != d * d) break;
        std::vector<double> v(d, 0.0);
        for (std::size_t j = 0; j < d; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < d; ++k) {
                s += Lcurr[k * d + j] * z[k];
            }
            v[j] = s;
        }

        // Rank-1 cov update: cov += (sigma_i / ||Z||^2) * v v^T.
        // PSD bound in EXACT arithmetic: with sigma_i in [-1, +inf) and
        // |Z|^2 = sum z_k^2, this preserves PSD when sigma_i / |Z|^2 *
        // v^T cov^-1 v >= -1.  Since v = L^T Z (with cov = L L^T), we
        // have v^T cov^-1 v = Z^T Z = |Z|^2, so the constraint reduces
        // to sigma_i >= -1.  This holds (eta_i <= 1, |alpha_i - target|
        // < 1 => sigma_i in (-1, +1)).
        //
        // But floating-point arithmetic can drive cov to the PSD boundary
        // (smallest eigenvalue ~ 0), and the next step's L would have
        // near-zero diag entries, causing v to be ill-conditioned.  Defensive
        // rollback: save pre-update cov + L, attempt update + chol, then
        // check chol-diag for non-finite or non-positive; if bad, revert.
        std::vector<double> cov_save = cov;
        std::vector<double> L_save   = *C.getCholeskyMatrix();

        const double scale = sigma_i / z2;
        for (std::size_t r = 0; r < d; ++r) {
            const double vr = v[r];
            for (std::size_t c2 = 0; c2 < d; ++c2) {
                cov[r * d + c2] += scale * vr * v[c2];
            }
        }

        // Refresh L from the updated cov so step i+1 sees the new shape.
        C.choleskyDecomposition();

        // Validate: every L diag entry must be finite and > 0.  In exact
        // arithmetic Vihola's rank-1 update preserves PSD (sigma_i in
        // (-1, +1) and v = L * Z so v^T * cov^-1 * v = |Z|^2; see PR
        // notes).  With the cholesky proposal-direction bug fixed
        // (fix/cholesky-proposal-index, 2026-05-22), this rollback should
        // essentially never fire -- if it does, that signals a
        // floating-point edge case worth investigating, not a routine
        // operating mode.  Keep the rollback as a safety net; do NOT add
        // a scalar-shrink fallback (the original workaround for the
        // cholesky bug, now obsolete).
        const std::vector<double>& Lcheck = *C.getCholeskyMatrix();
        bool bad = false;
        for (std::size_t k = 0; k < d; ++k) {
            double dk = Lcheck[k * d + k];
            if (!std::isfinite(dk) || dk <= 0.0) { bad = true; break; }
        }
        if (bad) {
            cov = cov_save;
            *C.getCholeskyMatrix() = L_save;
        }
    }

    // Advance per-AA step counter by the number of consumed steps so the
    // next fire's eta picks up where this one left off.
    aaStepCount[ctx.aaIndex] += static_cast<unsigned long>(nsteps);

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
