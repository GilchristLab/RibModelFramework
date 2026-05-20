/* ============================================================================
 * AndrieuThomsCSPAdapter.cpp -- implementation of Andrieu-Thoms 2008 A4
 *
 * See header for algorithm reference.  Range checks here mirror those in the
 * R-side AdaptiveScheme.AndrieuThoms() constructor; layer 4 of the dual
 * validation (R layer 2 + C++ layer 4) per docs/csp-adaptation-api.md.
 * ============================================================================ */

#include "include/base/AndrieuThomsCSPAdapter.h"
#include <cmath>
#include <stdexcept>
#include <sstream>

AndrieuThomsCSPAdapter::AndrieuThomsCSPAdapter(double target, double alpha,
                                               double c, double t0)
    : target(target), alpha(alpha), c(c), t0(t0)
{
    auto bad = [](const std::string& msg) {
        throw std::invalid_argument("AndrieuThomsCSPAdapter: " + msg);
    };

    if (!std::isfinite(target) || !(target > 0.0) || !(target < 1.0)) {
        std::ostringstream o; o << "target must be finite and in (0, 1); got " << target;
        bad(o.str());
    }
    if (!std::isfinite(alpha) || !(alpha > 0.5) || !(alpha <= 1.0)) {
        std::ostringstream o; o << "alpha must be finite and in (0.5, 1.0]; got " << alpha;
        bad(o.str());
    }
    if (!std::isfinite(c) || !(c > 0.0)) {
        std::ostringstream o; o << "c must be finite and > 0; got " << c;
        bad(o.str());
    }
    if (!std::isfinite(t0) || !(t0 >= 0.0)) {
        std::ostringstream o; o << "t0 must be finite and >= 0; got " << t0;
        bad(o.str());
    }
}

void AndrieuThomsCSPAdapter::update(const CSPAdaptContext& ctx) {
    if (aaFireCount.size() <= ctx.aaIndex)
        aaFireCount.resize(ctx.aaIndex + 1, 0u);

    double t      = static_cast<double>(aaFireCount[ctx.aaIndex]);
    double gamma  = c / std::pow(t + t0, alpha);
    double delta  = gamma * (ctx.acceptanceLevel - target);
    double scale  = std::exp(delta);

    for (unsigned k = ctx.aaStart; k < ctx.aaEnd; k++) {
        ctx.std_csp[k] *= scale;
    }
    // covarianceMatrix intentionally NOT scaled here; A-T tunes the scalar
    // proposal SD only.  Cholesky of the initial cov stays valid because
    // the proposal mechanism multiplies the decomposed factor by std_csp at
    // sample time.

    aaFireCount[ctx.aaIndex] += 1u;
}

std::unique_ptr<CSPAdaptationStrategy> AndrieuThomsCSPAdapter::clone() const {
    // Deep copy including the per-AA fire counter -- if a Parameter is
    // copied mid-fit, the clone resumes adaptation from the same schedule
    // position rather than restarting at t=0.
    auto p = std::unique_ptr<AndrieuThomsCSPAdapter>(
        new AndrieuThomsCSPAdapter(target, alpha, c, t0));
    p->aaFireCount = aaFireCount;
    return p;
}
