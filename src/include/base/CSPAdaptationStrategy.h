#ifndef CSP_ADAPTATION_STRATEGY_H
#define CSP_ADAPTATION_STRATEGY_H

/* ============================================================================
 * CSP adaptive proposal-width scheme: abstract interface.
 *
 * One scheme per fit, chosen at runtime from R or YAML.  Each concrete
 * subclass implements update(): per-AA, per-adapt-fire mutation of std_csp
 * and the per-AA proposal covariance matrix.  Parameter holds a
 * unique_ptr<CSPAdaptationStrategy> and delegates the inner-loop math
 * during adaptCodonSpecificParameterProposalWidth().
 *
 * Initial set: NativeCSPAdapter (the in-house 0.8/1.2 scheme, default;
 * bit-identical to pre-refactor behavior) and AndrieuThomsCSPAdapter
 * (Robbins-Monro update on log(std), 2008).
 *
 * See docs/csp-adaptation-api.md for the R-facing API and design
 * rationale (one scheme per fit, dual R+C++ validation, R-led naming).
 * ============================================================================ */

#include "../CovarianceMatrix.h"
#include "Trace.h"
#include <memory>
#include <string>
#include <vector>

/**
 * Per-AA context handed to a strategy on each adapt fire.
 *
 * Mutable refs (std_csp, covarianceMatrix) are the only state the strategy
 * may modify.  Strategy writes only the [aaStart, aaEnd) slice of std_csp.
 */
struct CSPAdaptContext {
    unsigned aaIndex;
    const std::string& aa;
    unsigned aaStart;
    unsigned aaEnd;
    double acceptanceLevel;
    unsigned adaptationWidth;
    unsigned lastIteration;
    unsigned samples;
    std::vector<double>& std_csp;
    CovarianceMatrix& covarianceMatrix;
    Trace& traces;
};

class CSPAdaptationStrategy {
public:
    virtual ~CSPAdaptationStrategy() = default;

    /**
     * One adapt fire, for one AA.  Called only when the outer adapt
     * flag is true (i.e. iteration <= stepsToAdapt).
     */
    virtual void update(const CSPAdaptContext& ctx) = 0;

    /** Lowercase snake_case canonical name (matches R schemes.available()). */
    virtual std::string name() const = 0;

    /**
     * Deep copy of this strategy (including any per-AA internal state).
     * Used by Parameter's copy ctor / operator= since unique_ptr is
     * non-copyable; the alternative (sharing across Parameter copies)
     * would alias mutable state, which is wrong.
     */
    virtual std::unique_ptr<CSPAdaptationStrategy> clone() const = 0;
};

#endif // CSP_ADAPTATION_STRATEGY_H
