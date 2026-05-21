#ifndef CSP_ADAPTATION_FACTORY_H
#define CSP_ADAPTATION_FACTORY_H

/* ============================================================================
 * Factory: construct a CSPAdaptationStrategy from a name + params map.
 *
 * Dispatch entry point used by Parameter::setCSPAdaptationScheme (Rcpp
 * setter) and by standalone C++ callers.  Validates the param-map shape
 * (layer 3 of the dual-validation scheme); the typed strategy ctor then
 * performs range checks (layer 4).
 *
 * Throws std::invalid_argument on:
 *   - unknown scheme name
 *   - unexpected params (e.g. unknown keys for the chosen scheme)
 *   - missing required params
 *   - range-invalid params (delegated to the strategy ctor)
 *
 * Recognized names (must match R schemes.available()):
 *   "native"        -- params must be empty
 *   "andrieu_thoms" -- params keys: target, alpha, c, t0
 * ============================================================================ */

#include "CSPAdaptationStrategy.h"
#include <map>
#include <string>

std::unique_ptr<CSPAdaptationStrategy> makeCSPAdapter(
    const std::string& name,
    const std::map<std::string, double>& params);

#endif // CSP_ADAPTATION_FACTORY_H
