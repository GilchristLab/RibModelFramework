/* ============================================================================
 * CSPAdaptationFactory.cpp -- name-driven CSPAdaptationStrategy constructor.
 *
 * See header for the dispatch table and contract.  All errors raised as
 * std::invalid_argument so layer-3 validation surfaces uniformly to both
 * the Rcpp seam (translated to R errors by Parameter::setCSPAdaptationScheme)
 * and direct C++ callers.
 * ============================================================================ */

#include "include/base/CSPAdaptationFactory.h"
#include "include/base/NativeCSPAdapter.h"
#include "include/base/AndrieuThomsCSPAdapter.h"
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>

static void requireExactKeys(const std::map<std::string, double>& params,
                             const std::vector<std::string>& expected,
                             const std::string& scheme_name)
{
    // Missing keys
    for (const std::string& k : expected) {
        if (params.find(k) == params.end()) {
            std::ostringstream o;
            o << "scheme '" << scheme_name << "' requires param '" << k
              << "' but it was not provided";
            throw std::invalid_argument(o.str());
        }
    }
    // Extra keys
    for (const auto& kv : params) {
        if (std::find(expected.begin(), expected.end(), kv.first) == expected.end()) {
            std::ostringstream o;
            o << "scheme '" << scheme_name << "' got unexpected param '" << kv.first
              << "'; allowed: ";
            for (size_t i = 0; i < expected.size(); ++i) {
                o << expected[i] << (i + 1 == expected.size() ? "" : ", ");
            }
            throw std::invalid_argument(o.str());
        }
    }
}

std::unique_ptr<CSPAdaptationStrategy> makeCSPAdapter(
    const std::string& name,
    const std::map<std::string, double>& params)
{
    if (name == "native") {
        if (!params.empty()) {
            std::ostringstream o;
            o << "scheme 'native' takes no params; got " << params.size()
              << " (e.g. '" << params.begin()->first << "')";
            throw std::invalid_argument(o.str());
        }
        return std::unique_ptr<CSPAdaptationStrategy>(new NativeCSPAdapter());
    }

    if (name == "andrieu_thoms") {
        requireExactKeys(params, {"target", "alpha", "c", "t0"}, "andrieu_thoms");
        return std::unique_ptr<CSPAdaptationStrategy>(
            new AndrieuThomsCSPAdapter(
                params.at("target"),
                params.at("alpha"),
                params.at("c"),
                params.at("t0")));
    }

    std::ostringstream o;
    o << "unknown CSP adaptation scheme: '" << name
      << "'; known: native, andrieu_thoms";
    throw std::invalid_argument(o.str());
}
