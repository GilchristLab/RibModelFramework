## ============================================================================
## adaptiveScheme.R -- R-side API for pluggable CSP adaptive proposal-width
## schemes.
##
## Public interface (see docs/csp-adaptation-api.md for the full design):
##
##   AdaptiveScheme.Native()
##       Returns an S3 object selecting the in-house scheme (the default).
##       Takes no parameters.
##
##   AdaptiveScheme.AndrieuThoms(target = 0.234, alpha = 0.7, c = 1.0, t0 = 10)
##       Returns an S3 object selecting Andrieu & Thoms 2008 Algorithm 4
##       (continuous Robbins-Monro update on log(SD)).  Parameters are
##       range-validated at the R layer (layer 2); the C++ ctor revalidates
##       (layer 4).  Defaults come from the literature: target = 0.234 is
##       Gelman et al's optimum-d acceptance rate; alpha = 0.7 satisfies the
##       diminishing-adaptation theorem (Andrieu & Thoms 2008, Theorem 2);
##       c = 1.0 is a standard initial step; t0 = 10 avoids huge early steps.
##
##   schemes.available()
##       Returns the canonical names suitable for YAML configs and for the
##       C++ factory dispatch.
##
##   is.AdaptiveScheme(x)
##       Predicate; returns TRUE if x was constructed by one of the
##       AdaptiveScheme.* functions.
##
##   print.AdaptiveScheme(x)
##       S3 print method.
##
##   adaptive.scheme.diagnostics(parameter, mcmc = NULL)
##       Accessor for the scheme currently bound to a Parameter; v1 returns
##       just the scheme name (the C++ side does not yet expose per-call
##       traces).  Trace-based diagnostics planned for v2.
##
## Naming: function names use dot-as-namespace (`AdaptiveScheme.Native`,
## `AdaptiveScheme.AndrieuThoms`).  This matches R's S3 dispatch convention
## (`print.AdaptiveScheme` etc.) and makes per-scheme constructors
## discoverable via tab-complete on `AdaptiveScheme.`.
## ============================================================================


#' Construct the native (in-house) CSP adaptive proposal-width scheme.
#'
#' @description The default scheme that has shipped with RibModelFramework
#'   since the early days.  Per-AA, on each adapt fire: if acceptance rate
#'   is below 0.225, multiply per-codon `std_csp` by 0.8 and blend the
#'   proposal covariance toward the sample covariance; if above 0.325,
#'   multiply by 1.2; otherwise no change.  Constants are not user-tunable
#'   in this scheme.
#'
#' @return An object of S3 class `c("AdaptiveScheme.Native", "AdaptiveScheme")`.
#' @seealso [AdaptiveScheme.AndrieuThoms()], [schemes.available()].
#' @examples
#' s <- AdaptiveScheme.Native()
#' print(s)
#' @export
AdaptiveScheme.Native <- function() {
    structure(
        list(scheme = "native", params = list()),
        class = c("AdaptiveScheme.Native", "AdaptiveScheme")
    )
}


#' Construct the Andrieu-Thoms 2008 CSP adaptive proposal-width scheme.
#'
#' @description Continuous Robbins-Monro update on `log(std_csp)`:
#'   `log(std)_{t+1} = log(std)_t + gamma_t * (acceptance_t - target)`
#'   with `gamma_t = c / (t + t0)^alpha`.  See Andrieu & Thoms (2008),
#'   "A tutorial on adaptive MCMC", Statistics and Computing 18:343-373,
#'   Algorithm 4.
#'
#' @param target  Target acceptance rate.  Must be in (0, 1).  Default
#'   0.234 (Gelman et al optimal-d).
#' @param alpha   Step-size decay exponent.  Must be in (0.5, 1.0] for
#'   the diminishing-adaptation theorem to hold.  Default 0.7.
#' @param c       Initial step size (positive).  Default 1.0.
#' @param t0      Step-size schedule offset (non-negative); larger values
#'   shrink early-iteration step magnitude.  Default 10.
#'
#' @return An object of S3 class `c("AdaptiveScheme.AndrieuThoms",
#'   "AdaptiveScheme")`.
#' @seealso [AdaptiveScheme.Native()], [schemes.available()].
#' @examples
#' s <- AdaptiveScheme.AndrieuThoms(target = 0.25, alpha = 0.7)
#' print(s)
#' @export
AdaptiveScheme.AndrieuThoms <- function(target = 0.234, alpha = 0.7,
                                        c = 1.0, t0 = 10) {
    stopifnot(
        is.numeric(target), length(target) == 1L, is.finite(target),
        target > 0, target < 1,
        is.numeric(alpha),  length(alpha)  == 1L, is.finite(alpha),
        alpha > 0.5, alpha <= 1.0,
        is.numeric(c),      length(c)      == 1L, is.finite(c),
        c > 0,
        is.numeric(t0),     length(t0)     == 1L, is.finite(t0),
        t0 >= 0
    )
    structure(
        list(
            scheme = "andrieu_thoms",
            params = list(target = target, alpha = alpha, c = c, t0 = t0)
        ),
        class = c("AdaptiveScheme.AndrieuThoms", "AdaptiveScheme")
    )
}


#' Canonical names of CSP adaptive proposal-width schemes available in
#' this build.
#'
#' @description The returned names are suitable for the `scheme:` field of
#'   the `fit.csp.adaptation` block in v.3 YAML configs and for the C++
#'   factory dispatch.  Each name maps to a constructor of the form
#'   `AdaptiveScheme.<PascalName>` via [.scheme.name.to.constructor()].
#'
#' @return Character vector.
#' @examples
#' schemes.available()
#' @export
schemes.available <- function() c("native", "andrieu_thoms")


#' Test whether an object was constructed by one of the AdaptiveScheme
#' constructors.
#'
#' @param x Any R object.
#' @return Logical.
#' @export
is.AdaptiveScheme <- function(x) inherits(x, "AdaptiveScheme")


# Internal helper: map a scheme's canonical lowercase snake_case name to its
# R constructor function.  Used by the v.3 YAML loader (lib/config.R) so a
# YAML `scheme: andrieu_thoms` block can be inflated to an
# AdaptiveScheme.AndrieuThoms(...) object before being handed to
# Parameter::setCSPAdaptationScheme(...).
.scheme.name.to.constructor <- function(name) {
    stopifnot(is.character(name), length(name) == 1L, !is.na(name))
    map <- list(
        native        = AdaptiveScheme.Native,
        andrieu_thoms = AdaptiveScheme.AndrieuThoms
    )
    fn <- map[[name]]
    if (is.null(fn)) {
        stop("unknown CSP adaptation scheme: '", name,
             "'; available: ", paste(schemes.available(), collapse = ", "),
             call. = FALSE)
    }
    fn
}


#' Format an AdaptiveScheme as a one-line summary string.
#'
#' @param x An AdaptiveScheme.
#' @param ... Unused.
#' @return Character vector of length 1 or 2.
#' @export
format.AdaptiveScheme <- function(x, ...) {
    header <- paste0("AdaptiveScheme: ", x$scheme)
    if (length(x$params) == 0L) return(header)
    kv <- paste0(names(x$params), " = ", unlist(x$params), collapse = ", ")
    c(header, paste0("  ", kv))
}


#' Print an AdaptiveScheme.
#'
#' @param x An AdaptiveScheme.
#' @param ... Unused.
#' @export
print.AdaptiveScheme <- function(x, ...) {
    cat(format(x), sep = "\n")
    cat("\n")
    invisible(x)
}


#' Diagnostics for the CSP adaptive scheme currently bound to a Parameter.
#'
#' @description v1 returns the scheme's canonical name only.  Future
#'   versions will populate `std_csp.trace`, `acceptance.trace`, and
#'   `scheme.specific` (e.g. the per-AA step-size history `gamma` for
#'   Andrieu-Thoms) from the underlying C++ Trace object.
#'
#' @param parameter A Parameter object returned by [initializeParameterObject()]
#'   or by a fit.
#' @param mcmc Optional MCMC object; reserved for v2 trace extraction.
#' @return A list with elements:
#'   - `scheme.name` (character)
#'   - `params` (list; NULL in v1 -- C++ does not yet expose them)
#'   - `std_csp.trace` (NULL in v1)
#'   - `acceptance.trace` (NULL in v1)
#'   - `scheme.specific` (empty list in v1)
#' @export
adaptive.scheme.diagnostics <- function(parameter, mcmc = NULL) {
    list(
        scheme.name      = parameter$getCSPAdapter()$name(),
        params           = NULL,
        std_csp.trace    = NULL,
        acceptance.trace = NULL,
        scheme.specific  = list()
    )
}
