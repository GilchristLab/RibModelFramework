# ============================================================================
# phiSpec.R -- phi prior specification constructors for initializeParameterObject
#
# Mode constructors:  constrained(), fixed(), estimated()
# Prior constructors: prior_uniform(), prior_normal(), prior_student_t(),
#                     prior_exponential()
#
# Usage:
#   phi.mphi = constrained(statistic = "median", value = 1)
#   phi.sphi = estimated(prior = prior_normal(mean = 1.4, sd = 0.125))
# ============================================================================


# ---- internal statistic code table ----------------------------------------

.PHI_STATISTIC_CODES <- c(
  mean     = 0L,
  median   = 1L,
  mode     = 2L,
  variance = 3L,
  sd       = 4L
)

.PHI_MU_MODE_CODES <- c(constrained = 0L, fixed = 1L)
.SPHI_PRIOR_CODES  <- c(flat = 0L, normal = 1L, uniform = 2L)


# ---- mode constructors -------------------------------------------------------

#' Specify a constrained mPhi or sphi
#'
#' @param statistic Character. Which distributional property of phi to pin.
#'   One of \code{"mean"}, \code{"median"}, \code{"mode"}, \code{"sd"},
#'   \code{"variance"}.  \code{"mean"} is the log-space mean (\code{meanlog}
#'   in \code{\link{dlnorm}}); \code{sd} / \code{variance} pin the
#'   scale of phi but allow \code{E[phi]} to vary with sphi (see Warning).
#' @param value Positive numeric. Target value for the statistic (default 1).
#' @return A \code{phi_spec} object of mode \code{"constrained"}.
#' @examples
#' # Recommended for fitting: anchor the median of phi to 1
#' constrained(statistic = "median", value = 1)
#' # Legacy convention: anchor the mean of phi to 1
#' constrained(statistic = "mean", value = 1)
#' @export
constrained <- function(statistic = "mean", value = 1) {
  if (!statistic %in% names(.PHI_STATISTIC_CODES))
    stop("statistic must be one of: ",
         paste(names(.PHI_STATISTIC_CODES), collapse = ", "), "\n")
  if (!is.numeric(value) || length(value) != 1L || value <= 0)
    stop("value must be a single positive number\n")

  if (statistic %in% c("sd", "variance"))
    warning("constrained(statistic='", statistic, "') pins the variance of phi ",
            "but not its mean. E[phi] shrinks toward 0 as sphi grows. ",
            "A prior on sphi (phi.sphi = estimated(..., prior = ...)) is strongly recommended.",
            call. = FALSE)
  if (statistic == "mode")
    warning("constrained(statistic='mode') causes E[phi] to grow rapidly with sphi ",
            "(E[phi] = value * exp(3*sphi^2/2)). ",
            "A prior on sphi is strongly recommended.", call. = FALSE)
  if (statistic == "median")
    warning("constrained(statistic='median') causes E[phi] to grow with sphi ",
            "(E[phi] = value * exp(sphi^2/2)). ",
            "Consider a prior on sphi.", call. = FALSE)

  structure(
    list(mode            = "constrained",
         statistic       = statistic,
         statistic_code  = .PHI_STATISTIC_CODES[[statistic]],
         value           = value),
    class = c("phi_spec_constrained", "phi_spec")
  )
}


#' Specify a fixed mPhi or sphi
#'
#' @param value Numeric. The fixed value. For \code{phi.mphi} this is the
#'   log-space mean (\code{meanlog}); for \code{phi.sphi} it is the log-space
#'   SD (\code{sdlog}) and must be positive.
#' @return A \code{phi_spec} object of mode \code{"fixed"}.
#' @examples
#' # Hold sphi fixed at 1.4
#' fixed(value = 1.4)
#' @export
fixed <- function(value) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value))
    stop("value must be a single finite number\n")
  structure(
    list(mode = "fixed", value = value),
    class = c("phi_spec_fixed", "phi_spec")
  )
}


#' Specify an estimated mPhi or sphi
#'
#' @param prior A \code{prior_dist} object (from \code{prior_uniform()},
#'   \code{prior_normal()}, etc.) or \code{NULL} for an improper flat prior
#'   (backward-compatible with the pre-phi-spec default).
#' @param init Numeric initial value, or \code{NULL} to derive from the
#'   constrained formula at the starting sphi. Only used for \code{phi.mphi}
#'   (deferred: not yet implemented for v1).
#' @return A \code{phi_spec} object of mode \code{"estimated"}.
#' @examples
#' # Estimate sphi with a weakly-informative normal prior
#' estimated(prior = prior_normal(mean = 1.4, sd = 0.125))
#' # Estimate sphi with a bounded-uniform prior
#' estimated(prior = prior_uniform(low = 0, high = 10))
#' @export
estimated <- function(prior = prior_uniform(low = 0, high = 10), init = NULL) {
  if (!is.null(prior) && !inherits(prior, "prior_dist"))
    stop("prior must be a prior_dist object (from prior_uniform, prior_normal, etc.) or NULL\n")
  structure(
    list(mode = "estimated", prior = prior, init = init),
    class = c("phi_spec_estimated", "phi_spec")
  )
}


# ---- prior distribution constructors ----------------------------------------

#' Uniform prior for MCMC parameter estimation
#'
#' @param low  Lower bound (default 0). Must be >= 0 for sphi.
#' @param high Upper bound (default 10).
#' @return A \code{prior_dist} object.
#' @examples
#' prior_uniform(low = 0, high = 10)
#' @export
prior_uniform <- function(low = 0, high = 10) {
  if (!is.numeric(low)  || !is.finite(low))  stop("low must be finite\n")
  if (!is.numeric(high) || !is.finite(high)) stop("high must be finite\n")
  if (low >= high) stop("low must be strictly less than high\n")
  structure(list(dist = "uniform", low = low, high = high), class = "prior_dist")
}


#' Normal prior for MCMC parameter estimation
#'
#' @param mean Prior mean (default 0).
#' @param sd   Prior standard deviation (default 1). Must be positive.
#' @return A \code{prior_dist} object.
#' @examples
#' prior_normal(mean = 1.4, sd = 0.125)
#' @export
prior_normal <- function(mean = 0, sd = 1) {
  if (!is.numeric(sd) || sd <= 0) stop("sd must be positive\n")
  structure(list(dist = "normal", mean = mean, sd = sd), class = "prior_dist")
}


#' Student-t prior for MCMC parameter estimation
#'
#' @param df   Degrees of freedom (default 3). Must be positive.
#' @param mean Location (default 0).
#' @param sd   Scale (default 1). Must be positive.
#' @return A \code{prior_dist} object.
#' @examples
#' prior_student_t(df = 3, mean = 0, sd = 1)
#' @export
prior_student_t <- function(df = 3, mean = 0, sd = 1) {
  if (!is.numeric(df) || df <= 0) stop("df must be positive\n")
  if (!is.numeric(sd) || sd <= 0) stop("sd must be positive\n")
  structure(list(dist = "student_t", df = df, mean = mean, sd = sd),
            class = "prior_dist")
}


#' Exponential prior for MCMC parameter estimation
#'
#' @param rate Rate parameter (default 1). Must be positive.
#' @return A \code{prior_dist} object.
#' @examples
#' prior_exponential(rate = 1)
#' @export
prior_exponential <- function(rate = 1) {
  if (!is.numeric(rate) || rate <= 0) stop("rate must be positive\n")
  structure(list(dist = "exponential", rate = rate), class = "prior_dist")
}


# ---- print methods ----------------------------------------------------------

#' Print a phi_spec object
#' @param x A \code{phi_spec} object.
#' @param ... Further arguments passed to or from other methods.
#' @return \code{x} invisibly.
#' @export
print.phi_spec <- function(x, ...) {
  cat("<phi_spec: mode =", x$mode)
  if (x$mode == "constrained")
    cat(", statistic =", x$statistic, ", value =", x$value)
  if (x$mode == "fixed")
    cat(", value =", x$value)
  if (x$mode == "estimated" && !is.null(x$prior))
    cat(", prior =", format(x$prior))
  cat(">\n")
  invisible(x)
}

#' Format a prior_dist object as a string
#' @param x A \code{prior_dist} object.
#' @param ... Further arguments passed to or from other methods.
#' @return A character string describing the distribution.
#' @export
format.prior_dist <- function(x, ...) {
  switch(x$dist,
    uniform   = paste0("Uniform(", x$low, ", ", x$high, ")"),
    normal    = paste0("Normal(mean=", x$mean, ", sd=", x$sd, ")"),
    student_t = paste0("StudentT(df=", x$df, ", mean=", x$mean, ", sd=", x$sd, ")"),
    exponential = paste0("Exponential(rate=", x$rate, ")"),
    paste0("prior_dist(", x$dist, ")")
  )
}

#' Print a prior_dist object
#' @param x A \code{prior_dist} object.
#' @param ... Further arguments passed to or from other methods.
#' @return \code{x} invisibly.
#' @export
print.prior_dist <- function(x, ...) {
  cat("<prior_dist:", format(x), ">\n")
  invisible(x)
}


# ---- internal helper: apply phi spec to a Parameter C++ object --------------

# Called from initializeParameterObject after the C++ object is constructed.
# phi.mphi: phi_spec or NULL (NULL -> use legacy phiPriorConstraint path)
# phi.sphi: phi_spec or NULL
.applyPhiSpec <- function(parameter, phi.mphi, phi.sphi) {

  # --- mPhi ---
  if (!is.null(phi.mphi)) {
    if (inherits(phi.mphi, "phi_spec_constrained")) {
      parameter$setPhiMuMode(0L)  # PHI_MU_CONSTRAINED
      parameter$setPhiPriorConstraint(phi.mphi$statistic_code)
      parameter$setPhiConstraintValue(phi.mphi$value)
    } else if (inherits(phi.mphi, "phi_spec_fixed")) {
      parameter$setPhiMuMode(1L)  # PHI_MU_FIXED
      parameter$setPhiMuFixed(phi.mphi$value)
    } else if (inherits(phi.mphi, "phi_spec_estimated")) {
      parameter$setPhiMuMode(2L)  # PHI_MU_ESTIMATED
      # Set initial mphi value for each synthesis rate category
      n_cat <- parameter$getNumSynthesisRateCategories()
      init_val <- if (!is.null(phi.mphi$init)) phi.mphi$init else 0.0
      for (i in seq_len(n_cat))
        parameter$setMuSynthesisRate(init_val, i - 1L)
      # Wire prior
      prior <- phi.mphi$prior
      if (is.null(prior)) {
        parameter$setPhiMuPriorType(0L)  # flat
      } else if (prior$dist == "normal") {
        parameter$setPhiMuPriorType(1L)
        parameter$setPhiMuPriorMu(prior$mean)
        parameter$setPhiMuPriorSd(prior$sd)
      } else if (prior$dist == "uniform") {
        parameter$setPhiMuPriorType(0L)  # treat as flat
      } else {
        stop("Only flat/uniform and normal priors are supported for phi.mphi = estimated()\n")
      }
    }
  }

  # --- sphi ---
  if (!is.null(phi.sphi)) {
    if (inherits(phi.sphi, "phi_spec_fixed")) {
      parameter$fixSphi()
    } else if (inherits(phi.sphi, "phi_spec_estimated")) {
      .applyPhiSphiPrior(parameter, phi.sphi$prior)
    }
  }
}


.applyPhiSphiPrior <- function(parameter, prior) {
  if (is.null(prior)) {
    # Improper flat: sphiPriorType = 0 (default), no setter call needed.
    parameter$setSphiPriorType(0L)
    return(invisible(NULL))
  }
  switch(prior$dist,
    uniform = {
      parameter$setSphiPriorType(2L)
      parameter$setSphiPriorBounds(prior$low, prior$high)
    },
    normal = {
      parameter$setSphiPriorType(1L)
      parameter$setSphiPrior(prior$mean, prior$sd)
    },
    student_t = stop("student_t prior for sphi not yet implemented\n"),
    exponential = stop("exponential prior for sphi not yet implemented in MCMC\n"),
    stop("Unknown prior distribution: ", prior$dist, "\n")
  )
}
