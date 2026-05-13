# =============================================================================
# Prototype: mixture-of-2 lognormal prior on phi for ROC
#
# Model: f(phi) = p * LN(mu1, sigma1) + (1 - p) * LN(mu2, sigma2)
# where mu2 is *derived* from one of two constraints on the mixture:
#   "mean":   p * exp(mu1 + s1^2/2) + (1-p) * exp(mu2 + s2^2/2) = 1
#   "median": p * pnorm(-mu1/s1)    + (1-p) * pnorm(-mu2/s2)    = 0.5
# Both have closed-form solutions for mu2; "median" via the inverse normal CDF.
#
# Identifiability: hard constraint mu2 >= mu1 (component 2 is the higher
# log-location mode). Plus soft priors that encode "p > 0.8" (Beta) and
# "sigma2 < sigma1" (half-Normal truncated by sigma1).
#
# Feasible region note: combining the median=1 constraint with mu2 >= mu1
# and a heavy lower component (large p) yields a narrow feasible region for
# mu1. Specifically, mu2 >= mu1 under median=1 requires
#   0.5 <= p * Phi(-mu1/sigma1) + (1-p) * Phi(-mu1/sigma2)
# When p is large (>= 0.9), the weighted average of the two CDFs at -mu1 must
# reach 0.5, which forces mu1 <= 0. Under the mean=1 constraint there is no
# analogous restriction. The sampler should explore mu1 in (-, 0] when p is
# large and the constraint is "median"; this is a feature, not a bug -- the
# data should drive p down if the lower component is well above phi=1.
#
# This file is exploratory R; ports to C++ once the math, the constraint
# enforcement, and the sampler all behave on simulated data. Not loaded by the
# package. See task #10 of the convergence/phi-mixture work plan.
#
# Author: Mike Gilchrist + Claude Opus 4.7
# =============================================================================

# ---------------- 1. Hyperparameters (settable) ------------------------------

default_hyperparams <- function() {
    list(
        # Beta on p; mode at (alpha-1)/(alpha+beta-2) = 7/8 = 0.875.
        p      = list(alpha = 8,    beta = 2),
        # Vague Normal on mu1.
        mu1    = list(mean  = 0,    sd   = 10),
        # Half-Normal(scale) on sigma1.
        sigma1 = list(scale = 1),
        # Half-Normal(scale) on sigma2, truncated above by sigma1.
        sigma2 = list(scale = 1)
    )
}

# Merge user overrides into defaults, element by element (so user can override
# just one sub-list and keep defaults for the rest).
merge_hyperparams <- function(user = list()) {
    defaults <- default_hyperparams()
    for (k in names(user)) {
        if (is.null(defaults[[k]])) {
            stop("unknown hyperparameter block: ", k,
                 " (valid: ", paste(names(defaults), collapse = ", "), ")")
        }
        for (nm in names(user[[k]])) defaults[[k]][[nm]] <- user[[k]][[nm]]
    }
    defaults
}


# ---------------- 2. Derive mu2 from constraint ------------------------------

# Returns NA if the constraint cannot be satisfied with mu2 finite (e.g.
# component 1 already accounts for more than the total mean/median mass).
derive_mu2 <- function(p, mu1, sigma1, sigma2,
                       constraint = c("mean", "median")) {
    constraint <- match.arg(constraint)

    # p in (0,1) strict: at the endpoints the mixture is degenerate.
    if (p <= 0 || p >= 1) return(NA_real_)
    if (sigma1 <= 0 || sigma2 <= 0) return(NA_real_)

    if (constraint == "mean") {
        # p*exp(mu1+s1^2/2) + (1-p)*exp(mu2+s2^2/2) = 1
        m1 <- exp(mu1 + sigma1^2 / 2)
        m2 <- (1 - p * m1) / (1 - p)
        if (m2 <= 0) return(NA_real_) # second component would need negative mean
        return(log(m2) - sigma2^2 / 2)
    } else {
        # p*pnorm(-mu1/s1) + (1-p)*pnorm(-mu2/s2) = 0.5
        # => pnorm(-mu2/s2) = (0.5 - p*pnorm(-mu1/s1)) / (1-p)
        q <- (0.5 - p * pnorm(-mu1 / sigma1)) / (1 - p)
        if (q <= 0 || q >= 1) return(NA_real_) # second-component CDF target out of range
        return(-sigma2 * qnorm(q))
    }
}


# ---------------- 3. Mixture density -----------------------------------------

# Pair-wise log-sum-exp: log(exp(a) + exp(b)), numerically stable, vectorized.
.lse2 <- function(a, b) {
    m <- pmax(a, b)
    res <- m + log1p(exp(pmin(a, b) - m))
    # both -Inf -> NaN; recover -Inf
    res[is.nan(res)] <- -Inf
    res
}

# Raw mixture log-density: caller supplies mu2 explicitly. No constraint check.
dmixture_lognormal_raw <- function(phi, p, mu1, sigma1, mu2, sigma2,
                                   log = TRUE) {
    log_f1 <- dlnorm(phi, meanlog = mu1, sdlog = sigma1, log = TRUE)
    log_f2 <- dlnorm(phi, meanlog = mu2, sdlog = sigma2, log = TRUE)
    res <- .lse2(log(p) + log_f1, log1p(-p) + log_f2)
    if (log) res else exp(res)
}

# Constrained density: mu2 is derived from the chosen constraint.
# Returns -Inf log-density for any phi if the constraint is infeasible or the
# label-switching guard (mu2 >= mu1) fires.
dmixture_lognormal <- function(phi, p, mu1, sigma1, sigma2,
                               constraint = c("mean", "median"),
                               log = TRUE) {
    constraint <- match.arg(constraint)
    mu2 <- derive_mu2(p, mu1, sigma1, sigma2, constraint)
    if (is.na(mu2) || mu2 < mu1) {
        return(rep(if (log) -Inf else 0, length(phi)))
    }
    dmixture_lognormal_raw(phi, p, mu1, sigma1, mu2, sigma2, log = log)
}


# ---------------- 4. Simulator -----------------------------------------------

# Draw n phi values from the mixture. Returns list(phi, component, mu2).
rmixture_lognormal <- function(n, p, mu1, sigma1, sigma2,
                               constraint = c("mean", "median")) {
    constraint <- match.arg(constraint)
    mu2 <- derive_mu2(p, mu1, sigma1, sigma2, constraint)
    if (is.na(mu2) || mu2 < mu1) {
        stop("Infeasible parameters: derive_mu2 returned ", mu2,
             " for mu1 = ", mu1, "; constraint = ", constraint)
    }
    component <- rbinom(n, size = 1, prob = p) # 1 with prob p (= component 1)
    phi <- numeric(n)
    n1 <- sum(component == 1)
    n2 <- n - n1
    if (n1 > 0) phi[component == 1] <- rlnorm(n1, meanlog = mu1, sdlog = sigma1)
    if (n2 > 0) phi[component == 0] <- rlnorm(n2, meanlog = mu2, sdlog = sigma2)
    list(phi = phi, component = component, mu2 = mu2)
}


# ---------------- 5. Log-prior on free parameters ----------------------------

.log_dhalfnorm <- function(x, scale) {
    if (x <= 0) return(-Inf)
    log(2) + dnorm(x, mean = 0, sd = scale, log = TRUE)
}

# Half-Normal truncated above by `upper`. f(x) = 2*dnorm(x,0,s) / Z, where
# Z = 2*pnorm(upper, 0, s) - 1 = Pr(|N(0,s)| <= upper).
.log_dhalfnorm_trunc_above <- function(x, scale, upper) {
    if (x <= 0 || x >= upper) return(-Inf)
    log_Z <- log(2 * pnorm(upper, mean = 0, sd = scale) - 1)
    log(2) + dnorm(x, mean = 0, sd = scale, log = TRUE) - log_Z
}

log_prior <- function(p, mu1, sigma1, sigma2,
                      hyperparams = default_hyperparams()) {
    if (p <= 0 || p >= 1) return(-Inf)
    if (sigma1 <= 0 || sigma2 <= 0) return(-Inf)
    if (sigma2 >= sigma1) return(-Inf) # truncation hard wall

    hp <- hyperparams
    dbeta(p, hp$p$alpha, hp$p$beta, log = TRUE) +
        dnorm(mu1, hp$mu1$mean, hp$mu1$sd, log = TRUE) +
        .log_dhalfnorm(sigma1, hp$sigma1$scale) +
        .log_dhalfnorm_trunc_above(sigma2, hp$sigma2$scale, sigma1)
}


# ---------------- 6. Joint log-posterior (used by the sampler in task #11) ---

log_posterior <- function(p, mu1, sigma1, sigma2, phi,
                          constraint = c("mean", "median"),
                          hyperparams = default_hyperparams()) {
    constraint <- match.arg(constraint)
    lp <- log_prior(p, mu1, sigma1, sigma2, hyperparams)
    if (!is.finite(lp)) return(-Inf)
    ll <- sum(dmixture_lognormal(phi, p, mu1, sigma1, sigma2,
                                  constraint = constraint, log = TRUE))
    if (!is.finite(ll)) return(-Inf) # label-switch or infeasibility
    lp + ll
}
