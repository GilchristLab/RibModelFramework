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


# ---------------- 7. Metropolis-Hastings sampler -----------------------------
#
# Free parameters (p, mu1, sigma1, sigma2). mu2 is derived from the constraint
# every time any of them changes. We update each parameter in turn via
# single-site random walk:
#   p      : random walk on logit(p) (handles (0,1) constraint, symmetric in
#            logit space; Jacobian p*(1-p) added to the M-H ratio)
#   mu1    : plain Gaussian random walk in mu1 space
#   sigma1 : random walk on log(sigma1) (handles s > 0, symmetric in log space;
#            Jacobian sigma added to the M-H ratio)
#   sigma2 : random walk on log(sigma2) (same as sigma1; sigma2 > sigma1 is
#            handled by log_prior returning -Inf, which rejects automatically)
#
# Proposal widths in `proposal` are on the *transformed* scale:
#   proposal$p       -> SD of the random walk on logit(p)
#   proposal$mu1     -> SD on mu1 directly
#   proposal$sigma1  -> SD on log(sigma1)
#   proposal$sigma2  -> SD on log(sigma2)

mh_phi_mixture <- function(phi,
                           constraint = c("mean", "median"),
                           hyperparams = default_hyperparams(),
                           init = list(p = 0.85, mu1 = -0.05,
                                       sigma1 = 0.5, sigma2 = 0.3),
                           proposal = list(p = 0.3, mu1 = 0.1,
                                           sigma1 = 0.15, sigma2 = 0.15),
                           n_iter = 5000,
                           n_burnin = 1000,
                           thin = 1,
                           seed = NULL,
                           verbose = TRUE) {
    constraint <- match.arg(constraint)
    if (!is.null(seed)) set.seed(seed)

    # State holds (p, mu1, sigma1, sigma2, mu2). mu2 is derived but cached so
    # the sampler outputs it alongside the free params.
    state <- list(p = init$p, mu1 = init$mu1,
                  sigma1 = init$sigma1, sigma2 = init$sigma2)
    state$mu2 <- derive_mu2(state$p, state$mu1, state$sigma1, state$sigma2,
                             constraint)
    if (is.na(state$mu2) || state$mu2 < state$mu1) {
        stop("Initial state is infeasible or violates ordering: ",
             "mu2 = ", state$mu2, ", mu1 = ", state$mu1)
    }
    lp <- log_posterior(state$p, state$mu1, state$sigma1, state$sigma2,
                         phi, constraint, hyperparams)
    if (!is.finite(lp)) {
        stop("Initial log-posterior is -Inf; check init values and hyperparams")
    }

    # Storage for post-burnin samples
    keep   <- (n_iter - n_burnin) %/% thin
    out    <- data.frame(
        iter    = integer(keep),
        p       = numeric(keep),
        mu1     = numeric(keep),
        sigma1  = numeric(keep),
        sigma2  = numeric(keep),
        mu2     = numeric(keep),
        logpost = numeric(keep)
    )
    accept_count <- c(p = 0L, mu1 = 0L, sigma1 = 0L, sigma2 = 0L)

    # M-H step for one parameter. Proposes on the appropriate transformed
    # scale and applies the corresponding Jacobian to the log-ratio.
    propose_one <- function(param, sd, state, lp) {
        proposed <- state
        log_jac <- 0
        if (param == "p") {
            logit_p_new <- log(state$p / (1 - state$p)) + rnorm(1, 0, sd)
            proposed$p <- 1 / (1 + exp(-logit_p_new))
            log_jac <- log(proposed$p * (1 - proposed$p)) -
                       log(state$p   * (1 - state$p))
        } else if (param == "mu1") {
            proposed$mu1 <- state$mu1 + rnorm(1, 0, sd)
        } else if (param == "sigma1") {
            proposed$sigma1 <- state$sigma1 * exp(rnorm(1, 0, sd))
            log_jac <- log(proposed$sigma1) - log(state$sigma1)
        } else if (param == "sigma2") {
            proposed$sigma2 <- state$sigma2 * exp(rnorm(1, 0, sd))
            log_jac <- log(proposed$sigma2) - log(state$sigma2)
        }
        proposed$mu2 <- derive_mu2(proposed$p, proposed$mu1,
                                    proposed$sigma1, proposed$sigma2,
                                    constraint)
        new_lp <- log_posterior(proposed$p, proposed$mu1,
                                 proposed$sigma1, proposed$sigma2,
                                 phi, constraint, hyperparams)
        log_ratio <- new_lp - lp + log_jac
        if (is.finite(log_ratio) && log(runif(1)) < log_ratio) {
            list(state = proposed, lp = new_lp, accepted = TRUE)
        } else {
            list(state = state, lp = lp, accepted = FALSE)
        }
    }

    # Main loop
    kept <- 0L
    for (i in seq_len(n_iter)) {
        for (param in c("p", "mu1", "sigma1", "sigma2")) {
            r  <- propose_one(param, proposal[[param]], state, lp)
            state <- r$state
            lp    <- r$lp
            if (r$accepted) accept_count[[param]] <- accept_count[[param]] + 1L
        }
        if (i > n_burnin && ((i - n_burnin) %% thin == 0L)) {
            kept <- kept + 1L
            out[kept, ] <- list(
                iter = i, p = state$p, mu1 = state$mu1,
                sigma1 = state$sigma1, sigma2 = state$sigma2,
                mu2 = state$mu2, logpost = lp
            )
        }
    }

    accept_rate <- accept_count / n_iter
    if (verbose) {
        cat("Acceptance rates:\n")
        for (param in names(accept_rate)) {
            cat(sprintf("  %-7s: %.3f\n", param, accept_rate[[param]]))
        }
    }
    list(samples = out[seq_len(kept), ],
         accept_rate = accept_rate,
         final_state = state,
         constraint = constraint)
}
