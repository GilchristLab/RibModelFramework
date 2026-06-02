#' Build Stan data and init lists from an AnaCoDa Genome object
#'
#' Converts an AnaCoDa Genome object and phi prior specifications to the data
#' and init lists required by the ROC Stan models (roc_sphi_est.stan and
#' roc_arcsine.stan).  Returns a two-slot list: \code{$data} (pass to
#' \code{mod$sample(data = ...)}) and \code{$init} (pass to
#' \code{mod$sample(init = list(...))}).
#'
#' @section phi.mphi modes:
#' \describe{
#'   \item{\code{constrained(statistic, value)}}{mphi is derived from sphi each
#'     MCMC iteration so that the specified distributional statistic of phi
#'     equals \code{value}.  Five statistics are supported: \code{"mean"},
#'     \code{"median"}, \code{"mode"}, \code{"variance"}, \code{"sd"}.}
#'   \item{\code{fixed(value)}}{mphi is a fixed data constant equal to
#'     \code{value} (log-space mean, i.e.\ \code{meanlog}).  sphi is still
#'     estimated.}
#' }
#' \code{estimated()} for phi.mphi is not yet implemented (deferred).
#'
#' @section phi.sphi modes:
#' \describe{
#'   \item{\code{estimated(prior_uniform(low, high))}}{sphi is sampled within
#'     \code{[low, high]}; no prior density statement is added (uniform on
#'     the parameter bounds).  Default: \code{prior_uniform(0, 10)}.}
#'   \item{\code{estimated(prior_normal(mean, sd))}}{sphi is sampled on
#'     \code{(0, Inf)} with prior \code{sphi ~ normal(mean, sd)} truncated
#'     to positive values.}
#'   \item{\code{estimated(NULL)}}{improper flat prior on positive reals.}
#' }
#' \code{fixed()} for phi.sphi is not yet implemented.
#'
#' @section dM prior from low-expression genes:
#' When \code{dM.prior = "scuo"} (default), the prior mean for each mutation
#' bias parameter is estimated from the codon frequencies of the
#' \code{scuo.low.frac} lowest-SCUO genes.  In low-expression genes, dEta*phi
#' contributes minimally, so observed codon frequencies approximate the
#' mutation-bias-only distribution:
#' \deqn{dM_k \approx \log(\hat{p}_{ref}) - \log(\hat{p}_k)}
#' where \eqn{\hat{p}} values use a 0.5 Laplace pseudocount per codon.
#' Pass \code{dM.prior = "flat"} or supply \code{dM_prior_mean} explicitly
#' to skip this estimation.
#'
#' @param genome An AnaCoDa Genome object.
#' @param phi.mphi A \code{phi_spec} object for the log-space mean of phi.
#'   Default: \code{constrained("mean", 1)} (mean(phi) = 1 convention).
#'   See Details.
#' @param phi.sphi A \code{phi_spec} object for the log-space SD of phi.
#'   Default: \code{estimated(prior_uniform(0, 10))}.  See Details.
#' @param approx_min_n Integer. Minimum total count per (gene, AA) for the
#'   arcsine branch in \code{roc_arcsine.stan} (default 20).  Included in
#'   \code{$data} for cross-model compatibility; ignored by
#'   \code{roc_sphi_est.stan}.
#' @param noncentered Integer 0/1.  0 = centered latent_phi (default, good
#'   for G >= 1000 full-genome fits), 1 = non-centered (better for data-sparse
#'   fits with poor sphi ESS).
#' @param dM.prior Character.  \code{"scuo"} (default): estimate
#'   \code{dM_prior_mean} from low-SCUO genes.  \code{"flat"}: use 0.
#'   Ignored if \code{dM_prior_mean} is supplied.
#' @param scuo Numeric vector of length G, or \code{NULL}.  Pre-computed SCUO
#'   (from \code{\link{calculateSCUO}}).  Computed internally if \code{NULL}.
#'   Pass a pre-computed vector to avoid redundant work.
#' @param scuo.low.frac Fraction of lowest-SCUO genes for dM prior estimation.
#'   Default 0.25; at least 5 genes are always used.
#' @param dM_prior_mean Prior mean for dM.  \code{NULL} (default): use
#'   \code{dM.prior}.  Scalar or length-K vector: used as-is.
#' @param dM_prior_sd Prior SD for dM.  Scalar or length-K vector.  Default 1.
#' @param dEta_prior_mean Prior mean for dEta.  Default 0.
#' @param dEta_prior_sd Prior SD for dEta.  Default 1.
#' @param deta_scale_anchor Integer 0/1.  Scale-anchor reparameterisation of
#'   dEta to reduce the dEta-phi ridge (default 0, disabled).
#' @param deta_anchor_ref Reference mphi level for the scale anchor
#'   (ignored when \code{deta_scale_anchor = 0}).
#' @param deta_phi_center Numeric.  Centering constant for the phi predictor
#'   (0 = disabled, default).
#' @param phi.init Character or numeric.  \code{"scuo"} (default): per-gene
#'   log-phi from SCUO.  \code{"uniform"}: phi = 1 everywhere.  Numeric
#'   vector of length G: used directly as starting phi values (log-transformed
#'   and scaled internally).
#' @param sphi.init Numeric.  Starting value for sphi and target SD for the
#'   initial log-phi distribution (default 1.0).  Must be within the sphi
#'   bounds implied by \code{phi.sphi}.
#' @param grainsize reduce_sum grain size (default 1; TBB auto-selects).
#'
#' @return A named list with two slots:
#'   \describe{
#'     \item{\code{$data}}{Stan data list.  Pass as \code{mod$sample(data = result$data)}.}
#'     \item{\code{$init}}{Stan init list.  Pass as \code{mod$sample(init = list(result$init))}.}
#'   }
#'
#' @seealso \code{\link{adviToWarmStart}} for warm-starting HMC from ADVI output.
#'   \code{\link{constrained}}, \code{\link{fixed}}, \code{\link{estimated}},
#'   \code{\link{prior_uniform}}, \code{\link{prior_normal}} for phi spec constructors.
#'
#' @examples
#' \dontrun{
#' genome <- initializeGenomeObject("mygenome.fasta")
#' scuo   <- calculateSCUO(genome)$SCUO   # compute once, reuse
#'
#' # default: mean(phi)=1 constraint, Uniform(0,10) prior on sphi
#' stan <- initializeStan(genome, scuo = scuo)
#'
#' # median(phi)=1, Normal(1.4, 0.125) prior on sphi
#' stan <- initializeStan(genome, scuo = scuo,
#'                        phi.mphi = constrained("median", 1),
#'                        phi.sphi = estimated(prior_normal(1.4, 0.125)))
#'
#' library(cmdstanr)
#' mod <- cmdstan_model("stan/roc_sphi_est.stan",
#'                      cpp_options = list(stan_threads = TRUE))
#' fit <- mod$variational(data = stan$data, init = list(stan$init), threads = 4L)
#' }
#'
#' @export
initializeStan <- function(genome,
                            phi.mphi          = constrained("mean", 1),
                            phi.sphi          = estimated(prior_uniform(0, 10)),
                            approx_min_n      = 20L,
                            noncentered       = 0L,
                            dM.prior          = "scuo",
                            scuo              = NULL,
                            scuo.low.frac     = 0.25,
                            dM_prior_mean     = NULL,
                            dM_prior_sd       = 1.0,
                            dEta_prior_mean   = 0.0,
                            dEta_prior_sd     = 1.0,
                            deta_scale_anchor = 0L,
                            deta_anchor_ref   = 0.0,
                            deta_phi_center   = 0.0,
                            phi.init          = "scuo",
                            sphi.init         = 1.0,
                            grainsize         = 1L) {

  if (!inherits(phi.mphi, "phi_spec"))
    stop("phi.mphi must be a phi_spec object (constrained(), fixed(), or estimated())\n")
  if (!inherits(phi.sphi, "phi_spec"))
    stop("phi.sphi must be a phi_spec object (estimated() or fixed())\n")

  # ---- genome-derived structural fields ------------------------------------
  group_list <- c("A","C","D","E","F","G","H","I","K","L",
                  "N","P","Q","R","S","T","V","Y","Z")
  A <- length(group_list)

  nonref_by_aa  <- lapply(group_list, function(aa) AAToCodon(aa, TRUE))
  nonref_codons <- unlist(nonref_by_aa)
  K             <- length(nonref_codons)

  aa_lengths <- vapply(nonref_by_aa, length, integer(1L))
  aa_end     <- cumsum(aa_lengths)
  aa_start   <- aa_end - aa_lengths + 1L

  counts64 <- getCodonCounts(genome)
  G        <- nrow(counts64)

  y_k <- as.matrix(counts64[, nonref_codons, drop = FALSE])
  storage.mode(y_k) <- "integer"
  dimnames(y_k)     <- NULL

  N_ga <- matrix(0L, nrow = G, ncol = A)
  for (a in seq_along(group_list)) {
    all_codons <- AAToCodon(group_list[a], FALSE)
    N_ga[, a]  <- as.integer(rowSums(counts64[, all_codons, drop = FALSE]))
  }

  expand_k <- function(x, name) {
    if (length(x) == 1L) return(rep(as.double(x), K))
    if (length(x) == K)  return(as.double(x))
    stop(sprintf("'%s' must be length 1 or K = %d, got length %d",
                 name, K, length(x)))
  }

  # ---- dM prior ------------------------------------------------------------
  if (!is.null(dM_prior_mean)) {
    dM_pm <- expand_k(dM_prior_mean, "dM_prior_mean")
  } else if (identical(dM.prior, "scuo")) {
    if (is.null(scuo))
      scuo <- calculateSCUO(genome)$SCUO
    dM_pm <- .dMPriorFromSCUO(scuo, counts64, group_list, nonref_codons,
                               aa_start, aa_end, scuo.low.frac)
  } else {
    dM_pm <- rep(0.0, K)
  }

  # ---- phi spec → Stan data fields -----------------------------------------
  phi_fields <- .phiSpecToStanData(phi.mphi, phi.sphi)

  # Clamp sphi.init into the declared bounds with a warning
  sphi_lo <- phi_fields$sphi_low + 1e-6
  sphi_hi <- phi_fields$sphi_high - 1e-6
  if (sphi.init <= phi_fields$sphi_low || sphi.init >= phi_fields$sphi_high)
    warning(sprintf(
      "sphi.init (%.4g) is outside (sphi_low=%.4g, sphi_high=%.4g). Clamping.",
      sphi.init, phi_fields$sphi_low, phi_fields$sphi_high),
      call. = FALSE)
  sphi.init <- max(sphi_lo, min(sphi_hi, sphi.init))

  # ---- assemble data list --------------------------------------------------
  stan_data <- c(
    list(
      G               = G,
      A               = A,
      K               = K,
      aa_start        = aa_start,
      aa_end          = aa_end,
      y_k             = y_k,
      N_ga            = N_ga,
      approx_min_n    = as.integer(approx_min_n),
      dM_prior_mean   = dM_pm,
      dM_prior_sd     = expand_k(dM_prior_sd,     "dM_prior_sd"),
      dEta_prior_mean = expand_k(dEta_prior_mean, "dEta_prior_mean"),
      dEta_prior_sd   = expand_k(dEta_prior_sd,   "dEta_prior_sd"),
      noncentered       = as.integer(noncentered),
      deta_scale_anchor = as.integer(deta_scale_anchor),
      deta_anchor_ref   = as.double(deta_anchor_ref),
      deta_phi_center   = as.double(deta_phi_center),
      grainsize         = as.integer(grainsize)
    ),
    phi_fields
  )

  # ---- build init list -----------------------------------------------------
  if (identical(phi.init, "scuo")) {
    if (is.null(scuo))
      scuo <- calculateSCUO(genome)$SCUO
    log_phi_init <- .scuoToLogPhi(scuo, sphi.init)
  } else if (is.numeric(phi.init) && length(phi.init) == G) {
    raw          <- log(pmax(phi.init, 1e-10))
    log_phi_init <- .scaleLogPhi(raw, sphi.init)
  } else {
    log_phi_init <- rep(0.0, G)
  }

  if (noncentered == 1L) {
    mphi_0      <- .computeMPhiR(phi.mphi, sphi.init)
    latent_init <- (log_phi_init - mphi_0) / sphi.init
  } else {
    latent_init <- log_phi_init
  }

  stan_init <- list(
    dM         = rep(0.0, K),
    dEta       = rep(0.0, K),
    latent_phi = latent_init,
    sphi       = sphi.init
  )

  list(data = stan_data, init = stan_init)
}


#' Build a warm-start init list and diagonal inv_metric from an ADVI fit
#'
#' Extracts posterior means (for the init point) and unconstrained marginal
#' variances (for the diagonal inverse mass matrix) from a Stan ADVI fit.
#' Pass the result to \code{mod$sample()} to give HMC a pre-estimated mass
#' matrix, reducing the warmup iterations needed for good scaling.
#'
#' @section Why inv_metric matters:
#' Stan's HMC adapts a mass matrix during warmup.  Starting from a good
#' diagonal estimate (ADVI marginal variances in unconstrained space) means
#' the sampler is already well-scaled on the first leapfrog step.  For
#' parameters with large variance, such as dEta under the dEta-phi ridge,
#' this can substantially reduce the number of max-treedepth proposals in
#' early warmup.
#'
#' @section Unconstrained transforms:
#' \code{sphi} has a \code{lower=sphi_low} constraint; Stan maps it to
#' unconstrained space via a logit transform, so its unconstrained variance
#' is approximated from \eqn{\text{Var}(\log(\text{sphi}_{\text{draws}}))}.
#' All other parameters (\code{dM}, \code{dEta}, \code{latent\_phi}) are
#' unconstrained and their ADVI variances are used directly.
#'
#' @section Stan parameter ordering:
#' The returned \code{inv_metric} vector follows Stan's internal unconstrained
#' parameter ordering:
#' \code{dM[1..K]}, \code{dEta[1..K]}, \code{latent\_phi[1..G]}, \code{sphi}.
#' Legacy models that include \code{mphi\_param} get an extra entry appended.
#'
#' @param fit A \code{CmdStanVB} object from \code{mod$variational()}.
#' @param data The Stan data list returned from \code{\link{initializeStan}}\code{$data}.
#'
#' @return A named list:
#'   \describe{
#'     \item{init}{Named list of parameter starting values (posterior means in
#'       constrained space).}
#'     \item{inv_metric}{Numeric vector of unconstrained parameter variances
#'       (diagonal inverse mass matrix).  Pass as
#'       \code{mod$sample(inv_metric = result$inv_metric)}.}
#'   }
#'
#' @seealso \code{\link{initializeStan}}
#'
#' @examples
#' \dontrun{
#' stan   <- initializeStan(genome, scuo = scuo)
#' fit_vi <- mod_arcsine$variational(data = stan$data, init = list(stan$init),
#'                                   threads = 4L)
#' ws     <- adviToWarmStart(fit_vi, stan$data)
#' fit_hmc <- mod_exact$sample(data       = stan$data,
#'                              init       = list(ws$init),
#'                              inv_metric = ws$inv_metric,
#'                              chains     = 4L, ...)
#' }
#'
#' @export
adviToWarmStart <- function(fit, data) {
  draws <- fit$draws(format = "df")
  G <- data$G
  K <- data$K

  dM_cols   <- paste0("dM[",         1:K, "]")
  dEta_cols <- paste0("dEta[",       1:K, "]")
  lphi_cols <- paste0("latent_phi[", 1:G, "]")
  has_mphi  <- "mphi_param" %in% names(draws)

  # ---- posterior means (constrained space) for init -----------------------
  dM_mean   <- as.numeric(colMeans(draws[, dM_cols,   drop = FALSE]))
  dEta_mean <- as.numeric(colMeans(draws[, dEta_cols, drop = FALSE]))
  lphi_mean <- as.numeric(colMeans(draws[, lphi_cols, drop = FALSE]))
  sphi_mean <- mean(draws[["sphi"]])

  init <- list(
    dM         = dM_mean,
    dEta       = dEta_mean,
    latent_phi = lphi_mean,
    sphi       = sphi_mean
  )
  if (has_mphi)
    init$mphi_param <- mean(draws[["mphi_param"]])

  # ---- unconstrained variances for inv_metric -----------------------------
  # sphi: lower-bounded -> approx unconstrained variance via log transform
  dM_var   <- as.numeric(apply(draws[, dM_cols,   drop = FALSE], 2L, var))
  dEta_var <- as.numeric(apply(draws[, dEta_cols, drop = FALSE], 2L, var))
  lphi_var <- as.numeric(apply(draws[, lphi_cols, drop = FALSE], 2L, var))
  sphi_var <- var(log(draws[["sphi"]]))

  # Stan parameter order: dM[1..K], dEta[1..K], latent_phi[1..G], sphi
  inv_metric <- c(dM_var, dEta_var, lphi_var, sphi_var)
  if (has_mphi)
    inv_metric <- c(inv_metric, var(draws[["mphi_param"]]))

  list(init = init, inv_metric = inv_metric)
}


# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

# Translate phi_spec objects into the scalar/integer Stan data fields that
# roc_sphi_est.stan expects for phi.mphi and phi.sphi.
.phiSpecToStanData <- function(phi.mphi, phi.sphi) {
  result <- list()

  # --- phi.mphi ---
  if (inherits(phi.mphi, "phi_spec_estimated")) {
    stop("phi.mphi = estimated() is not yet implemented for Stan\n")
  } else if (inherits(phi.mphi, "phi_spec_constrained")) {
    result$phi_mphi_mode      <- 0L
    result$phi_mphi_statistic <- phi.mphi$statistic_code
    result$phi_mphi_value     <- as.double(phi.mphi$value)
    result$phi_mphi_fixed     <- 0.0   # placeholder; Stan requires the field
  } else if (inherits(phi.mphi, "phi_spec_fixed")) {
    result$phi_mphi_mode      <- 1L
    result$phi_mphi_statistic <- 0L    # placeholder
    result$phi_mphi_value     <- 1.0   # placeholder; must be > 0 for Stan constraint
    result$phi_mphi_fixed     <- as.double(phi.mphi$value)
  } else {
    stop("Unrecognised phi.mphi class\n")
  }

  # --- phi.sphi ---
  if (inherits(phi.sphi, "phi_spec_fixed")) {
    stop("phi.sphi = fixed() is not yet implemented for Stan",
         " (requires removing sphi from the parameters block)\n")
  } else if (!inherits(phi.sphi, "phi_spec_estimated")) {
    stop("phi.sphi must be an estimated() phi_spec object\n")
  }

  prior <- phi.sphi$prior
  if (is.null(prior)) {
    # Improper flat on (0, Inf)
    result$sphi_low        <- 0.0
    result$sphi_high       <- 1e10
    result$sphi_prior_type <- 0L
    result$sphi_prior_mean <- 0.0
    result$sphi_prior_sd   <- 1.0
  } else if (prior$dist == "uniform") {
    result$sphi_low        <- as.double(prior$low)
    result$sphi_high       <- as.double(prior$high)
    result$sphi_prior_type <- 0L   # bounds enforce uniform; no prior statement
    result$sphi_prior_mean <- 0.0
    result$sphi_prior_sd   <- 1.0
  } else if (prior$dist == "normal") {
    result$sphi_low        <- 0.0
    result$sphi_high       <- 1e10
    result$sphi_prior_type <- 1L
    result$sphi_prior_mean <- as.double(prior$mean)
    result$sphi_prior_sd   <- as.double(prior$sd)
  } else {
    stop("Only uniform and normal priors are currently supported for phi.sphi in Stan\n")
  }

  result
}


# Compute mphi in R given a phi_spec and a sphi value.
# Used to correctly initialise latent_phi in non-centered parameterisation.
.computeMPhiR <- function(phi.mphi, sphi) {
  if (inherits(phi.mphi, "phi_spec_fixed"))
    return(phi.mphi$value)
  if (!inherits(phi.mphi, "phi_spec_constrained"))
    stop("phi.mphi = estimated() not yet implemented\n")
  s2 <- sphi^2
  v  <- phi.mphi$value
  switch(phi.mphi$statistic,
    mean     = log(v) - 0.5 * s2,
    median   = log(v),
    mode     = log(v) + s2,
    variance = 0.5 * (log(v / (exp(s2) - 1)) - s2),
    sd       = 0.5 * (log(v^2 / (exp(s2) - 1)) - s2),
    stop("Unknown statistic: ", phi.mphi$statistic, "\n")
  )
}


# Estimate dM_prior_mean from the codon frequencies of low-SCUO genes.
# In low-expression genes dEta*phi ≈ 0, so observed frequencies ≈ exp(-dM)/Z.
# dM_k = log(count_ref + 0.5) - log(count_k + 0.5)  [Laplace pseudocount]
.dMPriorFromSCUO <- function(scuo, counts64, group_list, nonref_codons,
                               aa_start, aa_end, scuo.low.frac) {
  G     <- nrow(counts64)
  K     <- length(nonref_codons)
  n_low <- max(5L, floor(G * scuo.low.frac))
  n_low <- min(n_low, G)

  low_idx    <- order(scuo)[seq_len(n_low)]
  counts_low <- counts64[low_idx, , drop = FALSE]

  dM_mean <- numeric(K)
  for (a in seq_along(group_list)) {
    aa         <- group_list[a]
    all_codons <- AAToCodon(aa, FALSE)
    nonref     <- AAToCodon(aa, TRUE)
    ref_codon  <- setdiff(all_codons, nonref)

    count_ref  <- sum(counts_low[, ref_codon]) + 0.5

    for (k in aa_start[a]:aa_end[a]) {
      count_k    <- sum(counts_low[, nonref_codons[k]]) + 0.5
      dM_mean[k] <- log(count_ref) - log(count_k)
    }
  }
  dM_mean
}


# Map SCUO values to centered, scaled log-phi initialisation.
.scuoToLogPhi <- function(scuo, sphi.init) {
  raw <- log(pmax(scuo, 1e-3))
  .scaleLogPhi(raw, sphi.init)
}

# Center and scale a log-phi vector to mean=0, sd=sphi.init.
.scaleLogPhi <- function(log_raw, sphi.init) {
  centered <- log_raw - mean(log_raw)
  s        <- sd(centered)
  if (s > 0) centered / s * sphi.init else centered
}
