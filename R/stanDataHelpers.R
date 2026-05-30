#' Build a Stan data list from an AnaCoDa Genome object
#'
#' Converts an AnaCoDa Genome object to the data list required by the ROC
#' Stan models (roc_arcsine.stan and roc_sphi_est.stan).
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
#' @section Output not in Stan data block:
#' The returned list contains only fields declared in the Stan model's
#' \code{data \{\}} block.  Use \code{\link{genomeToStanInit}} to build the
#' matching \code{init} list (with SCUO-based phi initialisation) for
#' \code{mod$sample()} / \code{mod$optimize()} / \code{mod$variational()}.
#'
#' @param genome An AnaCoDa Genome object.
#' @param approx_min_n Integer. Minimum N per (gene, AA) for the arcsine
#'   branch (default 20).  Ignored by roc_sphi_est.stan.
#' @param noncentered Integer 0/1.  0 = centered phi (default), 1 = non-centered.
#' @param anchor_phi Integer 0/1.  0 = mean(phi)=1 via mphi=-sphi^2/2 (default),
#'   1 = soft median anchor.
#' @param dM.prior Character.  \code{"scuo"} (default): estimate
#'   \code{dM_prior_mean} from low-SCUO genes.  \code{"flat"}: use 0 for all.
#'   Ignored if \code{dM_prior_mean} is supplied explicitly.
#' @param scuo Numeric vector of length G, or NULL.  Pre-computed SCUO values
#'   (from \code{\link{calculateSCUO}}).  Computed internally if NULL.
#'   Passing pre-computed values avoids redundant work when also calling
#'   \code{\link{genomeToStanInit}}.
#' @param scuo.low.frac Fraction of lowest-SCUO genes used to estimate the
#'   dM prior means.  Default 0.25.  At least 5 genes are always used.
#' @param dM_prior_mean Prior mean for dM.  NULL (default): use \code{dM.prior}.
#'   Scalar or length-K vector: use as-is, ignoring \code{dM.prior}.
#' @param dM_prior_sd Prior SD for dM.  Scalar or length-K vector.  Default 1.
#' @param dEta_prior_mean Prior mean for dEta.  Default 0.
#' @param dEta_prior_sd Prior SD for dEta.  Default 1.
#' @param sphi_prior_mean Prior mean for sphi.  Default 1 (half-normal away
#'   from 0 avoids exp(-750)=0 underflow warnings in HMC).
#' @param sphi_prior_sd Prior SD for sphi.  Default 2.
#' @param mphi_prior_sd Prior SD for the mphi soft anchor (anchor_phi=1 only).
#' @param grainsize reduce_sum grain size (default 1; TBB auto-selects).
#'
#' @return Named list suitable for \code{mod$sample(data = ...)}, etc.
#'
#' @seealso \code{\link{genomeToStanInit}} for the matching init list.
#'
#' @examples
#' \dontrun{
#' genome <- initializeGenomeObject("mygenome.fasta")
#' scuo   <- calculateSCUO(genome)$SCUO          # compute once, reuse
#' d      <- genomeToStanData(genome, scuo = scuo)
#' init   <- genomeToStanInit(genome, d, scuo = scuo)
#'
#' library(cmdstanr)
#' mod  <- cmdstan_model("stan/roc_arcsine.stan",
#'                       cpp_options = list(stan_threads = TRUE))
#' fit  <- mod$variational(data = d, init = list(init), threads = 4L)
#' }
#'
#' @export
genomeToStanData <- function(genome,
                              approx_min_n    = 20L,
                              noncentered     = 0L,
                              anchor_phi      = 0L,
                              dM.prior        = "scuo",
                              scuo            = NULL,
                              scuo.low.frac   = 0.25,
                              dM_prior_mean   = NULL,
                              dM_prior_sd     = 1.0,
                              dEta_prior_mean = 0.0,
                              dEta_prior_sd   = 1.0,
                              sphi_prior_mean = 1.0,
                              sphi_prior_sd   = 2.0,
                              mphi_prior_sd   = 0.5,
                              grainsize       = 1L) {

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

  # ---- dM prior -------------------------------------------------------
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
    sphi_prior_mean = as.double(sphi_prior_mean),
    sphi_prior_sd   = as.double(sphi_prior_sd),
    noncentered     = as.integer(noncentered),
    anchor_phi      = as.integer(anchor_phi),
    mphi_prior_sd   = as.double(mphi_prior_sd),
    grainsize       = as.integer(grainsize)
  )
}


#' Build a Stan init list from an AnaCoDa Genome object
#'
#' Constructs the per-chain \code{init} list for cmdstanr's \code{sample()},
#' \code{optimize()}, or \code{variational()}.  By default initializes
#' \code{latent_phi} from per-gene SCUO, which places each gene's starting
#' synthesis rate in the right ballpark from codon-bias alone — without any
#' MCMC required.
#'
#' @section SCUO to log-phi mapping:
#' SCUO in [0,1] is mapped to log-phi via:
#' \enumerate{
#'   \item \eqn{\ell_g = \log(\max(\text{SCUO}_g, \epsilon))} with \eqn{\epsilon = 10^{-3}};
#'   \item Center: \eqn{\ell_g \leftarrow \ell_g - \bar{\ell}};
#'   \item Scale to target spread: \eqn{\ell_g \leftarrow \ell_g / \text{sd}(\ell) \times \text{sphi.init}}.
#' }
#' This produces a phi distribution with geometric mean \eqn{\approx 1} and
#' log-scale SD \eqn{\approx \text{sphi.init}}, consistent with the
#' \code{anchor_phi=0} prior constraint \eqn{E[\phi]=1}.
#'
#' For \code{noncentered=1}, the returned \code{latent_phi} is the z-score
#' \eqn{(log\phi - m_\phi) / \text{sphi.init}} rather than \eqn{log\phi}.
#'
#' @param genome An AnaCoDa Genome object.
#' @param data The data list returned by \code{\link{genomeToStanData}}.
#'   Provides \code{G}, \code{K}, and \code{noncentered}.
#' @param phi.init Character or numeric.  \code{"scuo"} (default): per-gene
#'   log-phi from SCUO.  \code{"uniform"}: \eqn{\phi=1} for all genes
#'   (\code{latent_phi = 0}).  Numeric vector of length G: treated as
#'   phi values and log-transformed.
#' @param scuo Numeric vector of length G, or NULL.  Pre-computed SCUO.
#'   Computed internally if NULL.
#' @param sphi.init Numeric.  Target SD for the initial log-phi distribution
#'   (default 1.0).  Also used as the initial value for \code{sphi}.
#'
#' @return A named list with fields \code{dM}, \code{dEta}, \code{latent_phi},
#'   \code{sphi}, \code{mphi_param} — one element per Stan parameter.
#'   Pass as \code{mod$sample(init = list(result))} (wrap in a list to supply
#'   the same init to all chains; replicate the list for per-chain variation).
#'
#' @seealso \code{\link{genomeToStanData}}
#'
#' @examples
#' \dontrun{
#' scuo <- calculateSCUO(genome)$SCUO
#' d    <- genomeToStanData(genome, scuo = scuo)
#' init <- genomeToStanInit(genome, d, scuo = scuo)
#' fit  <- mod$sample(data = d, init = list(init), chains = 4L, ...)
#' }
#'
#' @export
genomeToStanInit <- function(genome, data,
                              phi.init  = "scuo",
                              scuo      = NULL,
                              sphi.init = 1.0) {

  G <- data$G
  K <- data$K

  if (identical(phi.init, "scuo")) {
    if (is.null(scuo))
      scuo <- calculateSCUO(genome)$SCUO
    log_phi_init <- .scuoToLogPhi(scuo, sphi.init)

  } else if (is.numeric(phi.init) && length(phi.init) == G) {
    raw          <- log(pmax(phi.init, 1e-10))
    log_phi_init <- .scaleLogPhi(raw, sphi.init)

  } else {
    log_phi_init <- rep(0.0, G)   # phi = 1 everywhere
  }

  # For non-centered, convert log_phi to z-score
  if (data$noncentered == 1L) {
    mphi_0       <- -0.5 * sphi.init^2
    latent_init  <- (log_phi_init - mphi_0) / sphi.init
  } else {
    latent_init  <- log_phi_init
  }

  list(
    dM         = rep(0.0, K),
    dEta       = rep(0.0, K),
    latent_phi = latent_init,
    sphi       = sphi.init,
    mphi_param = 0.0
  )
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
#' \code{sphi} has a \code{lower=0} constraint; Stan maps it to unconstrained
#' space as \eqn{x = \log(\text{sphi})}, so its unconstrained variance is
#' \eqn{\text{Var}(\log(\text{sphi}_{\text{draws}}))}.  All other parameters
#' (\code{dM}, \code{dEta}, \code{latent\_phi}, \code{mphi\_param}) are
#' unconstrained and their ADVI variances are used directly.
#'
#' @section Stan parameter ordering:
#' The returned \code{inv_metric} vector follows Stan's internal unconstrained
#' parameter ordering for \code{roc_arcsine.stan} and \code{roc_sphi_est.stan}:
#' \code{dM[1..K]}, \code{dEta[1..K]}, \code{latent\_phi[1..G]},
#' \code{sphi}, \code{mphi\_param}.
#'
#' @param fit A \code{CmdStanVB} object from \code{mod$variational()}.
#' @param data The Stan data list returned by \code{\link{genomeToStanData}}.
#'
#' @return A named list:
#'   \describe{
#'     \item{init}{Named list of parameter starting values (posterior means in
#'       constrained space).  Pass as \code{mod$sample(init = list(result$init))}.}
#'     \item{inv_metric}{Numeric vector of unconstrained parameter variances
#'       (diagonal inverse mass matrix).  Pass as
#'       \code{mod$sample(inv_metric = result$inv_metric)}.}
#'   }
#'
#' @seealso \code{\link{genomeToStanData}}, \code{\link{genomeToStanInit}}
#'
#' @examples
#' \dontrun{
#' d    <- genomeToStanData(genome, scuo = scuo)
#' init <- genomeToStanInit(genome, d, scuo = scuo)
#' fit_vi <- mod_arcsine$variational(data = d, init = list(init), threads = 4L)
#' ws     <- adviToWarmStart(fit_vi, d)
#' # Stage 2: exact HMC warm-started from ADVI
#' fit_hmc <- mod_exact$sample(data    = d,
#'                              init    = list(ws$init),
#'                              inv_metric = ws$inv_metric,
#'                              chains  = 4L, ...)
#' }
#'
#' @export
adviToWarmStart <- function(fit, data) {
  draws <- fit$draws(format = "df")
  G <- data$G
  K <- data$K

  dM_cols     <- paste0("dM[",         1:K, "]")
  dEta_cols   <- paste0("dEta[",       1:K, "]")
  lphi_cols   <- paste0("latent_phi[", 1:G, "]")

  # ---- posterior means (constrained space) for init ----------------------
  dM_mean     <- as.numeric(colMeans(draws[, dM_cols,   drop = FALSE]))
  dEta_mean   <- as.numeric(colMeans(draws[, dEta_cols, drop = FALSE]))
  lphi_mean   <- as.numeric(colMeans(draws[, lphi_cols, drop = FALSE]))
  sphi_mean   <- mean(draws[["sphi"]])
  mphi_mean   <- mean(draws[["mphi_param"]])

  init <- list(
    dM         = dM_mean,
    dEta       = dEta_mean,
    latent_phi = lphi_mean,
    sphi       = sphi_mean,
    mphi_param = mphi_mean
  )

  # ---- unconstrained variances for inv_metric ----------------------------
  # dM, dEta, latent_phi, mphi_param: already unconstrained -> var directly
  # sphi: lower=0 -> Stan unconstrained = log(sphi) -> var(log(draws))
  dM_var   <- as.numeric(apply(draws[, dM_cols,   drop = FALSE], 2L, var))
  dEta_var <- as.numeric(apply(draws[, dEta_cols, drop = FALSE], 2L, var))
  lphi_var <- as.numeric(apply(draws[, lphi_cols, drop = FALSE], 2L, var))
  sphi_var <- var(log(draws[["sphi"]]))
  mphi_var <- var(draws[["mphi_param"]])

  # Stan parameter order: dM[1..K], dEta[1..K], latent_phi[1..G], sphi, mphi_param
  inv_metric <- c(dM_var, dEta_var, lphi_var, sphi_var, mphi_var)

  list(init = init, inv_metric = inv_metric)
}


# -------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------

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
    ref_codon  <- setdiff(all_codons, nonref)   # single reference codon

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

# Center and scale a log-phi vector to have mean=0 and sd=sphi.init.
.scaleLogPhi <- function(log_raw, sphi.init) {
  centered <- log_raw - mean(log_raw)
  s        <- sd(centered)
  if (s > 0) centered / s * sphi.init else centered
}
