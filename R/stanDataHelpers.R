#' Build a Stan data list from an AnaCoDa Genome object
#'
#' Converts an AnaCoDa Genome object to the data list required by the ROC
#' Stan models (roc_arcsine.stan and roc_sphi_est.stan).  The output matches
#' the data block contract shared by both models; roc_arcsine.stan additionally
#' requires the \code{approx_min_n} field that is always included here.  To use
#' the output with roc_sphi_est.stan simply drop that field:
#' \code{d$approx_min_n <- NULL}.
#'
#' @param genome An AnaCoDa Genome object loaded via \code{initializeGenomeObject()}.
#' @param approx_min_n Integer.  Minimum total codon count per (gene, AA) pair
#'   for the arcsine approximation branch (default 20).  Has no effect when the
#'   data list is used with roc_sphi_est.stan.
#' @param noncentered Integer 0 or 1.  0 = centered phi parameterization
#'   (default; best when data strongly anchors phi).  1 = non-centered.
#' @param anchor_phi Integer 0 or 1.  0 = mean(phi) = 1 enforced via
#'   mphi = -sphi^2/2 (default).  1 = soft median(phi) ~ 1 via prior on
#'   mphi_param.
#' @param dM_prior_mean Prior mean for dM (mutation parameters).  Scalar
#'   (replicated to length K) or vector of length K.  Default 0.
#' @param dM_prior_sd Prior SD for dM.  Scalar or length-K vector.  Default 1.
#' @param dEta_prior_mean Prior mean for dEta (selection parameters).  Default 0.
#' @param dEta_prior_sd Prior SD for dEta.  Default 1.
#' @param sphi_prior_mean Prior mean for sphi (log-phi SD).  Default 1.
#' @param sphi_prior_sd Prior SD for sphi.  Default 2.
#' @param mphi_prior_sd Prior SD for the mphi soft anchor (used only when
#'   anchor_phi = 1).  Default 0.5.
#' @param grainsize reduce_sum grain size.  1 (default) lets TBB choose
#'   automatically; \code{ceiling(G / (2 * threads_per_chain))} is sometimes
#'   faster for small G.
#'
#' @return A named list suitable for passing directly to
#'   \code{mod$sample(data = ...)}, \code{mod$optimize(data = ...)}, or
#'   \code{mod$variational(data = ...)} via cmdstanr.
#'
#' @details
#' \strong{Data layout} (same as roc_sphi_est.stan / roc_arcsine.stan):
#' \itemize{
#'   \item \code{y_k[G, K]}: non-reference codon counts in codonArrayParameter
#'     order (K = 40 for the standard genetic code).
#'   \item \code{N_ga[G, A]}: total codon count per gene per AA group, summed
#'     over all synonymous codons including the reference.
#'   \item \code{aa_start[A]}, \code{aa_end[A]}: 1-indexed inclusive ranges
#'     into the flat K-vector for each AA group.
#' }
#'
#' The AA ordering follows \code{ROCParameter::groupList}:
#' A, C, D, E, F, G, H, I, K, L, N, P, Q, R, S, T, V, Y, Z  (A = 19 groups).
#' Met, Trp, and stop codons (M, W, X) are excluded (single-codon families).
#'
#' @examples
#' \dontrun{
#' genome <- initializeGenomeObject("mygenome.fasta")
#' d <- genomeToStanData(genome, approx_min_n = 20L)
#'
#' library(cmdstanr)
#' mod <- cmdstan_model("stan/roc_arcsine.stan",
#'                      cpp_options = list(stan_threads = TRUE))
#' fit <- mod$variational(data = d, threads = 4L)
#' }
#'
#' @export
genomeToStanData <- function(genome,
                              approx_min_n    = 20L,
                              noncentered     = 0L,
                              anchor_phi      = 0L,
                              dM_prior_mean   = 0.0,
                              dM_prior_sd     = 1.0,
                              dEta_prior_mean = 0.0,
                              dEta_prior_sd   = 1.0,
                              sphi_prior_mean = 1.0,
                              sphi_prior_sd   = 2.0,
                              mphi_prior_sd   = 0.5,
                              grainsize       = 1L) {

  # AA groups with >= 2 synonymous codons -- matches ROCParameter::groupList
  group_list <- c("A","C","D","E","F","G","H","I","K","L",
                  "N","P","Q","R","S","T","V","Y","Z")
  A <- length(group_list)  # 19

  # Non-reference codons per AA (forParamVector = TRUE excludes reference codon)
  nonref_by_aa  <- lapply(group_list, function(aa) AAToCodon(aa, TRUE))
  nonref_codons <- unlist(nonref_by_aa)  # flat K-vector (K = 40)
  K             <- length(nonref_codons)

  # 1-indexed inclusive index ranges into the K-vector, one range per AA
  aa_lengths <- vapply(nonref_by_aa, length, integer(1L))
  aa_end     <- cumsum(aa_lengths)
  aa_start   <- aa_end - aa_lengths + 1L

  # All 64-codon counts per gene (G x 64 data.frame, columns named by codon)
  counts64 <- getCodonCounts(genome)
  G        <- nrow(counts64)

  # y_k[G, K]: non-reference codon counts in parameter-vector order
  y_k <- as.matrix(counts64[, nonref_codons, drop = FALSE])
  storage.mode(y_k) <- "integer"
  dimnames(y_k)     <- NULL

  # N_ga[G, A]: total codon count per gene per AA (all K_a codons incl. reference)
  N_ga <- matrix(0L, nrow = G, ncol = A)
  for (a in seq_along(group_list)) {
    all_codons  <- AAToCodon(group_list[a], FALSE)
    aa_cols     <- counts64[, all_codons, drop = FALSE]
    N_ga[, a]   <- as.integer(rowSums(aa_cols))
  }

  # Helper: expand a scalar or K-vector prior parameter
  expand_k <- function(x, name) {
    if (length(x) == 1L) return(rep(as.double(x), K))
    if (length(x) == K)  return(as.double(x))
    stop(sprintf("'%s' must be length 1 or K = %d, got length %d",
                 name, K, length(x)))
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
    dM_prior_mean   = expand_k(dM_prior_mean,   "dM_prior_mean"),
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
