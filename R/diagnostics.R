#' Assert that MCMC chains are actually moving
#'
#' @description
#' Checks per-proposal-type acceptance rates (AR) in the trace object stored
#' inside a Parameter object after a run (or chunk) of MCMC.  If the maximum
#' AR across all groups for any checked proposal type falls below
#' \code{threshold}, the function fires a diagnostic: either a warning (the
#' default) or an error, naming the stuck type(s) so the operator can
#' diagnose root cause (misconfigured proposal covariance, likelihood too
#' steep at initialisation, etc.).
#'
#' The function is model-aware.  It detects the concrete parameter class
#' (\code{Rcpp_ROCParameter}, \code{Rcpp_FONSEParameter},
#' \code{Rcpp_PAParameter}, \code{Rcpp_PANSEParameter}) and dispatches to
#' the appropriate per-model collector.  Users call the single public
#' interface regardless of model.
#'
#' @details
#' **Proposal types checked per model**
#'
#' | Model  | Elongation (per-AA AR) | NSERate (per-codon AR) |
#' |--------|------------------------|------------------------|
#' | ROC    | Yes                    | No                     |
#' | FONSE  | Yes                    | No                     |
#' | PA     | Yes (per-codon)        | No                     |
#' | PANSE  | Yes (per-codon)        | Yes                    |
#'
#' For ROC and FONSE, "Elongation" means the per-amino-acid
#' CSP (mutation/selection) proposal group.  For PA and PANSE, it means
#' the per-codon Alpha/LambdaPrime (elongation shape/rate) proposal group.
#' The NSERate type is exclusive to PANSE.
#'
#' A type is declared stuck when it has at least one observable group AND
#' \code{max(AR)} across all groups is below \code{threshold}.  This means
#' one legitimately low-AR group does NOT trigger the alarm if other groups
#' in the same type are moving.
#'
#' **Why 0.001 as the default threshold?**
#' An AR of 0.001 means roughly 1 accepted proposal per 1000 iterations.
#' At typical thinning/adaptiveWidth settings this is effectively zero: the
#' chain is not exploring.  The value is conservative enough that it should
#' not produce false positives for legitimately sparse proposal groups (e.g.
#' high-constraint amino acids with few synonymous codons).
#'
#' **Validation history**
#' The PANSE implementation was first validated against a live silent-failure
#' case in the PANSE.data.analyses adapter stack (commit 940de04 in that
#' repo).  That catch motivated the generalisation here.
#'
#' @param parameter A live Rcpp parameter object returned by
#'   \code{initializeParameterObject()} or loaded from a restart file.
#'   Supported classes: \code{Rcpp_ROCParameter},
#'   \code{Rcpp_FONSEParameter}, \code{Rcpp_PAParameter},
#'   \code{Rcpp_PANSEParameter}.
#' @param threshold Numeric scalar (default \code{0.001}).  The maximum
#'   per-group AR below which a proposal type is declared stuck.
#' @param action Character, either \code{"warn"} (default) or
#'   \code{"stop"}.  Controls whether a stuck chain triggers
#'   \code{warning()} or \code{stop()}.
#' @param context Optional character string inserted into the diagnostic
#'   message for orientation (e.g. a chunk tag like \code{"chunk_003"} or
#'   a timestamp).  Ignored if \code{""} (the default).
#'
#' @return Invisibly returns a named list:
#' \describe{
#'   \item{stuck}{Logical.  \code{TRUE} if any proposal type was stuck.}
#'   \item{stuck.types}{Character vector of stuck proposal type names.}
#'   \item{ar}{Named list mapping each checked proposal type to a numeric
#'     vector of per-group acceptance rates (last recorded value per group).}
#' }
#'
#' @examples
#' \dontrun{
#' ## After a short ROC run ---------------------------------------------------
#' genome    <- initializeGenomeObject(file = "genome.fasta")
#' parameter <- initializeParameterObject(genome, sphi_init = 1,
#'                                        num.mixtures = 1,
#'                                        gene.assignment = rep(1, length(genome)))
#' mcmc      <- initializeMCMCObject(samples = 1000)
#' model     <- initializeModelObject(parameter, "ROC")
#' runMCMC(mcmc, genome, model, ncores = 1)
#'
#' result <- assertChainsMoving(parameter)
#' if (result$stuck) {
#'   message("Stuck types: ", paste(result$stuck.types, collapse = ", "))
#' }
#'
#' ## After a PANSE chunk, stop immediately on any stuck chain ----------------
#' assertChainsMoving(panse_parameter, action = "stop", context = "chunk_001")
#' }
#'
#' @export
assertChainsMoving <- function(parameter,
                               threshold = 0.001,
                               action    = c("warn", "stop"),
                               context   = "") {
  action <- match.arg(action)

  cls <- class(parameter)

  ## Dispatch to the per-model collector.  Each collector returns a named list
  ## mapping proposal-type names to numeric vectors of per-group ARs.
  ar.out <- if (cls %in% c("Rcpp_ROCParameter", "Rcpp_FONSEParameter")) {
    .collectAR.rocfonse(parameter)
  } else if (cls == "Rcpp_PAParameter") {
    .collectAR.pa(parameter)
  } else if (cls == "Rcpp_PANSEParameter") {
    .collectAR.panse(parameter)
  } else {
    warning("assertChainsMoving: unrecognised parameter class '", cls,
            "'.  No AR check performed.", call. = FALSE)
    return(invisible(list(stuck = FALSE, stuck.types = character(0L),
                          ar = list())))
  }

  ## Evaluate stuckness per type: must have data AND max(AR) < threshold.
  stuck.flags <- vapply(names(ar.out), function(nm) {
    v <- ar.out[[nm]]
    length(v) > 0L && max(v, na.rm = TRUE) < threshold
  }, logical(1L))

  any.stuck   <- any(stuck.flags)
  stuck.types <- names(stuck.flags)[stuck.flags]

  if (any.stuck) {
    ## Build a compact AR summary; mark stuck types with '*'.
    ar.summary <- paste(vapply(names(ar.out), function(nm) {
      v    <- ar.out[[nm]]
      mark <- if (stuck.flags[[nm]]) "*" else ""
      sprintf("%s=%.4f%s", nm, max(v, na.rm = TRUE), mark)
    }, character(1L)), collapse = ", ")

    msg <- paste0(
      "STUCK CHAIN(S) DETECTED",
      if (nchar(context) > 0L) paste0(" [", context, "]") else "",
      ": max per-group AR [", ar.summary, "]",
      " (*) < threshold=", threshold, ". ",
      "Stuck type(s): ", paste(stuck.types, collapse = ", "), ". ",
      "Chain is not exploring -- parameter(s) may be effectively frozen."
    )
    if (action == "stop") stop(msg) else warning(msg, immediate. = TRUE,
                                                  call. = FALSE)
  }

  invisible(list(stuck = any.stuck, stuck.types = stuck.types, ar = ar.out))
}


## ---------------------------------------------------------------------------
## Private per-model AR collectors
##
## Each returns a named list: list(TypeName = numeric_vector_of_per_group_AR).
## Vectors contain the last recorded AR value for each active proposal group.
## NA entries (missing trace data for a group) are silently dropped.
## ---------------------------------------------------------------------------

## ROC and FONSE: CSP proposals are per-amino-acid.
## The trace is indexed by AA; getCodonSpecificAcceptanceRateTraceForAA(aa)
## returns a vector of stored AR values (one per adaptiveWidth block).
## We take the last element as the current AR.
.collectAR.rocfonse <- function(parameter) {
  trace <- parameter$getTraceObject()
  aas   <- aminoAcids()  # all 22 one-letter AA codes inc. "X"

  elong.ar <- vapply(aas, function(aa) {
    v <- tryCatch(
      trace$getCodonSpecificAcceptanceRateTraceForAA(aa),
      error = function(e) numeric(0L)
    )
    if (length(v) == 0L) NA_real_ else as.numeric(tail(v, 1L))
  }, numeric(1L))

  elong.ar <- elong.ar[!is.na(elong.ar)]
  list(Elongation = elong.ar)
}

## PA: CSP proposals are per-codon (64 sense codons + stop codons mapped to
## 64 total entries, but PA's groupList is exactly the codons it estimates).
## Use parameter$getGroupList() to get the active codon set so we respect
## any custom groupList the user may have set.
.collectAR.pa <- function(parameter) {
  trace    <- parameter$getTraceObject()
  grp.list <- parameter$getGroupList()  # character vector of codon strings

  elong.ar <- vapply(grp.list, function(codon) {
    v <- tryCatch(
      trace$getCodonSpecificAcceptanceRateTraceForCodon(codon),
      error = function(e) numeric(0L)
    )
    if (length(v) == 0L) NA_real_ else as.numeric(tail(v, 1L))
  }, numeric(1L))

  elong.ar <- elong.ar[!is.na(elong.ar)]
  list(Elongation = elong.ar)
}

## PANSE: same Elongation check as PA (per-codon, 61 sense codons),
## plus a separate NSERate check via getNseRateSpecificAcceptanceRateTrace().
.collectAR.panse <- function(parameter) {
  trace    <- parameter$getTraceObject()
  grp.list <- parameter$getGroupList()

  ## Elongation (Alpha/LambdaPrime) proposals -- per-codon.
  elong.ar <- vapply(grp.list, function(codon) {
    v <- tryCatch(
      trace$getCodonSpecificAcceptanceRateTraceForCodon(codon),
      error = function(e) numeric(0L)
    )
    if (length(v) == 0L) NA_real_ else as.numeric(tail(v, 1L))
  }, numeric(1L))
  elong.ar <- elong.ar[!is.na(elong.ar)]

  ## NSERate proposals -- list of per-codon scalars (one per sense codon).
  nse.raw <- tryCatch(
    trace$getNseRateSpecificAcceptanceRateTrace(),
    error = function(e) list()
  )
  ## Each element is a vector of stored AR values; take the last.
  nse.ar <- if (length(nse.raw) > 0L) {
    vapply(nse.raw, function(v) {
      if (length(v) == 0L) NA_real_ else as.numeric(tail(v, 1L))
    }, numeric(1L))
  } else {
    numeric(0L)
  }
  nse.ar <- nse.ar[!is.na(nse.ar)]

  list(Elongation = elong.ar, NSERate = nse.ar)
}
