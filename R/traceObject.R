#### TODO, lets move it into parameterObject.R and use a parameter instead of trace. thats how it is done for the acf function

# Internal: extract a samples x nseries matrix (or samples vector) of trace
# values for one of the named groups. Used by as.mcmc.Rcpp_Trace and
# convergence.test.Rcpp_Trace.
.extract_trace_matrix <- function(object, what = "Mutation", mixture = 1)
{
  current.trace <- NULL
  if(what[1] == "Mutation" || what[1] == "Selection")
  {
    names.aa <- aminoAcids()
    numCodons <- 0
    for(aa in names.aa)
    {
      if (aa == "M" || aa == "W" || aa == "X") next
      codons <- AAToCodon(aa, T)
      numCodons <- numCodons + length(codons)
    }
    index <- 1
    cur.trace <- vector("list", numCodons)
    for(aa in names.aa)
    {
      if (aa == "M" || aa == "W" || aa == "X") next
      codons <- AAToCodon(aa, T)
      for(i in 1:length(codons))
      {
        if(what[1] == "Mutation"){
          cur.trace[[index]]<- object$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 0, T)
        }else{
          cur.trace[[index]] <- object$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codons[i], 1, T)
        }
        index <- index + 1
      }
    }
    ## Transpose matrix to get in correct format for coda::mcmc. Transposing results in same output from coda::geweke.test as performing the test separately on each codon specific parameter
    current.trace <- t(do.call("rbind", cur.trace))
  }
  else if(what[1] == "Alpha" || what[1] == "Lambda" || what[1] == "NSERate" || what[1] == "LambdaPrime")
  {
    codon.list <- codons()
    codon.list <- codon.list[1:(length(codon.list)-3)]
    cur.trace <- vector("list",length(codon.list))
    for (i in 1:length(codon.list))
    {
      if (what[1]=="Alpha")
      {
        cur.trace[[i]]<- object$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codon.list[i], 0, F)
      } else if (what[1]=="Lambda" || what[1]=="LambdaPrime"){
        cur.trace[[i]]<- object$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codon.list[i], 1, F)
      } else if (what[1]=="NSERate"){
        cur.trace[[i]]<- object$getCodonSpecificParameterTraceByMixtureElementForCodon(mixture, codon.list[i], 2, F)
      }
    }
    current.trace <- t(do.call("rbind", cur.trace))
  }
  else if(what[1] == "MixtureProbability")
  {
    numMixtures <- object$getNumberOfMixtures()
    cur.trace <- vector("list", numMixtures)
    for(i in 1:numMixtures)
    {
      cur.trace[[i]] <- object$getMixtureProbabilitiesTraceForMixture(i)
    }
    current.trace <- t(do.call("rbind", cur.trace))
  }
  else if(what[1] == "Sphi")
  {
    sphi <- object$getStdDevSynthesisRateTraces()
    current.trace <- t(do.call("rbind", sphi))
  }
  else if(what[1] == "Mphi")
  {
    sphi <- do.call("rbind", object$getStdDevSynthesisRateTraces())
    current.trace <- t(-(sphi * sphi) / 2)
  }
  else if(what[1] == "ExpectedPhi")
  {
    current.trace <- object$getExpectedSynthesisRateTrace()
  }
  else if(what[1] == "InitiationCost")
  {
    current.trace <- object$getInitiationCostTrace()
  }
  else if(what[1] == "AcceptanceCSP")
  {
    names.aa <- aminoAcids()
    index <- 1
    cur.trace <- vector("list", length(names.aa) - length(c("M","W","X")))
    for(aa in names.aa)
    {
      if (aa == "M" || aa == "W" || aa == "X") next
      cur.trace[[index]] <- object$getCodonSpecificAcceptanceRateTraceForAA(aa)
      index <- index + 1
    }
    current.trace <- t(do.call("rbind", cur.trace))
  }
  else if(what[1] %in% c("Aphi", "Sepsilon", "Expression"))
  {
    stop("convergence/as.mcmc for what=\"", what[1], "\" is not yet implemented")
  }
  else
  {
    stop("unknown `what`: ", what[1])
  }
  current.trace
}

#' Coerce an AnaCoDa Trace object to a \code{coda::mcmc} object
#'
#' @param x an \code{Rcpp_Trace} object (from \code{parameter$getTraceObject()}).
#' @param what which set of traces to extract. One of \code{"Mutation"} (default),
#'   \code{"Selection"}, \code{"Alpha"}, \code{"Lambda"}/\code{"LambdaPrime"},
#'   \code{"NSERate"}, \code{"MixtureProbability"}, \code{"Sphi"}, \code{"Mphi"},
#'   \code{"ExpectedPhi"}, \code{"InitiationCost"}, or \code{"AcceptanceCSP"}.
#' @param mixture mixture index for traces that are mixture-specific
#'   (Mutation/Selection/Alpha/Lambda/NSERate). Defaults to 1.
#' @param samples optional positive integer. If supplied, return only the last
#'   \code{samples} rows of the trace (clamped to the trace length). \code{NULL}
#'   (the default) returns the full trace.
#' @param thin thinning interval recorded as metadata on the returned mcmc object;
#'   does not subsample.
#' @param ... unused; present for S3 generic compatibility.
#'
#' @return a \code{coda::mcmc} object: vector for single-series traces
#'   (ExpectedPhi, InitiationCost), matrix with one column per codon/mixture for
#'   the others.
#'
#' @note Prefer the \code{samples} argument over post-hoc windowing with
#'   \code{utils::tail()}. \code{tail()} on an \code{mcmc} object dispatches to
#'   \code{coda::tail.mcmc}, which preserves the original iteration index in the
#'   \code{mcpar} attribute. Downstream diagnostics that read \code{start}/\code{end}
#'   metadata (e.g. \code{coda::geweke.diag}) then use different window boundaries
#'   than they would for a fresh-iteration matrix. Using \code{samples} here trims
#'   the raw trace before wrapping, so \code{mcpar} starts at iteration 1 either way.
#'
#' @export
as.mcmc.Rcpp_Trace <- function(x, what = "Mutation", mixture = 1, samples = NULL,
                               thin = 1, ...)
{
  if (!is.null(samples) && samples < 1)
    stop("`samples` must be a positive integer or NULL")
  trace <- .extract_trace_matrix(x, what, mixture)
  if (!is.null(samples))
    trace <- utils::tail(trace, n = min(samples, NROW(trace)))
  coda::mcmc(data = trace, thin = thin)
}

# see mcmc Object.R convergence.test function for documentation
convergence.test.Rcpp_Trace <- function(object, samples = NULL, frac1 = 0.1,
                                        frac2 = 0.5, thin = 1, plot = FALSE, what = "Mutation", mixture = 1)
{
  mcmcobj <- as.mcmc(object, what = what, mixture = mixture,
                     samples = samples, thin = thin)
  if(plot){
    coda::geweke.plot(mcmcobj, frac1=frac1, frac2=frac2)
  } else{
    coda::geweke.diag(mcmcobj, frac1=frac1, frac2=frac2)
  }
}
