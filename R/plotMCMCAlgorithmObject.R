#' Plot MCMC algorithm
#' 
#' @param x An Rcpp_MCMC object initialized with \code{initializeMCMCObject}.
#' 
#' @param zoom.window A vector describing the start and end of the zoom window.
#' 
#' @param what character defining if log(Posterior) (Default) or log(Likelihood) 
#' options are: LogPosterior or logLikelihood
#' 
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' 
#' @return This function has no return value.
#' 
#' @description This function will plot the logLikelihood trace, and if the Hmisc package is installed, it will 
#'  plot a subplot of the logLikelihood trace with the first few samples removed.
plot.Rcpp_MCMCAlgorithm <- function(x, what = "LogPosterior", zoom.window = NULL, ...)
{
  if(what[1] == "LogPosterior")
  {
    trace <- x$getLogPosteriorTrace()
    ylab = "log(Posterior Probability)"
  }else{
    trace <- x$getLogLikelihoodTrace()
    ylab = "log(Likelihood Probability)"
  }
  trace <- trace[-1]
  
  trace.length <- length(trace)
  
  zoomStart <- round(0.9*trace.length)
  zoomEnd <- trace.length
  logL <- mean(trace[zoomStart:trace.length])
  #TODO change main title
  plot(trace, type="l", main=paste0(ylab, ": ", logL), xlab="Sample", ylab=ylab)
  grid (NULL,NULL, lty = 6, col = "cornsilk2")
  trace[trace == -Inf] <- NA
  
  # TODO (Cedric): get rid of that line once problem with first element beeing 0 is solved
  trace <- trace[-1]
  
  if(!(is.null(zoom.window))) {
    zoomStart <- zoom.window[1]
    zoomEnd <- zoom.window[2]
  }
  else{
    warning("No window was given, zooming in at last 10% of trace")
  }
  
  Hmisc::subplot(
    plot(zoomStart:zoomEnd, trace[zoomStart:zoomEnd], type="l", xlab=NA, ylab=NA, las=2, cex.axis=0.55), 
    0.8*(round(0.9*trace.length)), (min(trace, na.rm = T)+max(trace, na.rm = T))/2, size=c(3,2))
}

#' Autocorrelation function for the likelihood or posterior trace
#'
#' @param mcmc object of class MCMC
#' @param what character vector of traces to include. Any combination of
#'   "LogPosterior" and "LogLikelihood". Defaults to both.
#' @param samples number of samples at the end of the trace used to calculate the acf.
#'   Defaults to 10*log10(N) where N is the trace length.
#' @param lag.max maximum lag for acf calculation. Default is 40.
#' @param plot logical. If TRUE (default) both traces are overlaid on a single plot,
#'   color-coded and marked with distinct point symbols.
#'
#' @return Invisibly returns a named list of acf objects, one per requested trace.
#'
#' @description Calculates and (by default) plots the ACF of the log-posterior and/or
#'   log-likelihood trace on a single combined panel.
#'
#' @seealso \code{\link{acfCSP}}
#'
acfMCMC <- function(mcmc, what = c("LogPosterior", "LogLikelihood"), samples = NULL, lag.max = 40, plot = TRUE)
{
  what <- match.arg(what, choices = c("LogPosterior", "LogLikelihood"), several.ok = TRUE)

  acf.list <- list()
  for (w in what) {
    tr <- if (w == "LogPosterior") mcmc$getLogPosteriorTrace() else mcmc$getLogLikelihoodTrace()
    if (is.null(samples)) { samples <- round(10 * log10(length(tr))) }
    acf.list[[w]] <- acf(x = tail(tr, samples), lag.max = lag.max, plot = FALSE)
  }

  if (plot) {
    lags    <- as.numeric(acf.list[[1]]$lag)
    acf.mat <- do.call("cbind", lapply(acf.list, function(a) as.numeric(a$acf)))
    ci      <- qnorm(0.975) / sqrt(samples)
    ylim    <- c(min(-ci * 2, min(acf.mat[-1, , drop = FALSE], na.rm = TRUE) - 0.02), 1.05)

    plot(NULL, NULL, xlim = range(lags), ylim = ylim,
         xlab = "Lag", ylab = "Autocorrelation",
         main = paste(paste(what, collapse = " & "), "Trace ACF"))
    abline(h = 0,          lty = 1, col = "gray70")
    abline(h = c(-ci, ci), lty = 2, col = "gray50")

    for (i in seq_along(what)) {
      lines(lags, acf.mat[, i], col = .mixtureColors[i], type = "o", pch = i, cex = 0.6)
    }
    legend("topright", legend = what, col = .mixtureColors[seq_along(what)],
           lty = 1, pch = seq_along(what), bty = "n")
  }

  invisible(acf.list)
}
