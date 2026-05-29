#'  Model Initialization
#'
#' @param parameter An object created with \code{initializeParameterObject}.
#'
#' @param model A string containing the model to run (ROC, FONSE, or PA), has to match parameter object.
#'
#' @param with.phi (ROC only) A boolean that determines whether or not to include empirical
#'    phi values (expression rates) for the calculations. Default value is FALSE
#'
#' @param fix.observation.noise (ROC only) Allows fixing the noise term sepsilon in the observed expression dataset to its initial condition.  This value should override the est.hyper=TRUE setting in \code{initializeMCMCObject()}
#'	The initial condition for the observed expression noise is set in the parameter object. Default value is FALSE.
#'
#' @param rfp.count.column (PA and PANSE only) A number representing the RFP count column to use. Default value is 1.
#'
#' @param approx (ROC only) Likelihood approximation method. \code{FALSE} (default) uses the
#'   exact multinomial log-likelihood. \code{TRUE} or \code{"hybrid.arcsine"} uses a
#'   variance-stabilising arcsine approximation for amino acid / gene combinations whose
#'   average expected codon count (\eqn{n / K}, where \eqn{K} is the number of synonymous
#'   codons) meets or exceeds \code{approx.min.expected}; exact multinomial is used as a
#'   fallback otherwise.  The arcsine approximation is at least as accurate as the normal
#'   approximation and is valid down to average expected counts of ~2--3.
#'
#' @param approx.min.expected (ROC only, used when \code{approx != FALSE}) Minimum average
#'   expected count per synonymous codon (\eqn{n / K}) required to apply the arcsine
#'   approximation.  Amino acid / gene combinations below this threshold fall back to the
#'   exact multinomial.  Default is \code{5.0}; lower values (e.g. \code{2.0}) extend
#'   coverage at the cost of somewhat larger approximation error.
#'
#' @return This function returns the model object created.
#'
#' @description initializes the model object.
#'
#' @details initializeModelObject initializes a model. The type of model is determined based on the string passed to the \code{model} argument.
#'  The Parameter object has to match the model that is initialized. E.g. to initialize a ROC model,
#'  it is required that a ROC parameter object is passed to the function.
#'
#' @examples
#'
#' #initializing a model object
#'
#' genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
#' expression_file <- system.file("extdata", "expression.csv", package = "AnaCoDa")
#'
#' genome <- initializeGenomeObject(file = genome_file,
#'                                  observed.expression.file = expression_file)
#' sphi_init <- c(1,1)
#' numMixtures <- 2
#' geneAssignment <- c(rep(1,floor(length(genome)/2)),rep(2,ceiling(length(genome)/2)))
#' parameter <- initializeParameterObject(genome = genome, sphi = sphi_init,
#'                                        num.mixtures = numMixtures,
#'                                        gene.assignment = geneAssignment,
#'                                        mixture.definition = "allUnique")
#'
#' # initializing a model object assuming we have observed expression (phi)
#' # values stored in the genome object.
#' initializeModelObject(parameter = parameter, model = "ROC", with.phi = TRUE)
#'
#' # initializing a model object ignoring observed expression (phi)
#' # values stored in the genome object.
#' initializeModelObject(parameter = parameter, model = "ROC", with.phi = FALSE)
#'
#' # initializing a ROC model with the hybrid arcsine approximation
#' initializeModelObject(parameter = parameter, model = "ROC", approx = TRUE)
#'
.approxToInt <- function(approx) {
  if (isFALSE(approx) || identical(approx, "exact"))          return(0L)
  if (isTRUE(approx)  || identical(approx, "hybrid.arcsine")) return(1L)
  stop("approx must be FALSE, TRUE, \"exact\", or \"hybrid.arcsine\"")
}

initializeModelObject <- function(parameter, model = "ROC", with.phi = FALSE,
                                  fix.observation.noise = FALSE, rfp.count.column = 1,
                                  approx = FALSE, approx.min.expected = 5.0) {
  if (model == "ROC") {
    approx_int <- .approxToInt(approx)
    c.model <- new(ROCModel, with.phi, fix.observation.noise, approx_int, approx.min.expected)
  } else if (model == "FONSE") {
    if (!isFALSE(approx))
      warning("approx is not yet implemented for FONSE; using exact multinomial.")
    c.model <- new(FONSEModel, with.phi, fix.observation.noise)
  } else if (model == "PA") {
    if (!isFALSE(approx))
      warning("approx is not applicable to the PA model.")
    c.model <- new(PAModel, rfp.count.column, with.phi, fix.observation.noise)
  } else if (model == "PANSE") {
    if (!isFALSE(approx))
      warning("approx is not applicable to the PANSE model.")
    c.model <- new(PANSEModel, rfp.count.column, with.phi, fix.observation.noise)
  } else {
    stop("Unknown model.")
  }
  c.model$setParameter(parameter)
  return(c.model)
}
