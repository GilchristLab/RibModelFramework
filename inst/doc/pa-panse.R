## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----echo = TRUE, eval = FALSE------------------------------------------------
# rfp_file <- system.file("extdata", "pa_rfpdata.csv", package = "AnaCoDa")


## ----echo = TRUE, eval = FALSE------------------------------------------------
# library(AnaCoDa)
# 
# rfp_file <- system.file("extdata", "pa_rfpdata.csv", package = "AnaCoDa")
# genome <- initializeGenomeObject(file = rfp_file, fasta = FALSE)


## ----echo = TRUE, eval = FALSE------------------------------------------------
# sphi_init    <- 2
# numMixtures  <- 1
# geneAssignment <- rep(1, length(genome))
# 
# parameter <- initializeParameterObject(
#   genome          = genome,
#   sphi            = sphi_init,
#   num.mixtures    = numMixtures,
#   gene.assignment = geneAssignment,
#   model           = "PA",
#   split.serine    = TRUE,
#   mixture.definition = "allUnique"
# )


## ----echo = TRUE, eval = FALSE------------------------------------------------
# model <- initializeModelObject(parameter, "PA")
# 
# mcmc <- initializeMCMCObject(
#   samples        = 1000,
#   thinning       = 10,
#   adaptive.width = 50,
#   est.expression = TRUE,
#   est.csp        = TRUE,
#   est.hyper      = TRUE
# )


## ----echo = TRUE, eval = FALSE------------------------------------------------
# runMCMC(mcmc, genome, model)


## ----echo = TRUE, eval = FALSE------------------------------------------------
# writeParameterObject(parameter, file = "pa_parameter.Rdat")
# writeMCMCObject(mcmc, file = "pa_mcmc.Rdat")
# 
# parameter <- loadParameterObject(file = "pa_parameter.Rdat")
# mcmc      <- loadMCMCObject(file = "pa_mcmc.Rdat")


## ----echo = TRUE, eval = FALSE------------------------------------------------
# logPost <- mcmc$getLogPosteriorTrace()
# burnin  <- ceiling(length(logPost) * 0.3)
# logPost_post_burnin <- logPost[burnin:length(logPost)]
# 
# plot(logPost_post_burnin, type = "l",
#      main = paste("log(Posterior), mean =", round(mean(logPost_post_burnin), 1)),
#      xlab = "Sample", ylab = "log(Posterior)")
# grid(NULL, NULL, lty = 6, col = "cornsilk2")


## ----echo = TRUE, eval = FALSE------------------------------------------------
# trace <- parameter$getTraceObject()
# plot(trace, what = "Sphi",             main = "s_phi hyperparameter")
# plot(trace, what = "MixtureProbability", main = "Mixture probability")
# plot(trace, what = "ExpectedPhi",      main = "Expected expression (psi)")


## ----echo = TRUE, eval = FALSE------------------------------------------------
# plot(trace, what = "Alpha",       mixture = 1)
# plot(trace, what = "LambdaPrime", mixture = 1)


## ----echo = TRUE, eval = FALSE------------------------------------------------
# samples    <- 1000
# codonList  <- codons()          # 61 sense codons
# cat        <- 1                 # mixture category index
# burnin_idx <- samples * 0.5    # use posterior half
# 
# alpha_mean <- numeric(61)
# alpha_ci   <- matrix(0, nrow = 61, ncol = 2)
# lp_mean    <- numeric(61)
# lp_ci      <- matrix(0, nrow = 61, ncol = 2)
# 
# for (i in seq_along(codonList)) {
#   codon <- codonList[i]
# 
#   alpha_mean[i] <- parameter$getCodonSpecificPosteriorMean(cat, burnin_idx, codon, 0, FALSE)
#   alpha_trace   <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 0, FALSE)
#   alpha_ci[i, ] <- quantile(alpha_trace[burnin_idx:samples], probs = c(0.025, 0.975))
# 
#   lp_mean[i]  <- parameter$getCodonSpecificPosteriorMean(cat, burnin_idx, codon, 1, FALSE)
#   lp_trace    <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 1, FALSE)
#   lp_ci[i, ]  <- quantile(lp_trace[burnin_idx:samples], probs = c(0.025, 0.975))
# }
# 
# # Plot alpha credible intervals
# plot(NULL, NULL,
#      xlim = c(1, 61), ylim = range(alpha_ci),
#      main = "Posterior mean and 95% CI: Alpha",
#      xlab = "Codon", ylab = "Alpha", axes = FALSE)
# confidenceInterval.plot(x = 1:61, y = alpha_mean, sd.y = alpha_ci)
# axis(2)
# axis(1, at = 1:61, labels = codonList, tck = 0.02, las = 2, cex.axis = 0.6)
# 
# # Plot lambda-prime credible intervals
# plot(NULL, NULL,
#      xlim = c(1, 61), ylim = range(lp_ci),
#      main = "Posterior mean and 95% CI: LambdaPrime",
#      xlab = "Codon", ylab = "LambdaPrime", axes = FALSE)
# confidenceInterval.plot(x = 1:61, y = lp_mean, sd.y = lp_ci)
# axis(2)
# axis(1, at = 1:61, labels = codonList, tck = 0.02, las = 2, cex.axis = 0.6)


## ----echo = TRUE, eval = FALSE------------------------------------------------
# pausing_times <- alpha_mean * lp_mean
# barplot(pausing_times, names.arg = codonList,
#         main = "Estimated pausing times per codon",
#         las = 2, cex.names = 0.6)


## ----echo = TRUE, eval = FALSE------------------------------------------------
# parameter_nse <- initializeParameterObject(
#   genome          = genome,
#   sphi            = 2,
#   num.mixtures    = 1,
#   gene.assignment = rep(1, length(genome)),
#   model           = "PANSE",
#   split.serine    = TRUE,
#   mixture.definition = "allUnique"
# )
# 
# model_nse <- initializeModelObject(parameter_nse, "PANSE")
# 
# mcmc_nse <- initializeMCMCObject(
#   samples        = 1000,
#   thinning       = 10,
#   adaptive.width = 50,
#   est.expression = TRUE,
#   est.csp        = TRUE,
#   est.hyper      = TRUE
# )
# 
# runMCMC(mcmc_nse, genome, model_nse)


## ----echo = TRUE, eval = FALSE------------------------------------------------
# setRestartSettings(mcmc, file = "pa_restart.rst",
#                    every = 50, initialize = TRUE)
# runMCMC(mcmc, genome, model)


## ----echo = TRUE, eval = FALSE------------------------------------------------
# parameter <- initializeParameterObject(model = "PA",
#                                         restart.file = "pa_restart.rst")

