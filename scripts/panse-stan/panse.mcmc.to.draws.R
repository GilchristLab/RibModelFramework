## panse.mcmc.to.draws.R -- convert AnaCoDa native MCMC PANSE traces to
## posterior::draws_array, using the same variable naming convention as the
## AnaCoDa Stan PANSE models so both backends feed the same evaluation and
## diagnostic tools.
##
## Public function:
##   panse_mcmc_to_draws(parameter, keep.fraction, codon.names, gene.names)
##
## Variable names (match AnaCoDa Stan panse_sphi_est_noncentered):
##   alpha[c]        per-codon alpha        (natural scale)
##   lambdaPrime[c]  per-codon lambda       (natural scale)
##   NSERate[c]      per-codon NSE rate     (or NSERate_shared when n_nse == 1)
##   log_phi[g]      log synthesis rate     (log of AnaCoDa's phi)
##   sphi            synthesis-rate SD

panse_mcmc_to_draws <- function(parameter,
                                keep.fraction = 0.5,
                                codon.names   = NULL,
                                gene.names    = NULL) {
    stopifnot(requireNamespace("posterior", quietly = TRUE))  # pre-existing cmdstanr dependency

    trace <- parameter$getTraceObject()

    alpha_tr  <- trace$getCodonSpecificParameterTrace(0L)  # [[mix]][[codon]]
    lambda_tr <- trace$getCodonSpecificParameterTrace(1L)
    nse_tr    <- trace$getCodonSpecificParameterTrace(2L)
    phi_tr    <- trace$getSynthesisRateTrace()              # [[mix]][[gene]]
    sphi_tr   <- trace$getStdDevSynthesisRateTraces()       # [[mix]] -> vector

    n_codons <- length(alpha_tr[[1L]])
    n_nse    <- length(nse_tr[[1L]])
    n_genes  <- length(phi_tr[[1L]])

    n_total <- length(alpha_tr[[1L]][[1L]])
    n_keep  <- ceiling(n_total * keep.fraction)
    idx     <- (n_total - n_keep + 1L):n_total

    if (is.null(codon.names)) codon.names <- seq_len(n_codons)
    if (is.null(gene.names))  gene.names  <- seq_len(n_genes)

    .extract <- function(tr_list) do.call(cbind, lapply(tr_list, `[`, idx))

    alpha_mat  <- .extract(alpha_tr[[1L]])
    lambda_mat <- .extract(lambda_tr[[1L]])
    nse_mat    <- .extract(nse_tr[[1L]])

    nse_names <- if (n_nse == 1L) "NSERate_shared" else paste0("NSERate[", codon.names, "]")

    blocks <- list(alpha_mat, lambda_mat, nse_mat)
    vnames <- c(paste0("alpha[",       codon.names, "]"),
                paste0("lambdaPrime[", codon.names, "]"),
                nse_names)

    has_phi  <- n_genes > 0L && length(phi_tr[[1L]][[1L]]) > 0L
    has_sphi <- length(sphi_tr) > 0L && length(sphi_tr[[1L]]) == n_total

    if (has_phi) {
        log_phi <- log(pmax(.extract(phi_tr[[1L]]), .Machine$double.eps))
        blocks  <- c(blocks, list(log_phi))
        vnames  <- c(vnames, paste0("log_phi[", gene.names, "]"))
    }
    if (has_sphi) {
        blocks <- c(blocks, list(matrix(sphi_tr[[1L]][idx], ncol = 1L)))
        vnames <- c(vnames, "sphi")
    }

    mat <- do.call(cbind, blocks)
    colnames(mat) <- vnames

    arr <- array(mat,
                 dim      = c(n_keep, 1L, ncol(mat)),
                 dimnames = list(NULL, "1", vnames))
    posterior::as_draws_array(arr)
}
