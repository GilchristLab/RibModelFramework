# fonseStan.R -- Stan/HMC backend for FONSE (fonse_sphi_est.stan).
#
# FONSE is an exact per-position categorical-logit (multinomial logistic
# regression with a latent per-gene phi):
#
#     eta_i = -dM_i - dEta_i*phi - dOmega_i*phi*beta(pos),  beta = a1 + a2*pos,
#     P(codon i | pos) = softmax(eta)_i,   reference codon eta = 0.
#
# Because beta(pos) sits inside the exponent, the per-AA normaliser depends on
# the codon's position and the likelihood CANNOT be collapsed to per-gene codon
# counts (as ROC/PANSE are).  The Stan data is therefore PER-POSITION: one row
# per sense codon occurrence.  This module turns a genome (or a long-format obs
# frame) into that layout, then compiles + samples the model with cmdstanr.
#
# Parameter indexing matches roc_sphi_est.stan: dM/dEta/dOmega are length-K
# vectors over the NON-reference codons; aa_start[a]..aa_end[a] is the param
# block for amino acid a, and the reference codon is implicit (eta = 0).
#
# This mirrors the PANSE dual-backend pattern (buildPanseStanData.R /
# fitPanseStan.R / panseStan.R).  The native MCMC backend is FONSEModel/
# FONSEParameter; this is the Stan backend for the same model.

# ---------------------------------------------------------------------------
# Stan model resolver (mirrors .panse_stan_file in panseStan.R)
# ---------------------------------------------------------------------------

.FONSE_STAN_MODELS <- c("fonse_sphi_est.stan")

.fonse_stan_file <- function(basename) {
    if (!basename %in% .FONSE_STAN_MODELS)
        stop(".fonse_stan_file: unknown model '", basename, "'; valid: ",
             paste(.FONSE_STAN_MODELS, collapse = ", "))

    # 1. Installed package path (post R CMD INSTALL).
    pkg_path <- system.file("stan", basename, package = "AnaCoDa")
    if (nzchar(pkg_path) && file.exists(pkg_path)) return(pkg_path)

    # 2. In-source fallback (worktree / devtools::load_all()).
    src_path <- system.file("../stan", basename, package = "AnaCoDa")
    if (nzchar(src_path) && file.exists(src_path)) return(src_path)

    stop(".fonse_stan_file: could not resolve '", basename, "' from the ",
         "installed package or in-source stan/. Run R CMD INSTALL or ",
         "devtools::load_all() from the package root.")
}

#' Return the path to the installed FONSE Stan model file
#'
#' @param basename Character; currently only \code{"fonse_sphi_est.stan"}.
#' @return Absolute path string.
#' @export
fonse_stan_model_path <- function(basename = "fonse_sphi_est.stan") {
    .fonse_stan_file(basename)
}

#' Compile the FONSE Stan model with cmdstanr
#'
#' @param basename FONSE Stan model basename (default \code{"fonse_sphi_est.stan"}).
#' @param build_dir Writable cache directory for the compiled binary.
#' @param threads Logical; compile with \code{stan_threads = TRUE} (default).
#' @param ... Forwarded to \code{cmdstanr::cmdstan_model()}.
#' @return A \code{CmdStanModel} object.
#' @export
fonse_compile_model <- function(basename  = "fonse_sphi_est.stan",
                                build_dir = tools::R_user_dir("AnaCoDa", "cache"),
                                threads   = TRUE,
                                ...) {
    if (!requireNamespace("cmdstanr", quietly = TRUE))
        stop("fonse_compile_model requires cmdstanr. ",
             "Install with: install.packages('cmdstanr', repos='https://mc-stan.org/r-packages/')")
    stan_file <- .fonse_stan_file(basename)
    if (!dir.exists(build_dir)) dir.create(build_dir, recursive = TRUE)
    cpp_opts <- if (threads) list(stan_threads = TRUE) else list()
    cmdstanr::cmdstan_model(stan_file, dir = build_dir, cpp_options = cpp_opts, ...)
}

# ---------------------------------------------------------------------------
# Codon bookkeeping
# ---------------------------------------------------------------------------

# Build the non-reference codon -> param-index map and aa_start/aa_end ranges.
# Returns:
#   aas        multi-codon amino acids (M/W/X excluded), AnaCoDa order
#   syn        list aa -> synonyms (AAToCodon focal=FALSE; LAST is reference)
#   ref_codon  named logical over all sense codons (TRUE = reference)
#   aa_index   named integer aa -> 1..A
#   aa_start   integer[A] first param index per AA (non-ref codons)
#   aa_end     integer[A] last  param index per AA
#   codon2param named integer: non-ref codon -> 1..K param index
#   codon2choice named integer: codon -> within-AA Stan choice index
#                (1 = reference, 2..n+1 = non-reference in synonym order)
#   nonref_codons character[K] param-index-ordered non-reference codons
.fonse_codon_map <- function() {
    aas <- setdiff(AnaCoDa::aminoAcids(), c("M", "W", "X"))
    A   <- length(aas)
    syn <- vector("list", A); names(syn) <- aas
    aa_start <- integer(A); aa_end <- integer(A)
    codon2param  <- integer(0)
    codon2choice <- integer(0)
    nonref_codons <- character(0)
    ref_codon <- logical(0)
    p <- 0L
    for (i in seq_along(aas)) {
        aa <- aas[i]
        cs <- AnaCoDa::AAToCodon(aa, focal = FALSE)   # last = reference
        syn[[aa]] <- cs
        k  <- length(cs)
        nref <- k - 1L
        aa_start[i] <- p + 1L
        aa_end[i]   <- p + nref
        # non-reference codons cs[1:nref] -> param indices, choice 2..k
        for (j in seq_len(nref)) {
            p <- p + 1L
            codon2param[cs[j]]  <- p
            codon2choice[cs[j]] <- j + 1L
            nonref_codons[p]    <- cs[j]
            ref_codon[cs[j]]    <- FALSE
        }
        # reference codon cs[k] -> choice 1
        codon2choice[cs[k]] <- 1L
        ref_codon[cs[k]]    <- TRUE
    }
    aa_index <- setNames(seq_along(aas), aas)
    list(aas = aas, A = A, K = p, syn = syn, ref_codon = ref_codon,
         aa_index = aa_index, aa_start = aa_start, aa_end = aa_end,
         codon2param = codon2param, codon2choice = codon2choice,
         nonref_codons = nonref_codons)
}

# Map a codon -> its amino acid (multi-codon AAs only; NA otherwise).
.fonse_codon2aa <- function(map) {
    out <- character(0)
    for (aa in map$aas) for (cod in map$syn[[aa]]) out[cod] <- aa
    out
}

# ---------------------------------------------------------------------------
# Genome -> long-format obs (gene, pos, codon)
# ---------------------------------------------------------------------------

#' Convert a FASTA genome to FONSE long-format codon observations
#'
#' One row per sense codon occurrence.  Position is the codon's ordinal index
#' along the CDS with the FIRST codon after the start (ATG) at \code{pos = 1}
#' (matching native FONSE, where \code{beta = a1 + a2*pos} and the start codon
#' is not modelled).  Stop codons and single-codon amino acids (Met/Trp) are
#' dropped.
#'
#' @param fasta Path to a FASTA file of coding sequences.
#' @param map Optional precomputed \code{.fonse_codon_map()} (for reuse).
#' @return \code{data.frame(gene, pos, codon)} with \code{gene} a 1..G integer.
#' @export
fonse_genome_to_obs <- function(fasta, map = NULL) {
    if (!file.exists(fasta)) stop("fonse_genome_to_obs: file not found: ", fasta)
    if (is.null(map)) map <- .fonse_codon_map()
    codon2aa <- .fonse_codon2aa(map)

    lines <- readLines(fasta, warn = FALSE)
    hdr   <- grepl("^>", lines)
    ids   <- sub("^>", "", lines[hdr])
    grp   <- cumsum(hdr)
    seqs  <- tapply(lines[!hdr], grp[!hdr], function(x) toupper(paste0(x, collapse = "")))
    seqs  <- seqs[order(as.integer(names(seqs)))]

    gene_v <- integer(0); pos_v <- integer(0); codon_v <- character(0)
    gi <- 0L
    for (s in seqs) {
        gi <- gi + 1L
        ncodon <- nchar(s) %/% 3L
        if (ncodon < 2L) next
        # codon ordinal index along CDS: 1 = ATG (start, skipped).
        # 2nd codon (CDS index 2) -> pos 1, etc.
        starts <- seq(4L, by = 3L, length.out = ncodon - 1L)  # skip codon 1 (ATG)
        cods   <- substring(s, starts, starts + 2L)
        keep   <- cods %in% names(codon2aa)
        if (!any(keep)) next
        np <- seq_len(ncodon - 1L)            # pos = 1,2,3,... for codons 2,3,4,...
        gene_v  <- c(gene_v,  rep.int(gi, sum(keep)))
        pos_v   <- c(pos_v,   np[keep])
        codon_v <- c(codon_v, cods[keep])
    }
    data.frame(gene = gene_v, pos = pos_v, codon = codon_v,
               stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# obs -> Stan data list
# ---------------------------------------------------------------------------

#' Assemble the Stan data list for fonse_sphi_est.stan
#'
#' @param obs Long-format \code{data.frame(gene, pos, codon)} (e.g. from
#'   \code{\link{fonse_genome_to_obs}}).  \code{gene} may be any vector; it is
#'   re-coded to a 1..G integer internally.
#' @param a1,a2 Initiation / elongation cost constants (default 4 / 4, the
#'   native fixed-by-default values).
#' @param phi.mphi,phi.sphi phi prior spec objects (see \code{phi.prior.spec};
#'   defaults: constrained mean 1, estimated sphi with a half-normal-ish flat).
#' @param dM.prior.mean,dM.prior.sd,dEta.prior.mean,dEta.prior.sd Scalars
#'   recycled to length K.  Defaults N(0,1).
#' @param dOmega.prior.mean,dOmega.prior.sd Scalars recycled to length K.
#'   dOmega is a small-scale parameter (odds ratio b/c ~ 5e-4) so the default
#'   prior sd is deliberately tight (0.01).
#' @param grainsize reduce_sum grainsize (default 1).
#' @return A named list ready for \code{cmdstanr}'s \code{$sample(data = ...)},
#'   plus attributes \code{codon_map} and \code{gene_levels} for downstream use.
#' @export
build_fonse_stan_data <- function(obs,
                                   a1 = 4, a2 = 4,
                                   phi.mphi = NULL, phi.sphi = NULL,
                                   dM.prior.mean = 0,    dM.prior.sd = 1,
                                   dEta.prior.mean = 0,  dEta.prior.sd = 1,
                                   dOmega.prior.mean = 0, dOmega.prior.sd = 0.01,
                                   grainsize = 1L,
                                   map = NULL) {
    stopifnot(all(c("gene", "pos", "codon") %in% colnames(obs)))
    if (is.null(map)) map <- .fonse_codon_map()
    codon2aa <- .fonse_codon2aa(map)

    # Drop non-sense / single-codon-AA codons.
    obs <- obs[obs$codon %in% names(codon2aa), , drop = FALSE]
    if (nrow(obs) == 0) stop("build_fonse_stan_data: no sense codons in obs")

    # Re-code genes to 1..G and sort rows by gene so each gene is contiguous.
    gene_levels <- sort(unique(obs$gene))
    obs$g <- match(obs$gene, gene_levels)
    obs <- obs[order(obs$g), , drop = FALSE]
    G   <- length(gene_levels)
    N   <- nrow(obs)

    obs_aa     <- map$aa_index[codon2aa[obs$codon]]
    obs_choice <- map$codon2choice[obs$codon]
    obs_pos    <- as.double(obs$pos)

    # gene_start / gene_end: contiguous row ranges. A gene with zero kept rows
    # gets start > end (gene_end = gene_start - 1), an empty reduce_sum range.
    gene_start <- integer(G); gene_end <- integer(G)
    rle_g <- rle(obs$g)
    ends  <- cumsum(rle_g$lengths)
    starts <- ends - rle_g$lengths + 1L
    present <- rle_g$values
    gene_start[present] <- starts
    gene_end[present]   <- ends
    # genes with no observations: empty range start=1,end=0-style
    missing <- setdiff(seq_len(G), present)
    if (length(missing)) { gene_start[missing] <- 1L; gene_end[missing] <- 0L }

    K <- map$K
    rec <- function(x, nm) {
        if (length(x) == 1L) return(rep(as.double(x), K))
        if (length(x) == K)  return(as.double(x))
        stop("build_fonse_stan_data: ", nm, " must be length 1 or K=", K)
    }

    # phi spec defaults: mean phi = 1 (constrained), sphi estimated (flat > 0).
    # constrained()/estimated() are the exported phi_spec constructors.
    if (is.null(phi.mphi)) phi.mphi <- constrained(statistic = "mean", value = 1)
    if (is.null(phi.sphi)) phi.sphi <- estimated()
    phi_fields <- .phiSpecToStanData(phi.mphi, phi.sphi)

    dat <- c(list(
        G = G, A = map$A, K = K, N = N,
        aa_start = as.array(map$aa_start), aa_end = as.array(map$aa_end),
        gene_start = as.array(gene_start), gene_end = as.array(gene_end),
        obs_aa = as.array(as.integer(obs_aa)),
        obs_pos = as.array(obs_pos),
        obs_choice = as.array(as.integer(obs_choice)),
        a1 = as.double(a1), a2 = as.double(a2),
        dM_prior_mean = rec(dM.prior.mean, "dM.prior.mean"),
        dM_prior_sd   = rec(dM.prior.sd,   "dM.prior.sd"),
        dEta_prior_mean = rec(dEta.prior.mean, "dEta.prior.mean"),
        dEta_prior_sd   = rec(dEta.prior.sd,   "dEta.prior.sd"),
        dOmega_prior_mean = rec(dOmega.prior.mean, "dOmega.prior.mean"),
        dOmega_prior_sd   = rec(dOmega.prior.sd,   "dOmega.prior.sd"),
        noncentered = 0L,
        grainsize = as.integer(grainsize)
    ), phi_fields)

    attr(dat, "codon_map")   <- map
    attr(dat, "gene_levels") <- gene_levels
    dat
}

# ---------------------------------------------------------------------------
# Fit driver
# ---------------------------------------------------------------------------

#' Fit FONSE via Stan/HMC (cmdstanr)
#'
#' Compiles \code{fonse_sphi_est.stan} and samples.  The data list must come
#' from \code{\link{build_fonse_stan_data}}.
#'
#' @param stan.data Data list from \code{build_fonse_stan_data}.
#' @param chains,parallel_chains,threads_per_chain,iter_warmup,iter_sampling,
#'   adapt_delta,max_treedepth,refresh,seed Standard cmdstanr sampling controls.
#' @param init Optional cmdstanr init (function/list/value).  If \code{NULL}, a
#'   mild jitter around 0 / log_phi=0 is used.
#' @param build_dir Writable compile cache dir.
#' @return The \code{CmdStanMCMC} fit object.
#' @export
fit_fonse_stan <- function(stan.data,
                           chains = 4, parallel_chains = chains,
                           threads_per_chain = 1,
                           iter_warmup = 1000, iter_sampling = 1000,
                           adapt_delta = 0.9, max_treedepth = 10,
                           refresh = 200, seed = 1234,
                           init = NULL,
                           build_dir = tools::R_user_dir("AnaCoDa", "cache")) {
    if (!requireNamespace("cmdstanr", quietly = TRUE))
        stop("fit_fonse_stan requires cmdstanr.")
    mod <- fonse_compile_model(build_dir = build_dir,
                               threads = threads_per_chain > 1)

    G <- stan.data$G; K <- stan.data$K
    if (is.null(init)) {
        init <- function(chain_id = 1) list(
            dM     = stats::rnorm(K, 0, 0.1),
            dEta   = stats::rnorm(K, 0, 0.1),
            dOmega = stats::rnorm(K, 0, 1e-3),
            latent_phi = stats::rnorm(G, 0, 0.1),
            sphi = 1.0,
            mphi_param = 0.0
        )
    }

    mod$sample(
        data = stan.data,
        chains = chains, parallel_chains = parallel_chains,
        threads_per_chain = threads_per_chain,
        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
        adapt_delta = adapt_delta, max_treedepth = max_treedepth,
        refresh = refresh, seed = seed, init = init
    )
}
