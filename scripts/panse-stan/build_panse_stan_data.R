# ============================================================================
# build_panse_stan_data.R -- assemble the Stan data list for panse_*.stan.
#
# Pure data preparation: no AnaCoDa required.  Reads the RFP CSV directly
# (GeneID, Position, Codon, [Mixture,] RFPCount), filters to the gene set
# defined by init.phi.file (which is canonically the paper's post-filter
# set), and builds a CSR-style flat layout:
#
#   gene_offset[G+1]   1-indexed offsets; gene g spans
#                      gene_offset[g] : gene_offset[g+1]-1
#   codon_at_pos[P]    1-indexed codon ID at each position
#   y[P]               RFP count at each position
#   like_mask[P]       1 = include in likelihood; 0 = sigma-only
#                      (RMF flag: positionMixture+1 < 0 -> mask=0)
#
# Codon ordering: taken from the alpha.csv `Codon` column.  AnaCoDa's
# canonical sense-codon ordering is alphabetical within amino acid,
# amino acids in some fixed order; alpha.csv from the paper's RMF run
# uses exactly that ordering, so matching against it gives index parity
# with the C++ posterior.
#
# Partition function:
#   U = Z / Y
#   Z = sum_g phi[g] * sum_c (alpha[c] / lambdaPrime[c]) * codon_count[g, c]
#   Y = sum(y) over likelihood-contributing positions
# Matches estimateStartingZ() in adapter.dev/lib/local.functions.R, which
# is itself lifted from Nonsense_error_rates/R_scripts/runPANSEMixModel.R.
#
# Prior translation:
#   YAML supplies natural-scale uniform bounds (alpha.prior.lower/upper,
#   lambda.prior.lower/upper) that approximate the RMF Uniform(0, 100)
#   prior.  We translate to a weak normal in log space such that ~99%
#   of the prior mass lands inside the natural-scale uniform range.
#
# Usage as a library:
#   source("scripts/lib/build_panse_stan_data.R")
#   stan_data <- build_panse_stan_data(config)
#
# Usage as a CLI debug tool:
#   Rscript scripts/lib/build_panse_stan_data.R <config.yaml> [--genes N]
#                                               [--out path.rds]
# ============================================================================

suppressPackageStartupMessages({
    library(yaml)
    library(data.table)
})


# --------------------------------------------------------------------------
# CSV readers

.read_codon_csv <- function(path, value_col = "Mean", codon_order = NULL) {
    if (!file.exists(path))
        stop(".read_codon_csv: file not found: ", path)
    tbl <- read.table(path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
    if (!"Codon" %in% colnames(tbl))
        stop(".read_codon_csv: expected `Codon` column in: ", path)
    if (!value_col %in% colnames(tbl))
        stop(".read_codon_csv: expected `", value_col, "` column in: ", path)
    if (is.null(codon_order)) return(setNames(tbl[[value_col]], tbl$Codon))
    out <- setNames(rep(NA_real_, length(codon_order)), codon_order)
    out[tbl$Codon] <- tbl[[value_col]]
    if (any(is.na(out)))
        stop(".read_codon_csv: missing codons in ", path, ": ",
             paste(codon_order[is.na(out)], collapse = ", "))
    out
}

.read_phi_csv <- function(path) {
    if (!file.exists(path))
        stop(".read_phi_csv: file not found: ", path)
    tbl <- read.table(path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
    if (ncol(tbl) < 2L)
        stop(".read_phi_csv: expected GeneID + value column in: ", path)
    setNames(tbl[, 2], tbl[, 1])
}


# --------------------------------------------------------------------------
# Prior translation: natural-scale uniform -> log-scale normal

.uniform_to_lognormal_prior <- function(lower, upper, tail_quantile = 0.995) {
    # Want ~99% (= 1 - 2*(1 - tail_quantile)) of normal mass in
    # (log(lower), log(upper)).  Default tail_quantile = 0.995 gives 99%.
    if (!is.numeric(lower) || lower <= 0)
        stop(".uniform_to_lognormal_prior: lower must be > 0; got ", lower)
    if (!is.numeric(upper) || upper <= lower)
        stop(".uniform_to_lognormal_prior: upper must be > lower; got ",
             upper, " vs ", lower)
    log_lo <- log(lower)
    log_hi <- log(upper)
    z      <- qnorm(tail_quantile)
    list(
        mean = 0.5 * (log_lo + log_hi),
        sd   = 0.5 * (log_hi - log_lo) / z
    )
}


# --------------------------------------------------------------------------
# Main builder

build_panse_stan_data <- function(config,
                                  gene_subset = NULL,
                                  verbose     = TRUE) {
    fit_cfg    <- config$fit
    genome_cfg <- config$genome

    # ------- Resolve input paths --------------------------------------------
    rfp_path <- file.path(genome_cfg$input.dir, genome_cfg$pattern)
    if (!file.exists(rfp_path))
        stop("build_panse_stan_data: RFP CSV not found: ", rfp_path)

    if (length(fit_cfg$init.alpha.files)   < 1L) stop("init.alpha.files missing")
    if (length(fit_cfg$init.lambda.files)  < 1L) stop("init.lambda.files missing")
    if (length(fit_cfg$init.nserate.files) < 1L) stop("init.nserate.files missing")
    if (is.null(fit_cfg$init.phi.file))         stop("init.phi.file missing")

    alpha_path  <- fit_cfg$init.alpha.files[[1]]
    lambda_path <- fit_cfg$init.lambda.files[[1]]
    nse_path    <- fit_cfg$init.nserate.files[[1]]
    phi_path    <- fit_cfg$init.phi.file

    # ------- Load init phi (defines the gene set) ---------------------------
    init_phi  <- .read_phi_csv(phi_path)
    gene_ids  <- names(init_phi)
    if (!is.null(gene_subset)) {
        if (is.numeric(gene_subset) && length(gene_subset) == 1L) {
            gene_ids <- head(gene_ids, gene_subset)
        } else if (is.character(gene_subset)) {
            missing <- setdiff(gene_subset, names(init_phi))
            if (length(missing) > 0L)
                stop("gene_subset has genes not in phi.csv: ",
                     paste(head(missing, 5), collapse = ", "))
            gene_ids <- gene_subset
        } else {
            stop("gene_subset must be integer N or character vector of GeneIDs")
        }
        init_phi <- init_phi[gene_ids]
    }
    G <- length(gene_ids)
    if (verbose) message("Gene set: ", G, " genes from ", basename(phi_path))

    # ------- Load alpha/lambda (defines codon ordering) ---------------------
    alpha_tbl  <- read.table(alpha_path, sep = ",", header = TRUE,
                             stringsAsFactors = FALSE)
    codon_order <- alpha_tbl$Codon
    C <- length(codon_order)
    if (verbose) message("Codon ordering: ", C, " codons from ",
                         basename(alpha_path))

    init_alpha  <- .read_codon_csv(alpha_path,  codon_order = codon_order)
    init_lambda <- .read_codon_csv(lambda_path, codon_order = codon_order)
    init_nse    <- .read_codon_csv(nse_path,    codon_order = codon_order)

    # ------- Load and filter the RFP CSV ------------------------------------
    if (verbose) message("Reading RFP CSV: ", basename(rfp_path))
    rfp <- fread(rfp_path)
    required_cols <- c("GeneID", "Position", "Codon", "RFPCount")
    missing_cols  <- setdiff(required_cols, colnames(rfp))
    if (length(missing_cols) > 0L)
        stop("RFP CSV missing columns: ", paste(missing_cols, collapse = ", "))
    has_mixture <- "Mixture" %in% colnames(rfp)
    if (!has_mixture) {
        if (verbose) message("RFP CSV has no Mixture column; like_mask = 1 for all positions")
        rfp[, Mixture := 1L]
    }

    # Filter to chosen gene set, preserve order from gene_ids
    rfp <- rfp[GeneID %in% gene_ids]
    rfp[, GeneID := factor(GeneID, levels = gene_ids)]
    setorder(rfp, GeneID, Position)
    rfp[, GeneID := as.character(GeneID)]

    # Sanity: every gene_id should be present
    actual_genes <- unique(rfp$GeneID)
    if (length(actual_genes) != G) {
        missing_in_rfp <- setdiff(gene_ids, actual_genes)
        if (length(missing_in_rfp) > 0L)
            stop("RFP CSV missing ", length(missing_in_rfp),
                 " genes from phi.csv; first few: ",
                 paste(head(missing_in_rfp, 5), collapse = ", "))
    }

    # Map codons to 1-indexed IDs
    codon_id <- match(rfp$Codon, codon_order)
    if (any(is.na(codon_id))) {
        bad <- unique(rfp$Codon[is.na(codon_id)])
        stop("RFP CSV has codons not in alpha.csv ordering: ",
             paste(head(bad, 5), collapse = ", "))
    }

    # like_mask: 1 if Mixture+1 > 0 (in RMF terms), 0 if Mixture+1 < 0.
    # Mixture == 0 is forbidden by RMF; we mirror that.
    if (any(rfp$Mixture == 0L))
        stop("RFP CSV has Mixture == 0 rows (forbidden by RMF: must be > 0 ",
             "to include in likelihood, < 0 for sigma-only).")
    like_mask <- as.integer(rfp$Mixture > 0L)

    # CSR offsets
    pos_per_gene <- rfp[, .N, by = GeneID]
    setkey(pos_per_gene, GeneID)
    pos_per_gene <- pos_per_gene[gene_ids]    # reorder to gene_ids order
    gene_offset  <- cumsum(c(1L, pos_per_gene$N))   # length G+1, 1-indexed
    P <- nrow(rfp)
    stopifnot(gene_offset[G + 1L] == P + 1L)

    # ------- Partition function U = Z / Y -----------------------------------
    # Z must match what the Stan model's likelihood implies for the SAME
    # init parameters:
    #
    #   Stan model implies mu[g,p] = alpha[c] * phi[g] * sigma_E[g,p] /
    #                                (U * lambdaPrime[c])
    #   summing E[counts] over all positions and setting = Y_obs:
    #     U = (sum_g phi[g] * sum_p (alpha[c]/lambdaPrime[c]) *
    #          sigma_E[g,p]) / Y_obs
    #
    # where sigma_E[g,p] = exp( sum_{p'<p} log_psuccess[c[p']] ) and
    # log_psuccess is the same 2nd-order Taylor of log E[v/(W+v)] that
    # PANSEModel.cpp::elongationUntilIndexApproximation2ProbabilityLog
    # computes (and that the Stan model's transformed parameters block
    # uses).
    #
    # Earlier versions of this builder used Z = sum phi * (alpha/lambda) *
    # codon_count (i.e. sigma == 1 everywhere), which inflated Z by
    # ~mean(1/sigma) - 1 (~4% at NSE=1e-5 with 500-codon ORFs).  That
    # over-inflated U fed into the Stan model, biasing recovered
    # lambdaPrime DOWN by the reciprocal factor (lambda_post ~ lambda_truth
    # / U_inflated, giving -4% log-bias and ~73% coverage on a nominal
    # 90% CI -- exactly what diagnostic comparison of nb-simulated vs
    # nb-2o-approx fits revealed 2026-05-24).
    v_per_codon  <- 1.0 / init_nse
    a_over_lv    <- init_alpha / (init_lambda * v_per_codon)
    log_psuccess <- -a_over_lv +
                    a_over_lv / (init_lambda * v_per_codon) +
                    0.5 * a_over_lv * a_over_lv
    wait_per_codon <- init_alpha / init_lambda
    # Position-aware sigma via grouped cumsum on the long-format RFP table.
    rfp[, .lp := log_psuccess[match(Codon, codon_order)]]
    rfp[, .log_survive := c(0, head(cumsum(.lp), .N - 1L)), by = GeneID]
    rfp[, .sigma := exp(.log_survive)]
    rfp[, .wait  := wait_per_codon[match(Codon, codon_order)]]
    phi_by_gene  <- setNames(as.numeric(init_phi), gene_ids)
    rfp[, .phi   := phi_by_gene[GeneID]]
    Z <- sum(rfp$.phi * rfp$.wait * rfp$.sigma)
    # Clean up the scratch columns so the rfp data.table stays consistent
    # with the downstream expectations (and doesn't bloat the saved stan-
    # data.rds).
    rfp[, c(".lp", ".log_survive", ".sigma", ".wait", ".phi") := NULL]
    Y <- sum(rfp$RFPCount)
    U <- Z / Y
    if (verbose) message(sprintf(
        "Partition function (sigma-aware): Z = %.6g, Y = %.6g, U = %.6g",
        Z, Y, U))

    # init.partition.function override (auto | numeric).  Mirror
    # .resolve.init.partition.function in lib/local.functions.R: if numeric,
    # use it directly; if "auto" (or NULL), use the computed U.
    init_pf <- fit_cfg$init.partition.function
    if (!is.null(init_pf) && is.numeric(init_pf) && length(init_pf) == 1L) {
        # User-supplied raw Z value; convert to U.
        if (verbose) message("Using YAML-supplied init.partition.function = ",
                             init_pf, " (Z, not U)")
        U <- init_pf / Y
    }

    # ------- Prior translation ----------------------------------------------
    alpha_lo  <- fit_cfg$alpha.prior.lower  %||% 1e-3
    alpha_hi  <- fit_cfg$alpha.prior.upper  %||% 100
    lambda_lo <- fit_cfg$lambda.prior.lower %||% 1e-3
    lambda_hi <- fit_cfg$lambda.prior.upper %||% 100
    alpha_prior  <- .uniform_to_lognormal_prior(alpha_lo,  alpha_hi)
    lambda_prior <- .uniform_to_lognormal_prior(lambda_lo, lambda_hi)

    nse_prior <- fit_cfg$nserate.prior %||% list()
    nse_type  <- nse_prior$type %||% "Log-Uniform"
    nse_lo    <- nse_prior$uniform.lower %||% 1e-7
    nse_hi    <- nse_prior$uniform.upper %||% 1e-3
    nse_log_uniform <- switch(
        nse_type,
        "Log-Uniform"     = 1L,
        "Natural-Uniform" = 0L,
        stop("nserate.prior.type must be 'Log-Uniform' or 'Natural-Uniform'; ",
             "got: ", nse_type)
    )

    # ------- Stan data list -------------------------------------------------
    # all_unmasked enables the vectorized NB2 fast path in panse_*.stan when
    # no sigma-only positions are present (common case: Weinberg, Wu, Mohammad
    # all have like_mask == 1 everywhere).
    all_unmasked <- as.integer(all(like_mask == 1L))

    # phi / sphi handoff:
    #   fit.with.phi: false (v0 csp-only)         -> phi is DATA (provide it)
    #   fit.with.phi: true  (v1+ basic / sphi-est) -> phi is a PARAMETER;
    #                                                provide sphi (fixed)
    #                                                for the log_phi prior
    with_phi_sampled <- isTRUE(fit_cfg$with.phi)
    estimates_sphi   <- grepl("sphi-est", fit_cfg$model %||% "")
    sphi_fixed       <- fit_cfg$sphi %||% sd(log(init_phi))
    sphi_prior_sd    <- as.numeric(fit_cfg$sphi.prior.sd %||% 2.5)

    stan_data <- list(
        G            = G,
        C            = C,
        P            = P,
        gene_offset  = as.integer(gene_offset),
        codon_at_pos = as.integer(codon_id),
        y            = as.integer(rfp$RFPCount),
        like_mask    = like_mask,
        all_unmasked = all_unmasked,
        U            = as.numeric(U),
        log_alpha_prior_mean  = alpha_prior$mean,
        log_alpha_prior_sd    = alpha_prior$sd,
        log_lambda_prior_mean = lambda_prior$mean,
        log_lambda_prior_sd   = lambda_prior$sd,
        # Hard bounds on log_alpha / log_lambdaPrime (Stan parameter
        # declarations).  Matches the natural-scale uniform YAML config:
        # alpha.prior.lower/upper -> log_alpha_lower/upper.
        log_alpha_lower       = log(alpha_lo),
        log_alpha_upper       = log(alpha_hi),
        log_lambda_lower      = log(lambda_lo),
        log_lambda_upper      = log(lambda_hi),
        log_nse_lower         = log(nse_lo),
        log_nse_upper         = log(nse_hi),
        nse_log_uniform       = nse_log_uniform,
        emit_log_lik          = as.integer(fit_cfg$emit.log.lik %||% 1L),
        grainsize             = as.integer(fit_cfg$grainsize %||% 1L)
    )
    if (with_phi_sampled) {
        if (estimates_sphi) {
            stan_data$sphi_prior_sd <- sphi_prior_sd
            if (verbose) message(sprintf(
                "with.phi = TRUE, sphi-est model: sphi is a parameter; half-normal(0, %.3f) prior",
                sphi_prior_sd))
        } else {
            stan_data$sphi <- as.numeric(sphi_fixed)
            if (verbose) message(sprintf(
                "with.phi = TRUE: phi is a parameter; sphi (fixed) = %.4f",
                sphi_fixed))
        }
    } else {
        stan_data$phi <- as.numeric(init_phi)
    }

    # Metadata for downstream tooling (not passed to Stan)
    attr(stan_data, "gene_ids")    <- gene_ids
    attr(stan_data, "codon_order") <- codon_order
    attr(stan_data, "init_alpha")  <- init_alpha
    attr(stan_data, "init_lambda") <- init_lambda
    attr(stan_data, "init_nse")    <- init_nse
    attr(stan_data, "init_phi")    <- init_phi
    attr(stan_data, "Z")           <- Z
    attr(stan_data, "Y")           <- Y

    invisible(stan_data)
}


# Minimal `%||%` shim (in case purrr / rlang not loaded)
`%||%` <- function(x, y) if (is.null(x)) y else x


# --------------------------------------------------------------------------
# CLI entry point

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    parse_args <- function(args) {
        out <- list(config = NULL, genes = NULL, out = NULL)
        i <- 1L
        while (i <= length(args)) {
            a <- args[[i]]
            if (a == "--genes")    { out$genes  <- as.integer(args[[i + 1L]]); i <- i + 2L }
            else if (a == "--out") { out$out    <- args[[i + 1L]];             i <- i + 2L }
            else if (a == "--help" || a == "-h") {
                cat("Usage: build_panse_stan_data.R <config.yaml> [--genes N] [--out path.rds]\n")
                quit(save = "no", status = 0)
            }
            else if (is.null(out$config))        { out$config <- a; i <- i + 1L }
            else stop("Unknown argument: ", a)
            }
        if (is.null(out$config)) stop("Usage: build_panse_stan_data.R <config.yaml> [...]")
        out
    }
    opts <- parse_args(args)

    cfg_path <- opts$config
    if (!file.exists(cfg_path)) stop("Config not found: ", cfg_path)
    cfg_abs  <- normalizePath(cfg_path)
    # Native runner convention (see adapter.dev/fit.R): setwd to the dir
    # containing runs/, so YAML paths like "../data/..." resolve to the
    # species-root data tree.  Mirror that here.
    setwd(dirname(dirname(cfg_abs)))
    config <- yaml.load_file(cfg_abs)

    stan_data <- build_panse_stan_data(config, gene_subset = opts$genes)

    cat("\n--- Stan data summary ---\n")
    cat("G              =", stan_data$G, "\n")
    cat("C              =", stan_data$C, "\n")
    cat("P              =", stan_data$P, "\n")
    cat("y total        =", sum(stan_data$y), "\n")
    cat("y > 0 (count)  =", sum(stan_data$y > 0), "\n")
    cat("like_mask sum  =", sum(stan_data$like_mask), "\n")
    cat("phi range      =", range(stan_data$phi), "\n")
    cat("U              =", stan_data$U, "\n")
    cat("log_alpha_prior= N(", round(stan_data$log_alpha_prior_mean, 3),
        ",", round(stan_data$log_alpha_prior_sd, 3), ")\n")
    cat("log_nse bounds = [", round(stan_data$log_nse_lower, 3),
        ",", round(stan_data$log_nse_upper, 3), "]  (log_uniform =",
        stan_data$nse_log_uniform, ")\n")

    if (!is.null(opts$out)) {
        saveRDS(stan_data, opts$out)
        cat("\nSaved Stan data list to:", opts$out, "\n")
    }
}
