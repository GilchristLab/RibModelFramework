# buildPanseStanData.R -- assemble the Stan data list for panse_*.stan.
#
# Pure data preparation: no AnaCoDa C++ required.  Reads the RFP CSV directly
# (GeneID, Position, Codon, [Mixture,] RFPCount), filters to the gene set
# defined by init.phi.file, and builds a CSR-style flat layout:
#
#   gene_offset[G+1]   1-indexed offsets; gene g spans
#                      gene_offset[g] : gene_offset[g+1]-1
#   codon_at_pos[P]    1-indexed codon ID at each position
#   y[P]               RFP count at each position
#   like_mask[P]       1 = include in likelihood; 0 = sigma-only
#
# Partition function (U = Z / Y) is sigma-aware to avoid the ~4% inflation
# that results from ignoring elongation survival (see comment in function body).

utils::globalVariables(c(
    ".N", ":=",
    "GeneID", "Position", "Codon", "RFPCount", "Mixture",
    ".lp", ".log_survive", ".sigma", ".wait", ".phi"
))


# --------------------------------------------------------------------------
# CSV readers (internal)

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
# Phi prior helpers (internal)

# Parse a phi prior column spec "data(path, col=name)" and return a G-length
# numeric vector aligned to gene_ids.
.load_phi_prior_col <- function(spec, gene_ids) {
    m <- regexec("^data\\((.+),\\s*col=([^)]+)\\)$", trimws(spec))[[1]]
    if (m[1] == -1L)
        stop(".load_phi_prior_col: spec must be 'population' or ",
             "'data(file, col=name)'; got: ", spec)
    sub_m <- regmatches(spec, regexec("^data\\((.+),\\s*col=([^)]+)\\)$",
                                      trimws(spec)))[[1]]
    path <- trimws(sub_m[2])
    col  <- trimws(sub_m[3])
    if (!file.exists(path))
        stop(".load_phi_prior_col: file not found: ", path)
    tbl <- read.csv(path, stringsAsFactors = FALSE)
    if (!"gene_id" %in% names(tbl))
        stop(".load_phi_prior_col: file must have a 'gene_id' column: ", path)
    if (!col %in% names(tbl))
        stop(".load_phi_prior_col: column '", col, "' not found in ", path,
             "; available: ", paste(names(tbl), collapse = ", "))
    val_map <- setNames(tbl[[col]], tbl$gene_id)
    missing <- setdiff(gene_ids, names(val_map))
    if (length(missing) > 0L)
        stop(".load_phi_prior_col: file ", path, " missing genes: ",
             paste(head(missing, 5L), collapse = ", "))
    unname(val_map[gene_ids])
}

# Parse the config$phi block and return a list:
#   phi_use_data    0L (population) or 1L (any data-sourced component)
#   phi_prior_mu    G-length numeric vector (placeholder rep(0,G) when population)
#   phi_prior_sigma G-length numeric vector (placeholder rep(1,G) when population)
.parse_phi_prior_spec <- function(phi_cfg, gene_ids) {
    G <- length(gene_ids)
    # Default: no phi block -> population/population
    if (is.null(phi_cfg)) {
        return(list(phi_use_data    = 0L,
                    phi_prior_mu    = rep(0.0, G),
                    phi_prior_sigma = rep(1.0, G)))
    }
    gene_cfg   <- phi_cfg$gene %||% list()
    mu_spec    <- gene_cfg$prior.mu    %||% "population"
    sigma_spec <- gene_cfg$prior.sigma %||% "population"

    use_data_mu    <- !identical(mu_spec,    "population")
    use_data_sigma <- !identical(sigma_spec, "population")
    phi_use_data   <- as.integer(use_data_mu || use_data_sigma)

    phi_prior_mu <- if (use_data_mu) {
        .load_phi_prior_col(mu_spec, gene_ids)
    } else {
        rep(0.0, G)   # placeholder; Stan ignores when phi_use_data=0
    }
    phi_prior_sigma <- if (use_data_sigma) {
        .load_phi_prior_col(sigma_spec, gene_ids)
    } else {
        rep(1.0, G)   # placeholder; Stan ignores when phi_use_data=0
    }
    if (any(phi_prior_sigma <= 0, na.rm = TRUE))
        stop(".parse_phi_prior_spec: phi prior sigma must be > 0 for all genes")

    list(phi_use_data    = phi_use_data,
         phi_prior_mu    = phi_prior_mu,
         phi_prior_sigma = phi_prior_sigma)
}


# --------------------------------------------------------------------------
# Prior translation: natural-scale uniform -> log-scale normal (internal)

.uniform_to_lognormal_prior <- function(lower, upper, tail_quantile = 0.995) {
    # ~99% (= 1 - 2*(1 - tail_quantile)) of normal mass inside (log(lower), log(upper)).
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

#' Assemble the Stan data list for panse_*.stan
#'
#' Reads RFP CSV + init parameter CSVs and assembles the flat CSR-layout Stan
#' data list consumed by all \code{panse_*.stan} models.  Returns the list
#' invisibly; parameter init values and partition-function intermediates are
#' attached as attributes for downstream tooling.
#'
#' @param config Parsed YAML config list (as returned by \code{yaml::read_yaml}).
#' @param gene_subset Integer N (first N genes) or character vector of GeneIDs.
#'   \code{NULL} (default) uses all genes defined by \code{init.phi.file}.
#' @param verbose Logical; if \code{TRUE} (default) print progress messages.
#' @return Invisibly, a named list suitable for \code{cmdstanr::CmdStanModel$sample()}.
#'   Attributes: \code{gene_ids}, \code{codon_order}, \code{init_alpha},
#'   \code{init_lambda}, \code{init_nse}, \code{init_phi}, \code{Z}, \code{Y}.
#' @export
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

    # ------- Partition function U = Z / Y (sigma-aware) ---------------------
    # Earlier versions used Z = sum phi*(alpha/lambda)*codon_count (sigma==1
    # everywhere), inflating Z by ~mean(1/sigma)-1 (~4% at NSE=1e-5 with
    # 500-codon ORFs) and biasing recovered lambdaPrime DOWN.  Here we compute
    # sigma_E via the same 2nd-order Taylor that the Stan model uses.
    v_per_codon  <- 1.0 / init_nse
    a_over_lv    <- init_alpha / (init_lambda * v_per_codon)
    log_psuccess <- -a_over_lv +
                    a_over_lv / (init_lambda * v_per_codon) +
                    0.5 * a_over_lv * a_over_lv
    wait_per_codon <- init_alpha / init_lambda
    rfp[, .lp := log_psuccess[match(Codon, codon_order)]]
    rfp[, .log_survive := c(0, head(cumsum(.lp), .N - 1L)), by = GeneID]
    rfp[, .sigma := exp(.log_survive)]
    rfp[, .wait  := wait_per_codon[match(Codon, codon_order)]]
    phi_by_gene  <- setNames(as.numeric(init_phi), gene_ids)
    rfp[, .phi   := phi_by_gene[GeneID]]
    Z <- sum(rfp$.phi * rfp$.wait * rfp$.sigma)
    rfp[, c(".lp", ".log_survive", ".sigma", ".wait", ".phi") := NULL]
    Y <- sum(rfp$RFPCount)
    U <- Z / Y
    if (verbose) message(sprintf(
        "Partition function (sigma-aware): Z = %.6g, Y = %.6g, U = %.6g",
        Z, Y, U))

    # init.partition.function override (auto | numeric)
    init_pf <- fit_cfg$init.partition.function
    if (!is.null(init_pf) && is.numeric(init_pf) && length(init_pf) == 1L) {
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
    # Gamma prior params (mean = shape/rate; default Gamma(2, 2e5) -> mean 1e-5).
    # Ignored by the Stan model unless nse_gamma == 1, but always emitted.
    nse_gamma       <- 0L
    nse_gamma_shape <- as.numeric(nse_prior$gamma.shape %||% 2.0)
    nse_gamma_rate  <- as.numeric(nse_prior$gamma.rate  %||% 2e5)
    nse_log_uniform <- switch(
        nse_type,
        "Log-Uniform"     = 1L,
        "Natural-Uniform" = 0L,
        "Gamma"           = 1L,   # placeholder; nse_gamma below activates Gamma
        stop("nserate.prior.type must be 'Log-Uniform', 'Natural-Uniform', or ",
             "'Gamma'; got: ", nse_type)
    )
    if (identical(nse_type, "Gamma")) nse_gamma <- 1L

    # ------- Stan data list -------------------------------------------------
    all_unmasked <- as.integer(all(like_mask == 1L))

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
        log_alpha_lower       = log(alpha_lo),
        log_alpha_upper       = log(alpha_hi),
        log_lambda_lower      = log(lambda_lo),
        log_lambda_upper      = log(lambda_hi),
        log_nse_lower         = log(nse_lo),
        log_nse_upper         = log(nse_hi),
        nse_log_uniform       = nse_log_uniform,
        nse_gamma             = nse_gamma,
        nse_gamma_shape       = nse_gamma_shape,
        nse_gamma_rate        = nse_gamma_rate,
        emit_log_lik          = as.integer(fit_cfg$emit.log.lik %||% 1L),
        grainsize             = as.integer(fit_cfg$grainsize %||% 1L)
    )
    if (with_phi_sampled) {
        if (estimates_sphi) {
            stan_data$sphi_prior_sd <- sphi_prior_sd
            if (verbose) message(sprintf(
                "with.phi = TRUE, sphi-est model: sphi is a parameter; half-normal(0, %.3f) prior",
                sphi_prior_sd))
            # Generalized phi prior spec (phi_use_data / phi_prior_mu / phi_prior_sigma)
            phi_spec <- .parse_phi_prior_spec(config$phi, gene_ids)
            stan_data$phi_use_data    <- phi_spec$phi_use_data
            stan_data$phi_prior_mu    <- phi_spec$phi_prior_mu
            stan_data$phi_prior_sigma <- phi_spec$phi_prior_sigma
            if (verbose && phi_spec$phi_use_data == 1L)
                message("phi prior: gene-specific mu/sigma from data")
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
