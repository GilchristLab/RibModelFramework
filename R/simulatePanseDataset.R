# simulatePanseDataset.R -- forward simulator for PANSE under controlled truth.
#
# Supports three generative processes:
#   "nb-simulated"  -- wait_times ~ Gamma; sigma from realized waits (matches
#                      AnaCoDa::simulateGenome exactly).
#   "nb-2o-approx"  -- deterministic E[sigma] via 2nd-order Taylor (matches
#                      the Stan likelihood exactly; fitting Stan to this data
#                      gives near-perfect recovery).
#   "nb-1o-approx"  -- 1st-order Taylor only; useful for sanity checks.
#
# See scripts/panse-stan/simulate_panse_dataset.R header for the math.


# --------------------------------------------------------------------------
# Default 22-codon set: 4 x 2-codon AAs, 2 x 4-codon, 1 x 6-codon (internal)

.default_codon_table <- function() {
    data.frame(
        AA    = c(rep("F", 2), rep("Y", 2), rep("H", 2), rep("N", 2),
                  rep("V", 4), rep("A", 4),
                  rep("L", 6)),
        Codon = c("TTT", "TTC", "TAT", "TAC", "CAT", "CAC", "AAT", "AAC",
                  "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
                  "TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
        stringsAsFactors = FALSE
    )
}


# --------------------------------------------------------------------------
# Truth-value loaders (internal)

.load_truth_codon_values <- function(path, value_col, codon_order) {
    if (!file.exists(path))
        stop(".load_truth_codon_values: file not found: ", path)
    tbl <- read.csv(path, stringsAsFactors = FALSE)
    if (!"Codon" %in% colnames(tbl))
        stop(".load_truth_codon_values: missing `Codon` column in: ", path)
    if (!value_col %in% colnames(tbl))
        stop(".load_truth_codon_values: missing `", value_col, "` column in: ", path)
    # Subset to our codon_order first to handle full 60/61-codon RMF CSVs.
    tbl <- tbl[tbl$Codon %in% codon_order, , drop = FALSE]
    missing_codons <- setdiff(codon_order, tbl$Codon)
    if (length(missing_codons) > 0L)
        stop(".load_truth_codon_values: missing codons in ", path, ": ",
             paste(missing_codons, collapse = ", "))
    out <- setNames(tbl[[value_col]], tbl$Codon)
    out[codon_order]
}


# --------------------------------------------------------------------------

#' Simulate a PANSE dataset from controlled truth values
#'
#' Generates a synthetic ribosome footprinting dataset from specified truth
#' parameters.  Supports three generative processes; see the package vignette
#' for the math.
#'
#' @param config Parsed YAML config list with a \code{sim:} block.
#' @param seed Integer seed override (overrides \code{config$sim$seed}).
#' @param verbose Logical; if \code{TRUE} (default) print progress messages.
#' @return A named list: \code{rfp} (data.table), \code{codon_table},
#'   \code{truth_alpha}, \code{truth_lambda}, \code{truth_nse},
#'   \code{nse_shared}, \code{phi_truth}, \code{Z}, \code{Y_target},
#'   \code{Y_obs}, \code{U_truth}, \code{seed}, \code{n_genes},
#'   \code{gene_length}, \code{process}.
#' @export
simulate_panse_dataset <- function(config, seed = NULL, verbose = TRUE) {
    sim_cfg <- config$sim %||% stop("config missing top-level `sim:` block")

    # ---- Resolve config -----------------------------------------------------
    n_genes        <- as.integer(sim_cfg$n.genes        %||% 1000L)
    gene_length    <- as.integer(sim_cfg$gene.length    %||% 500L)
    nse_shared     <- as.numeric(sim_cfg$nse.shared     %||% 1.0e-5)
    sigma_phi      <- as.numeric(sim_cfg$phi.sdlog      %||% 1.5)
    phi_scale      <- as.numeric(sim_cfg$phi.scale      %||% 1.0)
    y_target_mode  <- sim_cfg$y.target.mode %||% "U-one"
    y_target_fixed <- as.numeric(sim_cfg$y.target.value  %||% NA_real_)
    process        <- sim_cfg$process       %||% "nb-simulated"
    base_seed      <- as.integer(seed %||% sim_cfg$seed  %||% 20260524L)

    if (!(process %in% c("nb-simulated", "nb-2o-approx", "nb-1o-approx")))
        stop("sim.process must be one of 'nb-simulated', 'nb-2o-approx', ",
             "'nb-1o-approx'; got: ", process)

    truth_alpha_path  <- sim_cfg$truth.alpha.file  %||% stop("sim.truth.alpha.file missing")
    truth_lambda_path <- sim_cfg$truth.lambda.file %||% stop("sim.truth.lambda.file missing")

    codon_table <- if (!is.null(sim_cfg$codon.table)) {
        # sim_cfg$codon.table may arrive as a list of per-codon maps (YAML
        # "- {AA: A, Codon: GCA}") -> each row is a named list, so rbind would
        # leave list-valued columns.  Coerce explicitly to character columns.
        raw <- sim_cfg$codon.table
        data.frame(
            AA    = vapply(raw, function(r) as.character(r[["AA"]]),    character(1)),
            Codon = vapply(raw, function(r) as.character(r[["Codon"]]), character(1)),
            stringsAsFactors = FALSE)
    } else {
        .default_codon_table()
    }
    codon_order <- codon_table$Codon
    aa_of_codon <- setNames(codon_table$AA, codon_order)
    C <- length(codon_order)

    # ---- Truth values for alpha[c], lambda[c] ------------------------------
    truth_alpha  <- .load_truth_codon_values(truth_alpha_path,  "Mean", codon_order)
    truth_lambda <- .load_truth_codon_values(truth_lambda_path, "Mean", codon_order)
    truth_nse    <- setNames(rep(nse_shared, C), codon_order)

    if (verbose) {
        message("Codons: ", C, " across AAs: ",
                paste(sort(unique(aa_of_codon)), collapse = ","))
        message(sprintf("alpha range:  [%.3g, %.3g]", min(truth_alpha), max(truth_alpha)))
        message(sprintf("lambda range: [%.3g, %.3g]", min(truth_lambda), max(truth_lambda)))
        message(sprintf("NSE shared:   %.3g", nse_shared))
    }

    # ---- Synthesize genome (codon-at-position) -----------------------------
    set.seed(base_seed)
    gene_ids <- sprintf("gene%04d", seq_len(n_genes))

    codon_at_pos <- vector("list", n_genes)
    for (g in seq_len(n_genes)) {
        codon_at_pos[[g]] <- sample.int(C, size = gene_length, replace = TRUE)
    }

    # ---- phi[g] ~ LN(-sigma^2/2, sigma) * phi_scale ------------------------
    log_phi   <- rnorm(n_genes, mean = -0.5 * sigma_phi * sigma_phi, sd = sigma_phi)
    phi_truth <- exp(log_phi) * phi_scale
    if (verbose)
        message(sprintf("phi: mean=%.3g, median=%.3g, range=[%.3g, %.3g] (sigma=%.2f, scale=%.2g)",
                        mean(phi_truth), median(phi_truth),
                        min(phi_truth), max(phi_truth), sigma_phi, phi_scale))

    # ---- Pass 1: compute Z (process-specific) ------------------------------
    if (process == "nb-simulated") {
        v_per_codon <- 1.0 / truth_nse
        wait_times   <- vector("list", n_genes)
        sigma_at_pos <- vector("list", n_genes)
        Z <- 0.0
        for (g in seq_len(n_genes)) {
            cidx <- codon_at_pos[[g]]
            wt   <- rgamma(gene_length, shape = truth_alpha[cidx],
                           rate = truth_lambda[cidx])
            v_c   <- v_per_codon[cidx]
            ratio <- v_c / (wt + v_c)
            sigma_pos <- c(1.0, head(cumprod(ratio), gene_length - 1))
            wait_times[[g]]   <- wt
            sigma_at_pos[[g]] <- sigma_pos
            Z <- Z + phi_truth[[g]] * sum(wt * sigma_pos)
        }
        if (verbose)
            message(sprintf("Pass 1 (nb-simulated, wait_times + Z): Z = %.6g", Z))
    } else {
        v_per_codon <- 1.0 / truth_nse
        a_over_lv <- truth_alpha / (truth_lambda * v_per_codon)
        log_psuccess <- if (process == "nb-1o-approx") {
            -a_over_lv
        } else {
            -a_over_lv +
                a_over_lv / (truth_lambda * v_per_codon) +
                0.5 * a_over_lv * a_over_lv
        }
        a_over_l <- truth_alpha / truth_lambda
        sigma_at_pos <- vector("list", n_genes)
        wait_times   <- NULL
        Z <- 0.0
        for (g in seq_len(n_genes)) {
            cidx <- codon_at_pos[[g]]
            lp   <- log_psuccess[cidx]
            log_survive <- c(0, head(cumsum(lp), gene_length - 1))
            sg          <- exp(log_survive)
            sigma_at_pos[[g]] <- sg
            Z <- Z + phi_truth[[g]] * sum(a_over_l[cidx] * sg)
        }
        if (verbose)
            message(sprintf("Pass 1 (%s, E[sigma] via Taylor): Z = %.6g", process, Z))
    }

    # ---- Y_target -> U ------------------------------------------------------
    Y_target <- switch(
        y_target_mode,
        "U-one"   = Z,
        "fixed-Y" = {
            if (is.na(y_target_fixed))
                stop("y.target.mode = 'fixed-Y' requires y.target.value")
            y_target_fixed
        },
        stop("y.target.mode must be 'U-one' or 'fixed-Y'; got ", y_target_mode)
    )
    U_truth <- Z / Y_target
    if (verbose)
        message(sprintf("Y_target = %.6g (mode = %s); U_truth = %.6g",
                        Y_target, y_target_mode, U_truth))

    # ---- Pass 2: draw counts -----------------------------------------------
    rows <- vector("list", n_genes)
    if (process == "nb-simulated") {
        for (g in seq_len(n_genes)) {
            cidx <- codon_at_pos[[g]]
            wt   <- wait_times[[g]]
            sg   <- sigma_at_pos[[g]]
            mu   <- phi_truth[[g]] * wt * sg / U_truth
            y    <- rpois(gene_length, mu)
            rows[[g]] <- data.table(
                GeneID   = gene_ids[[g]],
                Position = seq_len(gene_length) - 1L,
                Codon    = codon_order[cidx],
                RFPCount = as.integer(y)
            )
        }
    } else {
        a_over_l <- truth_alpha / truth_lambda
        for (g in seq_len(n_genes)) {
            cidx <- codon_at_pos[[g]]
            sg   <- sigma_at_pos[[g]]
            mu   <- phi_truth[[g]] * a_over_l[cidx] * sg / U_truth
            y    <- rnbinom(gene_length, size = truth_alpha[cidx], mu = mu)
            rows[[g]] <- data.table(
                GeneID   = gene_ids[[g]],
                Position = seq_len(gene_length) - 1L,
                Codon    = codon_order[cidx],
                RFPCount = as.integer(y)
            )
        }
    }
    rfp <- rbindlist(rows)
    Y_obs <- sum(rfp$RFPCount)
    if (verbose)
        message(sprintf("Pass 2 (counts, %s): Y_obs = %d (target %.6g, U_obs = %.6g)",
                        process, Y_obs, Y_target, Z / Y_obs))

    list(
        rfp           = rfp,
        codon_table   = codon_table,
        truth_alpha   = truth_alpha,
        truth_lambda  = truth_lambda,
        truth_nse     = truth_nse,
        nse_shared    = nse_shared,
        phi_truth     = setNames(phi_truth, gene_ids),
        Z             = Z,
        Y_target      = Y_target,
        Y_obs         = Y_obs,
        U_truth       = U_truth,
        seed          = base_seed,
        n_genes       = n_genes,
        gene_length   = gene_length,
        process       = process
    )
}


# --------------------------------------------------------------------------
# Write outputs (internal helper)

.write_codon_csv <- function(df, value, out_path) {
    out <- data.frame(
        AA    = df$AA,
        Codon = df$Codon,
        Mean  = as.numeric(value[df$Codon]),
        stringsAsFactors = FALSE
    )
    write.csv(out, out_path, row.names = FALSE, quote = FALSE)
}


#' Write simulate_panse_dataset outputs to disk
#'
#' Writes the full set of output files expected by \code{build_panse_stan_data}
#' and analysis scripts: \code{rfp_counts.csv}, per-codon truth CSVs,
#' \code{truth_meta.rds}, \code{sim_config.yaml}, and \code{sim.log}.
#'
#' @param sim Return value from \code{\link{simulate_panse_dataset}}.
#' @param out_dir Output directory (created if it does not exist).
#' @param config The config list passed to \code{simulate_panse_dataset}.
#' @param verbose Logical; if \code{TRUE} (default) print paths written.
#' @return Invisibly, \code{out_dir}.
#' @export
write_sim_outputs <- function(sim, out_dir, config, verbose = TRUE) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    rfp_path <- file.path(out_dir, "rfp_counts.csv")
    fwrite(sim$rfp, rfp_path)

    .write_codon_csv(sim$codon_table, sim$truth_alpha,
                     file.path(out_dir, "truth_alpha.csv"))
    .write_codon_csv(sim$codon_table, sim$truth_lambda,
                     file.path(out_dir, "truth_lambda.csv"))
    .write_codon_csv(sim$codon_table, sim$truth_nse,
                     file.path(out_dir, "truth_nse.csv"))

    write.csv(
        data.frame(GeneID = names(sim$phi_truth),
                   phi    = as.numeric(sim$phi_truth)),
        file.path(out_dir, "truth_phi.csv"),
        row.names = FALSE, quote = FALSE
    )

    saveRDS(list(
        Z           = sim$Z,
        Y_target    = sim$Y_target,
        Y_obs       = sim$Y_obs,
        U_truth     = sim$U_truth,
        nse_shared  = sim$nse_shared,
        seed        = sim$seed,
        n_genes     = sim$n_genes,
        gene_length = sim$gene_length,
        codon_table = sim$codon_table,
        config      = config
    ), file.path(out_dir, "truth_meta.rds"))

    yaml::write_yaml(config, file.path(out_dir, "sim_config.yaml"))

    log_path <- file.path(out_dir, "sim.log")
    cat(sprintf(
        paste0("PANSE simulator log\n",
               "Timestamp:    %s\n",
               "Seed:         %d\n",
               "n_genes:      %d\n",
               "gene_length:  %d\n",
               "C (codons):   %d\n",
               "NSE_shared:   %.6g\n",
               "Z:            %.6g\n",
               "Y_target:     %.6g\n",
               "Y_obs:        %d\n",
               "U_truth:      %.6g\n",
               "phi range:    [%.4g, %.4g]\n",
               "phi mean:     %.4g (median %.4g)\n"),
        format(Sys.time()),
        sim$seed, sim$n_genes, sim$gene_length, nrow(sim$codon_table),
        sim$nse_shared, sim$Z, sim$Y_target, sim$Y_obs, sim$U_truth,
        min(sim$phi_truth), max(sim$phi_truth),
        mean(sim$phi_truth), median(sim$phi_truth)
    ), file = log_path)

    if (verbose) {
        cat("[written]\n")
        cat("  ", rfp_path, "\n")
        cat("  ", file.path(out_dir, "truth_alpha.csv"), "\n")
        cat("  ", file.path(out_dir, "truth_lambda.csv"), "\n")
        cat("  ", file.path(out_dir, "truth_nse.csv"), "\n")
        cat("  ", file.path(out_dir, "truth_phi.csv"), "\n")
        cat("  ", file.path(out_dir, "truth_meta.rds"), "\n")
        cat("  ", log_path, "\n")
    }

    invisible(out_dir)
}
