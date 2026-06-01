#!/usr/bin/env Rscript
# ============================================================================
# simulate_panse_dataset.R -- forward simulator for PANSE under controlled
# truth.  Pure R; supports THREE generative processes:
#
#   1. process = "nb-simulated"   (default; matches AnaCoDa::simulateGenome
#      exactly; cross-check against PANSEModel.cpp:1271-1399)
#
#      Pass 1 (wait_times + Z):
#        wait_time[g,p] ~ Gamma(shape = alpha[c], rate = lambda[c])
#        sigma[g,p]      = product_{p'<p} v[c'] / (wait_time[g,p'] + v[c'])
#                          (sigma[g,1] = 1; v[c] = 1 / NSERate[c])
#        Z += phi[g] * wait_time[g,p] * sigma[g,p]
#      Pass 2 (counts):
#        count[g,p] ~ Poisson(phi[g] * wait_time[g,p] * sigma[g,p] / U)
#        (same realized wait_times and sigma from pass 1)
#
#      Marginally (integrating wait_time), at a fixed position this is
#      Gamma-Poisson = NB2.  Across positions, sigma[g,p] is random in
#      upstream wait_times so the joint is NOT simply NB2 -- there is
#      extra variance from realized sigma that the nb-*-approx variants
#      lack.  See notes/panse.upstream-integration-approximations.md.
#
#   2. process = "nb-2o-approx"    (Stan-likelihood-matching, 2nd-order)
#
#      Replace random sigma[g,p] with E[sigma[g,p]] computed by the same
#      2nd-order Taylor expansion of log E[v/(W+v)] that the Stan model
#      uses (PANSEModel.cpp::elongationUntilIndexApproximation2-
#      ProbabilityLog):
#        log_psuccess_2o[c] = -a/(lv) + a/(l^2 v^2) + 0.5 (a/lv)^2
#        log_survive[g,p]   = sum_{p'<p} log_psuccess_2o[c[p']]
#        sigma_E[g,p]       = exp(log_survive[g,p])
#        count[g,p] ~ NB2(mu = phi * (alpha/lambda) * sigma_E / U,
#                         size = alpha[c])
#
#      This matches the Stan likelihood exactly; fitting Stan to this
#      data should give near-perfect recovery.
#
#   3. process = "nb-1o-approx"    (cheapest approximation; for testing)
#
#      Identical to nb-2o-approx but drop the higher-order Taylor terms,
#      keeping only the leading 1st-order in 1/(lv):
#        log_psuccess_1o[c] = -a/(lv)
#
#      Quantifies how much we lose by truncating to 1st order.  Useful
#      as a stepping-stone "is the answer roughly right" check before
#      paying for the 2nd-order term, and as a sanity-check that the
#      2nd-order improvement actually matters.
#
# All three processes share the same Z bookkeeping (Y_target controls U).
# Per-codon and per-gene mean counts agree across processes up to small
# bias from the approximation; per-position variance differs slightly
# (nb-simulated > nb-2o-approx ~ nb-1o-approx).
#
# For better-than-2nd-order options (exact closed form via incomplete
# gamma, generalized Gauss-Laguerre quadrature, etc.), see
# notes/panse.upstream-integration-approximations.md.
#
# Codon set (22 codons across 7 amino acids):
#   4 x 2-codon AAs:  Phe (TTT/TTC), Tyr (TAT/TAC), His (CAT/CAC),
#                     Asn (AAT/AAC)
#   2 x 4-codon AAs:  Val (GTT/GTC/GTA/GTG), Ala (GCT/GCC/GCA/GCG)
#   1 x 6-codon AA:   Leu (TTA/TTG/CTT/CTC/CTA/CTG)
#
# Truth values:
#   alpha[c], lambda[c]: read from RMF Parameter_est CSVs (config$truth.*)
#                        and subset to the chosen 22 codons.
#   NSE_shared:          single scalar specified in config (default 1e-5).
#   phi[g]:              drawn from LN(-sigma^2/2, sigma); mean(phi)=1.
#                        Scale by config$phi.scale (default 1).
#
# Outputs (under <out.dir>/):
#   rfp_counts.csv      RFP CSV in Wu/Weinberg format
#                       (GeneID, Position, Codon, RFPCount)
#   truth_alpha.csv     22-codon Alpha CSV  (AA, Codon, Mean)
#   truth_lambda.csv    22-codon Lambda_Prime CSV
#   truth_nse.csv       22-codon NSERate CSV (Codon, Mean = NSE_shared)
#   truth_phi.csv       per-gene phi (GeneID, phi)
#   truth_meta.rds      list(U, Z, Y, NSE_shared, codon_table, seed, ...)
#   sim_config.yaml     resolved config snapshot
#   sim.log             plain-text simulation log
#
# Usage:
#   Rscript scripts/sim/simulate_panse_dataset.R <config.yaml>
#   Rscript scripts/sim/simulate_panse_dataset.R <config.yaml> --seed 42
#
# Options:
#   --seed N        Override config$seed
#   --out DIR       Override config$out.dir
#   --dry-run       Print spec, skip simulation
#   -h, --help      Show this message
# ============================================================================

suppressPackageStartupMessages({
    library(yaml)
    library(data.table)
})

`%||%` <- function(x, y) if (is.null(x)) y else x


# --------------------------------------------------------------------------
# Default 22-codon set: 4 x 2-codon AAs, 2 x 4-codon, 1 x 6-codon.
# Edit here (or override via YAML config$codon.table) to change the set.

.default_codon_table <- function() {
    data.frame(
        AA    = c(rep("F", 2), rep("Y", 2), rep("H", 2), rep("N", 2),    # 2-codon
                  rep("V", 4), rep("A", 4),                              # 4-codon
                  rep("L", 6)),                                          # 6-codon
        Codon = c("TTT", "TTC", "TAT", "TAC", "CAT", "CAC", "AAT", "AAC",
                  "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG",
                  "TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
        stringsAsFactors = FALSE
    )
}


# --------------------------------------------------------------------------
# CLI parser

.parse_cli <- function(argv) {
    opts <- list(config = NULL, seed = NULL, out = NULL, dry_run = FALSE)
    i <- 1L
    while (i <= length(argv)) {
        a <- argv[[i]]
        if      (a == "--seed")      { opts$seed     <- as.integer(argv[[i + 1L]]); i <- i + 2L }
        else if (a == "--out")       { opts$out      <- argv[[i + 1L]];             i <- i + 2L }
        else if (a == "--dry-run")   { opts$dry_run  <- TRUE;                       i <- i + 1L }
        else if (a %in% c("-h", "--help")) {
            cat("Usage: simulate_panse_dataset.R <config.yaml> [--seed N] [--out DIR] [--dry-run]\n")
            quit(save = "no", status = 0)
        }
        else if (is.null(opts$config)) { opts$config <- a; i <- i + 1L }
        else stop("Unknown argument: ", a)
    }
    if (is.null(opts$config)) stop("Usage: simulate_panse_dataset.R <config.yaml> [...]")
    opts
}


# --------------------------------------------------------------------------
# Truth-value loaders

.load_truth_codon_values <- function(path, value_col, codon_order) {
    if (!file.exists(path))
        stop(".load_truth_codon_values: file not found: ", path)
    tbl <- read.csv(path, stringsAsFactors = FALSE)
    if (!"Codon" %in% colnames(tbl))
        stop(".load_truth_codon_values: missing `Codon` column in: ", path)
    if (!value_col %in% colnames(tbl))
        stop(".load_truth_codon_values: missing `", value_col, "` column in: ", path)
    # Subset tbl to only our codon_order codons FIRST, then assign.
    # (Previous version used `out[tbl$Codon] <- tbl[[value_col]]` which
    # grew the named vector beyond length(codon_order) when tbl had codons
    # not in our subset -- exactly what happens when we point at the
    # full 60/61-codon RMF Parameter_est CSV.)
    tbl <- tbl[tbl$Codon %in% codon_order, , drop = FALSE]
    missing_codons <- setdiff(codon_order, tbl$Codon)
    if (length(missing_codons) > 0L)
        stop(".load_truth_codon_values: missing codons in ", path, ": ",
             paste(missing_codons, collapse = ", "))
    out <- setNames(tbl[[value_col]], tbl$Codon)
    out[codon_order]    # reorder to match codon_order exactly
}


# --------------------------------------------------------------------------
# Main simulator

simulate_panse_dataset <- function(config, seed = NULL, verbose = TRUE) {
    sim_cfg <- config$sim %||% stop("config missing top-level `sim:` block")

    # ---- Resolve config -----------------------------------------------------
    n_genes        <- as.integer(sim_cfg$n.genes        %||% 1000L)
    gene_length    <- as.integer(sim_cfg$gene.length    %||% 500L)
    nse_shared     <- as.numeric(sim_cfg$nse.shared     %||% 1.0e-5)
    sigma_phi      <- as.numeric(sim_cfg$phi.sdlog      %||% 1.5)
    phi_scale      <- as.numeric(sim_cfg$phi.scale      %||% 1.0)
    y_target_mode  <- sim_cfg$y.target.mode %||% "U-one"  # 'U-one' or 'fixed-Y'
    y_target_fixed <- as.numeric(sim_cfg$y.target.value  %||% NA_real_)
    process        <- sim_cfg$process       %||% "nb-simulated"
    base_seed      <- as.integer(seed %||% sim_cfg$seed  %||% 20260524L)

    if (!(process %in% c("nb-simulated", "nb-2o-approx", "nb-1o-approx")))
        stop("sim.process must be one of 'nb-simulated', 'nb-2o-approx', ",
             "'nb-1o-approx'; got: ", process)

    truth_alpha_path  <- sim_cfg$truth.alpha.file  %||% stop("sim.truth.alpha.file missing")
    truth_lambda_path <- sim_cfg$truth.lambda.file %||% stop("sim.truth.lambda.file missing")

    # Codon table (default 22 codons; override via YAML if needed)
    codon_table <- if (!is.null(sim_cfg$codon.table)) {
        as.data.frame(do.call(rbind, sim_cfg$codon.table), stringsAsFactors = FALSE)
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

    # Uniform usage across the 22 codons (each appears with prob 1/22).
    # AA-level proportions implied: equal weight per codon = 2/22 for 2-codon
    # AAs (=0.091), 4/22 for 4-codon AAs (=0.18), 6/22 for the 6-codon AA
    # (=0.27).  Each codon expected to appear gene_length/C times per gene
    # (e.g., 500/22 ~= 23 positions/codon/gene).
    codon_at_pos <- vector("list", n_genes)
    for (g in seq_len(n_genes)) {
        codon_at_pos[[g]] <- sample.int(C, size = gene_length, replace = TRUE)
    }

    # ---- phi[g] ~ LN(-sigma^2/2, sigma) then * phi_scale --------------------
    log_phi   <- rnorm(n_genes, mean = -0.5 * sigma_phi * sigma_phi, sd = sigma_phi)
    phi_truth <- exp(log_phi) * phi_scale
    if (verbose)
        message(sprintf("phi: mean=%.3g, median=%.3g, range=[%.3g, %.3g] (sigma=%.2f, scale=%.2g)",
                        mean(phi_truth), median(phi_truth),
                        min(phi_truth), max(phi_truth), sigma_phi, phi_scale))

    # ---- Pass 1: compute Z (depends on process) ----------------------------
    if (process == "nb-simulated") {
        # Draw wait_times ~ Gamma per position; sigma from realized waits.
        v_per_codon <- 1.0 / truth_nse   # constant across codons here
        wait_times   <- vector("list", n_genes)
        sigma_at_pos <- vector("list", n_genes)
        Z <- 0.0
        for (g in seq_len(n_genes)) {
            cidx <- codon_at_pos[[g]]
            wt   <- rgamma(gene_length, shape = truth_alpha[cidx],
                           rate = truth_lambda[cidx])
            # Running sigma: sigma[1] = 1;
            # sigma[p+1] = sigma[p] * v[c_p]/(wt[p]+v[c_p])
            # (matches PANSEModel.cpp:1296-1300)
            v_c   <- v_per_codon[cidx]
            ratio <- v_c / (wt + v_c)
            sigma_pos <- c(1.0, head(cumprod(ratio), gene_length - 1))
            wait_times[[g]]   <- wt
            sigma_at_pos[[g]] <- sigma_pos
            Z <- Z + phi_truth[[g]] * sum(wt * sigma_pos)
        }
        if (verbose)
            message(sprintf(
                "Pass 1 (nb-simulated, wait_times + Z): Z = %.6g", Z))
    } else {
        # nb-{1o,2o}-approx: deterministic E[sigma] via Taylor of
        # log E[v/(W+v)] (matches PANSEModel.cpp's
        # elongationUntilIndexApproximation2ProbabilityLog and the Stan
        # model's transformed parameters block when order == 2).
        v_per_codon <- 1.0 / truth_nse
        a_over_lv <- truth_alpha / (truth_lambda * v_per_codon)
        log_psuccess <- if (process == "nb-1o-approx") {
            # 1st-order: keep only the leading -a/(lv) term.
            -a_over_lv
        } else {
            # 2nd-order: include the +a/(lv)^2 and +0.5*(a/(lv))^2 terms.
            -a_over_lv +
                a_over_lv / (truth_lambda * v_per_codon) +
                0.5 * a_over_lv * a_over_lv
        }
        a_over_l <- truth_alpha / truth_lambda  # E[wait_time]
        sigma_at_pos <- vector("list", n_genes)
        wait_times   <- NULL   # not used in this process
        Z <- 0.0
        for (g in seq_len(n_genes)) {
            cidx <- codon_at_pos[[g]]
            lp   <- log_psuccess[cidx]
            # log_survive[1] = 0; log_survive[p+1] = log_survive[p] + lp[p]
            log_survive <- c(0, head(cumsum(lp), gene_length - 1))
            sg          <- exp(log_survive)
            sigma_at_pos[[g]] <- sg
            # Z = sum phi * E[wait_time] * sigma_E
            #   = sum phi * (alpha/lambda) * sigma_E
            Z <- Z + phi_truth[[g]] * sum(a_over_l[cidx] * sg)
        }
        if (verbose)
            message(sprintf(
                "Pass 1 (%s, E[sigma] via Taylor): Z = %.6g", process, Z))
    }

    # ---- Y_target -> U ------------------------------------------------------
    Y_target <- switch(
        y_target_mode,
        "U-one"   = Z,                           # U = 1 by construction
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
        # nb-{1o,2o}-approx: y ~ NB2(mu, size = alpha[c]).  rnbinom with
        # parameterization (size, mu): Var = mu + mu^2/size.  Same draw
        # form regardless of which Taylor order produced sigma_at_pos.
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
# Write outputs

.write_codon_csv <- function(df, value, out_path) {
    out <- data.frame(
        AA    = df$AA,
        Codon = df$Codon,
        Mean  = as.numeric(value[df$Codon]),
        stringsAsFactors = FALSE
    )
    write.csv(out, out_path, row.names = FALSE, quote = FALSE)
}

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


# --------------------------------------------------------------------------
# CLI entry

if (sys.nframe() == 0L) {
    opts <- .parse_cli(commandArgs(trailingOnly = TRUE))
    cfg_abs <- normalizePath(opts$config)
    if (!file.exists(cfg_abs)) stop("Config not found: ", cfg_abs)

    # Match adapter.dev convention: setwd to adapter.dev/ so relative paths
    # in YAML resolve consistently.
    setwd(dirname(dirname(dirname(cfg_abs))))   # sim/<cfg>.yaml -> adapter.dev/
    cat("[setwd] ", getwd(), "\n", sep = "")

    config <- yaml.load_file(cfg_abs)
    if (!is.null(opts$out))  config$sim$out.dir <- opts$out
    out_dir <- config$sim$out.dir %||% stop("config$sim$out.dir missing")

    if (opts$dry_run) {
        cat("[--dry-run] config:\n")
        str(config$sim)
        quit(save = "no", status = 0)
    }

    sim <- simulate_panse_dataset(config, seed = opts$seed)
    write_sim_outputs(sim, out_dir, config)
    cat("\n[ok] sim outputs -> ", normalizePath(out_dir), "\n", sep = "")
}
