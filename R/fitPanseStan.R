# fitPanseStan.R -- exported fit_panse_stan() function (PANSE Stan/HMC driver).
#
# Logic lifted from scripts/panse-stan/fit.stan.R; CLI arg parsing moved to
# inst/scripts/fit_panse_stan.R shim.  Key changes vs the script:
#   - Stan model resolved via .panse_stan_file() (inst/stan/) instead of
#     file.path(rmf_root, "stan", ...).
#   - build_dir defaults to tools::R_user_dir("AnaCoDa","cache").
#   - .load_init_util / RMF_INIT_UTILS_DIR removed; use package exports directly.
#   - source() calls removed; build_panse_stan_data / panse_stan_finalize are
#     now in-package.
#   - cmdstanr guarded by requireNamespace; stays in Suggests.


# Internal: read per-codon Std.Dev column for rmf-posterior init.
.read_codon_sd <- function(path, codon_order) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    tbl <- read.csv(path, stringsAsFactors = FALSE)
    if (!("Codon" %in% colnames(tbl)) || !("Std.Dev" %in% colnames(tbl)))
        return(NULL)
    setNames(tbl$Std.Dev, tbl$Codon)[codon_order]
}

# Internal: map config$fit$model + parameterization -> stan_basename + flags.
.panse_select_model <- function(config) {
    model_key       <- config$fit$model          %||% "csp-only"
    parameterization <- config$fit$parameterization %||% "centered"
    stan_basename <- switch(
        paste0(model_key, ":", parameterization),
        "csp-only:centered"              = "panse_csp_only",
        "csp-only:noncentered"           = "panse_csp_only",
        "csp-only-sharednse:centered"    = "panse_csp_only_sharednse",
        "basic:centered"                 = "panse_basic",
        "basic-sharednse:centered"       = "panse_basic_sharednse",
        "sphi-est:centered"              = "panse_sphi_est_centered",
        "sphi-est-sharednse:centered"    = "panse_sphi_est_centered_sharednse",
        "sphi-est:noncentered"           = "panse_sphi_est_noncentered",
        "sphi-est:sumzero"               = "panse_sphi_est_sumzero",
        "sphi-est:sumzero-gm"            = "panse_sphi_est_sumzero_geomean",
        "sphi-est-sharednse:noncentered" = "panse_sphi_est_noncentered_sharednse",
        "sphi-est-sharednse:sumzero"     = "panse_sphi_est_sumzero_sharednse",
        stop("Unsupported fit.model + parameterization combo: ", model_key,
             " / ", parameterization)
    )
    list(
        stan_basename    = stan_basename,
        shared_nse       = model_key %in% c("csp-only-sharednse", "basic-sharednse",
                                            "sphi-est-sharednse"),
        samples_phi      = model_key %in% c("basic", "basic-sharednse",
                                            "sphi-est", "sphi-est-sharednse"),
        estimates_sphi   = model_key %in% c("sphi-est", "sphi-est-sharednse"),
        noncentered_phi  = (model_key %in% c("sphi-est", "sphi-est-sharednse") &&
                            parameterization %in% c("noncentered", "sumzero", "sumzero-gm")),
        sumzero_phi      = (model_key %in% c("sphi-est", "sphi-est-sharednse") &&
                            parameterization %in% c("sumzero", "sumzero-gm"))
    )
}


#' Fit a PANSE Stan/HMC model from a YAML config
#'
#' Assembles the Stan data list, compiles the appropriate \code{panse_*.stan}
#' model, samples with \code{cmdstanr}, and calls
#' \code{\link{panse_stan_finalize}} to write outputs and print diagnostics.
#' Returns the \code{CmdStanMCMC} fit object invisibly.
#'
#' @param config Parsed YAML config list (as returned by
#'   \code{yaml::read_yaml}).  Must have \code{name}, \code{fit}, and
#'   \code{genome} blocks.
#' @param out_dir Output directory.  If \code{NULL} (default), derived as
#'   \code{file.path("Output", paste0(config$name, "-stan"))} relative to the
#'   calling process's working directory.
#' @param gene_subset Integer N (first N genes) or character gene-ID vector
#'   passed to \code{\link{build_panse_stan_data}}.
#' @param init_mode Init strategy: \code{"fixed"}, \code{"rmf-posterior"},
#'   \code{"scuo"}, \code{"enc_prime"}, \code{"mixed_scuo_encp"}, or
#'   \code{"warm-start"}.  \code{NULL} (default) falls back to
#'   \code{config$fit$init.mode}, then \code{"fixed"}.
#' @param warm_start_from Path to a prior fit output directory (overrides
#'   \code{config$fit$warm_start_from} when non-\code{NULL}).
#' @param no_log_lik Logical; if \code{TRUE}, override
#'   \code{stan_data$emit_log_lik <- 0L}.
#' @param no_compile Logical; if \code{TRUE}, skip recompilation and require a
#'   precompiled exe in \code{build_dir}.
#' @param dry_run Logical; if \code{TRUE}, build data + config but skip
#'   sampling.
#' @param build_dir Writable cache for compiled Stan binaries.  Default:
#'   \code{tools::R_user_dir("AnaCoDa", "cache")}.
#' @param chains,warmup,sampling,threads,adapt_delta,seed Sampler settings;
#'   \code{NULL} (default) defers to \code{config$stan}.
#' @param git_sha Short git SHA string recorded in the fit RDS.  If
#'   \code{NULL} (default), attempts \code{git rev-parse --short HEAD}; falls
#'   back to \code{"unknown"}.
#' @param verbose Logical; if \code{TRUE} (default) print progress.
#' @return Invisibly, the \code{CmdStanMCMC} fit object.
#' @export
fit_panse_stan <- function(config,
                           out_dir       = NULL,
                           gene_subset   = NULL,
                           init_mode     = NULL,
                           warm_start_from = NULL,
                           no_log_lik    = FALSE,
                           no_compile    = FALSE,
                           dry_run       = FALSE,
                           build_dir     = tools::R_user_dir("AnaCoDa", "cache"),
                           chains        = NULL,
                           warmup        = NULL,
                           sampling      = NULL,
                           threads       = NULL,
                           adapt_delta   = NULL,
                           seed          = NULL,
                           git_sha       = NULL,
                           verbose       = TRUE) {

    if (!requireNamespace("cmdstanr", quietly = TRUE))
        stop("fit_panse_stan requires cmdstanr. ",
             "Install with: install.packages('cmdstanr', repos='https://mc-stan.org/r-packages/')")

    if (is.null(config$name))
        stop("Config missing top-level `name`")

    # ---- Model selection ----------------------------------------------------
    ms <- .panse_select_model(config)
    stan_basename  <- ms$stan_basename
    shared_nse     <- ms$shared_nse
    samples_phi    <- ms$samples_phi
    estimates_sphi <- ms$estimates_sphi
    noncentered_phi <- ms$noncentered_phi
    sumzero_phi    <- ms$sumzero_phi

    # Resolve stan file + build dir
    stan_file <- .panse_stan_file(paste0(stan_basename, ".stan"))
    if (!dir.exists(build_dir)) dir.create(build_dir, recursive = TRUE)

    if (is.null(git_sha)) {
        git_sha <- tryCatch(
            system2("git", c("rev-parse", "--short", "HEAD"),
                    stdout = TRUE, stderr = FALSE),
            error = function(e) "unknown"
        )
        if (length(git_sha) == 0 || nzchar(git_sha[[1]]) == FALSE)
            git_sha <- "unknown"
        git_sha <- git_sha[[1]]
    }
    cat("[model] ", stan_basename, "\n",
        "[build] ", build_dir, "\n", sep = "")

    # ---- Build Stan data ----------------------------------------------------
    stan_data <- build_panse_stan_data(config, gene_subset = gene_subset,
                                       verbose = verbose)
    if (no_log_lik) {
        stan_data$emit_log_lik <- 0L
        cat("[data] --no-log-lik: skipping log_lik in generated quantities\n")
    }
    cat("[data] G =", stan_data$G, " C =", stan_data$C, " P =", stan_data$P,
        " Y =", sum(stan_data$y), " U =", round(stan_data$U, 6),
        " emit_log_lik =", stan_data$emit_log_lik, "\n")

    # ---- Compile (or load existing exe) -------------------------------------
    exe_path <- file.path(build_dir, stan_basename)
    if (!no_compile || !file.exists(exe_path)) {
        cat("[compile] ", stan_file, "\n", sep = "")
        mod <- cmdstanr::cmdstan_model(
            stan_file   = stan_file,
            dir         = build_dir,
            cpp_options = list(stan_threads = TRUE),
            quiet       = TRUE
        )
    } else {
        cat("[load] precompiled ", exe_path, "\n", sep = "")
        mod <- cmdstanr::cmdstan_model(exe_file = exe_path,
                                       stan_file = stan_file,
                                       compile  = FALSE)
    }

    # ---- Init function ------------------------------------------------------
    init_mode <- init_mode %||% config$fit$init.mode %||% "fixed"
    warm_start_metrics <- NULL

    init_alpha  <- attr(stan_data, "init_alpha")
    init_lambda <- attr(stan_data, "init_lambda")
    init_nse    <- attr(stan_data, "init_nse")
    codon_order <- attr(stan_data, "codon_order")

    if (init_mode == "rmf-posterior") {
        alpha_path  <- config$fit$init.alpha.files[[1]]
        lambda_path <- config$fit$init.lambda.files[[1]]
        a_sd <- .read_codon_sd(alpha_path,  codon_order)
        l_sd <- .read_codon_sd(lambda_path, codon_order)
        if (is.null(a_sd) || is.null(l_sd))
            stop("init.mode = 'rmf-posterior' requires init.alpha.files / ",
                 "init.lambda.files with a `Std.Dev` column (e.g. RMF ",
                 "Parameter_est/*.csv).  Got: ", alpha_path, " / ", lambda_path)

        log_a_mean <- log(init_alpha);  log_a_jsd <- a_sd / init_alpha
        log_l_mean <- log(init_lambda); log_l_jsd <- l_sd / init_lambda
        log_nse_mean <- log(init_nse)
        log_nse_jsd  <- rep((stan_data$log_nse_upper - stan_data$log_nse_lower) / 6,
                            stan_data$C)

        cat(sprintf("[init] mode = rmf-posterior (per-chain jitter)\n"))
        cat(sprintf("[init]   log_alpha       jitter SD: [%.3f, %.3f]\n",
                    min(log_a_jsd), max(log_a_jsd)))
        cat(sprintf("[init]   log_lambdaPrime jitter SD: [%.3f, %.3f]\n",
                    min(log_l_jsd), max(log_l_jsd)))
        cat(sprintf("[init]   log_NSERate     jitter SD: %.3f (fallback; no RMF SD)\n",
                    log_nse_jsd[[1]]))

        base_seed_init <- (config$stan %||% list())$seed %||% 20260523L
        init_fn <- function(chain_id = 1) {
            set.seed(base_seed_init + 1000L * chain_id)
            log_nse <- rnorm(stan_data$C, log_nse_mean, log_nse_jsd)
            log_nse <- pmin(pmax(log_nse,
                                 stan_data$log_nse_lower + 1e-6),
                            stan_data$log_nse_upper - 1e-6)
            list(
                log_alpha       = rnorm(stan_data$C, log_a_mean, log_a_jsd),
                log_lambdaPrime = rnorm(stan_data$C, log_l_mean, log_l_jsd),
                log_NSERate     = log_nse
            )
        }
    } else if (init_mode == "fixed") {
        cat("[init] mode = fixed (all chains start from init CSV values, no jitter)\n")
        init_fn <- function(chain_id = 1) {
            log_nse <- pmin(pmax(log(init_nse),
                                 stan_data$log_nse_lower + 1e-6),
                            stan_data$log_nse_upper - 1e-6)
            list(
                log_alpha       = log(init_alpha),
                log_lambdaPrime = log(init_lambda),
                log_NSERate     = log_nse
            )
        }
    } else if (init_mode %in% c("scuo", "enc_prime", "mixed_scuo_encp")) {
        cat(sprintf("[init] mode = %s (CSP from init CSV files; phi overridden below)\n", init_mode))
        init_fn <- function(chain_id = 1) {
            log_nse <- pmin(pmax(log(init_nse),
                                 stan_data$log_nse_lower + 1e-6),
                            stan_data$log_nse_upper - 1e-6)
            list(
                log_alpha       = log(init_alpha),
                log_lambdaPrime = log(init_lambda),
                log_NSERate     = log_nse
            )
        }
    } else if (init_mode == "warm-start") {
        ws_dir <- warm_start_from %||% config$fit$warm_start_from
        if (is.null(ws_dir)) stop("warm-start requires fit.warm_start_from in YAML or warm_start_from arg")
        ws_rds <- file.path(ws_dir, "panse-stan-fit.rds")
        if (!file.exists(ws_rds)) stop("warm-start fit not found: ", ws_rds)

        cat(sprintf("[init] mode = warm-start from %s\n", ws_dir))
        ws_obj  <- readRDS(ws_rds)
        ws_fit  <- ws_obj$fit
        ws_data <- ws_obj$stan_data

        ws_summ_path <- file.path(ws_dir, "stan-summary.rds")
        if (file.exists(ws_summ_path)) {
            ws_summ <- readRDS(ws_summ_path)
            ws_med  <- setNames(ws_summ$mean, ws_summ$variable)
            cat("[init] warm-start: loaded prior summary (using posterior means)\n")
        } else {
            ws_draws <- ws_fit$draws(format = "matrix")
            ws_med   <- apply(ws_draws, 2, median)
            names(ws_med) <- colnames(ws_draws)
            cat("[init] warm-start: computed medians from draws\n")
        }

        .ws_vec <- function(log_stem, nat_stem, n) {
            v <- ws_med[grep(paste0("^", log_stem, "\\["), names(ws_med))]
            if (length(v) == n) return(as.numeric(v))
            v <- ws_med[grep(paste0("^", nat_stem, "\\["), names(ws_med))]
            if (length(v) == n) return(log(as.numeric(v)))
            stop("warm-start: neither ", log_stem, " nor ", nat_stem,
                 " found with length ", n)
        }
        ws_log_alpha  <- .ws_vec("log_alpha", "alpha", stan_data$C)
        ws_log_lambda <- .ws_vec("log_lambdaPrime", "lambdaPrime", stan_data$C)
        ws_sphi       <- ws_med[["sphi"]]

        ws_z_phi <- ws_med[grep("^z_phi\\[", names(ws_med))]
        if (length(ws_z_phi) == stan_data$G) {
            ws_z_phi <- as.numeric(ws_z_phi)
        } else {
            ws_log_phi <- ws_med[grep("^log_phi\\[", names(ws_med))]
            ws_z_phi   <- (as.numeric(ws_log_phi) + 0.5 * ws_sphi^2) / ws_sphi
            cat(sprintf("[init] warm-start: derived z_phi from log_phi (sphi=%.4f)\n", ws_sphi))
        }

        .ws_get <- function(nm) if (nm %in% names(ws_med)) ws_med[[nm]] else NULL
        ws_nse_shared <- .ws_get("log_NSERate_shared")
        if (is.null(ws_nse_shared)) {
            nat <- .ws_get("NSERate_shared")
            if (!is.null(nat)) ws_nse_shared <- log(nat)
        }
        ws_nse_vec <- ws_med[grep("^log_NSERate\\[", names(ws_med))]
        if (length(ws_nse_vec) == 0) {
            nat_vec <- ws_med[grep("^NSERate\\[", names(ws_med))]
            if (length(nat_vec) > 0) ws_nse_vec <- log(as.numeric(nat_vec))
        }

        if (!is.null(ws_nse_shared) && length(ws_nse_vec) == 0) {
            ws_log_nse <- rep(ws_nse_shared, stan_data$C)
            cat(sprintf("[init] warm-start: expanding NSERate_shared (%.3e) -> %d per-codon\n",
                        exp(ws_nse_shared), stan_data$C))
        } else if (length(ws_nse_vec) == stan_data$C) {
            ws_log_nse <- as.numeric(ws_nse_vec)
        } else {
            stop("warm-start: can't resolve NSE from prior fit")
        }
        ws_log_nse <- pmin(pmax(ws_log_nse,
                                stan_data$log_nse_lower + 1e-6),
                           stan_data$log_nse_upper - 1e-6)

        cat(sprintf("[init] warm-start: alpha[%d] lambda[%d] nse[%d] sphi=%.4f z_phi[%d]\n",
                    length(ws_log_alpha), length(ws_log_lambda),
                    length(ws_log_nse), ws_sphi, length(ws_z_phi)))

        init_fn <- function(chain_id = 1) {
            list(
                log_alpha       = as.numeric(ws_log_alpha),
                log_lambdaPrime = as.numeric(ws_log_lambda),
                log_NSERate     = ws_log_nse
            )
        }
        gene_ids <- attr(stan_data, "gene_ids")
        ws_log_phi_vals <- as.numeric(ws_z_phi) * ws_sphi - 0.5 * ws_sphi^2
        init_phi_attr <- setNames(exp(ws_log_phi_vals), gene_ids)
        sphi_init_override <- ws_sphi
        z_phi_init_override <- as.numeric(ws_z_phi)

        ws_metrics_raw <- tryCatch(ws_fit$inv_metric(matrix = TRUE),
                                   error = function(e) NULL)
        if (!is.null(ws_metrics_raw)) {
            n_chains_ws <- length(ws_metrics_raw)
            ws_C <- ws_data$C
            n_csp <- 2L * ws_C
            ws_has_shared <- (!is.null(ws_nse_shared) &&
                              (is.null(ws_nse_vec) || length(ws_nse_vec) == 0))

            n_chains_new <- chains %||% (config$stan %||% list())$chains %||% 4L
            ws_metric_type <- (config$stan %||% list())$metric %||% "diag_e"
            warm_start_metrics <- vector("list", n_chains_new)
            for (wi in seq_len(n_chains_new)) {
                wm <- ws_metrics_raw[[(wi - 1L) %% n_chains_ws + 1L]]
                wm_diag <- if (is.matrix(wm)) diag(wm) else wm
                if (ws_has_shared) {
                    csp_block  <- wm_diag[1:n_csp]
                    nse_metric <- wm_diag[n_csp + 1L]
                    tail_block <- wm_diag[(n_csp + 2L):length(wm_diag)]
                    new_diag   <- c(csp_block, rep(nse_metric, stan_data$C), tail_block)
                } else {
                    new_diag <- wm_diag
                }
                warm_start_metrics[[wi]] <- if (ws_metric_type == "dense_e") {
                    diag(new_diag)
                } else {
                    new_diag
                }
            }
            n_orig <- length(if (is.matrix(ws_metrics_raw[[1]])) {
                diag(ws_metrics_raw[[1]])
            } else {
                ws_metrics_raw[[1]]
            })
            n_new  <- length(new_diag)
            cat(sprintf("[init] warm-start: inv_metric %s diagonal (%d -> %d params)\n",
                        ws_metric_type, n_orig, n_new))
        } else {
            cat("[warn] warm-start: inv_metric extraction failed; will re-adapt from scratch\n")
        }
    } else {
        stop("init.mode must be one of: fixed, rmf-posterior, scuo, enc_prime, ",
             "mixed_scuo_encp, warm-start; got: ", init_mode)
    }

    # Shared-NSE collapse: log_NSERate[C] -> log_NSERate_shared scalar
    if (shared_nse) {
        cat("[init] shared-NSE model: collapsing log_NSERate -> log_NSERate_shared\n")
        inner_init_fn <- init_fn
        init_fn <- function(chain_id = 1) {
            ll <- inner_init_fn(chain_id)
            ll$log_NSERate_shared <- mean(ll$log_NSERate)
            ll$log_NSERate <- NULL
            ll
        }
    }

    # Phi init for phi-sampling models
    init_phi_attr <- attr(stan_data, "init_phi")

    if (init_mode %in% c("scuo", "enc_prime", "mixed_scuo_encp")) {
        gene_ids <- attr(stan_data, "gene_ids")
        rfp_csv  <- Sys.glob(file.path(config$genome$input.dir,
                                       config$genome$pattern))[[1]]
        cat(sprintf("[init] loading genome from %s\n", basename(rfp_csv)))
        genome <- initializeGenomeObject(rfp_csv, fasta = FALSE, positional = TRUE)
    }

    if (init_mode %in% c("scuo", "mixed_scuo_encp")) {
        cat("[init] computing SCUO\n")
        scuo_df      <- calculateSCUO(genome)
        scuo_all     <- setNames(scuo_df$SCUO, scuo_df$ORF)
        log_phi_scuo <- scuo_to_log_phi(scuo_all[gene_ids])
        cat(sprintf("[init] scuo: G=%d  phi in [%.3f, %.3f]  sd(log phi)=%.4f\n",
                    length(gene_ids), min(exp(log_phi_scuo)), max(exp(log_phi_scuo)),
                    sd(log_phi_scuo)))
    }

    if (init_mode %in% c("enc_prime", "mixed_scuo_encp")) {
        cat("[init] computing ENC' (Novembre 2002)\n")
        cc       <- as.matrix(getCodonCounts(genome))
        mode(cc) <- "integer"
        cc       <- cc[gene_ids, , drop = FALSE]
        aa_cmap  <- build_aa_codon_map()
        null_f   <- derive_null_from_genome(cc, aa_cmap)
        encp_raw <- calc_enc_prime(cc, aa_cmap, null_freqs = null_f)
        log_phi_encp <- scuo_to_log_phi(-encp_raw[gene_ids])
        cat(sprintf("[init] enc_prime: G=%d  phi in [%.3f, %.3f]  sd(log phi)=%.4f\n",
                    length(gene_ids), min(exp(log_phi_encp)), max(exp(log_phi_encp)),
                    sd(log_phi_encp)))
    }

    if (init_mode == "scuo")
        init_phi_attr <- setNames(exp(log_phi_scuo), gene_ids)
    if (init_mode == "enc_prime")
        init_phi_attr <- setNames(exp(log_phi_encp), gene_ids)

    if (samples_phi && !noncentered_phi) {
        cat("[init] phi-sampling model (centered): adding log_phi vector\n")
        inner_phi_init_fn <- init_fn
        init_fn <- function(chain_id = 1) {
            ll <- inner_phi_init_fn(chain_id)
            ll$log_phi <- log(as.numeric(init_phi_attr))
            ll
        }
    }

    if (estimates_sphi) {
        sphi_init <- if (exists("sphi_init_override")) sphi_init_override
                     else sd(log(as.numeric(init_phi_attr)))
        cat(sprintf("[init] sphi-est model: adding sphi init = %.4f\n", sphi_init))
        inner_sphi_init_fn <- init_fn
        init_fn <- function(chain_id = 1) {
            ll <- inner_sphi_init_fn(chain_id)
            ll$sphi <- sphi_init
            ll
        }
    }

    if (noncentered_phi) {
        if (exists("z_phi_init_override")) {
            z_phi_init <- z_phi_init_override
        } else {
            log_phi_truth <- log(as.numeric(init_phi_attr))
            z_phi_init    <- (log_phi_truth + 0.5 * sphi_init * sphi_init) / sphi_init
        }
        if (sumzero_phi) {
            z_phi_init <- z_phi_init - mean(z_phi_init)
            cat(sprintf("[init] sumzero phi: centered z_phi to mean 0 (residual = %.2e)\n",
                        mean(z_phi_init)))
        }
        cat(sprintf("[init] noncentered phi: z_phi init range [%.3f, %.3f] (G=%d)\n",
                    min(z_phi_init), max(z_phi_init), length(z_phi_init)))
        inner_z_init_fn <- init_fn
        init_fn <- function(chain_id = 1) {
            ll <- inner_z_init_fn(chain_id)
            ll$z_phi <- z_phi_init
            ll
        }
    }

    if (init_mode == "mixed_scuo_encp" && noncentered_phi) {
        z_scuo <- (log_phi_scuo + 0.5 * sphi_init^2) / sphi_init
        z_encp <- (log_phi_encp + 0.5 * sphi_init^2) / sphi_init
        cat(sprintf(paste0("[init] mixed_scuo_encp: chains 1-2 SCUO z_phi [%.3f, %.3f];",
                           " chains 3-4 ENC' z_phi [%.3f, %.3f]\n"),
                    min(z_scuo), max(z_scuo), min(z_encp), max(z_encp)))
        inner_mixed_fn <- init_fn
        init_fn <- function(chain_id = 1) {
            ll       <- inner_mixed_fn(chain_id)
            ll$z_phi <- if (chain_id <= 2L) z_scuo else z_encp
            ll$sphi  <- sphi_init
            ll
        }
    }

    # ---- Sampler settings ---------------------------------------------------
    stan_cfg  <- config$stan %||% list()
    n_chains  <- chains      %||% stan_cfg$chains             %||% 4L
    n_threads <- threads     %||% stan_cfg$threads_per_chain  %||% 1L
    n_warmup  <- warmup      %||% stan_cfg$warmup             %||% 1000L
    n_sample  <- sampling    %||% stan_cfg$sampling           %||% 1000L
    adapt_dl  <- adapt_delta %||% stan_cfg$adapt_delta        %||% 0.9
    max_td    <- stan_cfg$max_treedepth    %||% 10L
    parallel  <- min(n_chains, stan_cfg$parallel_chains %||% n_chains)
    run_seed  <- seed %||% stan_cfg$seed %||% 20260523L

    cat("[sample] chains=", n_chains, " parallel=", parallel,
        " warmup=", n_warmup, " sampling=", n_sample,
        " threads_per_chain=", n_threads, " adapt_delta=", adapt_dl,
        " max_treedepth=", max_td, "\n", sep = "")

    # ---- Output dir ---------------------------------------------------------
    if (is.null(out_dir))
        out_dir <- file.path("Output", paste0(config$name, "-stan"))
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    cat("[out] ", normalizePath(out_dir), "\n", sep = "")

    saveRDS(stan_data, file.path(out_dir, "stan-data.rds"))
    yaml::write_yaml(config, file.path(out_dir, "config.yaml"))

    if (dry_run) {
        cat("[--dry-run] skipping sampling\n")
        return(invisible(NULL))
    }

    # ---- Sample -------------------------------------------------------------
    stan_csv_dir <- file.path(out_dir, "stan-csv")
    if (!dir.exists(stan_csv_dir)) dir.create(stan_csv_dir, recursive = TRUE)

    metric_choice <- stan_cfg$metric %||% "diag_e"
    if (!metric_choice %in% c("diag_e", "dense_e"))
        stop("stan.metric must be 'diag_e' or 'dense_e'; got: ", metric_choice)
    cat("[sample] metric=", metric_choice, "\n", sep = "")

    sample_args <- list(
        data              = stan_data,
        init              = init_fn,
        chains            = n_chains,
        parallel_chains   = parallel,
        threads_per_chain = n_threads,
        iter_warmup       = n_warmup,
        iter_sampling     = n_sample,
        adapt_delta       = adapt_dl,
        max_treedepth     = max_td,
        metric            = metric_choice,
        seed              = run_seed,
        refresh           = max(1L, (n_warmup + n_sample) %/% 20L),
        output_dir        = stan_csv_dir
    )
    if (!is.null(warm_start_metrics)) {
        sample_args$inv_metric <- warm_start_metrics
        cat(sprintf("[sample] warm-start: seeding inv_metric (%d entries/chain)\n",
                    length(warm_start_metrics[[1]])))
    }

    t0 <- proc.time()[["elapsed"]]
    fit <- do.call(mod$sample, sample_args)
    wall_sec <- proc.time()[["elapsed"]] - t0
    cat(sprintf("\n[done] sample wall = %.1f sec (%.2f min)\n",
                wall_sec, wall_sec / 60))

    # ---- Diagnostics + save -------------------------------------------------
    panse_stan_finalize(fit, config, stan_data, out_dir, wall_sec, git_sha,
                        diagnose_fn = function(fit) fit$cmdstan_diagnose())

    invisible(fit)
}
