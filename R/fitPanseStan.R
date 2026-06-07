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
    # Orthogonal dwell-time-distribution axis (model comparison via LOO).
    # gamma    : W ~ Gamma(alpha, lambda)  -> NB2 counts + exact survival (current)
    # singular : W = E[W] deterministic    -> Poisson counts + point-value survival
    # invgamma : W ~ InvGamma(alpha,lambda)-> Sichel counts (NB2 surrogate) + exact survival
    # See panse/notes/dwell-time-distributions.md (Analyses-RibModelFramework).
    dwell_dist <- config$fit$dwell.dist %||% "gamma"
    if (!dwell_dist %in% c("gamma", "singular", "invgamma"))
        stop("Unsupported fit.dwell.dist '", dwell_dist,
             "'; valid: gamma, singular, invgamma")
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
        "sphi-est-pointval:noncentered"  = "panse_sphi_est_noncentered_pointval",
        "sphi-est-sharednse:noncentered" = "panse_sphi_est_noncentered_sharednse",
        "sphi-est-sharednse:sumzero"     = "panse_sphi_est_sumzero_sharednse",
        "sphi-est-aanse:noncentered"     = "panse_sphi_est_noncentered_aanse",
        stop("Unsupported fit.model + parameterization combo: ", model_key,
             " / ", parameterization)
    )
    # Apply the dwell.dist axis.  Only the per-codon-NSE sphi-est:noncentered
    # base has singular/invgamma variants; gamma is the unchanged default for
    # every combo.  Reject non-gamma dwell on combos that lack a variant so a
    # misplaced dwell.dist fails loudly instead of silently fitting gamma.
    if (dwell_dist != "gamma") {
        if (stan_basename != "panse_sphi_est_noncentered")
            stop("fit.dwell.dist '", dwell_dist, "' is only supported with ",
                 "fit.model 'sphi-est' + parameterization 'noncentered' ",
                 "(resolved base '", stan_basename, "')")
        stan_basename <- switch(
            dwell_dist,
            "singular" = "panse_sphi_est_noncentered_singular",
            "invgamma" = "panse_sphi_est_noncentered_invgamma"
        )
    }
    list(
        dwell_dist       = dwell_dist,
        stan_basename    = stan_basename,
        shared_nse       = model_key %in% c("csp-only-sharednse", "basic-sharednse",
                                            "sphi-est-sharednse"),
        aa_nse           = model_key == "sphi-est-aanse",
        samples_phi      = model_key %in% c("basic", "basic-sharednse",
                                            "sphi-est", "sphi-est-pointval",
                                            "sphi-est-sharednse", "sphi-est-aanse"),
        estimates_sphi   = model_key %in% c("sphi-est", "sphi-est-pointval",
                                            "sphi-est-sharednse", "sphi-est-aanse"),
        noncentered_phi  = (model_key %in% c("sphi-est", "sphi-est-pointval",
                                             "sphi-est-sharednse", "sphi-est-aanse") &&
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
#'   \code{"scuo"}, \code{"enc_prime"}, \code{"mixed_scuo_encp"},
#'   \code{"warm-start"}, or \code{"advi"} (ADVI variational warm-start).
#'   \code{NULL} (default) falls back to \code{config$fit$init.mode}, then
#'   \code{"fixed"}.  With \code{"advi"}, the model is first fitted by
#'   mean-field variational inference; the posterior means seed the HMC init
#'   and the marginal variances (transformed to unconstrained space) seed the
#'   diagonal inverse mass matrix, reducing warmup adaptation cost.
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
    aa_nse         <- ms$aa_nse
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

    # ---- AA-codon mapping (aanse model) ----------------------------------------
    # Computed here (before init block) only because aanse stan_data needs N_AA
    # and aa_of_codon.  For ADVI with aa-nse source on a non-aanse target, the
    # mapping is built lazily in the ADVI block below.
    aa_of_codon_vec <- NULL
    n_aa_groups     <- NULL
    if (aa_nse) {
        aa_map <- .build_nse_aa_map()   # 21 groups: 19 bias + M + W singletons
        aa_map          <- aa_map[sort(names(aa_map))]   # alphabetical order
        aa_names_ord    <- names(aa_map)
        codon_to_aa_idx <- setNames(
            rep(seq_along(aa_names_ord), lengths(aa_map)),
            unlist(aa_map, use.names = FALSE)
        )
        codon_order_local <- attr(stan_data, "codon_order")
        aa_of_codon_vec   <- unname(codon_to_aa_idx[codon_order_local])
        n_aa_groups       <- length(aa_names_ord)
        if (any(is.na(aa_of_codon_vec)))
            stop("[data] aa_of_codon_vec has NAs for codons: ",
                 paste(codon_order_local[is.na(aa_of_codon_vec)], collapse = ", "))
        stan_data$N_AA        <- n_aa_groups
        stan_data$aa_of_codon <- aa_of_codon_vec
        cat("[data] AA-NSE model: N_AA =", n_aa_groups, "\n")
    }

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
    } else if (init_mode %in% c("advi", "advi-cross", "pathfinder")) {
        # ADVI / Pathfinder warm-start: fixed CSP init here; inference runs after
        # phi/sphi init is fully resolved (see warm-start blocks below).
        cat(sprintf("[init] mode = %s (CSP init = fixed; %s will run before sampling)\n",
                    init_mode,
                    if (init_mode == "pathfinder") "Pathfinder" else "ADVI"))
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
    } else {
        stop("init.mode must be one of: fixed, rmf-posterior, scuo, enc_prime, ",
             "mixed_scuo_encp, warm-start, advi, advi-cross, pathfinder; got: ", init_mode)
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

    # AA-NSE collapse: log_NSERate[C] -> log_NSERate_aa[N_AA] via aa_of_codon
    if (aa_nse) {
        cat("[init] AA-NSE model: grouping log_NSERate by AA family (N_AA =",
            n_aa_groups, ")\n")
        inner_init_fn <- init_fn
        aoc <- aa_of_codon_vec   # close over local copy
        naa <- n_aa_groups
        init_fn <- function(chain_id = 1) {
            ll <- inner_init_fn(chain_id)
            ll$log_NSERate_aa <- vapply(seq_len(naa), function(i)
                mean(ll$log_NSERate[aoc == i]), 0.0)
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

    # ---- ADVI warm-start (runs after phi/sphi init is fully resolved) --------
    if (init_mode %in% c("advi", "advi-cross")) {
        advi_n_chains  <- chains  %||% (config$stan %||% list())$chains            %||% 4L
        advi_n_threads <- threads %||% (config$stan %||% list())$threads_per_chain %||% 1L
        advi_threads   <- advi_n_chains * advi_n_threads
        advi_seed      <- seed %||% (config$stan %||% list())$seed %||% 20260523L
        advi_metric    <- (config$stan %||% list())$metric %||% "diag_e"

        # init.advi.source: which model to run ADVI on (may differ from HMC target).
        # advi-cross (legacy) == init.mode:advi + init.advi.source:sharednse on per-codon.
        # NULL means ADVI on the same model as HMC.
        advi_source <- config$fit$init.advi.source %||%
            if (init_mode == "advi-cross") "sharednse" else NULL

        if (!is.null(advi_source)) {
            # Cross-model ADVI: compile source model, strip/augment stan_data as needed.
            advi_src_stan <- switch(advi_source,
                "sharednse" = "panse_sphi_est_noncentered_sharednse.stan",
                "aa-nse"    = "panse_sphi_est_noncentered_aanse.stan",
                stop("Unknown init.advi.source: ", advi_source)
            )
            advi_model_file <- .panse_stan_file(advi_src_stan)
            advi_build_dir  <- tools::R_user_dir("AnaCoDa", "cache")
            if (!dir.exists(advi_build_dir)) dir.create(advi_build_dir, recursive = TRUE)
            cat("[advi-cross] compiling source model (", advi_source,
                ") for ADVI...\n", sep = "")
            mod_vi <- cmdstanr::cmdstan_model(
                advi_model_file, dir = advi_build_dir,
                cpp_options = list(stan_threads = TRUE)
            )
            # Build advi_stan_data: strip or add AA fields to match source model.
            advi_stan_data <- stan_data
            if (advi_source == "sharednse") {
                advi_stan_data$N_AA        <- NULL   # sharednse has no N_AA/aa_of_codon
                advi_stan_data$aa_of_codon <- NULL
            } else if (advi_source == "aa-nse") {
                if (!aa_nse) {
                    # Target is per-codon; build AA mapping on the fly for the source model.
                    if (is.null(aa_of_codon_vec)) {
                        aa_map_loc  <- .build_nse_aa_map()
                        aa_map_loc  <- aa_map_loc[sort(names(aa_map_loc))]
                        aa_names_loc <- names(aa_map_loc)
                        c2a          <- setNames(
                            rep(seq_along(aa_names_loc), lengths(aa_map_loc)),
                            unlist(aa_map_loc, use.names = FALSE)
                        )
                        aa_of_codon_vec <- unname(c2a[attr(stan_data, "codon_order")])
                        n_aa_groups     <- length(aa_names_loc)
                    }
                    advi_stan_data$N_AA        <- n_aa_groups
                    advi_stan_data$aa_of_codon <- aa_of_codon_vec
                }
                # If target is aanse, N_AA/aa_of_codon already present in stan_data.
            }
            raw <- init_fn(chain_id = 1)
            if (advi_source == "sharednse") {
                nse_mean <- if (!is.null(raw$log_NSERate_aa)) mean(raw$log_NSERate_aa)
                            else mean(raw$log_NSERate)
                advi_start <- list(
                    log_alpha          = raw$log_alpha,
                    log_lambdaPrime    = raw$log_lambdaPrime,
                    log_NSERate_shared = nse_mean,
                    sphi               = if (!is.null(raw$sphi)) raw$sphi else 1.0,
                    z_phi              = if (!is.null(raw$z_phi)) raw$z_phi
                                         else rep(0.0, stan_data$G)
                )
            } else {   # aa-nse source
                nse_aa <- if (!is.null(raw$log_NSERate_aa)) raw$log_NSERate_aa
                          else vapply(seq_len(n_aa_groups), function(i)
                                  mean(raw$log_NSERate[aa_of_codon_vec == i]), 0.0)
                advi_start <- list(
                    log_alpha       = raw$log_alpha,
                    log_lambdaPrime = raw$log_lambdaPrime,
                    log_NSERate_aa  = nse_aa,
                    sphi            = if (!is.null(raw$sphi)) raw$sphi else 1.0,
                    z_phi           = if (!is.null(raw$z_phi)) raw$z_phi
                                      else rep(0.0, stan_data$G)
                )
            }
        } else {
            mod_vi         <- mod
            advi_stan_data <- stan_data
            advi_start     <- init_fn(chain_id = 1)
        }

        cat(sprintf("[advi] running variational inference (threads=%d, seed=%d)...\n",
                    advi_threads, advi_seed))
        t0_vi <- proc.time()[["elapsed"]]
        fit_vi <- tryCatch(
            mod_vi$variational(data    = advi_stan_data,
                               init    = list(advi_start),
                               threads = advi_threads,
                               seed    = advi_seed,
                               refresh = 200L),
            error = function(e) {
                cat("[warn] ADVI failed: ", conditionMessage(e), "\n", sep = "")
                NULL
            }
        )
        t_vi <- proc.time()[["elapsed"]] - t0_vi
        cat(sprintf("[advi] wall = %.1f sec\n", t_vi))

        if (!is.null(fit_vi)) {
            advi_src_type <- if (is.null(advi_source)) {
                if (shared_nse) "sharednse" else if (aa_nse) "aa-nse" else "percodon"
            } else advi_source
            ws <- .panse_advi_warm_start(fit_vi, advi_stan_data, advi_src_type)

            # Expand metric from source model dims to HMC target model dims
            if (!is.null(advi_source) && advi_source == "sharednse" && !shared_nse) {
                if (aa_nse) {
                    ws <- .expand_shared_to_aa_metric(ws, stan_data)
                } else {
                    ws <- .expand_cross_advi_metric(ws, stan_data)  # shared -> per-codon
                }
            } else if (!is.null(advi_source) && advi_source == "aa-nse" && !aa_nse) {
                ws <- .expand_aa_to_percodon_metric(ws, advi_stan_data)
            }
            inv_m <- ws$inv_metric
            inv_m[!is.finite(inv_m) | inv_m <= 0] <- 1.0
            cat(sprintf("[advi] inv_metric diagonal: n=%d  range=[%.3e, %.3e]\n",
                        length(inv_m), min(inv_m), max(inv_m)))

            advi_means  <- ws$init
            init_fn     <- function(chain_id = 1) advi_means
            n_rep       <- advi_n_chains
            warm_start_metrics <- lapply(seq_len(n_rep), function(i) {
                if (advi_metric == "dense_e") diag(inv_m) else inv_m
            })
        } else {
            cat("[warn] ADVI failed; sampling will use fixed init with default metric\n")
        }
    }

    # ---- Pathfinder warm-start (same model as HMC) ---------------------------
    if (init_mode == "pathfinder") {
        pf_cfg        <- config$stan %||% list()
        pf_n_chains   <- chains   %||% pf_cfg$chains             %||% 4L
        pf_n_threads  <- threads  %||% pf_cfg$threads_per_chain  %||% 1L
        pf_num_paths  <- pf_cfg$pathfinder_num_paths       %||% 4L
        pf_max_iter   <- pf_cfg$pathfinder_max_lbfgs_iters %||% 100L
        pf_seed       <- seed %||% pf_cfg$seed %||% 20260523L
        pf_metric     <- pf_cfg$metric %||% "diag_e"
        # Run paths in parallel: one thread per path (up to available n_threads)
        pf_threads    <- pf_n_chains * pf_n_threads

        cat(sprintf("[pathfinder] num_paths=%d max_lbfgs_iters=%d threads=%d seed=%d...\n",
                    pf_num_paths, pf_max_iter, pf_threads, pf_seed))
        t0_pf <- proc.time()[["elapsed"]]
        pf_init <- lapply(seq_len(pf_num_paths), function(i) init_fn(chain_id = i))
        fit_pf <- tryCatch(
            mod$pathfinder(
                data            = stan_data,
                init            = pf_init,
                num_paths       = as.integer(pf_num_paths),
                max_lbfgs_iters = as.integer(pf_max_iter),
                psis_resample   = TRUE,
                calculate_lp    = TRUE,
                num_threads     = as.integer(pf_threads),
                seed            = as.integer(pf_seed),
                refresh         = 100L
            ),
            error = function(e) {
                cat("[warn] Pathfinder failed: ", conditionMessage(e), "\n", sep = "")
                NULL
            }
        )
        t_pf <- proc.time()[["elapsed"]] - t0_pf
        cat(sprintf("[pathfinder] wall = %.1f sec\n", t_pf))

        if (!is.null(fit_pf)) {
            pf_src_type <- if (shared_nse) "sharednse" else if (aa_nse) "aa-nse" else "percodon"
            ws <- .panse_pathfinder_warm_start(fit_pf, stan_data, pf_src_type)

            if (isTRUE(ws$is_dense)) {
                cd <- diag(ws$inv_metric)
                cat(sprintf("[pathfinder] inv_metric (dense %dx%d): diag range=[%.3e, %.3e]\n",
                            nrow(ws$inv_metric), ncol(ws$inv_metric), min(cd), max(cd)))
                warm_start_metrics <- lapply(seq_len(pf_n_chains), function(i) ws$inv_metric)
            } else {
                inv_m <- ws$inv_metric
                cat(sprintf("[pathfinder] inv_metric (diag fallback): n=%d  range=[%.3e, %.3e]\n",
                            length(inv_m), min(inv_m), max(inv_m)))
                warm_start_metrics <- lapply(seq_len(pf_n_chains), function(i) {
                    if (pf_metric == "dense_e") diag(inv_m) else inv_m
                })
            }
            pf_means <- ws$init
            init_fn  <- function(chain_id = 1) pf_means
        } else {
            cat("[warn] Pathfinder failed; sampling will use fixed init with default metric\n")
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


# ---- Internal helpers -------------------------------------------------------

# Extract init point and diagonal inv_metric from a PANSE ADVI fit.
#
# Parameters are box-constrained on the log scale (log_alpha, log_lambdaPrime,
# log_NSERate).  Their unconstrained Stan transforms are logit-scaled:
#   y = log((x - lower) / (upper - x))
# Variances are computed in unconstrained space.
#
# sphi (lower=0) uses var(log(sphi)) as the unconstrained variance proxy.
# z_phi is unconstrained; raw variance used directly.
#
# Stan parameter declaration order:
#   log_alpha[1..C], log_lambdaPrime[1..C],
#   log_NSERate[1..C] or log_NSERate_shared(1) or log_NSERate_aa[1..N_AA],
#   sphi, z_phi[1..G]
#
# src_type: "sharednse" | "aa-nse" | "percodon"
.panse_advi_warm_start <- function(fit_vi, stan_data, src_type) {
    draws <- fit_vi$draws(format = "df")
    C <- stan_data$C
    G <- stan_data$G

    a_cols <- paste0("log_alpha[",       seq_len(C), "]")
    l_cols <- paste0("log_lambdaPrime[", seq_len(C), "]")
    z_cols <- paste0("z_phi[",           seq_len(G), "]")

    # Posterior means in constrained space (used as HMC init)
    init <- list(
        log_alpha       = as.numeric(colMeans(draws[, a_cols, drop = FALSE])),
        log_lambdaPrime = as.numeric(colMeans(draws[, l_cols, drop = FALSE])),
        sphi            = mean(draws[["sphi"]]),
        z_phi           = as.numeric(colMeans(draws[, z_cols, drop = FALSE]))
    )
    if (src_type == "sharednse") {
        init$log_NSERate_shared <- mean(draws[["log_NSERate_shared"]])
    } else if (src_type == "aa-nse") {
        N_AA <- stan_data$N_AA
        nse_aa_cols <- paste0("log_NSERate_aa[", seq_len(N_AA), "]")
        init$log_NSERate_aa <- as.numeric(colMeans(draws[, nse_aa_cols, drop = FALSE]))
    } else {   # percodon
        nse_cols <- paste0("log_NSERate[", seq_len(C), "]")
        init$log_NSERate <- as.numeric(colMeans(draws[, nse_cols, drop = FALSE]))
    }

    # Unconstrained variances for diagonal inv_metric
    a_lo <- stan_data$log_alpha_lower;  a_hi <- stan_data$log_alpha_upper
    l_lo <- stan_data$log_lambda_lower; l_hi <- stan_data$log_lambda_upper
    n_lo <- stan_data$log_nse_lower;    n_hi <- stan_data$log_nse_upper

    .logit_unc_var <- function(x, lower, upper) {
        y <- log((x - lower) / (upper - x))
        v <- var(y[is.finite(y)])
        if (!is.finite(v) || v <= 0) 1.0 else v
    }

    a_var <- vapply(a_cols, function(col)
                        .logit_unc_var(draws[[col]], a_lo, a_hi), 0.0)
    l_var <- vapply(l_cols, function(col)
                        .logit_unc_var(draws[[col]], l_lo, l_hi), 0.0)

    if (src_type == "sharednse") {
        nse_var <- .logit_unc_var(draws[["log_NSERate_shared"]], n_lo, n_hi)
    } else if (src_type == "aa-nse") {
        N_AA <- stan_data$N_AA
        nse_aa_cols <- paste0("log_NSERate_aa[", seq_len(N_AA), "]")
        nse_var  <- vapply(nse_aa_cols, function(col)
                               .logit_unc_var(draws[[col]], n_lo, n_hi), 0.0)
    } else {
        nse_cols <- paste0("log_NSERate[", seq_len(C), "]")
        nse_var  <- vapply(nse_cols, function(col)
                               .logit_unc_var(draws[[col]], n_lo, n_hi), 0.0)
    }

    sphi_draws <- draws[["sphi"]]
    sphi_var   <- var(log(sphi_draws[sphi_draws > 0]))
    if (!is.finite(sphi_var) || sphi_var <= 0) sphi_var <- 1.0

    z_var <- as.numeric(apply(draws[, z_cols, drop = FALSE], 2L, var))
    z_var[!is.finite(z_var) | z_var <= 0] <- 1.0
    # ADVI mean-field underestimates z_phi variance (noncentered prior is N(0,1));
    # floor at 1.0 so the seeded inv_metric does not produce pathologically tiny steps.
    z_var <- pmax(z_var, 1.0)

    inv_metric <- c(a_var, l_var, nse_var, sphi_var, z_var)

    list(init = init, inv_metric = inv_metric)
}


# Extract warm-start from a Pathfinder fit.
#
# Unlike ADVI (mean-field, diagonal), Pathfinder draws capture correlations via
# the L-BFGS Hessian approximation.  Returns a full covariance matrix in
# unconstrained space as inv_metric, enabling dense_e adaptation from iter 1.
# Falls back to diagonal variance if fewer than n_params+5 finite draws survive.
#
# Returns: list(init, inv_metric, is_dense)
#   init:       constrained-space parameter means (HMC starting point)
#   inv_metric: covariance matrix (is_dense=TRUE) or variance vector (is_dense=FALSE)
#   is_dense:   logical; TRUE when inv_metric is a matrix
.panse_pathfinder_warm_start <- function(fit_pf, stan_data, src_type) {
    draws <- fit_pf$draws(format = "df")
    C <- stan_data$C
    G <- stan_data$G

    a_cols <- paste0("log_alpha[",       seq_len(C), "]")
    l_cols <- paste0("log_lambdaPrime[", seq_len(C), "]")
    z_cols <- paste0("z_phi[",           seq_len(G), "]")

    # Constrained-space means -> HMC init (same logic as .panse_advi_warm_start)
    init <- list(
        log_alpha       = as.numeric(colMeans(draws[, a_cols, drop = FALSE])),
        log_lambdaPrime = as.numeric(colMeans(draws[, l_cols, drop = FALSE])),
        sphi            = mean(draws[["sphi"]]),
        z_phi           = as.numeric(colMeans(draws[, z_cols, drop = FALSE]))
    )
    if (src_type == "sharednse") {
        init$log_NSERate_shared <- mean(draws[["log_NSERate_shared"]])
    } else if (src_type == "aa-nse") {
        N_AA <- stan_data$N_AA
        nse_aa_cols <- paste0("log_NSERate_aa[", seq_len(N_AA), "]")
        init$log_NSERate_aa <- as.numeric(colMeans(draws[, nse_aa_cols, drop = FALSE]))
    } else {
        nse_cols <- paste0("log_NSERate[", seq_len(C), "]")
        init$log_NSERate <- as.numeric(colMeans(draws[, nse_cols, drop = FALSE]))
    }

    # Transform constrained draws to unconstrained space
    a_lo <- stan_data$log_alpha_lower;  a_hi <- stan_data$log_alpha_upper
    l_lo <- stan_data$log_lambda_lower; l_hi <- stan_data$log_lambda_upper
    n_lo <- stan_data$log_nse_lower;    n_hi <- stan_data$log_nse_upper

    # Box-constrained [lo, hi] -> unconstrained via logit transform
    .logit_unc <- function(x, lo, hi) log((x - lo) / (hi - x))

    a_unc <- .logit_unc(as.matrix(draws[, a_cols, drop = FALSE]), a_lo, a_hi)
    l_unc <- .logit_unc(as.matrix(draws[, l_cols, drop = FALSE]), l_lo, l_hi)

    if (src_type == "sharednse") {
        nse_unc <- matrix(.logit_unc(as.numeric(draws[["log_NSERate_shared"]]),
                                     n_lo, n_hi), ncol = 1L)
    } else if (src_type == "aa-nse") {
        N_AA <- stan_data$N_AA
        nse_aa_cols <- paste0("log_NSERate_aa[", seq_len(N_AA), "]")
        nse_unc <- .logit_unc(as.matrix(draws[, nse_aa_cols, drop = FALSE]), n_lo, n_hi)
    } else {
        nse_cols <- paste0("log_NSERate[", seq_len(C), "]")
        nse_unc <- .logit_unc(as.matrix(draws[, nse_cols, drop = FALSE]), n_lo, n_hi)
    }

    # Lower-bounded [0, inf) -> unconstrained via log
    sphi_unc <- matrix(log(as.numeric(draws[["sphi"]])), ncol = 1L)
    # z_phi is unconstrained (identity)
    z_unc <- as.matrix(draws[, z_cols, drop = FALSE])

    unc_mat <- cbind(a_unc, l_unc, nse_unc, sphi_unc, z_unc)
    finite_mask <- rowSums(!is.finite(unc_mat)) == 0L
    unc_mat <- unc_mat[finite_mask, , drop = FALSE]
    n_draws  <- nrow(unc_mat)
    n_params <- ncol(unc_mat)

    if (n_draws < n_params + 5L) {
        warning(".panse_pathfinder_warm_start: only ", n_draws, " finite draws for ",
                n_params, " params; falling back to diagonal inv_metric")
        dv <- apply(unc_mat, 2L, var)
        dv[!is.finite(dv) | dv <= 0] <- 1.0
        z_idx <- seq.int(n_params - G + 1L, n_params)
        dv[z_idx] <- pmax(dv[z_idx], 1.0)
        return(list(init = init, inv_metric = dv, is_dense = FALSE))
    }

    # Full sample covariance matrix = optimal inv_metric for dense_e
    cov_mat <- cov(unc_mat)

    # Floor z_phi diagonal at 1.0 (noncentered N(0,1) prior; same as ADVI fix)
    z_idx <- seq.int(n_params - G + 1L, n_params)
    diag(cov_mat)[z_idx] <- pmax(diag(cov_mat)[z_idx], 1.0)

    # Ridge regularization for numerical stability (1e-6 * mean diagonal)
    diag(cov_mat) <- diag(cov_mat) + 1e-6 * mean(diag(cov_mat))

    list(init = init, inv_metric = cov_mat, is_dense = TRUE)
}


# Expand a shared-NSE warm-start (424 dims) to per-codon NSE layout (484 dims).
#
# Shared-NSE inv_metric order: a[C], l[C], nse_shared(1), sphi, z[G]
# Per-codon   inv_metric order: a[C], l[C], nse[C],        sphi, z[G]
#
# Transformation: replicate the single nse_shared variance C times.
# Also updates the init list: log_NSERate_shared -> log_NSERate[C].
.expand_cross_advi_metric <- function(ws_shared, stan_data) {
    C     <- stan_data$C
    inv_m <- ws_shared$inv_metric

    # Locate the nse_shared scalar at index 2C+1
    idx_nse   <- 2L * C + 1L
    nse_var   <- inv_m[idx_nse]
    after_nse <- inv_m[(idx_nse + 1L):length(inv_m)]   # sphi_var + z_var[G]

    inv_metric_expanded <- c(
        inv_m[seq_len(2L * C)],   # a_var[C] + l_var[C]
        rep(nse_var, C),           # nse_var replicated to C codons
        after_nse                   # sphi_var, z_var[G]
    )

    # Expand init: log_NSERate_shared scalar -> log_NSERate[C] vector
    init <- ws_shared$init
    nse_mean               <- init$log_NSERate_shared
    init$log_NSERate        <- rep(nse_mean, C)
    init$log_NSERate_shared <- NULL

    list(init = init, inv_metric = inv_metric_expanded)
}


# Expand a shared-NSE warm-start (2C+2+G dims) to AA-NSE layout (2C+N_AA+1+G dims).
#
# Shared-NSE inv_metric order: a[C], l[C], nse_shared(1), sphi, z[G]
# AA-NSE      inv_metric order: a[C], l[C], nse_aa[N_AA],  sphi, z[G]
#
# Transformation: replicate the single nse_shared variance N_AA times.
# Also updates init: log_NSERate_shared -> log_NSERate_aa[N_AA].
.expand_shared_to_aa_metric <- function(ws_shared, stan_data) {
    C    <- stan_data$C
    N_AA <- stan_data$N_AA
    inv_m <- ws_shared$inv_metric

    idx_nse   <- 2L * C + 1L
    nse_var   <- inv_m[idx_nse]
    after_nse <- inv_m[(idx_nse + 1L):length(inv_m)]   # sphi_var + z_var[G]

    inv_metric_expanded <- c(
        inv_m[seq_len(2L * C)],
        rep(nse_var, N_AA),
        after_nse
    )

    init <- ws_shared$init
    nse_mean              <- init$log_NSERate_shared
    init$log_NSERate_aa   <- rep(nse_mean, N_AA)
    init$log_NSERate_shared <- NULL

    list(init = init, inv_metric = inv_metric_expanded)
}


# Expand an AA-NSE warm-start (2C+N_AA+1+G dims) to per-codon NSE (2C+C+1+G dims).
#
# AA-NSE    inv_metric order: a[C], l[C], nse_aa[N_AA], sphi, z[G]
# Per-codon inv_metric order: a[C], l[C], nse[C],       sphi, z[G]
#
# Transformation: replicate each AA-group variance to all codons in that group.
# Also updates init: log_NSERate_aa[N_AA] -> log_NSERate[C].
.expand_aa_to_percodon_metric <- function(ws_aa, stan_data) {
    C           <- stan_data$C
    N_AA        <- stan_data$N_AA
    aa_of_codon <- stan_data$aa_of_codon
    inv_m <- ws_aa$inv_metric

    idx_nse   <- 2L * C + 1L
    nse_vars  <- inv_m[idx_nse:(idx_nse + N_AA - 1L)]   # N_AA values
    after_nse <- inv_m[(idx_nse + N_AA):length(inv_m)]   # sphi_var + z_var[G]

    # Expand N_AA variances to C codons via aa_of_codon index
    nse_var_percodon <- nse_vars[aa_of_codon]

    inv_metric_expanded <- c(
        inv_m[seq_len(2L * C)],
        nse_var_percodon,
        after_nse
    )

    init <- ws_aa$init
    init$log_NSERate    <- init$log_NSERate_aa[aa_of_codon]
    init$log_NSERate_aa <- NULL

    list(init = init, inv_metric = inv_metric_expanded)
}


# Extended AA-to-codon map for NSE grouping.
#
# build_aa_codon_map() excludes Met (M) and Trp (W) because they have d=1
# (no synonymous codons, no codon-bias selection).  For NSE rates, every codon
# still has a rate, so we add singleton groups for M (ATG) and W (TGG).
# This gives N_AA = 21 groups (19 bias + 2 singletons) in alphabetical order:
# A C D E F G H I K L M N P Q R S T V W Y Z.
.build_nse_aa_map <- function() {
    aa_map        <- build_aa_codon_map()   # 19 groups; M and W excluded
    aa_map[["M"]] <- "ATG"                  # Met singleton
    aa_map[["W"]] <- "TGG"                  # Trp singleton
    aa_map
}
