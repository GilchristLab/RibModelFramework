#!/usr/bin/env Rscript
# ============================================================================
# fit_panse_stan.R -- thin CLI shim: parse args, call AnaCoDa::fit_panse_stan()
#
# Usage:
#   Rscript fit_panse_stan.R runs/<config>.yaml [OPTIONS]
#
# Options:
#   --genes N         Subset to first N genes (for smoke tests).
#   --warmup N        Override stan.warmup from YAML.
#   --sampling N      Override stan.sampling from YAML.
#   --chains N        Override stan.chains from YAML.
#   --threads N       Override stan.threads_per_chain from YAML.
#   --adapt-delta X   Override stan.adapt_delta from YAML.
#   --out DIR         Override output directory.
#   --init-mode MODE  Init strategy (fixed, rmf-posterior, scuo, enc_prime,
#                     mixed_scuo_encp, warm-start).
#   --lib-loc DIR     Prepend DIR to R .libPaths().
#   --no-log-lik      Override emit_log_lik <- 0L.
#   --no-compile      Require precompiled exe; skip recompilation.
#   --dry-run         Build data + config; skip sampling.
#   -h, --help        Show this message.
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || any(args %in% c("-h", "--help"))) {
    cat("Usage: Rscript fit_panse_stan.R runs/<config>.yaml [OPTIONS]\n")
    cat("Run with --help for the full option list.\n")
    quit(status = if (length(args) == 0L) 1L else 0L)
}

opts <- list(
    config_path = NULL,
    lib_loc     = NULL,
    genes       = NULL,
    warmup      = NULL,
    sampling    = NULL,
    chains      = NULL,
    threads     = NULL,
    adapt_delta = NULL,
    out_dir     = NULL,
    init_mode   = NULL,
    no_log_lik  = FALSE,
    no_compile  = FALSE,
    dry_run     = FALSE
)
i <- 1L
while (i <= length(args)) {
    a <- args[[i]]
    if      (a == "--lib-loc")     { opts$lib_loc     <- args[[i + 1L]];             i <- i + 2L }
    else if (a == "--genes")       { opts$genes       <- as.integer(args[[i + 1L]]); i <- i + 2L }
    else if (a == "--warmup")      { opts$warmup      <- as.integer(args[[i + 1L]]); i <- i + 2L }
    else if (a == "--sampling")    { opts$sampling    <- as.integer(args[[i + 1L]]); i <- i + 2L }
    else if (a == "--chains")      { opts$chains      <- as.integer(args[[i + 1L]]); i <- i + 2L }
    else if (a == "--threads")     { opts$threads     <- as.integer(args[[i + 1L]]); i <- i + 2L }
    else if (a == "--adapt-delta") { opts$adapt_delta <- as.numeric(args[[i + 1L]]); i <- i + 2L }
    else if (a == "--out")         { opts$out_dir     <- args[[i + 1L]];             i <- i + 2L }
    else if (a == "--init-mode")   { opts$init_mode   <- args[[i + 1L]];             i <- i + 2L }
    else if (a == "--no-log-lik")  { opts$no_log_lik  <- TRUE;                       i <- i + 1L }
    else if (a == "--no-compile")  { opts$no_compile  <- TRUE;                       i <- i + 1L }
    else if (a == "--dry-run")     { opts$dry_run     <- TRUE;                       i <- i + 1L }
    else if (startsWith(a, "--"))  { stop("unknown option: ", a, call. = FALSE) }
    else                           { opts$config_path <- a;                          i <- i + 1L }
}

if (!is.null(opts$lib_loc))
    .libPaths(c(opts$lib_loc, .libPaths()))

suppressPackageStartupMessages(library(AnaCoDa))

if (is.null(opts$config_path))
    stop("missing positional arg: path to config.yaml")
cfg_abs <- normalizePath(opts$config_path)
if (!file.exists(cfg_abs)) stop("Config not found: ", cfg_abs)

config <- yaml::read_yaml(cfg_abs)

AnaCoDa::fit_panse_stan(
    config       = config,
    out_dir      = opts$out_dir,
    gene_subset  = opts$genes,
    init_mode    = opts$init_mode,
    no_log_lik   = opts$no_log_lik,
    no_compile   = opts$no_compile,
    dry_run      = opts$dry_run,
    chains       = opts$chains,
    warmup       = opts$warmup,
    sampling     = opts$sampling,
    threads      = opts$threads,
    adapt_delta  = opts$adapt_delta
)
