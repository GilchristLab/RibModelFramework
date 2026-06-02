#!/usr/bin/env Rscript
# ============================================================================
# simulate_panse_dataset.R -- thin CLI shim: parse args, call
#   AnaCoDa::simulate_panse_dataset() + AnaCoDa::write_sim_outputs()
#
# Usage:
#   Rscript simulate_panse_dataset.R <config.yaml> [--seed N] [--out DIR]
#                                    [--dry-run] [--lib-loc DIR] [-h]
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L || any(args %in% c("-h", "--help"))) {
    cat("Usage: Rscript simulate_panse_dataset.R <config.yaml> [--seed N] [--out DIR] [--dry-run]\n")
    quit(status = if (length(args) == 0L) 1L else 0L)
}

opts <- list(config = NULL, lib_loc = NULL, seed = NULL, out = NULL,
             dry_run = FALSE)
i <- 1L
while (i <= length(args)) {
    a <- args[[i]]
    if      (a == "--lib-loc")           { opts$lib_loc <- args[[i + 1L]];             i <- i + 2L }
    else if (a == "--seed")              { opts$seed    <- as.integer(args[[i + 1L]]); i <- i + 2L }
    else if (a == "--out")               { opts$out     <- args[[i + 1L]];             i <- i + 2L }
    else if (a == "--dry-run")           { opts$dry_run <- TRUE;                       i <- i + 1L }
    else if (startsWith(a, "--"))        { stop("unknown option: ", a, call. = FALSE) }
    else if (is.null(opts$config))       { opts$config  <- a;                          i <- i + 1L }
    else stop("unexpected positional arg: ", a)
}

if (!is.null(opts$lib_loc))
    .libPaths(c(opts$lib_loc, .libPaths()))

suppressPackageStartupMessages(library(AnaCoDa))

if (is.null(opts$config)) stop("Usage: simulate_panse_dataset.R <config.yaml> [...]")
cfg_abs <- normalizePath(opts$config)
if (!file.exists(cfg_abs)) stop("Config not found: ", cfg_abs)

config <- yaml::read_yaml(cfg_abs)
if (!is.null(opts$out)) config$sim$out.dir <- opts$out
out_dir <- config$sim$out.dir %||% stop("config$sim$out.dir missing")

if (opts$dry_run) {
    cat("[--dry-run] config:\n")
    str(config$sim)
    quit(save = "no", status = 0)
}

sim <- AnaCoDa::simulate_panse_dataset(config, seed = opts$seed)
AnaCoDa::write_sim_outputs(sim, out_dir, config)
cat("\n[ok] sim outputs -> ", normalizePath(out_dir), "\n", sep = "")
