# ============================================================================
# helper-panse-stan-harness.R -- helpers for testing the MIGRATED PANSE-Stan
# R harness (data builder, simulator, shared init/codon-bias utils).  Kept
# separate from helper-panse-stan.R (the .stan MODEL tests) to avoid colliding
# with concurrent edits there.  These tests are pure R -- no cmdstan/AnaCoDa.
# ============================================================================

# Locate repo root by walking up to scripts/lib/init_helpers.R (self-contained;
# does not depend on the model-test helper).
panse_harness_root <- function() {
    d <- normalizePath(testthat::test_path("."), mustWork = FALSE)
    for (i in 1:6) {
        if (file.exists(file.path(d, "scripts", "lib", "init_helpers.R")))
            return(normalizePath(d))
        d <- dirname(d)
    }
    testthat::skip("cannot locate repo root (scripts/lib/init_helpers.R)")
}

# Source a migrated harness file by repo-relative path (skips if absent).
source_panse_harness <- function(relpath) {
    f <- file.path(panse_harness_root(), relpath)
    if (!file.exists(f)) testthat::skip(paste("migrated harness file not found:", relpath))
    sys.source(f, envir = globalenv())
}
