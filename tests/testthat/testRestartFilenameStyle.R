library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Regression tests for setRestartFileSettings' leaf-aware filename
# handling.  Pins the two branches of the filename-format contract:
#
#   A. Caller's filename leaf has an extension (".rst" or any "*.X"):
#      strip-and-canonicalise to ".rst", then per-sample writes append
#      "_<sample>" AFTER the extension ("<base>.rst_<sample>"; "_final"
#      at end of run).  This is the 2016-era convention used by Loki
#      and other downstream tooling.
#
#   B. Caller's filename leaf has NO extension (e.g. bare prefix
#      "checkpoint", or any path ending in a directory component):
#      use the filename verbatim as the base, then per-sample writes
#      insert "_<sample>" BEFORE ".rst" ("<base>_<sample>.rst";
#      "_final.rst" at end of run).
#
# Also serves as a regression test against the path-stripping bug:
# absolute paths whose directory components contain dots (e.g.
# "/.../s.cerevisiae/.../<bare-leaf>") used to fool find_last_of(".")
# into stripping at the directory dot, mis-routing all output to the
# wrong directory.  The leaf-aware logic must not strip in that case.
# ======================================================================

genome_file <- system.file("extdata", "genome.fasta", package = "AnaCoDa")
genome      <- initializeGenomeObject(file = genome_file)
gene_assign <- rep(1, length(genome))

# Helper: run a tiny MCMC with the given restart filename and return
# the list of files produced in dirname(prefix).
runTinyMCMC <- function(prefix, write.multiple = TRUE, samples = 4,
                        interval.samples = 2) {
    set.seed(42)
    parameter <- initializeParameterObject(genome        = genome,
                                            sphi          = 1,
                                            num.mixtures  = 1,
                                            gene.assignment = gene_assign,
                                            model         = "ROC")
    mcmc <- initializeMCMCObject(samples = samples, thinning = 1,
                                  adaptive.width = 2,
                                  est.expression = TRUE,
                                  est.csp        = TRUE,
                                  est.hyper      = TRUE)
    setRestartSettings(mcmc, filename = prefix,
                       samples = interval.samples,
                       write.multiple = write.multiple)
    model <- initializeModelObject(parameter, "ROC", with.phi = FALSE)
    sink(tempfile())
    runMCMC(mcmc = mcmc, genome = genome, model = model,
            ncores = 1, divergence.iteration = 0)
    sink()
}

test_that("filename leaf with extension -> <base>.rst_<sample> + _final (legacy)", {
    d <- tempfile(pattern = "rstFilenameStyle_ext_")
    dir.create(d)
    prefix <- file.path(d, "myrun.rst")
    runTinyMCMC(prefix, write.multiple = TRUE)

    produced <- list.files(d)
    expect_true("myrun.rst_2"     %in% produced)
    expect_true("myrun.rst_4"     %in% produced)
    expect_true("myrun.rst_final" %in% produced)
    # And nothing of the "new" form
    expect_false(any(grepl("^myrun_[0-9]+\\.rst$", produced)))
})

test_that("filename leaf without extension -> <base>_<sample>.rst + _final.rst (new)", {
    d <- tempfile(pattern = "rstFilenameStyle_bare_")
    dir.create(d)
    prefix <- file.path(d, "myrun")  # bare leaf, no extension
    runTinyMCMC(prefix, write.multiple = TRUE)

    produced <- list.files(d)
    expect_true("myrun_2.rst"     %in% produced)
    expect_true("myrun_4.rst"     %in% produced)
    expect_true("myrun_final.rst" %in% produced)
    # And nothing of the legacy "_<n>" / "_final" (without trailing .rst)
    expect_false(any(grepl("^myrun\\.rst_", produced)))
    expect_false(any(grepl("^myrun_(final|[0-9]+)$", produced)))
})

test_that("path with dotted directory components is not mis-stripped (regression)", {
    # Construct a directory tree that mimics the production layout where
    # the parent contains a dot (e.g. "s.cerevisiae") but the LEAF has no
    # extension. Pre-fix this caused find_last_of(".") to strip inside
    # the parent name and write files into the grandparent.
    base   <- tempfile(pattern = "rstFilenameStyle_dottedParent_")
    dir.create(base)
    dotted <- file.path(base, "s.cerevisiae")
    leaf   <- file.path(dotted, "round-1")
    dir.create(leaf, recursive = TRUE)
    prefix <- file.path(leaf, "checkpoint")   # bare leaf, dotted ancestor

    runTinyMCMC(prefix, write.multiple = TRUE)

    # All output must land in `leaf`, none in `dotted` or `base`.
    expect_true(file.exists(file.path(leaf, "checkpoint_2.rst")))
    expect_true(file.exists(file.path(leaf, "checkpoint_4.rst")))
    expect_true(file.exists(file.path(leaf, "checkpoint_final.rst")))

    # No files should leak up to `dotted` or `base`.
    stray.dotted <- list.files(dotted, pattern = "^s\\.rst",   full.names = FALSE)
    stray.base   <- list.files(base,   pattern = "^s\\.rst",   full.names = FALSE)
    expect_length(stray.dotted, 0)
    expect_length(stray.base,   0)
})

test_that("write.multiple = FALSE writes a single rolling file regardless of style", {
    d <- tempfile(pattern = "rstFilenameStyle_single_")
    dir.create(d)
    prefix <- file.path(d, "rolling.rst")  # leaf with ext
    runTinyMCMC(prefix, write.multiple = FALSE)

    produced <- list.files(d)
    # Loki-style: single rolling file at the base prefix, plus the
    # end-of-run _final file.
    expect_true("rolling.rst"        %in% produced)
    expect_true("rolling.rst_final"  %in% produced)
    # No per-sample suffixed files for write.multiple = FALSE.
    expect_false(any(grepl("^rolling\\.rst_[0-9]+$", produced)))
})
