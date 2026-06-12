library(testthat)
library(AnaCoDa)

context("ACF diagnostic functions: acfMCMC, acfCSP, helpers")

# ============================================================
# .acf.ess.per.1k  -- pure function, no MCMC needed
# ============================================================

# Build a minimal acf-class object with prescribed lag-1..n autocorrelations.
.fake_acf <- function(rhos) {
    n <- length(rhos)
    structure(
        list(acf = array(c(1.0, rhos), dim = c(n + 1L, 1L, 1L)),
             lag = array(0L:n,         dim = c(n + 1L, 1L, 1L)),
             n.used = 100L),
        class = "acf"
    )
}

test_that(".acf.ess.per.1k is 1000 when all autocorrelations are zero", {
    a <- .fake_acf(c(0, 0, 0))
    expect_equal(.acf.ess.per.1k(a), 1000.0)
})

test_that(".acf.ess.per.1k is 1000 when lag-1 is immediately non-positive", {
    a <- .fake_acf(c(-0.3, 0.5, 0.5))
    expect_equal(.acf.ess.per.1k(a), 1000.0)
})

test_that(".acf.ess.per.1k correct for single positive lag", {
    # rho_1=0.5, then negative -> sum = 0.5, ESS/1k = 1000/2 = 500
    a <- .fake_acf(c(0.5, -0.1))
    expect_equal(.acf.ess.per.1k(a), 500.0)
})

test_that(".acf.ess.per.1k sums all positive lags before first non-positive", {
    # rho = 0.4, 0.2, -0.1, 0.3  -> sum positive prefix = 0.6
    # ESS/1k = round(1000 / 2.2, 1)
    a <- .fake_acf(c(0.4, 0.2, -0.1, 0.3))
    expect_equal(.acf.ess.per.1k(a),
                 round(1000.0 / (1 + 2 * 0.6), 1L))
})

test_that(".acf.ess.per.1k is independent of number of samples (chunk-size)", {
    # Same autocorrelation structure => same ESS/1k regardless of n.used
    a1 <- .fake_acf(c(0.3, 0.1, -0.2))
    a2 <- .fake_acf(c(0.3, 0.1, -0.2))
    a2$n.used <- 1000L
    expect_equal(.acf.ess.per.1k(a1),
                 .acf.ess.per.1k(a2))
})

# ============================================================
# .buildACFLayout  -- pure function, no MCMC needed
# ============================================================

test_that(".buildACFLayout returns list with mat, widths, heights", {
    lay <- .buildACFLayout(17L)
    expect_named(lay, c("mat", "widths", "heights"))
})

test_that(".buildACFLayout dimensions correct for n.slots=17", {
    lay <- .buildACFLayout(17L)
    n.rows.data <- ceiling(17L / 4L)  # 5
    expect_equal(nrow(lay$mat), n.rows.data + 2L)  # title row + data rows + x.lbl row
    expect_equal(ncol(lay$mat), 5L)                # y.lbl col + 4 data cols
    expect_equal(length(lay$widths),  5L)
    expect_equal(length(lay$heights), n.rows.data + 2L)
})

test_that(".buildACFLayout y.lbl fills entire left column", {
    for (n in c(1L, 4L, 5L, 17L)) {
        lay <- .buildACFLayout(n)
        expect_true(all(lay$mat[, 1L] == lay$mat[1L, 1L]),
                    info = paste("n.slots =", n))
    }
})

test_that(".buildACFLayout all data panels 2..(n+1) present exactly once in data block", {
    n <- 17L
    lay <- .buildACFLayout(n)
    nr  <- nrow(lay$mat)
    data.block <- lay$mat[2L:(nr - 1L), 2L:5L]
    present    <- sort(data.block[data.block > 0L])
    expect_equal(present, 2L:(n + 1L))
})

test_that(".buildACFLayout title (panel 1) fills non-y.lbl cells of first row", {
    lay <- .buildACFLayout(8L)
    expect_true(all(lay$mat[1L, 2L:5L] == 1L))
})

test_that(".buildACFLayout x.lbl fills non-y.lbl cells of last row", {
    lay <- .buildACFLayout(8L)
    nr    <- nrow(lay$mat)
    x.lbl <- 8L + 2L  # n.slots + 2
    expect_true(all(lay$mat[nr, 2L:5L] == x.lbl))
})

test_that(".buildACFLayout row count scales with n.slots", {
    for (n in c(1L, 4L, 5L, 8L, 17L)) {
        lay <- .buildACFLayout(n)
        expect_equal(nrow(lay$mat), ceiling(n / 4L) + 2L,
                     info = paste("n.slots =", n))
    }
})

test_that(".buildACFLayout n.slots=4 fits in one data row", {
    lay <- .buildACFLayout(4L)
    expect_equal(nrow(lay$mat), 3L)
})

test_that(".buildACFLayout n.slots=5 spans two data rows", {
    lay <- .buildACFLayout(5L)
    expect_equal(nrow(lay$mat), 4L)
})

# ============================================================
# acfMCMC and acfCSP -- require a live MCMC trace
# ============================================================

data.dir <- file.path("UnitTestingData", "testMCMCROCFiles")
fasta    <- file.path(data.dir, "simulatedAllUniqueR.fasta")
phi.csv  <- file.path(data.dir, "simulatedAllUniqueR_phi_withPhiSet.csv")
sel1     <- file.path(data.dir, "selection_1.csv")
sel2     <- file.path(data.dir, "selection_2.csv")
mut1     <- file.path(data.dir, "mutation_1.csv")
mut2     <- file.path(data.dir, "mutation_2.csv")

skip_if_not(file.exists(fasta), "ROC test FASTA not found -- skipping ACF integration tests")

set.seed(42L)
genome <- initializeGenomeObject(file = fasta,
                                 observed.expression.file = phi.csv,
                                 match.expression.by.id = FALSE)
gene.assignment <- c(rep(1L, 250L), rep(2L, length(genome) - 250L))
parameter <- initializeParameterObject(genome, c(1, 1), 2L, gene.assignment,
                                       split.serine = TRUE,
                                       mixture.definition = "allUnique")
parameter$initSelectionCategories(c(sel1, sel2), 2L, FALSE)
parameter$initMutationCategories(c(mut1, mut2), 2L, FALSE)
model <- initializeModelObject(parameter, "ROC", with.phi = TRUE)
mcmc  <- initializeMCMCObject(samples = 20L, thinning = 10L,
                              adaptive.width = 10L,
                              est.expression = TRUE, est.csp = TRUE,
                              est.hyper = TRUE)
runMCMC(mcmc = mcmc, genome = genome, model = model, ncores = 1L,
        divergence.iteration = 0L)

# ---- acfMCMC return structure ----

test_that("acfMCMC returns named list with both traces by default", {
    result <- acfMCMC(mcmc, plot = FALSE)
    expect_type(result, "list")
    expect_named(result, c("LogPosterior", "LogLikelihood"))
})

test_that("acfMCMC each element is an acf-class object", {
    result <- acfMCMC(mcmc, plot = FALSE)
    expect_s3_class(result$LogPosterior,   "acf")
    expect_s3_class(result$LogLikelihood,  "acf")
})

test_that("acfMCMC what='LogPosterior' returns single-element list", {
    result <- acfMCMC(mcmc, what = "LogPosterior", plot = FALSE)
    expect_length(result, 1L)
    expect_named(result, "LogPosterior")
})

test_that("acfMCMC what='LogLikelihood' returns single-element list", {
    result <- acfMCMC(mcmc, what = "LogLikelihood", plot = FALSE)
    expect_length(result, 1L)
    expect_named(result, "LogLikelihood")
})

test_that("acfMCMC lag.max is respected", {
    result <- acfMCMC(mcmc, samples = 20L, lag.max = 10L, plot = FALSE)
    expect_equal(max(as.integer(result$LogPosterior$lag)), 10L)
})

test_that("acfMCMC lag-0 autocorrelation is 1", {
    result <- acfMCMC(mcmc, plot = FALSE)
    expect_equal(as.numeric(result$LogPosterior$acf)[1L], 1.0)
})

test_that("acfMCMC plots without error (both traces)", {
    pdf(NULL)
    on.exit(dev.off())
    expect_no_error(acfMCMC(mcmc, lag.max = 10L))
})

test_that("acfMCMC plots without error (single trace)", {
    pdf(NULL)
    on.exit(dev.off())
    expect_no_error(acfMCMC(mcmc, what = "LogPosterior", lag.max = 10L))
})

# ---- acfCSP return structure ----

test_that("acfCSP returns named list keyed by valid amino acids", {
    result <- acfCSP(parameter, what = "Mutation", mixture = 1L, plot = FALSE)
    expect_type(result, "list")
    valid <- Filter(function(aa) !(aa %in% c("M", "W", "X")) &&
                        length(AAToCodon(aa, TRUE)) > 0L,
                    aminoAcids())
    expect_setequal(names(result), valid)
})

test_that("acfCSP each AA element is a list of acf objects keyed by codon", {
    result <- acfCSP(parameter, what = "Mutation", mixture = 1L, plot = FALSE)
    aa.acfs <- result[["A"]]
    expect_type(aa.acfs, "list")
    expect_true(all(vapply(aa.acfs, inherits, logical(1L), "acf")))
    expect_setequal(names(aa.acfs), AAToCodon("A", TRUE))
})

test_that("acfCSP lag-0 autocorrelation is 1 for all codons", {
    result <- acfCSP(parameter, what = "Mutation", mixture = 1L, plot = FALSE)
    lag0.vals <- unlist(lapply(result, function(aa.acfs)
        lapply(aa.acfs, function(a) as.numeric(a$acf)[1L])))
    expect_true(all(abs(lag0.vals - 1.0) < 1e-10))
})

test_that("acfCSP what='Selection' works and returns correct codons", {
    result <- acfCSP(parameter, what = "Selection", mixture = 1L, plot = FALSE)
    expect_type(result, "list")
    expect_setequal(names(result[["A"]]), AAToCodon("A", TRUE))
})

test_that("acfCSP lag.max is respected", {
    result <- acfCSP(parameter, what = "Mutation", mixture = 1L,
                     lag.max = 10L, plot = FALSE)
    aa    <- names(result)[1L]
    codon <- names(result[[aa]])[1L]
    expect_equal(max(as.integer(result[[aa]][[codon]]$lag)), 10L)
})

test_that("acfCSP scales='fixed' plots without error", {
    pdf(NULL)
    on.exit(dev.off())
    expect_no_error(acfCSP(parameter, what = "Mutation", mixture = 1L,
                           lag.max = 10L, scales = "fixed"))
})

test_that("acfCSP scales='free_y' plots without error", {
    pdf(NULL)
    on.exit(dev.off())
    expect_no_error(acfCSP(parameter, what = "Mutation", mixture = 1L,
                           lag.max = 10L, scales = "free_y"))
})

test_that("acfCSP mixture=2 works", {
    pdf(NULL)
    on.exit(dev.off())
    expect_no_error(acfCSP(parameter, what = "Mutation", mixture = 2L,
                           lag.max = 10L, scales = "fixed"))
})

test_that("acfCSP invalid what= raises error", {
    expect_error(acfCSP(parameter, what = "Bogus", plot = FALSE))
})
