
library(testthat)
library(AnaCoDa)
rm(list = ls(all.names = TRUE))

# ======================================================================
# Tests for simulateGenome() -- task #6
#
# simulateGenome(genome) modifies the genome in-place: for ROC/FONSE it
# re-draws codon sequences from the current parameter distribution; for
# PA/PANSE it re-draws RFP counts. The method returns void (invisible NULL
# in R). Tests verify:
#   - The call completes without error
#   - The number of genes is unchanged
#   - Simulated sequences are valid DNA and have the same length
#   - Simulated sequences differ from the originals (probabilistic; may
#     occasionally fail for a tiny genome; seed is set for reproducibility)
# ======================================================================

context("simulateGenome -- ROC and FONSE (task #6)")

fastaFile <- file.path("UnitTestingData", "testMCMCROCFiles", "simulatedAllUniqueR.fasta")

# ----------------------------------------------------------------------
# Helper: given a genome, return a character vector of all gene sequences
# ----------------------------------------------------------------------
genome_seqs <- function(genome) {
  genes <- genome$getGenes(simulated = FALSE)
  vapply(genes, function(g) g$seq, character(1))
}

sim_seqs <- function(genome) {
  genes <- genome$getGenes(simulated = TRUE)
  vapply(genes, function(g) g$seq, character(1))
}

# ----------------------------------------------------------------------
# ROC simulateGenome
# ----------------------------------------------------------------------
set.seed(42)
genome_roc    <- initializeGenomeObject(file = fastaFile)
original_seqs <- genome_seqs(genome_roc)
n_genes       <- length(genome_roc)
original_lens <- nchar(original_seqs)

parameter_roc <- initializeParameterObject(
  genome         = genome_roc,
  sphi           = 1,
  num.mixtures   = 1,
  gene.assignment = rep(1, n_genes)
)
model_roc <- initializeModelObject(parameter_roc, "ROC")

test_that("ROC simulateGenome returns NULL invisibly", {
  result <- model_roc$simulateGenome(genome_roc)
  expect_null(result)
})

test_that("ROC simulateGenome preserves gene count", {
  expect_equal(length(genome_roc), n_genes)
})

test_that("ROC simulated sequences have the same length as originals", {
  sim <- sim_seqs(genome_roc)
  expect_equal(nchar(sim), original_lens)
})

test_that("ROC simulated sequences contain only valid DNA bases", {
  sim <- sim_seqs(genome_roc)
  invalid <- grep("[^ACGTacgt]", sim, value = TRUE)
  expect_equal(length(invalid), 0L)
})

test_that("ROC simulated sequences differ from originals (seed 42, 1000 genes)", {
  sim <- sim_seqs(genome_roc)
  # With 1000 genes and random parameters, at least one sequence must change.
  expect_false(all(sim == original_seqs))
})

test_that("ROC simulated sequences are multiples of 3 (valid codon frames)", {
  sim <- sim_seqs(genome_roc)
  expect_true(all(nchar(sim) %% 3 == 0))
})

# ----------------------------------------------------------------------
# FONSE simulateGenome
# ----------------------------------------------------------------------
set.seed(42)
genome_fonse    <- initializeGenomeObject(file = fastaFile)
original_seqs_f <- genome_seqs(genome_fonse)

parameter_fonse <- initializeParameterObject(
  genome          = genome_fonse,
  sphi            = 1,
  num.mixtures    = 1,
  gene.assignment = rep(1, length(genome_fonse)),
  model           = "FONSE",
  init.initiation.cost = 4
)
model_fonse <- initializeModelObject(parameter_fonse, "FONSE")

test_that("FONSE simulateGenome returns NULL invisibly", {
  result <- model_fonse$simulateGenome(genome_fonse)
  expect_null(result)
})

test_that("FONSE simulateGenome preserves gene count", {
  expect_equal(length(genome_fonse), length(genome_roc))
})

test_that("FONSE simulated sequences have the same length as originals", {
  sim <- sim_seqs(genome_fonse)
  expect_equal(nchar(sim), original_lens)
})

test_that("FONSE simulated sequences contain only valid DNA bases", {
  sim <- sim_seqs(genome_fonse)
  invalid <- grep("[^ACGTacgt]", sim, value = TRUE)
  expect_equal(length(invalid), 0L)
})

test_that("FONSE simulated sequences are multiples of 3", {
  sim <- sim_seqs(genome_fonse)
  expect_true(all(nchar(sim) %% 3 == 0))
})


# ======================================================================
# Tests for ID/description preservation across simulateGenome + the new
# writeFasta(simulated=TRUE, includeDescription=TRUE) option.
#
# Prior to this branch, ROCModel::simulateGenome and FONSEModel::simulateGenome
# constructed simulated Gene objects with the wrong constructor arg order,
# leaving id = "Simulated Gene" and description = source-gene-id.  This
# made writeFasta(simulated=TRUE) emit bare ">ID" headers and discarded any
# original SGD-style annotation.  The fix preserves id and description on
# the simulated gene; writeFasta(includeDescription=TRUE) splices a
# descriptionTag (default "[simulated]") between the ID and the rest of
# the original description so downstream tools can tell sim from real.
# ======================================================================

context("simulateGenome + writeFasta ID / description preservation")

test_that("simulated genes retain the source gene's id (not 'Simulated Gene')", {
  sim.genes <- genome_roc$getGenes(simulated = TRUE)
  src.genes <- genome_roc$getGenes(simulated = FALSE)
  expect_equal(length(sim.genes), length(src.genes))
  sim.ids <- vapply(sim.genes, function(g) g$id, character(1))
  src.ids <- vapply(src.genes, function(g) g$id, character(1))
  expect_equal(sim.ids, src.ids)
  expect_false(any(sim.ids == "Simulated Gene"))
})

test_that("simulated genes retain the source gene's description", {
  sim.descs <- vapply(genome_roc$getGenes(simulated = TRUE),
                       function(g) g$description, character(1))
  src.descs <- vapply(genome_roc$getGenes(simulated = FALSE),
                       function(g) g$description, character(1))
  expect_equal(sim.descs, src.descs)
})

test_that("genome$writeFasta(file, TRUE) writes bare-ID headers (backward compat)", {
  out <- tempfile(fileext = ".fasta")
  on.exit(unlink(out), add = TRUE)
  genome_roc$writeFasta(out, TRUE)
  hdrs <- grep("^>", readLines(out), value = TRUE)
  src.ids <- vapply(genome_roc$getGenes(simulated = FALSE),
                    function(g) g$id, character(1))
  expect_equal(length(hdrs), length(src.ids))
  expect_equal(hdrs, paste0(">", src.ids))       # bare ID, no description
})

test_that("genome$writeFasta(file, FALSE) writes full descriptions (unchanged)", {
  out <- tempfile(fileext = ".fasta")
  on.exit(unlink(out), add = TRUE)
  genome_roc$writeFasta(out, FALSE)
  hdrs <- grep("^>", readLines(out), value = TRUE)
  src.descs <- vapply(genome_roc$getGenes(simulated = FALSE),
                      function(g) g$description, character(1))
  expect_equal(hdrs, paste0(">", src.descs))
})

test_that("genome$writeFastaWithDescription splices the tag between ID and original description", {
  out <- tempfile(fileext = ".fasta")
  on.exit(unlink(out), add = TRUE)
  genome_roc$writeFastaWithDescription(out, TRUE, "[simulated]")
  hdrs <- grep("^>", readLines(out), value = TRUE)

  src.genes <- genome_roc$getGenes(simulated = FALSE)
  for (i in seq_along(hdrs)) {
    src.id   <- src.genes[[i]]$id
    src.desc <- src.genes[[i]]$description
    after.id <- if (nchar(src.desc) > nchar(src.id) &&
                    startsWith(src.desc, paste0(src.id, " "))) {
                    substring(src.desc, nchar(src.id) + 2L)
                } else ""
    expected <- if (nzchar(after.id)) {
        paste0(">", src.id, " [simulated] ", after.id)
    } else {
        paste0(">", src.id, " [simulated]")
    }
    expect_equal(hdrs[i], expected,
                 info = paste("header mismatch at gene", src.id))
  }
})

test_that("genome$writeFastaWithDescription accepts a custom tag", {
  out <- tempfile(fileext = ".fasta")
  on.exit(unlink(out), add = TRUE)
  genome_roc$writeFastaWithDescription(out, TRUE, "[sim-from-foo]")
  hdrs <- grep("^>", readLines(out), value = TRUE)
  expect_true(all(grepl("\\[sim-from-foo\\]", hdrs)))
  expect_false(any(grepl("\\[simulated\\]", hdrs)))
})

test_that("genome$writeFastaWithDescription with empty tag writes description with no marker", {
  out <- tempfile(fileext = ".fasta")
  on.exit(unlink(out), add = TRUE)
  genome_roc$writeFastaWithDescription(out, TRUE, "")
  hdrs <- grep("^>", readLines(out), value = TRUE)
  src.descs <- vapply(genome_roc$getGenes(simulated = FALSE),
                      function(g) g$description, character(1))
  expect_equal(hdrs, paste0(">", src.descs))
})

test_that("R wrapper writeFasta() dispatches based on includeDescription", {
  out <- tempfile(fileext = ".fasta")
  on.exit(unlink(out), add = TRUE)

  ## Default: simulated=FALSE -> full description
  writeFasta(genome_roc, out)
  hdrs <- grep("^>", readLines(out), value = TRUE)
  src.descs <- vapply(genome_roc$getGenes(simulated = FALSE),
                      function(g) g$description, character(1))
  expect_equal(hdrs, paste0(">", src.descs))

  ## simulated=TRUE, includeDescription=FALSE -> bare ID
  writeFasta(genome_roc, out, simulated = TRUE)
  hdrs <- grep("^>", readLines(out), value = TRUE)
  src.ids <- vapply(genome_roc$getGenes(simulated = FALSE),
                    function(g) g$id, character(1))
  expect_equal(hdrs, paste0(">", src.ids))

  ## simulated=TRUE, includeDescription=TRUE -> ID + tag + description
  writeFasta(genome_roc, out, simulated = TRUE, includeDescription = TRUE)
  hdrs <- grep("^>", readLines(out), value = TRUE)
  expect_true(all(grepl("\\[simulated\\]", hdrs)))

  ## simulated=TRUE, includeDescription=TRUE with custom tag
  writeFasta(genome_roc, out, simulated = TRUE,
             includeDescription = TRUE,
             descriptionTag = "[custom-tag]")
  hdrs <- grep("^>", readLines(out), value = TRUE)
  expect_true(all(grepl("\\[custom-tag\\]", hdrs)))
})
