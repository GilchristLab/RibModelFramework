library(testthat)
library(AnaCoDa)

context("compact layout helpers")

# .buildModelLayout: 5-col matrix (y.lbl | 4x data panels)
# 20 data slots: 19 AA codon-freq panels + 1 phi histogram (20th/bottom-right slot).
# Skipped AAs (M/W/X) mean the phi histogram naturally lands in slot 20 after
# 19 AA plot.new() calls -- no blank fill needed.
# Panel numbering:
#   1=title, 2..(n+1)=data slots, n+2=x.lbl, n+3=y.lbl
#
# .buildTraceLayout: 9-col matrix (y.lbl | 4x(wide+narrow trace))
# Panels: 1=title, 2..(2n+1)=interleaved trace+ECDF pairs, x.lbl, y.lbl

# --- .buildModelLayout ---

test_that(".buildModelLayout returns correct structure for n.slots=20", {
    lay <- .buildModelLayout(20L)
    n.rows.data <- ceiling(20L / 4L)  # 5

    expect_equal(nrow(lay$mat), n.rows.data + 2L)  # title + 5 data + x.lbl
    expect_equal(ncol(lay$mat), 5L)
    expect_equal(length(lay$widths),  5L)
    expect_equal(length(lay$heights), n.rows.data + 2L)

    # special panel indices
    expect_equal(lay$x.lbl, 20L + 2L)  # 22
    expect_equal(lay$y.lbl, 20L + 3L)  # 23
})

test_that(".buildModelLayout y.lbl always fills left column", {
    for (n in c(1L, 4L, 5L, 20L)) {
        lay <- .buildModelLayout(n)
        expect_true(all(lay$mat[, 1L] == lay$y.lbl),
                    info = paste("n.slots =", n))
    }
})

test_that(".buildModelLayout all data panels present and in order", {
    lay <- .buildModelLayout(20L)
    nr  <- nrow(lay$mat)
    # rows 2..(nr-1) cols 2..5 are the data block
    data.block <- lay$mat[2L:(nr - 1L), 2L:5L]
    present    <- sort(data.block[data.block > 0L])
    expect_equal(present, 2L:(20L + 1L))
})

test_that(".buildModelLayout n.slots=1 edge case", {
    lay <- .buildModelLayout(1L)
    expect_equal(nrow(lay$mat), 3L)   # title + 1 data row + x.lbl
    expect_equal(ncol(lay$mat), 5L)
    expect_equal(lay$x.lbl, 3L)
    expect_equal(lay$y.lbl, 4L)
})

test_that(".buildModelLayout n.slots=4 fits in one data row", {
    lay <- .buildModelLayout(4L)
    expect_equal(nrow(lay$mat), 3L)
})

test_that(".buildModelLayout n.slots=5 spans two data rows", {
    lay <- .buildModelLayout(5L)
    expect_equal(nrow(lay$mat), 4L)
})

# --- .buildTraceLayout ---

test_that(".buildTraceLayout returns correct structure for n.slots=20", {
    lay <- .buildTraceLayout(20L)
    n.rows.data <- ceiling(20L / 4L)  # 5

    expect_equal(nrow(lay$mat), n.rows.data + 2L)  # title + 5 data + x.lbl
    expect_equal(ncol(lay$mat), 9L)
    expect_equal(length(lay$widths),  9L)
    expect_equal(length(lay$heights), n.rows.data + 2L)

    expect_equal(lay$x.lbl, 2L * 20L + 2L)  # 42
    expect_equal(lay$y.lbl, 2L * 20L + 3L)  # 43
})

test_that(".buildTraceLayout y.lbl always fills left column", {
    for (n in c(1L, 4L, 5L, 20L)) {
        lay <- .buildTraceLayout(n)
        expect_true(all(lay$mat[, 1L] == lay$y.lbl),
                    info = paste("n.slots =", n))
    }
})

test_that(".buildTraceLayout trace+marginal panels all present in data block", {
    lay <- .buildTraceLayout(20L)
    nr  <- nrow(lay$mat)
    data.block <- lay$mat[2L:(nr - 1L), 2L:9L]
    present    <- sort(data.block[data.block > 0L])
    expect_equal(present, 2L:(2L * 20L + 1L))
})

test_that(".buildTraceLayout x.lbl fills bottom row (excluding y.lbl column)", {
    lay <- .buildTraceLayout(20L)
    nr  <- nrow(lay$mat)
    expect_true(all(lay$mat[nr, 2L:9L] == lay$x.lbl))
})

test_that(".buildTraceLayout n.slots=1 edge case", {
    lay <- .buildTraceLayout(1L)
    expect_equal(nrow(lay$mat), 3L)
    expect_equal(ncol(lay$mat), 9L)
    expect_equal(lay$x.lbl, 4L)
    expect_equal(lay$y.lbl, 5L)
})

# --- options functions ---

test_that(".plotModelOptions.common returns correct defaults", {
    opts <- .plotModelOptions.common()
    expect_equal(opts$layout, "original")
    expect_false(opts$show.gene.hist)
    expect_true(opts$show.date)
    expect_named(opts, c("layout", "show.gene.hist", "show.date"))
})

test_that("plotROCOptions returns common defaults", {
    opts <- plotROCOptions()
    expect_equal(opts$layout, "original")
    expect_false(opts$show.gene.hist)
    expect_true(opts$show.date)
})

test_that("plotROCOptions overrides work", {
    opts <- plotROCOptions(layout = "compact-v1", show.gene.hist = TRUE, show.date = FALSE)
    expect_equal(opts$layout, "compact-v1")
    expect_true(opts$show.gene.hist)
    expect_false(opts$show.date)
})

test_that("plotROCOptions accepts compact as layout value", {
    opts <- plotROCOptions(layout = "compact")
    expect_equal(opts$layout, "compact")
})

test_that("plotFONSEOptions includes codon.window with NULL default", {
    opts <- plotFONSEOptions()
    expect_null(opts$codon.window)
    expect_equal(opts$layout, "original")
    expect_named(opts, c("layout", "show.gene.hist", "show.date", "codon.window"))
})

test_that("plotFONSEOptions codon.window override works", {
    opts <- plotFONSEOptions(codon.window = c(1, 100))
    expect_equal(opts$codon.window, c(1, 100))
})

test_that("plotROCOptions and plotFONSEOptions share common fields", {
    roc  <- plotROCOptions()
    fons <- plotFONSEOptions()
    common <- c("layout", "show.gene.hist", "show.date")
    expect_equal(roc[common], fons[common])
})

# --- both helpers ---

test_that("row count scales with n.slots for both layout helpers", {
    for (n in c(1L, 4L, 5L, 8L, 20L)) {
        ml <- .buildModelLayout(n)
        tl <- .buildTraceLayout(n)
        expected.data.rows <- ceiling(n / 4L)
        expect_equal(nrow(ml$mat), expected.data.rows + 2L,
                     info = paste("model n.slots =", n))
        expect_equal(nrow(tl$mat), expected.data.rows + 2L,
                     info = paste("trace n.slots =", n))
    }
})
