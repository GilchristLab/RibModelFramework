library(testthat)
library(AnaCoDa)

context("compact layout helpers")

# .buildModelLayout: 10-col matrix (y.lbl | 4x(wide+narrow AA) | phi.hist)
# Panels: 1=title, 2..(2n+1)=interleaved AA+ECDF pairs, phi.hist, n.strip, x.lbl, y.lbl
#
# .buildTraceLayout: 9-col matrix (y.lbl | 4x(wide+narrow trace) )
# Panels: 1=title, 2..(2n+1)=interleaved trace+ECDF pairs, x.lbl, y.lbl

# --- .buildModelLayout ---

test_that(".buildModelLayout returns correct structure for n.slots=20", {
    lay <- .buildModelLayout(20L)
    n.rows.data <- ceiling(20L / 4L)  # 5

    expect_equal(nrow(lay$mat), n.rows.data + 3L)  # title + 5 data + n.strip + x.lbl
    expect_equal(ncol(lay$mat), 10L)
    expect_equal(length(lay$widths),  10L)
    expect_equal(length(lay$heights), n.rows.data + 3L)

    # special panel indices
    expect_equal(lay$phi.hist, 2L * 20L + 2L)  # 42
    expect_equal(lay$n.strip,  2L * 20L + 3L)  # 43
    expect_equal(lay$x.lbl,   2L * 20L + 4L)  # 44
    expect_equal(lay$y.lbl,   2L * 20L + 5L)  # 45
})

test_that(".buildModelLayout y.lbl always fills left column", {
    for (n in c(1L, 4L, 5L, 20L)) {
        lay <- .buildModelLayout(n)
        expect_true(all(lay$mat[, 1L] == lay$y.lbl),
                    info = paste("n.slots =", n))
    }
})

test_that(".buildModelLayout phi.hist fills right column except x-label row", {
    lay <- .buildModelLayout(20L)
    nr  <- nrow(lay$mat)
    expect_true(all(lay$mat[-nr, 10L] == lay$phi.hist))
    expect_equal(lay$mat[nr, 10L], 0L)  # x.lbl row has 0 in phi.hist column
})

test_that(".buildModelLayout AA+marginal panels are all present and in order", {
    lay <- .buildModelLayout(20L)
    nr  <- nrow(lay$mat)
    # rows 2..(nr-2) cols 2..9 are the data block
    data.block <- lay$mat[2L:(nr - 2L), 2L:9L]
    present    <- sort(data.block[data.block > 0L])
    expect_equal(present, 2L:(2L * 20L + 1L))
})

test_that(".buildModelLayout n.slots=1 edge case", {
    lay <- .buildModelLayout(1L)
    expect_equal(nrow(lay$mat), 4L)   # title + 1 data row + n.strip + x.lbl
    expect_equal(ncol(lay$mat), 10L)
    expect_equal(lay$phi.hist, 4L)
    expect_equal(lay$n.strip,  5L)
    expect_equal(lay$x.lbl,   6L)
    expect_equal(lay$y.lbl,   7L)
})

test_that(".buildModelLayout n.slots=4 fits in one data row", {
    lay <- .buildModelLayout(4L)
    expect_equal(ceiling(4L / 4L), 1L)
    expect_equal(nrow(lay$mat), 4L)
})

test_that(".buildModelLayout n.slots=5 spans two data rows", {
    lay <- .buildModelLayout(5L)
    expect_equal(ceiling(5L / 4L), 2L)
    expect_equal(nrow(lay$mat), 5L)
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
    expect_equal(nrow(lay$mat), 3L)  # title + 1 data row + x.lbl
    expect_equal(ncol(lay$mat), 9L)
    expect_equal(lay$x.lbl, 4L)
    expect_equal(lay$y.lbl, 5L)
})

# --- both helpers ---

test_that("row count scales with n.slots for both layout helpers", {
    for (n in c(1L, 4L, 5L, 8L, 20L)) {
        ml <- .buildModelLayout(n)
        tl <- .buildTraceLayout(n)
        expected.data.rows <- ceiling(n / 4L)
        expect_equal(nrow(ml$mat), expected.data.rows + 3L,
                     info = paste("model n.slots =", n))
        expect_equal(nrow(tl$mat), expected.data.rows + 2L,
                     info = paste("trace n.slots =", n))
    }
})
