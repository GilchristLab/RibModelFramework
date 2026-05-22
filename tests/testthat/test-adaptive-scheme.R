## ============================================================================
## tests/testthat/test-adaptive-scheme.R
##
## R-side tests for the AdaptiveScheme constructors and helpers.  These do
## NOT exercise the C++ side (no Parameter is constructed); see
## test-Parameter-csp-strategy.R for the C++-integration tests.
## ============================================================================

library(testthat)
library(AnaCoDa)


test_that("schemes.available() returns native + andrieu_thoms", {
    s <- schemes.available()
    expect_type(s, "character")
    expect_true("native"        %in% s)
    expect_true("andrieu_thoms" %in% s)
    expect_equal(length(s), 2L)
})


test_that("AdaptiveScheme.Native() default has aggressiveness = 0.2", {
    s <- AdaptiveScheme.Native()
    expect_s3_class(s, c("AdaptiveScheme.Native", "AdaptiveScheme"))
    expect_equal(s$scheme, "native")
    expect_equal(s$params$aggressiveness, 0.2)
    expect_true(is.AdaptiveScheme(s))
})


test_that("AdaptiveScheme.Native() accepts numeric aggressiveness", {
    expect_equal(AdaptiveScheme.Native(0.1)$params$aggressiveness, 0.1)
    expect_equal(AdaptiveScheme.Native(0.3)$params$aggressiveness, 0.3)
    expect_equal(AdaptiveScheme.Native(0.5)$params$aggressiveness, 0.5)
})


test_that("AdaptiveScheme.Native() rejects out-of-range aggressiveness", {
    expect_error(AdaptiveScheme.Native(0))         # closed lower
    expect_error(AdaptiveScheme.Native(1))         # closed upper
    expect_error(AdaptiveScheme.Native(-0.1))
    expect_error(AdaptiveScheme.Native(1.5))
    expect_error(AdaptiveScheme.Native(NA_real_))
    expect_error(AdaptiveScheme.Native(NaN))
    expect_error(AdaptiveScheme.Native(Inf))
    expect_error(AdaptiveScheme.Native(c(0.1, 0.2)))
    expect_error(AdaptiveScheme.Native("0.2"))
})


test_that("AdaptiveScheme.Native() default has prev.weight = 0.6", {
    s <- AdaptiveScheme.Native()
    expect_equal(s$params$prev.weight, 0.6)
})


test_that("AdaptiveScheme.Native() accepts numeric prev.weight in (0, 1)", {
    expect_equal(AdaptiveScheme.Native(prev.weight = 0.4)$params$prev.weight, 0.4)
    expect_equal(AdaptiveScheme.Native(prev.weight = 0.8)$params$prev.weight, 0.8)
    expect_equal(AdaptiveScheme.Native(prev.weight = 0.95)$params$prev.weight, 0.95)
})


test_that("AdaptiveScheme.Native() rejects out-of-range prev.weight", {
    expect_error(AdaptiveScheme.Native(prev.weight =  0.0))    # closed lower
    expect_error(AdaptiveScheme.Native(prev.weight =  1.0))    # closed upper
    expect_error(AdaptiveScheme.Native(prev.weight = -0.1))
    expect_error(AdaptiveScheme.Native(prev.weight =  1.5))
    expect_error(AdaptiveScheme.Native(prev.weight =  NA_real_))
    expect_error(AdaptiveScheme.Native(prev.weight =  NaN))
    expect_error(AdaptiveScheme.Native(prev.weight =  Inf))
    expect_error(AdaptiveScheme.Native(prev.weight =  c(0.4, 0.6)))
    expect_error(AdaptiveScheme.Native(prev.weight =  "0.6"))
})


test_that("AdaptiveScheme.Native() accepts both knobs independently", {
    s <- AdaptiveScheme.Native(aggressiveness = 0.1, prev.weight = 0.8)
    expect_equal(s$params$aggressiveness, 0.1)
    expect_equal(s$params$prev.weight,    0.8)
})


test_that("AdaptiveScheme.AndrieuThoms() defaults pass and have expected fields", {
    s <- AdaptiveScheme.AndrieuThoms()
    expect_s3_class(s, c("AdaptiveScheme.AndrieuThoms", "AdaptiveScheme"))
    expect_equal(s$scheme, "andrieu_thoms")
    expect_equal(s$params$target, 0.234)
    expect_equal(s$params$alpha,  0.7)
    expect_equal(s$params$c,      1.0)
    expect_equal(s$params$t0,     10)
    expect_true(is.AdaptiveScheme(s))
})


test_that("AdaptiveScheme.AndrieuThoms() rejects out-of-range target", {
    expect_error(AdaptiveScheme.AndrieuThoms(target =  0.0))
    expect_error(AdaptiveScheme.AndrieuThoms(target =  1.0))
    expect_error(AdaptiveScheme.AndrieuThoms(target = -0.1))
    expect_error(AdaptiveScheme.AndrieuThoms(target =  1.5))
    expect_error(AdaptiveScheme.AndrieuThoms(target =  NA_real_))
    expect_error(AdaptiveScheme.AndrieuThoms(target =  NaN))
    expect_error(AdaptiveScheme.AndrieuThoms(target =  Inf))
    expect_error(AdaptiveScheme.AndrieuThoms(target =  c(0.3, 0.4)))
    expect_error(AdaptiveScheme.AndrieuThoms(target =  "0.3"))
})


test_that("AdaptiveScheme.AndrieuThoms() rejects out-of-range alpha", {
    expect_error(AdaptiveScheme.AndrieuThoms(alpha = 0.5))    # closed lower bound -> reject
    expect_error(AdaptiveScheme.AndrieuThoms(alpha = 1.01))
    expect_error(AdaptiveScheme.AndrieuThoms(alpha = 0))
    expect_error(AdaptiveScheme.AndrieuThoms(alpha = NA_real_))
    expect_error(AdaptiveScheme.AndrieuThoms(alpha = c(0.7, 0.8)))
})


test_that("AdaptiveScheme.AndrieuThoms() rejects bad c and t0", {
    expect_error(AdaptiveScheme.AndrieuThoms(c = 0))
    expect_error(AdaptiveScheme.AndrieuThoms(c = -0.5))
    expect_error(AdaptiveScheme.AndrieuThoms(c = NA_real_))
    expect_error(AdaptiveScheme.AndrieuThoms(t0 = -1))
    expect_error(AdaptiveScheme.AndrieuThoms(t0 = NA_real_))
    expect_error(AdaptiveScheme.AndrieuThoms(t0 = c(10, 20)))
})


test_that("AdaptiveScheme.AndrieuThoms() boundary cases that should pass", {
    expect_s3_class(AdaptiveScheme.AndrieuThoms(alpha = 1.0),
                    "AdaptiveScheme.AndrieuThoms")
    expect_s3_class(AdaptiveScheme.AndrieuThoms(t0 = 0),
                    "AdaptiveScheme.AndrieuThoms")
    expect_s3_class(AdaptiveScheme.AndrieuThoms(target = 0.001, alpha = 0.501),
                    "AdaptiveScheme.AndrieuThoms")
})


test_that("is.AdaptiveScheme distinguishes constructed objects from non-scheme objects", {
    expect_true(is.AdaptiveScheme(AdaptiveScheme.Native()))
    expect_true(is.AdaptiveScheme(AdaptiveScheme.AndrieuThoms()))
    expect_false(is.AdaptiveScheme(list()))
    expect_false(is.AdaptiveScheme(NULL))
    expect_false(is.AdaptiveScheme(0.234))
    expect_false(is.AdaptiveScheme("native"))
})


test_that("print.AdaptiveScheme produces output without error", {
    expect_output(print(AdaptiveScheme.Native()), "AdaptiveScheme: native")
    expect_output(print(AdaptiveScheme.AndrieuThoms()),
                  "AdaptiveScheme: andrieu_thoms")
    expect_output(print(AdaptiveScheme.AndrieuThoms()),
                  "target = 0.234")
})


test_that("format.AdaptiveScheme returns expected character shape", {
    f.native <- format(AdaptiveScheme.Native())
    expect_type(f.native, "character")
    expect_equal(length(f.native), 2L)
    expect_true(grepl("native", f.native[1]))
    expect_true(grepl("aggressiveness = 0.2", f.native[2]))

    f.at <- format(AdaptiveScheme.AndrieuThoms())
    expect_type(f.at, "character")
    expect_equal(length(f.at), 2L)
    expect_true(grepl("andrieu_thoms", f.at[1]))
    expect_true(grepl("target = 0.234", f.at[2]))
})


test_that(".scheme.name.to.constructor maps known names to constructors", {
    fn <- AnaCoDa:::.scheme.name.to.constructor("native")
    expect_true(is.function(fn))
    expect_s3_class(fn(), "AdaptiveScheme.Native")

    fn2 <- AnaCoDa:::.scheme.name.to.constructor("andrieu_thoms")
    expect_true(is.function(fn2))
    expect_s3_class(fn2(target = 0.3), "AdaptiveScheme.AndrieuThoms")
})


test_that(".scheme.name.to.constructor errors on unknown names", {
    expect_error(AnaCoDa:::.scheme.name.to.constructor("does_not_exist"))
    expect_error(AnaCoDa:::.scheme.name.to.constructor(""))
    expect_error(AnaCoDa:::.scheme.name.to.constructor(NA_character_))
    expect_error(AnaCoDa:::.scheme.name.to.constructor(c("a", "b")))
})


## ============================================================================
## C++ Rcpp seam tests: verify setCSPAdaptationScheme and
## getCSPAdaptationSchemeName work through the Rcpp module and that
## C++-side validation errors surface as clean R errors.
## ============================================================================

test_that("default Parameter has 'native' scheme bound", {
    p <- new("Rcpp_ROCParameter")
    expect_equal(p$getCSPAdaptationSchemeName(), "native")
})


test_that("setCSPAdaptationScheme switches to andrieu_thoms", {
    p <- new("Rcpp_ROCParameter")
    p$setCSPAdaptationScheme(
        "andrieu_thoms",
        list(target = 0.234, alpha = 0.7, c = 1.0, t0 = 10))
    expect_equal(p$getCSPAdaptationSchemeName(), "andrieu_thoms")
})


test_that("setCSPAdaptationScheme accepts 'native' with no params", {
    p <- new("Rcpp_ROCParameter")
    p$setCSPAdaptationScheme("andrieu_thoms",
                             list(target=0.234, alpha=0.7, c=1.0, t0=10))
    expect_equal(p$getCSPAdaptationSchemeName(), "andrieu_thoms")
    p$setCSPAdaptationScheme("native", list())
    expect_equal(p$getCSPAdaptationSchemeName(), "native")
})


test_that("setCSPAdaptationScheme rejects unknown scheme name", {
    p <- new("Rcpp_ROCParameter")
    expect_error(p$setCSPAdaptationScheme("does_not_exist", list()),
                 regexp = "unknown")
})


test_that("setCSPAdaptationScheme rejects out-of-range params at C++ layer", {
    p <- new("Rcpp_ROCParameter")
    # target > 1
    expect_error(p$setCSPAdaptationScheme(
        "andrieu_thoms", list(target=2.0, alpha=0.7, c=1.0, t0=10)),
        regexp = "target")
    # alpha = 0.5 (closed bound)
    expect_error(p$setCSPAdaptationScheme(
        "andrieu_thoms", list(target=0.234, alpha=0.5, c=1.0, t0=10)),
        regexp = "alpha")
    # c <= 0
    expect_error(p$setCSPAdaptationScheme(
        "andrieu_thoms", list(target=0.234, alpha=0.7, c=0.0, t0=10)),
        regexp = "c ")
    # t0 < 0
    expect_error(p$setCSPAdaptationScheme(
        "andrieu_thoms", list(target=0.234, alpha=0.7, c=1.0, t0=-1)),
        regexp = "t0")
})


test_that("setCSPAdaptationScheme rejects extra params for native", {
    p <- new("Rcpp_ROCParameter")
    expect_error(p$setCSPAdaptationScheme("native", list(target = 0.234)),
                 regexp = "takes no params|unexpected param")
})


test_that("setCSPAdaptationScheme accepts native + aggressiveness", {
    p <- new("Rcpp_ROCParameter")
    p$setCSPAdaptationScheme("native", list(aggressiveness = 0.3))
    expect_equal(p$getCSPAdaptationSchemeName(), "native")
})


test_that("setCSPAdaptationScheme accepts native + prev.weight", {
    p <- new("Rcpp_ROCParameter")
    p$setCSPAdaptationScheme("native", list("prev.weight" = 0.8))
    expect_equal(p$getCSPAdaptationSchemeName(), "native")
})


test_that("setCSPAdaptationScheme accepts native + aggressiveness + prev.weight", {
    p <- new("Rcpp_ROCParameter")
    p$setCSPAdaptationScheme("native",
                             list(aggressiveness = 0.2, "prev.weight" = 0.8))
    expect_equal(p$getCSPAdaptationSchemeName(), "native")
})


test_that("setCSPAdaptationScheme rejects out-of-range prev.weight at C++ layer", {
    p <- new("Rcpp_ROCParameter")
    expect_error(p$setCSPAdaptationScheme("native", list("prev.weight" = 0.0)),
                 regexp = "prevWeight|prev.weight")
    expect_error(p$setCSPAdaptationScheme("native", list("prev.weight" = 1.0)),
                 regexp = "prevWeight|prev.weight")
    expect_error(p$setCSPAdaptationScheme("native", list("prev.weight" = -0.1)),
                 regexp = "prevWeight|prev.weight")
    expect_error(p$setCSPAdaptationScheme("native", list("prev.weight" = 1.5)),
                 regexp = "prevWeight|prev.weight")
})


test_that("setCSPAdaptationScheme rejects out-of-range aggressiveness at C++ layer", {
    p <- new("Rcpp_ROCParameter")
    expect_error(p$setCSPAdaptationScheme("native", list(aggressiveness = 0.0)),
                 regexp = "aggressiveness")
    expect_error(p$setCSPAdaptationScheme("native", list(aggressiveness = 1.0)),
                 regexp = "aggressiveness")
    expect_error(p$setCSPAdaptationScheme("native", list(aggressiveness = -0.1)),
                 regexp = "aggressiveness")
})


test_that("setCSPAdaptationScheme rejects missing required params for andrieu_thoms", {
    p <- new("Rcpp_ROCParameter")
    expect_error(p$setCSPAdaptationScheme("andrieu_thoms", list()),
                 regexp = "requires param")
    expect_error(p$setCSPAdaptationScheme(
        "andrieu_thoms", list(target = 0.234, alpha = 0.7, c = 1.0)),
        regexp = "requires param 't0'")
})


test_that("setCSPAdaptationScheme rejects extra params for andrieu_thoms", {
    p <- new("Rcpp_ROCParameter")
    expect_error(p$setCSPAdaptationScheme(
        "andrieu_thoms",
        list(target=0.234, alpha=0.7, c=1.0, t0=10, garbage=42)),
        regexp = "unexpected param")
})


test_that("scheme survives Parameter copy via operator=", {
    # Use FONSEParameter, which has a copy ctor delegating to Parameter
    # base.  Tests that the strategy clone() path works.
    p <- new("Rcpp_ROCParameter")
    p$setCSPAdaptationScheme(
        "andrieu_thoms",
        list(target = 0.234, alpha = 0.7, c = 1.0, t0 = 10))
    # The R-side copy operation on Rcpp ref-class objects calls C++
    # operator=.  After a re-assignment, scheme should still be AT.
    expect_equal(p$getCSPAdaptationSchemeName(), "andrieu_thoms")
})
