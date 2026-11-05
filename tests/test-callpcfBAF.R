# Tests for Haplarithmisis_callpcfBAF.
# Run: Rscript -e 'testthat::test_file("tests/test-callpcfBAF.R")'
# Needs only testthat and zoo (no Bioconductor).
library(testthat)

root <- normalizePath(if (dir.exists("Haplarithmisis")) "." else "..")
scratch <- tempfile("script")
dir.create(scratch)
file.symlink(file.path(root, "Haplarithmisis"), file.path(scratch, "analyses"))
source(file.path(root, "Haplarithmisis", "Haplarithmisis", "functions", "Haplarithmisis_callpcfBAF.R"))

# Each call runs under a time limit so a hang counts as a failure
run <- function(df, gamma = 1, secs = 20) {
  setTimeLimit(elapsed = secs, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf))
  invisible(capture.output(res <- Haplarithmisis_callpcfBAF(scratch, df, gamma, 25, "E1")))
  res
}

# Two chromosomes of 150 SNPs with a step so MAD is non-zero
make <- function(chr1 = NULL, chr2 = NULL) {
  set.seed(1)
  b1 <- c(rnorm(75, 0.3, 0.05), rnorm(75, 0.7, 0.05))
  b2 <- c(rnorm(75, 0.5, 0.05), rnorm(75, 0.2, 0.05))
  if (!is.null(chr1)) b1 <- chr1(b1)
  if (!is.null(chr2)) b2 <- chr2(b2)
  data.frame(Names = paste0("s", 1:300), Chr = rep(c("1", "2"), each = 150),
             Position = rep(1:150, 2) * 1000, E1 = c(b1, b2), stringsAsFactors = FALSE)
}

test_that("no NA input still works", {
  res <- run(make())
  expect_equal(nrow(res), 300)
  expect_false(anyNA(res$E1))
})

test_that("one leading NA is filled", {
  res <- run(make(function(x) replace(x, 1, NA)))
  expect_equal(nrow(res), 300)
  expect_false(anyNA(res$E1))
})

test_that("several leading NAs do not hang", {
  res <- run(make(function(x) replace(x, 1:5, NA)))
  expect_equal(nrow(res), 300)
  expect_false(anyNA(res$E1))
})

test_that("internal NAs are filled with the preceding value", {
  res <- run(make(function(x) replace(x, c(10, 11, 80), NA)))
  expect_false(anyNA(res$E1))
})

test_that("an all-NA chromosome is passed through as NA, other chromosomes segmented", {
  res <- run(make(chr2 = function(x) rep(NA_real_, length(x))))
  expect_equal(nrow(res), 300)
  expect_true(all(is.na(res$E1[res$Chr == "2"])))
  expect_false(anyNA(res$E1[res$Chr == "1"]))
})

test_that("all-NA genome (sdevGenome NA, no fallback) returns raw values without error", {
  df <- make()
  df$E1 <- NA_real_
  res <- run(df)
  expect_equal(nrow(res), 300)
  expect_true(all(is.na(res$E1)))
})

test_that("a constant chromosome (MAD 0) falls back to the genome scale", {
  res <- run(make(chr1 = function(x) rep(0.5, length(x))))
  expect_equal(nrow(res), 300)
  expect_false(anyNA(res$E1))
})

test_that("sdevGenome of 0 everywhere uses the per-chromosome fallback without error", {
  res <- run(make(chr1 = function(x) rep(0.5, length(x)), chr2 = function(x) rep(0.5, length(x))))
  expect_equal(nrow(res), 300)
  expect_false(anyNA(res$E1))
})

test_that("Y chromosome is dropped and short chromosomes are returned raw", {
  df <- make()
  df <- rbind(df, data.frame(Names = "y1", Chr = "Y", Position = 1, E1 = 0.4))
  df <- rbind(df, data.frame(Names = paste0("x", 1:10), Chr = "X", Position = 1:10,
                             E1 = c(NA, NA, runif(8))))
  res <- run(df)
  expect_false("Y" %in% res$Chr)
  expect_equal(sum(res$Chr == "X"), 10)
  expect_false(anyNA(res$E1))
})
