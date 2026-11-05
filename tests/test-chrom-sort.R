library(testthat)

# Snippet test: re-evaluates the ordering/intersect expressions used in Haplarithmisis.R
# (the script itself needs the full pipeline, so it is only checked by parsing and grep).
test_that("Gtypes are ordered 1, 2, 10, X and Chroms follows that order", {
  Gtypes <- data.frame(Chr = c("X", "10", "2", "1", "10", "1"),
                       Position = c(5, 9, 3, 8, 2, 1), stringsAsFactors = FALSE)
  Gtypes <- Gtypes[order(factor(Gtypes$Chr, levels = c(1:22, "X")), Gtypes$Position), ]
  expect_equal(Gtypes$Chr, c("1", "1", "2", "10", "10", "X"))
  expect_equal(Gtypes$Position, c(1, 8, 3, 2, 9, 5))
  Chroms <- intersect(c(1:22, "X"), unique(Gtypes$Chr))
  expect_equal(Chroms, c("1", "2", "10", "X"))
})

test_that("Haplarithmisis.R parses and uses the new expressions", {
  f <- file.path("Haplarithmisis", "Haplarithmisis", "Haplarithmisis.R")
  if (!file.exists(f)) f <- file.path("..", f)
  expect_silent(invisible(parse(f)))
  src <- readLines(f)
  expect_true(any(grepl('Chroms <- intersect(c(1:22, "X"), unique(Gtypes$Chr))', src, fixed = TRUE)))
  expect_equal(sum(grepl('order(factor(', src, fixed = TRUE)), 2)
})
