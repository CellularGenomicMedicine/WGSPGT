library(testthat)

root <- file.path(getwd(), "Haplarithmisis", "Haplarithmisis")
if (!dir.exists(root)) root <- file.path(getwd(), "..", "Haplarithmisis", "Haplarithmisis")
root <- normalizePath(root)
fn_file <- file.path(root, "functions", "Haplarithmisis_Htyping_Grandparents.R")
main_file <- file.path(root, "Haplarithmisis.R")

test_that("Htyping_Grandparents takes Window and Int", {
  e <- new.env()
  sys.source(fn_file, envir = e)
  expect_true(all(c("Window", "Int") %in% names(formals(e$Haplarithmisis_Htyping_Grandparents))))
})

test_that("main script passes Window and Int to Htyping_Grandparents", {
  src <- paste(readLines(main_file), collapse = "\n")
  m <- regmatches(src, gregexpr("Haplarithmisis_Htyping_Grandparents\\([^)]*\\)", src))[[1]]
  expect_length(m, 1)
  args <- trimws(strsplit(sub("^[^(]*\\(", "", sub("\\)$", "", m)), ",")[[1]])
  expect_true(all(c("Window", "Int") %in% args))
  f <- new.env()
  sys.source(fn_file, envir = f)
  expect_equal(length(args), length(formals(f$Haplarithmisis_Htyping_Grandparents)))
})
