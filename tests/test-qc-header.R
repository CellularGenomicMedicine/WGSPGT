library(testthat)

root <- if (file.exists("Haplarithmisis/Haplarithmisis/functions")) "." else ".."
fn <- file.path(root, "Haplarithmisis/Haplarithmisis/functions")

test_that("qcbyparents writes header with same field count as data rows", {
  source(file.path(fn, "Haplarithmisis_qcbyparents.R"))
  set.seed(1)
  g <- c("AA", "AB", "BB")
  Gtypes <- data.frame(Chr = rep(c("1", "2"), each = 30), Position = 1:60,
                       Fa = sample(g, 60, TRUE), Mo = sample(g, 60, TRUE),
                       E1 = sample(g, 60, TRUE), stringsAsFactors = FALSE)
  out <- tempfile(); dir.create(out)
  invisible(capture.output(
    Haplarithmisis_qcbyparents("Fa", "Mo", Gtypes, "E1", c("1", "2"), out)))
  for (f in c("E1_ChrSpec_ADO.txt", "E1_ChrSpec_ADI.txt")) {
    lines <- readLines(file.path(out, f))
    nf <- lengths(strsplit(lines, "\t", fixed = TRUE))
    expect_equal(length(unique(nf)), 1L, info = f)
    expect_equal(nf[1], 4L)
    expect_equal(strsplit(lines[1], "\t")[[1]][1], "")
  }
})

test_that("every row.names = TRUE write.table uses col.names = NA", {
  for (f in c("Haplarithmisis_qcgtype.R", "Haplarithmisis_qcbyparents.R")) {
    txt <- paste(readLines(file.path(fn, f), warn = FALSE), collapse = "\n")
    expect_false(grepl("col.names = TRUE, row.names = TRUE", txt, fixed = TRUE), info = f)
    expect_equal(lengths(regmatches(txt, gregexpr("row.names = TRUE", txt))),
                 lengths(regmatches(txt, gregexpr("col.names = NA, row.names = TRUE", txt))), info = f)
  }
})
