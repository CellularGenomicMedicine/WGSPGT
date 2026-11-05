library(testthat)
source("../Haplarithmisis/Haplarithmisis/functions/Haplarithmisis_disthap.R")

test_that("disthap writes a prefixed file per parent", {
  Intp <- list(dataHap = data.frame(Chr = rep("1", 6), Position = c(10, 20, 30, 40, 50, 60),
                                    E1_Mat = c(1, 1, 2, 2, 2, 1), E1_Pat = c(1, 2, 2, 1, 1, 1),
                                    stringsAsFactors = FALSE))
  Int <- data.frame(Chr = "1", Start = 30, End = 40, stringsAsFactors = FALSE)
  d <- tempfile(); dir.create(d)
  capture.output({
    Haplarithmisis_disthap(Intp, Int, d, "E1", "1", "Father", "Fam", 0.5)
    Haplarithmisis_disthap(Intp, Int, d, "E1", "1", "Mother", "Fam", 0.5)
  })
  expect_true(file.exists(file.path(d, "Father_E1_Fam_0.5_DistanceToHR.txt")))
  expect_true(file.exists(file.path(d, "Mother_E1_Fam_0.5_DistanceToHR.txt")))
  expect_false(file.exists(file.path(d, "DistanceToHR.txt")))
})
