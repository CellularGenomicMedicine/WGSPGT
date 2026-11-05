library(testthat)
library(stats)

root <- if (file.exists("Haplarithmisis/Haplarithmisis/functions")) "Haplarithmisis/Haplarithmisis" else "../Haplarithmisis/Haplarithmisis"
source(file.path(root, "functions/Haplarithmisis_testmedfilt.R"))

test_that("window is scaled by chromosome size relative to chr 1", {
  set.seed(1)
  # chr 1: 200 rows, chr 2: 100 rows (half)
  Hap <- data.frame(Names = paste0("s", 1:300),
                    Chr = c(rep("1", 200), rep("2", 100)),
                    Position = 1:300,
                    V = sample(c(0, 1), 300, replace = TRUE),
                    stringsAsFactors = FALSE)
  Window <- 41
  res <- Haplarithmisis_testmedfilt(Hap, Window, c("1", "2"))
  expect_length(res, 300)
  # chr 1 unchanged: full window (41)
  expect_equal(as.numeric(res[1:200]), as.numeric(runmed(Hap$V[1:200], 41, "median")))
  # chr 2: scaled window round(41*100/200)=21 (already odd)
  expect_equal(as.numeric(res[201:300]), as.numeric(runmed(Hap$V[201:300], 21, "median")))
  # and differs from the unscaled (buggy) window
  expect_false(isTRUE(all.equal(as.numeric(res[201:300]),
                                as.numeric(runmed(Hap$V[201:300], 41, "median")))))
})
