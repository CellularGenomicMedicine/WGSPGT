# Tests for Haplarithmisis_intphappropser3
# Run from the repo root: Rscript -e 'testthat::test_file("tests/test-intphappropser3.R")'

library(testthat)
library(MASS)

source(test_path("..", "Haplarithmisis", "Haplarithmisis", "functions", "Haplarithmisis_intphappropser3.R"))

StdCols <- c("Chr", "Start", "Stop", "Length", "CumLength", "Value", "#InfSNPs", "#InfSNPsDisc", "Frac")
ResCols <- c(StdCols, "Fitted", "Resid", "Prob")

# Build Haps for embryo E1 on chromosomes 1-3.
# Blocks: named list per chromosome of list(len, val, zero), where zero is the number of
# SNPs at the start of the block that are 0 in dataHapRaw (lowers #InfSNPs, keeps Frac at 1)
makeHaps <- function(Blocks) {
  dataHap <- NULL
  dataHapRaw <- NULL
  for (chr in names(Blocks)) {
    Hap <- unlist(lapply(Blocks[[chr]], function(b) rep(b$val, b$len)))
    Raw <- unlist(lapply(Blocks[[chr]], function(b) c(rep(0, b$zero), rep(b$val, b$len - b$zero))))
    n <- length(Hap)
    Meta <- data.frame(Names = paste0("rs", chr, "_", seq_len(n)), Chr = chr, Position = seq_len(n) * 100)
    dataHap <- rbind(dataHap, cbind(Meta, E1_Pat = Hap, E1_Mat = Hap))
    dataHapRaw <- rbind(dataHapRaw, cbind(Meta, E1_Pat = Raw, E1_Mat = Raw))
  }
  list(dataHap = dataHap, dataHapRaw = dataHapRaw)
}

makeParScore <- function(Par, Chroms = c("1", "2", "3")) {
  m <- matrix(1, length(Chroms), 6, dimnames = list(Chroms, c("Pat", "Mat", "PatMat", "MatPat", "CallRate", "Par")))
  m[, "Par"] <- Par
  list(E1 = m)
}

blk <- function(len, val, zero = 0) list(len = len, val = val, zero = zero)

# Paternal chr 1 and 3, maternal chr 2; blocks of distinct length, #InfSNPs not proportional to Length
HappyBlocks <- list(
  "1" = list(blk(400, 1, 10), blk(600, 2, 120)),
  "2" = list(blk(500, 1, 60), blk(700, 2, 5)),
  "3" = list(blk(350, 1, 40), blk(650, 2, 30))
)
HappyPar <- c(11, 22, 11)

run <- function(Haps, ParScore, Chroms) {
  out <- NULL
  capture.output(out <- Haplarithmisis_intphappropser3(Haps, ParScore, "E1", Chroms))
  out
}

sortRes <- function(d) {
  d <- d[order(d$Chr, d$Start), ]
  rownames(d) <- NULL
  d
}

test_that("happy path gives finite fit and probabilities", {
  Res <- run(makeHaps(HappyBlocks), makeParScore(HappyPar), c("1", "2", "3"))$DataCompl
  expect_true(is.data.frame(Res))
  expect_equal(colnames(Res), ResCols)
  expect_equal(nrow(Res), 6)
  expect_true(all(is.finite(Res$Fitted)))
  expect_true(all(is.finite(Res$Resid)))
  expect_true(all(Res$Prob >= 0.5 & Res$Prob <= 1))
})

test_that("Chroms in a different order than the ParScore rows gives the same result", {
  Haps <- makeHaps(HappyBlocks)
  Ref <- run(Haps, makeParScore(HappyPar), c("1", "2", "3"))
  Rev <- run(Haps, makeParScore(HappyPar), c("3", "1", "2"))
  expect_equal(sortRes(Rev$DataCompl), sortRes(Ref$DataCompl))
  expect_equal(Rev$dataHap, Ref$dataHap)
})

test_that("a CumLengths object in the global environment is not picked up", {
  Haps <- makeHaps(HappyBlocks)
  Ref <- run(Haps, makeParScore(HappyPar), c("1", "2", "3"))
  assign("CumLengths", data.frame(junk = 1), envir = globalenv())
  on.exit(rm("CumLengths", envir = globalenv()))
  Res <- run(Haps, makeParScore(HappyPar), c("3", "1", "2"))
  expect_equal(sortRes(Res$DataCompl), sortRes(Ref$DataCompl))
  expect_true(exists("CumLengths", envir = globalenv()))
})

test_that("every Par NA gives a zero-row result with the standard columns", {
  Res <- NULL
  expect_no_error(Res <- run(makeHaps(HappyBlocks), makeParScore(c(NA, NA, NA)), c("1", "2", "3")))
  expect_true(is.data.frame(Res$DataCompl))
  expect_equal(nrow(Res$DataCompl), 0)
  expect_equal(colnames(Res$DataCompl), ResCols)
  expect_true(all(StdCols %in% colnames(Res$DataCompl)))
})

test_that("no blocks passing the Length filter gives an empty result", {
  # Blocks of 100 SNPs are all shorter than MinLength (300)
  Short <- lapply(c("1" = 1, "2" = 2, "3" = 3), function(x) rep(list(blk(100, 1), blk(100, 2)), 4))
  Res <- NULL
  expect_no_warning(expect_no_error(Res <- run(makeHaps(Short), makeParScore(HappyPar), c("1", "2", "3"))))
  expect_equal(nrow(Res$DataCompl), 0)
  expect_equal(colnames(Res$DataCompl), ResCols)
})

test_that("no blocks passing the Frac filter gives an empty result", {
  # Every SNP in dataHapRaw is discordant with the haplotype, so Frac is 0
  Haps <- makeHaps(HappyBlocks)
  for (s in c("E1_Pat", "E1_Mat")) {
    Haps$dataHapRaw[, s] <- ifelse(Haps$dataHap[, s] == 1, 2, 1)
  }
  Res <- NULL
  expect_no_error(Res <- run(Haps, makeParScore(HappyPar), c("1", "2", "3")))
  expect_equal(nrow(Res$DataCompl), 0)
  expect_equal(colnames(Res$DataCompl), ResCols)
})

test_that("a singular rlm fit gives NaN with a warning", {
  # All qualifying blocks have the same Length
  Same <- lapply(c("1" = 1, "2" = 2, "3" = 3), function(x) list(blk(500, 1, 10 * x), blk(500, 2, 20 * x)))
  Res <- NULL
  expect_warning(Res <- run(makeHaps(Same), makeParScore(HappyPar), c("1", "2", "3")), "singular")
  expect_equal(nrow(Res$DataCompl), 6)
  expect_equal(colnames(Res$DataCompl), ResCols)
  expect_true(all(is.nan(Res$DataCompl$Fitted)))
  expect_true(all(is.nan(Res$DataCompl$Resid)))
  expect_true(all(is.nan(Res$DataCompl$Prob)))
})

test_that("a single qualifying block gives NaN with a warning", {
  One <- list("1" = list(blk(400, 1, 10), blk(100, 2)), "2" = list(blk(100, 1), blk(100, 2)), "3" = list(blk(100, 1), blk(100, 2)))
  Res <- NULL
  expect_warning(Res <- run(makeHaps(One), makeParScore(HappyPar), c("1", "2", "3")))
  expect_equal(nrow(Res$DataCompl), 1)
  expect_true(is.nan(Res$DataCompl$Prob))
})
