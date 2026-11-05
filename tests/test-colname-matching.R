# Tests for exact embryo column-name matching
# Run from the repo root: Rscript -e 'testthat::test_file("tests/test-colname-matching.R")'

library(testthat)
library(MASS)

source(test_path("..", "Haplarithmisis", "Haplarithmisis", "functions", "Haplarithmisis_intphappropser3.R"))

# Embryo IDs that act as regular expressions or as substrings of other columns
Ids <- c("E.1", "E(1", "Pat1", "Pat", "Names", "Position")

# Meta columns plus, for each ID, its own _Pat/_Mat pair and decoys that a regex match would also hit
makeDecoyHaps <- function(Id) {
  Blocks <- list("1" = list(c(400, 1, 10), c(600, 2, 120)),
                 "2" = list(c(500, 1, 60), c(700, 2, 5)),
                 "3" = list(c(350, 1, 40), c(650, 2, 30)))
  Decoys <- c("E11", "E.10", "E(10", "Pat1_Mat", "Pat12", "Xnames", "Xposition")
  Cols <- c(paste0(Id, c("_Pat", "_Mat")), as.vector(outer(Decoys, c("_Pat", "_Mat"), paste0)))
  Raw <- NULL
  for (chr in names(Blocks)) {
    Hap <- unlist(lapply(Blocks[[chr]], function(b) rep(b[2], b[1])))
    Zer <- unlist(lapply(Blocks[[chr]], function(b) c(rep(0, b[3]), rep(b[2], b[1] - b[3]))))
    n <- length(Hap)
    Meta <- data.frame(Names = paste0("rs", chr, "_", seq_len(n)), Chr = chr, Position = seq_len(n) * 100)
    Dat <- as.data.frame(setNames(rep(list(Zer), length(Cols)), Cols), check.names = FALSE)
    Raw <- rbind(Raw, cbind(Meta, Dat))
  }
  Hp <- Raw
  for (cn in Cols) Hp[[cn]] <- ifelse(Raw[[cn]] == 0, 1, Raw[[cn]])
  list(dataHap = Hp, dataHapRaw = Raw)
}

test_that("intphappropser3 selects only the embryo's own _Pat/_Mat columns", {
  for (Id in Ids) {
    Haps <- makeDecoyHaps(Id)
    Par <- matrix(1, 3, 6, dimnames = list(c("1", "2", "3"), c("Pat", "Mat", "PatMat", "MatPat", "CallRate", "Par")))
    Par[, "Par"] <- c(11, 22, 11)
    ParScore <- setNames(list(Par), Id)
    out <- NULL
    capture.output(out <- suppressWarnings(
      Haplarithmisis_intphappropser3(Haps, ParScore, Id, c("1", "2", "3"))))
    expect_gt(nrow(out$DataCompl), 0)
    Samples <- unique(sub("\\.[0-9]+$", "", rownames(out$DataCompl)))
    expect_true(all(Samples %in% paste0(Id, c("_Pat", "_Mat"))), info = Id)
  }
})

# Same selection expression as Haplarithmisis_inthapnew1 L14 on a tiny table
test_that("inthapnew1 zeroing hits only the embryo's own columns", {
  for (Id in Ids) {
    Dat <- makeDecoyHaps(Id)$dataHapRaw
    Int <- data.frame(Chr = "1", S = 0, E = 1e9, ID = Id, stringsAsFactors = FALSE)
    Dat[Dat$Chr == Int[1, 1] & as.numeric(Dat$Position) >= as.numeric(Int[1, 2]) & as.numeric(Dat$Position) <= as.numeric(Int[1, 3]),
        colnames(Dat) %in% paste0(Int[1, 4], c("_Pat", "_Mat"))] <- 0
    Zeroed <- colnames(Dat)[sapply(Dat, function(x) all(x[Dat$Chr == "1"] == 0))]
    expect_equal(sort(Zeroed), sort(paste0(Id, c("_Pat", "_Mat"))), info = Id)
  }
})
