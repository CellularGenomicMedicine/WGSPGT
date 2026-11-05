# Tests for ParScore (zero informative SNPs) and sexdetermination (no chrX row).
# Run from repo root: Rscript tests/test-parscore-chrx.R   (base R only)
fdir <- file.path("Haplarithmisis", "Haplarithmisis", "functions")
source(file.path(fdir, "Haplarithmisis_ParScore.R"))
source(file.path(fdir, "Haplarithmisis_sexdetermination.R"))

fails <- 0
check <- function(name, ok) {
  if (isTRUE(ok)) cat("PASS:", name, "\n") else { cat("FAIL:", name, "\n"); fails <<- fails + 1 }
}

# Chromosomes "1" and "2" have informative SNPs; "10" and "X" have no rows at all
Gtypes <- data.frame(
  Names = paste0("s", 1:4), Chr = c("1", "1", "2", "2"), Position = 1:4,
  Father = c("AA", "AA", "AB", "AB"), Mother = c("BB", "BB", "AB", "AB"),
  stringsAsFactors = FALSE)
Chroms <- c("1", "2", "10", "X")
dataPo <- data.frame(Names = Gtypes$Names, Chr = Gtypes$Chr, Position = Gtypes$Position,
                     E1 = c(-1, 1, 0, 0), stringsAsFactors = FALSE)
# more than one sample row: ParScore indexes CallRates[EmbryoID, chr] on a matrix
CR <- matrix(90, 2, 4, dimnames = list(c("E1", "E2"), Chroms))
QC <- list(CallRateChrsInd = CR)

res <- tryCatch(suppressWarnings(invisible(capture.output(
         out <- Haplarithmisis_ParScore("Father", "Mother", dataPo, QC, Chroms, Gtypes, "E1")))),
       error = function(e) e)
check("ParScore does not error with zero-informative chromosomes", !inherits(res, "error"))
if (!inherits(res, "error")) {
  s <- out[["E1"]]
  check("zero-informative chr 10 has NA scores", all(is.na(s["10", c("Pat", "Mat", "PatMat", "MatPat")])))
  check("zero-informative chr 10 has NA Par", is.na(s["10", "Par"]))
  check("zero-informative chr X has NA Par", is.na(s["X", "Par"]))
  check("informative chr 1 still scored", !is.na(s["1", "Pat"]))
  sx <- matrix(NA, 1, 2, dimnames = list("E1", NULL))
  invisible(capture.output(r <- Haplarithmisis_sexdetermination(sx, out)))
  check("sexdetermination gives NA sex for NA chrX Par", is.na(r["E1", 2]))
}

# No "X" row in the ParScore at all
cn <- c("Pat", "Mat", "PatMat", "MatPat", "CallRate", "Par")
noX <- list(E1 = matrix(c(1, 2, 3, 4, 90, 12), 1, 6, dimnames = list("1", cn)))
sx <- matrix(NA, 1, 2, dimnames = list("E1", NULL))
r <- tryCatch({ invisible(capture.output(o <- Haplarithmisis_sexdetermination(sx, noX))); o },
              error = function(e) e)
check("sexdetermination with no X row does not error", !inherits(r, "error"))
check("sexdetermination with no X row gives NA sex", !inherits(r, "error") && is.na(r["E1", 2]))

# Known X row still works
withX <- list(E1 = rbind(noX$E1, X = c(1, 2, 3, 4, 90, 22)))
invisible(capture.output(r2 <- Haplarithmisis_sexdetermination(sx, withX)))
check("sexdetermination still reports male for Par 22", identical(unname(r2["E1", 2]), "male"))

cat(fails, "failure(s)\n")
quit(status = if (fails) 1 else 0)
