Haplarithmisis_qcgtype <- function(script,colID,Father,Mother, Gtypes, ChrPos, EmbryoID, outPath,Chroms) {
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_qcgtype_functions.R",sep="/"))
  print(paste("Chr-specific QC analysis for family"))
  QC <- vector("list", 4)

  print("#1 ==> Call-rate computation...")
  for (ind in colID) {
    Gtype <- as.character(Gtypes[, ind])
    CallRate <- Haplarithmisis_callrate(Gtypes[, ind])

    for (chr in Chroms) {

      CallRateChr <- Haplarithmisis_callrate(Gtypes[ChrPos$Chr == chr, ind])

      if (chr == Chroms[1]) {
        CallRateChrs <- CallRateChr[1]
        CallRateChrsHet <- CallRateChr[4]
        CallRateChrsHom <- CallRateChr[5]
      } else {
        CallRateChrs <- cbind(CallRateChrs, CallRateChr[1])
        CallRateChrsHet <- cbind(CallRateChrsHet, CallRateChr[4])
        CallRateChrsHom <- cbind(CallRateChrsHom, CallRateChr[5])
      }

    } #end chr loop

    CallRateChrs <- cbind(CallRateChrs, CallRate[1])
    CallRateChrsHet <- cbind(CallRateChrsHet, CallRate[4])
    CallRateChrsHom <- cbind(CallRateChrsHom, CallRate[5])

    if (ind == colID[1]) {
      CallRates <- CallRate
      CallRateChrsInd <- CallRateChrs
      CallRateChrsHetInd <- CallRateChrsHet
      CallRateChrsHomInd <- CallRateChrsHom
    } else {
      CallRates <- cbind(CallRates, CallRate)
      CallRateChrsInd <- rbind(CallRateChrsInd, CallRateChrs)
      CallRateChrsHetInd <- rbind(CallRateChrsHetInd, CallRateChrsHet)
      CallRateChrsHomInd <- rbind(CallRateChrsHomInd, CallRateChrsHom)

    }
    print(ind)
  } #end ind loop

  rownames(CallRateChrsInd) <- colID
  rownames(CallRateChrsHetInd) <- colID
  rownames(CallRateChrsHomInd) <- colID

  colnames(CallRateChrsInd) <- c(Chroms, "Genome")
  colnames(CallRateChrsHetInd) <- c(Chroms, "Genome")
  colnames(CallRateChrsHomInd) <- c(Chroms, "Genome")

  print("#2 ==> Mendelian inconsistency computation...")

  Child <- Gtypes[, EmbryoID]
  #Genome-wide mendelina inconsistency rate for autosomes
  MendIncRateAut <- Haplarithmisis_mendinc(Gtypes[,Father][Gtypes$Chr != "X" |
                                                         Gtypes$Chr == "XY" |
                                                         Gtypes$Chr == "Y"], Gtypes[,Mother][Gtypes$Chr != "X" |
                                                                                             Gtypes$Chr == "XY" |
                                                                                             Gtypes$Chr == "Y"], Child[Gtypes$Chr != "X" |
                                                                                                                         Gtypes$Chr == "XY" |
                                                                                                                         Gtypes$Chr == "Y"])

  #Chr-specific mendelina inconsistency rate for autosomes
  for (chr in Chroms) {
    MendIncRateChr <- Haplarithmisis_mendinc(Gtypes[,Father][Gtypes$Chr == chr], Gtypes[,Mother][Gtypes$Chr == chr], Child[Gtypes$Chr == chr])
    if (chr == Chroms[1]) {
        MendIncRateChrs <- MendIncRateChr
    } else {
      MendIncRateChrs <- cbind(MendIncRateChrs, MendIncRateChr)
    }
  } #end chr loop

  MendIncRateChrs <- cbind(MendIncRateChrs, MendIncRateAut)
  MendIncRateChrsInd <- MendIncRateChrs
  print(EmbryoID)
  MendIncRateChrsInd <- as.matrix(MendIncRateChrsInd)

  rownames(MendIncRateChrsInd) <- EmbryoID
  colnames(MendIncRateChrsInd) <- c(Chroms, "GenomeAut")

  QC[[1]] <- CallRateChrsInd
  QC[[2]] <- CallRateChrsHetInd
  QC[[3]] <- CallRateChrsHomInd
  QC[[4]] <- MendIncRateChrsInd
  names(QC) <- c("CallRateChrsInd", "CallRateChrsHetInd", "CallRateChrsHomInd", "MendIncRateChrsInd")
  write.table(MendIncRateChrsInd, paste(outPath, paste(EmbryoID,"ChrSpec_MendInc.txt",sep="_"), sep = "/"), sep = "\t", quote = F, col.names = T, row.names = T)
  if (any(MendIncRateChrsInd > 15)) {
    write.table(MendIncRateChrsInd, paste(outPath, paste(EmbryoID,"Chr_with_MendalianInconsistency_over_15_percent.txt",sep="_"), sep = "/"), sep = "\t", quote = F, col.names = T, row.names = T)
  }

  return(QC)
}