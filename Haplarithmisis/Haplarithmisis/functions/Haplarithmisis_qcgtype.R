# Haplarithmisis_qcgtype - function to perform quality control on genotype data

Haplarithmisis_qcgtype <- function(script, colID, Father, Mother, Gtypes, ChrPos, EmbryoID, outPath,Chroms) {
  #' Source additional function required for quality control 
  #' Haplarithmisis_callrate and Haplarithmisis_mendinc
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_qcgtype_functions.R",sep="/"))
  
  print(paste("Chr-specific QC analysis for family"))
  
  # Initialize a list to store QC results
  QC <- vector("list", 4)

  print("#1 ==> Call-rate computation...")
  
  # Perform QC for each individual sample in column IDs
  for (ind in colID) {
    Gtype <- as.character(Gtypes[, ind])
    
    # Calculate call rates for individual 
    CallRate <- Haplarithmisis_callrate(Gtypes[, ind])

    # Loop though each chromosome
    for (chr in Chroms) {

      CallRateChr <- Haplarithmisis_callrate(Gtypes[ChrPos$Chr == chr, ind])

      # If it is the first chromosome, initialize matrices
      if (chr == Chroms[1]) {
        CallRateChrs <- CallRateChr[1]
        CallRateChrsHet <- CallRateChr[4]
        CallRateChrsHom <- CallRateChr[5]
      } else {
        # If not, concatenate to existing matrices
        CallRateChrs <- cbind(CallRateChrs, CallRateChr[1])
        CallRateChrsHet <- cbind(CallRateChrsHet, CallRateChr[4])
        CallRateChrsHom <- cbind(CallRateChrsHom, CallRateChr[5])
      }

    } # End of chromosome loop

    # Concatenate overall call rates to the matrices
    CallRateChrs <- cbind(CallRateChrs, CallRate[1])
    CallRateChrsHet <- cbind(CallRateChrsHet, CallRate[4])
    CallRateChrsHom <- cbind(CallRateChrsHom, CallRate[5])

    # If it is the first individual, initialize matrices
    if (ind == colID[1]) {
      CallRates <- CallRate
      CallRateChrsInd <- CallRateChrs
      CallRateChrsHetInd <- CallRateChrsHet
      CallRateChrsHomInd <- CallRateChrsHom
    } else {
      # If not, concatenate to existing matrices
      CallRates <- cbind(CallRates, CallRate)
      CallRateChrsInd <- rbind(CallRateChrsInd, CallRateChrs)
      CallRateChrsHetInd <- rbind(CallRateChrsHetInd, CallRateChrsHet)
      CallRateChrsHomInd <- rbind(CallRateChrsHomInd, CallRateChrsHom)

    }
    print(ind)
  } #end of individual sample loop

  # Set row names for call rate matrices
  rownames(CallRateChrsInd) <- colID
  rownames(CallRateChrsHetInd) <- colID
  rownames(CallRateChrsHomInd) <- colID
  
  # Set column names for call rate matrices
  colnames(CallRateChrsInd) <- c(Chroms, "Genome")
  colnames(CallRateChrsHetInd) <- c(Chroms, "Genome")
  colnames(CallRateChrsHomInd) <- c(Chroms, "Genome")

  # Perform Mendelian inconsistency computation
  print("#2 ==> Mendelian inconsistency computation...")

  # Get child genotype data
  Child <- Gtypes[, EmbryoID]
  
  # Compute genome-wide Mendelian inconsistency rate for autosomes
  MendIncRateAut <- Haplarithmisis_mendinc(Gtypes[,Father][Gtypes$Chr != "X" | Gtypes$Chr == "XY" | Gtypes$Chr == "Y"], 
                                           Gtypes[,Mother][Gtypes$Chr != "X" | Gtypes$Chr == "XY" | Gtypes$Chr == "Y"], 
                                           Child[Gtypes$Chr != "X" | Gtypes$Chr == "XY" | Gtypes$Chr == "Y"])

  # Compute chromosome-specific Mendelian inconsistency rate for autosomes
  for (chr in Chroms) {
    MendIncRateChr <- Haplarithmisis_mendinc(Gtypes[,Father][Gtypes$Chr == chr], 
                                             Gtypes[,Mother][Gtypes$Chr == chr], 
                                             Child[Gtypes$Chr == chr])
    # If it is the first chromosome, initialize matrices
    if (chr == Chroms[1]) {
        MendIncRateChrs <- MendIncRateChr
    } else {
      # If not, concatenate to existing matrix
      MendIncRateChrs <- cbind(MendIncRateChrs, MendIncRateChr)
    }
  } # End of chromosome loop

  # Concatenate genome-wide Mendelian inconsistency rate to chromosome-specific rates
  MendIncRateChrs <- cbind(MendIncRateChrs, MendIncRateAut)
  MendIncRateChrsInd <- MendIncRateChrs
  print(EmbryoID)
  MendIncRateChrsInd <- as.matrix(MendIncRateChrsInd)

  # Set row and column names for Mendelian inconsistency matrix
  rownames(MendIncRateChrsInd) <- EmbryoID
  colnames(MendIncRateChrsInd) <- c(Chroms, "GenomeAut")

  # Store QC results in the list
  QC[[1]] <- CallRateChrsInd
  QC[[2]] <- CallRateChrsHetInd
  QC[[3]] <- CallRateChrsHomInd
  QC[[4]] <- MendIncRateChrsInd
  names(QC) <- c("CallRateChrsInd", "CallRateChrsHetInd", "CallRateChrsHomInd", "MendIncRateChrsInd")
  
  # Write Call rate per chromosome matrix to file
  write.table(CallRateChrsInd, paste(outPath, paste(EmbryoID, "ChrSpec_CallRate.txt", sep = "_"), sep = "/"), 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # Write Heterozygous Call rate per chromosome matrix to file
  write.table(CallRateChrsHetInd, paste(outPath, paste(EmbryoID, "ChrSpec_HetCallRate.txt", sep = "_"), sep = "/"), 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # Write Homozygous Call rate per chromosome matrix to file
  write.table(CallRateChrsHomInd, paste(outPath, paste(EmbryoID, "ChrSpec_HomCallRate.txt", sep = "_"), sep = "/"), 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # Write Mendelian inconsistency matrix to file
  write.table(MendIncRateChrsInd, paste(outPath, paste(EmbryoID, "ChrSpec_MendInc.txt", sep = "_"), sep = "/"), 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # Check for chromosomes with Mendelian inconsistency rate over 15% and write to file
  if (any(MendIncRateChrsInd > 15)) {
    write.table(MendIncRateChrsInd, paste(outPath, paste(EmbryoID, "Chr_with_MendalianInconsistency_over_15_percent.txt", sep = "_"), sep = "/"), 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
  }

  # Return the QC results list
  return(QC)
}