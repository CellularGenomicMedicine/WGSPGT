# Haplarithmisis_callpcfBAF: A function to perform PCF segmentation on phased B-allele frequencies
# Inputs:
# - script: THe path to the R script directory
# - PhBAFseg: Phased B-allele frequencies to be segmented
# - Gamma_value: Parameter for discontinuity adjustment
# - plateau: Parameter for the MAD calculation
# - EmbryoID: Identifier for the embryo

Haplarithmisis_callpcfBAF <- function(script,PhBAFseg,Gamma_value,plateau,EmbryoID){
  # Load R script to calculate the Median Absolute Deviation (MAD)
  source(paste(script, "analyses", "functions", "getMad.R", sep = "/"))
  # Load R script to perform Piecewise Constant Fitting (PCF)
  source(paste(script, "analyses",  "functions","selectFastPcf.R", sep = "/"))
  
  # Filter out the Y-chromosome from the Phased B-allele-frequencies (BAF)
  PhBAFseg <- PhBAFseg[PhBAFseg[,"Chr"] != "Y",]
  
  print("PCF segmentation is applying...")
  
  # Initialize variables SegGenomes and SegGenome
  SegGenomes <- NULL
  SegGenome <- NULL
  
  # Calculate the standard deviation of the BAF data
  sdevGenome <- getMad(na.omit(PhBAFseg[,colnames(PhBAFseg) %in% EmbryoID]), k = plateau)
  
  # if the calculated sd is zero, it computes sdevGenome as the mean of the standard deviation for each chromosome separately
  if (sdevGenome == 0) {
    sdevCat <- NULL;
    for(chr in unique(PhBAFseg[,"Chr"])){
      # For each chromosome in the BAF data assign to PhBAFChr variable
      PhBAFChr <- PhBAFseg[as.character(PhBAFseg$Chr) == chr, EmbryoID];
      # while there are more than 1 NAs in the BAF data per chromosome, assign the preceding BAF value to the missing value.
      while(sum(is.na(PhBAFChr)) >= 1){
        PhBAFChr[which(is.na(PhBAFChr))] <- PhBAFChr[which(is.na(PhBAFChr)) - 1]
      }
      
      # Calculate the standard deviation of the adjusted BAF data
      sdev <- getMad(PhBAFChr, k = plateau);
      sdevCat <- c(sdevCat, sdev)
    }
    sdevGenome <- mean(sdevCat)
  } # end if statement sdevGenome == 0
  
  print(paste("sdevGenome", sdevGenome, sep = "_"))
  print(paste("gamma", Gamma_value, sep = "_"))
  
  for(chr in unique(PhBAFseg[, "Chr"])){
    # For each chromosome in the BAF data assign to PhBAFChr
    PhBAFChr <- PhBAFseg[as.character(PhBAFseg$Chr) == chr, EmbryoID]
    # while there are more than 1 NAs in the BAF data per chromosome, assign the preceding PhBAF value to the missing value.
    while(sum(is.na(PhBAFChr)) >= 1){ 
      PhBAFChr[which(is.na(PhBAFChr))] <- PhBAFChr[which(is.na(PhBAFChr)) - 1]
    }
    
    # Calculate the MAD for the current chromosome
    sdev <- getMad(PhBAFChr, k = plateau)
    
    # Missing values in the BAF data will be replaced by the mean of the standard deviation for each chromosome separately
    if (is.na(sdev) | sdev == 0) {
      print(paste(chr, "sdev", sdev, sep = "_"))
      # assign the sdev to scale the piece-wise constant filtering
      sdev <- sdevGenome
    }
    # Depending on the number of data points per chromosome for the BAF data, the raw values are directly assigned to the segmented chromosome
    if (length(PhBAFChr) < 100){
      SegChr <- PhBAFChr; 
      print(paste(chr, "PhBAFChr", length(PhBAFChr), sep = "_"))
    } else { 
      # Perform piecewise constant fillting (PCF) using selectFastPcf function
      # Gamma_value is a parameter for discontinuity adjusted by the standard deviation (sdev)
      res <- selectFastPcf(PhBAFChr, 3, Gamma_value * sdev, TRUE)
      SegChr <- res$yhat
    }
    SegGenome <- c(SegGenome, SegChr)
  }#end chr loop

  SegGenomes <- cbind(SegGenomes, SegGenome)
  print(paste(EmbryoID, "==> gamma", Gamma_value, "is applied"))

  PhBAFSeg <- data.frame(PhBAFseg[,c("Names","Chr","Position")], SegGenomes, stringsAsFactors = FALSE)
  colnames(PhBAFSeg) <- colnames(PhBAFseg)

  return(PhBAFSeg)
}#end function
