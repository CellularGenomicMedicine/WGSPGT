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
  
  # Fill NAs with the preceding value (leading NAs take the first non-NA value); an all-NA vector is returned unchanged
  fillNA <- function(x){
    if (all(is.na(x))) return(x)
    x <- zoo::na.locf(x, na.rm = FALSE)
    zoo::na.locf(x, fromLast = TRUE)
  }
  
  # Filter out the Y-chromosome from the Phased B-allele-frequencies (BAF)
  PhBAFseg <- PhBAFseg[PhBAFseg[,"Chr"] != "Y",]
  
  print("PCF segmentation is applying...")
  
  # Initialize variables SegGenomes and SegGenome
  SegGenomes <- NULL
  SegGenome <- NULL
  
  # Calculate the standard deviation of the BAF data
  PhBAFAll <- na.omit(PhBAFseg[,colnames(PhBAFseg) %in% EmbryoID])
  sdevGenome <- if (length(PhBAFAll) > 0) getMad(PhBAFAll, k = plateau) else NA_real_
  
  # if the calculated sd is zero, it computes sdevGenome as the mean of the standard deviation for each chromosome separately
  if (is.na(sdevGenome) || sdevGenome == 0) {
    sdevCat <- NULL;
    for(chr in unique(PhBAFseg[,"Chr"])){
      # For each chromosome in the BAF data assign to PhBAFChr variable
      PhBAFChr <- PhBAFseg[as.character(PhBAFseg$Chr) == chr, EmbryoID];
      # replace the NAs in the BAF data per chromosome with the preceding BAF value.
      PhBAFChr <- fillNA(PhBAFChr)
      # an all-NA chromosome has no usable scale
      if (all(is.na(PhBAFChr))) next
      
      # Calculate the standard deviation of the adjusted BAF data
      sdev <- getMad(PhBAFChr, k = plateau);
      sdevCat <- c(sdevCat, sdev)
    }
    # no chromosome with a usable scale gives NA (not NaN)
    sdevGenome <- if (any(!is.na(sdevCat))) mean(sdevCat, na.rm = TRUE) else NA_real_
  } # end if statement sdevGenome == 0
  
  print(paste("sdevGenome", sdevGenome, sep = "_"))
  print(paste("gamma", Gamma_value, sep = "_"))
  
  for(chr in unique(PhBAFseg[, "Chr"])){
    # For each chromosome in the BAF data assign to PhBAFChr
    PhBAFChr <- PhBAFseg[as.character(PhBAFseg$Chr) == chr, EmbryoID]
    # replace the NAs in the BAF data per chromosome with the preceding PhBAF value.
    PhBAFChr <- fillNA(PhBAFChr)
    
    # an all-NA chromosome cannot be segmented: raw values are passed through
    if (all(is.na(PhBAFChr))) {
      print(paste(chr, "PhBAFChr_allNA", length(PhBAFChr), sep = "_"))
      SegGenome <- c(SegGenome, PhBAFChr)
      next
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
    # (also when there is no usable scale, i.e. sdev is still NA)
    if (length(PhBAFChr) < 100 || is.na(sdev)){
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
