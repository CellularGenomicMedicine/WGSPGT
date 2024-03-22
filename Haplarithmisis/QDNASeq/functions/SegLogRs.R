#' SegLogRs - Function to 
#' 
#' Input:
#' - script: path to script directory
#' - logRs: dataframe containing logR values
#' - gammaMC: gamma value to be used for non-embryos 
#' - gammaSC: gamma value to be used for embryos
#' - plateau: window value to calculate the median absolute deviation (MAD)
#' - outPath: path to output segmented logR values
#' - fam_members_Embryo: dataframe containing information about family members / samples 

#' Output:
#' - Segmented LogR dataframe

SegLogRs <- function(script, logRs, gammaMC, gammaSC, plateau, outPath, fam_members_Embryo){

  # Print message indicating the start of PCF segmentation
  print("PCF segmentation is applying...")

  # Initialize variables
	SegGenomes <- NULL
	segLogRs <- NULL
	
	# Load required functions
	source(paste(script, "analyses", "functions", "getMad.R", sep = "/"))
	source(paste(script, "analyses", "functions", "selectFastPcf.R", sep = "/"))
	
	# Iterate over each column starting from the column with the logR data
	for(i in 4:ncol(logRs)){
	  # Determine the segmentation value based on whether the sample is Embryo (gammaSC) or Bulk: parent or referent (gammaMC)
		if(any(colnames(logRs)[i]) %in% fam_members_Embryo[, "SampleID"]) {
		  Gamma = gammaSC 
		} else { 
		  Gamma = gammaMC
		}
	  
	  SegGenome <- NULL
	  
		for(chr in unique(logRs[, "Chr"])){
			logRChr <- logRs[as.character(logRs$Chr) == chr, i]
			# Replace NA values in logRChr with previous non-NA values
			while(sum(is.na(logRChr)) >= 1){
			  logRChr[which(is.na(logRChr))] <- logRChr[which(is.na(logRChr)) - 1]
			}
			
			# Calculate Median Absolute Deviation (MAD)
			sdev <- getMad(logRChr, k = plateau)
			
			# Check conditions for segmenting chromosome
			if ( is.na(sdev) | sdev == 0 | length(logRChr) < 100 ) {
				SegChr <- logRChr
				# Save segment with less than 100 logRseg values to a file
				write.table(logRChr, paste(paste(outPath, "", sep = "/"), "Chr", chr, "_for_", colnames(logRs)[i], "_has_less_then_100_logRseg_values", ".txt", sep = ""), 
				            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
			} else {
			  # Apply PCF segmentation algorithm
        res <- selectFastPcf(logRChr, 3, Gamma * sdev, T)
        SegChr <- res$yhat
			}
			
			SegGenome <- c(SegGenome, SegChr)
		} # end chromosome loop
	  
	  # Append SegGenome to SegGenomes matrix
		SegGenomes <- cbind(SegGenomes, SegGenome)
		
		# Print message indicating completion of segmentation for a column
		print(paste(colnames(logRs)[i], "==> gamma", Gamma, "is applied"))
	} # end column loop

	# Create a data frame containing segmented logRs data
	segLogRs <- data.frame(logRs[, c("Names", "Chr", "Position")], SegGenomes, stringsAsFactors = FALSE, check.names = FALSE)
	colnames(segLogRs) <- colnames(logRs)

	# Return the segmented logRs data frame
	return(segLogRs)
	
} # end function

