# Haplarithmisis_intphappropser3 - Function for calculating the probabilty scores for the phased haplotypes compared to raw haplotypes to generate interpreted haplotypes

Haplarithmisis_intphappropser3 <- function(Haps, ParScore, EmbryoID, Chroms){
  # Extract required data from Haps list
  dataHap <- Haps[["dataHap"]]
	dataHapRaw <- Haps[["dataHapRaw"]]
	
	# Define metadata columns
	MetaInfoCols <- c("Names", "Chr", "Position")
	MinLength = 300 # Minimum length for considering a block
	MinFrac = 0.6 # Minimum fraction for informative SNPs
	
	dataHap1 <- dataHap # Create a copy of dataHap for processing
	
	# Extract sample names from dataHap1
	Samps <- vector("list", ncol(dataHap1) - length(MetaInfoCols))
	names(Samps) <- colnames(dataHap1[, !names(dataHap1) %in% MetaInfoCols])

	ParScore2 <- ParScore # Backup ParScore data
	
	# Select the chromosomes for which there is a ParScore
	ParScore[[EmbryoID]] <- ParScore[[EmbryoID]][!is.na(ParScore[[EmbryoID]][, "Par"]), ]
	
	# Process ParScore data
	if (sum(nrow(ParScore[[EmbryoID]])) > 1){
		# Determine whether chromsomes originate from paternal and maternal origin
	  PatChr <- rownames(ParScore[[EmbryoID]])[ParScore[[EmbryoID]][, "Par"] == 11]
		MatChr <- rownames(ParScore[[EmbryoID]])[ParScore[[EmbryoID]][, "Par"] == 22]
		NullChr <- rownames(ParScore[[EmbryoID]])[ParScore[[EmbryoID]][, "Par"] == 0 ]
		
		# Update haplotype data based on chromosome origins
		
		# when for the specific chromosome there is paternal origin, the maternal haplotypes are set to 0
		if (length(PatChr) != 0) {
		  for(chr in PatChr) {
		    dataHap1[dataHap1[,"Chr"] == chr, paste(EmbryoID, "_Mat", sep = "")] <- 0
		  }
		}
		
		# when for the specific chromosome there is maternal origin, the paternal haplotypes are set to 0
		if (length(MatChr) != 0) {
		  for(chr in MatChr) {
		    dataHap1[dataHap1[,"Chr"] == chr, paste(EmbryoID, "_Pat", sep = "")] <- 0 
		  }
		}
		
		# when for the specific chromosome there is no origin determined, the paternal and maternal haplotypes are set to 0
		if (length(NullChr) != 0) {
		  for(chr in NullChr) {
		    dataHap1[dataHap1[, "Chr"] == chr, paste(EmbryoID, "_Pat", sep = "")] <- 0
		    dataHap1[dataHap1[, "Chr"] == chr, paste(EmbryoID, "_Mat", sep = "")] <- 0
		  }
		}
		
		print(EmbryoID)
	} else {
		print(paste(EmbryoID,"sample sbould be excluded for diagnosis, as the origin of the chromosomes could not be determined"))
	}

	# Begin probability computation
	ParScore <- ParScore2
	print("--------------------------------------------------------------")
	print("Probability computation...")

	# Extract columns for the current embryo from dataHap1
	Blasts <- colnames(dataHap1)[grepl(EmbryoID, colnames(dataHap1))]
	
	for (samp in Blasts){
	  
		for(chr in Chroms){
			if(!is.na(ParScore[[EmbryoID]][chr, "Par"])){
			  # Filter haplotype to the current chromosome
				dataHapChr <- dataHap1[dataHap1$Chr == chr, ]

				# Compute the run-length encoding of haplotype data for the current chromosome and assign to IntBlock
				IntBlock <- rle(dataHap1[dataHap1$Chr == chr, samp])
				IntBlock <- cbind(IntBlock$values, IntBlock$lengths)
				colnames(IntBlock) <- c("values", "lengths")
				IntBlock <- data.frame(rbind(IntBlock,c(0, 0)))

				# Calculate the cumulative sum of lengths to track haplotype block positions in CumLength
				CumLength <- rbind(cumsum(IntBlock$lengths), IntBlock$values)

				# If IntBlock is empty, initialize Cumlength with 0s
				if(nrow(IntBlock) == 0){
					CumLength <- cbind(0, 0, 0, 0)
				} else {
				  # Initilaize InfSNPs (Informative SNPs), InfSNPsDisc (Informative SNPs discordant), Frac (Fraction), Annot (Annotation)
					InfSNPs <- rep(NA, ncol(CumLength))
					InfSNPsDisc <- rep(NA, ncol(CumLength))
					Frac <- rep(NA, ncol(CumLength))
					Annot <- matrix(NA, ncol(CumLength), 3)

					for (i in 1:ncol(CumLength)){
					  # Determine starting position for each haplotype block
						if(i == 1){
						  start = 1
						} else {
						  start = CumLength[1, i-1] + 1
						}
					  
						if (start > CumLength[1, i]) {
						  start = CumLength[1, i]
						}
					  # Compute the number of informative SNPs, informative SNP discrepancies, and their fraction within each haplotype block
						InfSNPs[i] <- sum(dataHapRaw[dataHap1$Chr == chr, samp][start:CumLength[1, i]] == CumLength[2, i])
						InfSNPsDisc[i] <- sum(dataHapRaw[dataHap1$Chr == chr, samp][start:CumLength[1,i]] != CumLength[2, i] & dataHapRaw[dataHap1$Chr == chr, samp][start:CumLength[1, i]] != 0)
						Frac[i] <- round(InfSNPs[i] / (InfSNPs[i] + InfSNPsDisc[i]), digits = 3)
						# Store annotation data related to the haplotype blocks and informative SNPs
						Annot[i,] <- cbind(chr, dataHapChr$Position[start], dataHapChr$Position[CumLength[1, i]])
						if (Frac[i] == "NaN"){
						  Frac[i] <- 0
						  } # Update dataHap1 to set haplotype values to zero for block with fraction below the minimum threshold (MinFrac)
						if (Frac[i] < MinFrac){
						  dataHap1[dataHap1$Chr == chr, samp][start:CumLength[1, i]] <- 0
						  }
					} # end i loop
					
					# Combine annotation data, block lengths, SNP counts, and fraction information into Cumulative length
					CumLength <- cbind.data.frame(Annot,IntBlock$lengths, t(CumLength), InfSNPs, InfSNPsDisc, Frac)
				}
				
				# Checks chromosome origins based on phased haplotype probabilities in ParScore and updates CumLengths accordingly
				Chroms2 <- rownames(ParScore[[EmbryoID]])
				if (chr == Chroms2[which(ParScore[[EmbryoID]][, "Par"] != "NA")[1]]){
				  CumLengths <- CumLength 
				} else {
				    CumLengths <- rbind(CumLengths, CumLength)
				}
				
			} else {
			  print(paste("The origin of chr.", chr, "of", EmbryoID, "could not be determined"))
			}
		  
		}#end chr loop

	  # If Cumlengths does not exists, create an empty matrix with 2 rows and 9 columns
		if (exists("CumLengths") == FALSE){
		  CumLengths <- matrix(0, 2, 9)
		}
		
	  # Assign column names to the CumLengths
		colnames(CumLengths) <- c("Chr", "Start", "Stop", "Length", "CumLength", "Value", "#InfSNPs", "#InfSNPsDisc", "Frac")
		# Set the Vlaue column to 0 when the Fraction column is less than minimum fraction (MinFrac)
		CumLengths[CumLengths[, "Frac"] < MinFrac, "Value"] <- 0
		# Remove rows from CumLengths where the value is equal to 0.
		CumLengths <- CumLengths[CumLengths[,"Value"] !=0 , ]
		Samps[[samp]] <- CumLengths
		print(samp)
		rm(list = "CumLengths")
	}#end samp loop

	# Process and filter computed data
	MDAs <- do.call(rbind, Samps[grepl(EmbryoID, names(Samps))])

	# Separate the rows with chromsome X into a separate matrix and add them back if they contain Mat in the row names
	X <- MDAs[which(MDAs$Chr == "X"), ]
	MDAs <- rbind(MDAs[MDAs$Chr != "X", ], X[grep("Mat", rownames(X)), ])
	MDAs <- MDAs[MDAs[, "Value"] != 0, ]
	MDAs <- MDAs[MDAs[, "Length"] >= MinLength, ]

	# Convert Chr to character type and Start and Stop columns to numeric 
	MDAs[, "Chr"] <- as.character(MDAs[, "Chr"])
	MDAs[, "Start"] <- as.numeric(as.character(MDAs[, "Start"]))
	MDAs[, "Stop"] <- as.numeric(as.character(MDAs[, "Stop"]))

	# Perform linear regression with Informative SNPs and length as predictor
	LmFit <- rlm(MDAs[,"#InfSNPs"] ~ MDAs[, "Length"], data = MDAs)
	# Combine the original data, fitted values from linear model, residuals and fills the probability column with NA values.
	DataCompl <- cbind(MDAs, LmFit$fitted.values, LmFit$resid, rep(NA, nrow(MDAs)))
	colnames(DataCompl)[(ncol(DataCompl) - 2) : ncol(DataCompl)] <- cbind("Fitted", "Resid", "Prob")
	# Create a subset of DataCompl containing rows where residual is lower than 0.
	NegRes <- subset(DataCompl, Resid < 0)

	# Computes Z-scores for the residuals in DataCompl
	zscores <- scale(DataCompl$Resid, center = TRUE, scale = TRUE)
	datawithzscores <- cbind(DataCompl, zscores)
	DataCompl[,"Prob"] <- round(pnorm(-abs(datawithzscores$zscores)) + 0.5, digits = 3)
	DataCompl[DataCompl[, "Resid"] > 0, "Prob"] <- 1

	# Organize and return interpreted data
	Intp <- vector("list", 2)
	names(Intp) <- c("dataHap", "DataCompl")
	Intp[["dataHap"]] <- dataHap1
	Intp[["DataCompl"]] <- DataCompl
	
	return(Intp)

}

