#' Haplarithmisis_ParScore - Function to compute parental scores

Haplarithmisis_ParScore <- function(Father, Mother, dataPo, QC, Chroms, Gtypes, EmbryoID){
	 
  # Threshold parameters for computing parental scores
	AvgMDA <- 2.2862648
	SdMDA <- 0.4721314
	MinCallRate <- 16.38

	print("Computing parental scores...")
	
	# Extract call rates from QC data
	dataCallRate <- QC[["CallRateChrsInd"]]
	
	CallRates <- dataCallRate[, Chroms]
	
	# Extract parental origin data
	Pos <- dataPo[, EmbryoID]
	Pos <- as.matrix(Pos)
	colnames(Pos) <- EmbryoID
	SamplesPo <- vector("list", ncol(Pos))
	names(SamplesPo) <- EmbryoID
	
	# Initialize matrix to store parental scores	
	PoScore <- matrix(NA, length(Chroms), 6)
	rownames(PoScore) <- Chroms
	
	# Compute parental scores for each chromosome
	for(chr in Chroms){
	  # Extract genotypes of parents for the current chromosome
		PatGtype <- Gtypes[Gtypes[,"Chr"] == chr, Father]
		MatGtype <- Gtypes[Gtypes[,"Chr"] == chr, Mother]
		
		# Count informative parent of origin markers
		InfPo <- sum((MatGtype == "AA" & PatGtype == "BB") |
							    (MatGtype == "BB" & PatGtype == "AA") | 
							    (MatGtype == "AB" & PatGtype == "AA") | 
							    (MatGtype == "AA" & PatGtype == "AB") | 
							    (MatGtype == "BB" & PatGtype == "AB") | 
							    (MatGtype == "AB" & PatGtype == "BB"))
		
		# Extract parent of origin data for the current chromosome
		dataPoChr <- Pos[dataPo[, "Chr"] == chr, EmbryoID]
		
		# Calculate paternal and maternal scores
		# - If dataPo is lower than 0, the parent of origin was paternal
		# - If dataPo is higher than 0, the parent of origin was maternal
		Pat <- (abs(sum(dataPoChr[dataPoChr < 0])) / InfPo) * 100
		Mat <- (sum(dataPoChr[dataPoChr > 0]) / InfPo) * 100
		
		# Avoid division by zero
		if (Pat == 0) { Pat = 0.0001}
		if (Mat == 0) { Mat = 0.0001}
	
		# Calculate paternal/maternal and maternal/paternal ratios
		PatMat <- Pat / (Mat + Pat)
		MatPat <- Mat / (Mat + Pat)
		
		# Adjust parental scores based on thresholds
		if(Pat <= (AvgMDA - (3 * SdMDA)) & Mat <= (AvgMDA - (3 * SdMDA))) {
		  PatMat <- 0
		  MatPat <- 0
		  }
		
		# Store computed scores in PoScore matrix
		PoScore[chr, ] <- cbind(Pat, Mat, PatMat, MatPat, CallRates[EmbryoID, chr], NA)
	}
	
	# Set column names of PoScore matrix
	colnames(PoScore) <- c("Pat", "Mat", "PatMat", "MatPat", "CallRate", "Par")
	
	# Store PoScore matrix in SamplesPo list
	SamplesPo[[EmbryoID]] <- PoScore
	
	# Print the ID of the embryo being processed
	print(EmbryoID)
	
	# initialize list to store modified parental scores
	SamplesPo2 <- vector("list", ncol(Pos))
	names(SamplesPo2) <- colnames(Pos)

	# Modify parental scores based on thresholds for each individual embryo
	for(ind in names(SamplesPo)){

	  # Extract scores for the current embryo
		Sample <- SamplesPo[[ind]]
		
		# Adjust parental scores based on thresholds
		# 1: paternal haplotype
		# 2: maternal haplotype
		
		#' Condition 1: Call rate is greater than the MinCallRate and either both parental scores are below their average plus standard deviation,
		#' or both Paternal/Maternal and Maternal/Paternal score ratios are below 0.9.
		Sample[Sample[, "CallRate"] > MinCallRate & 
		         ((Sample[, "Pat"] < (AvgMDA + SdMDA) & Sample[, "Mat"] < (AvgMDA + SdMDA)) | 
		            (Sample[, "PatMat"] < 0.9 & Sample[, "MatPat"] < 0.9)),"Par"] <- 12 # 12: paternal and maternal haplotypes
		
		
		#' Condition 2: Paternal score is above average plus standard deviation,
		#' call rate is above minimum, and paternal/maternal ratio is at least 0.9
		Sample[Sample[, "Pat"] > (AvgMDA + SdMDA) & 
		         Sample[, "CallRate"] > MinCallRate & 
		         Sample[, "PatMat"] >= 0.9, "Par"] <- 11 # 11: two paternal haplotypes
		
		#' Condition 3: Maternal score is above average plus standard deviation,
		#' call rate is above minimum, and maternal/paternal ratio is at least 0.9
		Sample[Sample[, "Mat"] > (AvgMDA + SdMDA) & 
		         Sample[, "CallRate"] > MinCallRate & 
		         Sample[, "MatPat"] >= 0.9, "Par"] <- 22 # two maternal haplotypes
		
		#' Condition 4: Call rate is below minimum and both paternal/maternal 
		#' and maternal/paternal ratios are below 0.5
		Sample[Sample[,"CallRate"] < MinCallRate & 
		         (Sample[,"PatMat"] < 0.5 & Sample[, "MatPat"] < 0.5), "Par"] <- 00 # the ParScore could not be determined
 		
		# Store modified scores in SamplesPo2 list per individual embryo
		SamplesPo2[[ind]] <- Sample

	}# end ind loop
	
# Return the list of modified parental scores
return(SamplesPo2)

}#end function