# Haplarithmisis_chrxhtyping_Grandparents - Function for haplotyping chromosome X using grandparents as reference

Haplarithmisis_chrxhtyping_Grandparents <- function(Father, Mother, Grandfather, Grandmother, Gtypes, parent, ScSexes, EmbryoID){
  # Subset genetic data for chromosome X
  GtypesChrX <- Gtypes[Gtypes$Chr == "X", ]
	
  # Extract relevant columns for each family member on chromosome X
  MotherChrX <- GtypesChrX[, Mother]
	FatherChrX <- GtypesChrX[, Father]
	GrandfatherChrX <- GtypesChrX[, Grandfather]
	GrandmotherChrX <- GtypesChrX[, Grandmother]

	# Initialize matrix to store maternal haplotypes for chromosome X
	HapXMat <- matrix(0, nrow(GtypesChrX), ncol(as.matrix(GtypesChrX[, which(colnames(Gtypes) %in% EmbryoID)])))
	colnames(HapXMat) <- EmbryoID

	# Handle haplotype assignment based on the indication parent (Mother)
	if(parent == "Mother"){
	  # Assign "BA" haplotype when the Grandfather's genotype is BB and the Mother's genotype is "AB"
		MotherChrX[(GtypesChrX[, Grandfather] == "BB" & GtypesChrX[, Mother] == "AB")] <- "BA"
		#  Assign "NC" (NoCall) genotype in case of errorenous conditions
		MotherChrX[(GtypesChrX[, Grandfather] == "NC" & GtypesChrX[, Mother] == "AB") | 
		           (GtypesChrX[, Grandfather] == "BB" & GtypesChrX[, Grandmother] == "BB" & GtypesChrX[, Mother] == "AB") |
					     (GtypesChrX[, Grandfather] == "AA" & GtypesChrX[, Grandmother] == "AA" & GtypesChrX[, Mother] == "AB") ] <- "NC"
	}	   

	na.omit(ScSexes[ScSexes[, 2] == "male", 1])

	
	# Handle haplotype assignment based on sex (male)
	if (ScSexes[, 2] == "male") {
		print(paste(EmbryoID, "is a male sample but has", sum(GtypesChrX[, EmbryoID] == "AB"), "heterozygous SNP-calls; these are treated as NoCalls"))
	  #  Assign "NC" (NoCall) haplotype to male embryo with heterozygous SNP calls on chromosome X.
	  GtypesChrX[GtypesChrX[, EmbryoID] == "AB", EmbryoID] <- "NC"
		
	  # Assign haplotypes based on specific conditions for males
		HapXMat[(MotherChrX == "AB" & GtypesChrX[, EmbryoID] == "BB") | (MotherChrX == "BA" & GtypesChrX[, EmbryoID] == "AA"), EmbryoID] <- 2
		HapXMat[(MotherChrX == "AB" & GtypesChrX[, EmbryoID] == "AA") | (MotherChrX == "BA" & GtypesChrX[, EmbryoID] == "BB"), EmbryoID] <- 1
		# For positions with NoCalls set haplotype to 0
		HapXMat[GrandfatherChrX == "NC" | GrandmotherChrX == "NC" , EmbryoID] <- 0
		print(EmbryoID)
	}

	# Handle haplotype assignment based on sex (female)
	if(ScSexes[, 2] == "female") {
	  
	  # Assign haplotypes based on specific conditions for females
		HapXMat[(MotherChrX == "AB" & GtypesChrX[, EmbryoID] == "BB") | (MotherChrX == "BA" & GtypesChrX[, EmbryoID] == "AA") |
				    (MotherChrX == "AB" & FatherChrX == "AA" & GtypesChrX[, EmbryoID] == "AB") |
				    (MotherChrX == "BA" & FatherChrX == "BB" & GtypesChrX[, EmbryoID] == "AB"), EmbryoID] <- 2
		
		HapXMat[(MotherChrX == "AB" & GtypesChrX[, EmbryoID] == "AA") |
				    (MotherChrX == "BA" & GtypesChrX[, EmbryoID] == "BB") |
					  (MotherChrX == "BA" & FatherChrX == "AA" & GtypesChrX[, EmbryoID] == "AB") |
					  (MotherChrX == "AB" & FatherChrX == "BB" & GtypesChrX[, EmbryoID] == "AB"), EmbryoID] <- 1
		
		# For positions with NoCalls set haplotype to 0
		HapXMat[GrandfatherChrX == "NC" | GrandmotherChrX == "NC" | FatherChrX == "NC", EmbryoID] <- 0
		print(EmbryoID)
	}
	
	# Return the haplotype matrix for chromosome X
	return(HapXMat)
}#end chrxhtypingOpt1 function