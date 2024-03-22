# Haplarithmisis_chrxhtyping_Grandmother - Function for haplotyping chromosome X using grandmother as reference

Haplarithmisis_chrxhtyping_Grandmother <- function(Father, Mother, RefSampleID, Gtypes, parent, ScSexes, EmbryoID){
  
  # Subset genetic data for chromosome X
	GtypesChrX <- Gtypes[Gtypes$Chr == "X", ]
	
	# Extract relevant columns for each family member on chromosome X
	MotherChrX <- GtypesChrX[, Mother]
	FatherChrX <- GtypesChrX[, Father]
	GrandmotherChrX <- GtypesChrX[, RefSampleID]

	# Initialize matrix to store maternal haplotypes for chromosome X
	HapXMat <- matrix(0, nrow(GtypesChrX), ncol(as.matrix(GtypesChrX[, colnames(GtypesChrX) %in% EmbryoID])))
	colnames(HapXMat) <- paste(colnames(Gtypes)[colnames(Gtypes) %in% EmbryoID])
	
	# Handle haplotype assignment based on the indication parent (Mother)
	if(parent == "Mother"){
	  # Assign "BA" haplotype when the Grandmother's genotype is BB and the Mother's genotype is "AB"
		MotherChrX[(GtypesChrX[, RefSampleID] == "BB" & GtypesChrX[, Mother] == "AB")] <- "BA"
		MotherChrX[(GtypesChrX[, RefSampleID] == "AA" & GtypesChrX[, Mother] == "AB")] <- "AB"
		#  Assign "NC" (NoCall) genotype in case of errorenous conditions
		MotherChrX[(GtypesChrX[, RefSampleID] == "NC" & GtypesChrX[, Mother] == "AB") |
               (GtypesChrX[, RefSampleID] == "AB" & GtypesChrX[, Mother] == "AB") ] <- "NC"
	}
	
	# Handle haplotype assignment based on sex (male)
	if(ScSexes[,2] == "male"){
		print(paste(EmbryoID,"is a male sample but has",sum(GtypesChrX[,EmbryoID] == "AB"),"heterozygous SNP-calls; these are treated as NoCalls"))
		GtypesChrX[GtypesChrX[, EmbryoID] == "AB", EmbryoID] <-"NC"
	
		# Assign haplotypes based on specific conditions for males
		HapXMat[(MotherChrX=="AB" & GtypesChrX[, EmbryoID] == "BB") | (MotherChrX=="BA" & GtypesChrX[, EmbryoID] == "AA"), EmbryoID] <- 2
		HapXMat[(MotherChrX=="AB" & GtypesChrX[, EmbryoID] == "AA") | (MotherChrX=="BA" & GtypesChrX[, EmbryoID] == "BB"), EmbryoID] <- 1
		
		# For positions with NoCalls set haplotype to 0
		HapXMat[GrandmotherChrX=="NC", EmbryoID] <- 0
		print(EmbryoID)
	}
	
	# Handle haplotype assignment based on sex (female)
	if(ScSexes[,2] == "female"){
	  
	  # Assign haplotypes based on specific conditions for females
		HapXMat[(MotherChrX=="AB" & GtypesChrX[, EmbryoID] == "BB") |
				    (MotherChrX=="BA" & GtypesChrX[, EmbryoID] == "AA") |
				    (MotherChrX=="AB" & FatherChrX == "AA" & GtypesChrX[, EmbryoID] == "AB") |
				    (MotherChrX=="BA" & FatherChrX == "BB" & GtypesChrX[, EmbryoID] == "AB"), EmbryoID] <- 2
		HapXMat[(MotherChrX=="AB" & GtypesChrX[, EmbryoID] == "AA") |
				    (MotherChrX=="BA" & GtypesChrX[, EmbryoID] == "BB") |
					  (MotherChrX=="BA" & FatherChrX == "AA" & GtypesChrX[, EmbryoID] == "AB")|
					  (MotherChrX=="AB" & FatherChrX == "BB" & GtypesChrX[, EmbryoID] == "AB"), EmbryoID] <- 1
		
		# For positions with NoCalls set haplotype to 0
		HapXMat[GrandmotherChrX=="NC" | FatherChrX=="NC", EmbryoID] <- 0
		print(EmbryoID)
 	}#end s loop
	
	# Return the haplotype matrix for chromosome X
	return(HapXMat)

}#end chrxhtypingOpt1 function
