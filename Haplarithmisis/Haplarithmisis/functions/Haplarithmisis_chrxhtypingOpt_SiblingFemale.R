# Haplarithmisis_chrxhtypingOpt_SiblingFemale - Function for haplotyping of chromosome X using female sibling as reference

Haplarithmisis_chrxhtypingOpt_SiblingFemale <- function(script, Father, Mother, RefSampleID, Gtypes, ParScore, EmbryoID){
	
  # Source the functions used in this function
  
  # Source the function used for sex determination of the embryo
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_sexdetermination.R", sep = "/"))
  # Source the function for handling AB SNPs on chromosome X of Father
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_ChrXFatherAB.R", sep = "/"))

  # Matrix to store sex determination results
  ScSexes <- matrix(NA, length(ParScore), 2)
  rownames(ScSexes) <- names(ParScore)
  
  # Perform sex determination using the ParScore
  ScSexes <- Haplarithmisis_sexdetermination(ScSexes, ParScore)
  
  # Print a message indicating the start of haplotyping of chromosome X
  print("Haplotype reconstrucion of the sex chromosome...")
	print("Ammended chrxhtyping")

	# Change "AB" SNPs on chromosome X of father to NoCalls "NC"
	Gtypes[Gtypes[,"Chr"] == "X", Father] <- Haplarithmisis_ChrXFatherAB(Gtypes, Father)

	# Subset genetic data for chromosome X
	GTFather <- Gtypes[, Father][Gtypes$Chr == "X"]
	GTMother <- Gtypes[, Mother][Gtypes$Chr == "X"]
	GTRef <- Gtypes[Gtypes$Chr == "X", RefSampleID]
	
	# Adjust genotypes based on specific genetic patterns
	GTRef[GTRef == "AB" & GTFather == "BB"]="BA"#There could not be htz SNPs on chromosome X of male affected
	GTMother[GTRef == "AB" & GTFather== "AA" & GTMother == "AB"] = "BA" #M1 is the carrier chromosome
	GTMother[GTRef == "BB" & GTMother == "AB"] = "BA" #M1 is the carrier chromosome
	GTMother[GTRef == "BB" & GTFather== "AA" & GTMother == "AB"] = "NC"
	GTMother[GTRef == "AA" & GTFather== "BB" & GTMother == "AB"] = "NC"

	MatHaps <- do.call("rbind", strsplit(GTMother, split = ""))
	PatHap <-	as.matrix(do.call("rbind", strsplit(GTFather, split = ""))[, 1])

	M1 <- as.matrix(MatHaps[, 1])
	M2 <- as.matrix(MatHaps[, 2])

	Bls <- Gtypes[Gtypes$Chr == "X", EmbryoID]#Blastomeres
	Bls <- as.matrix(Bls)

	UnInf <- which(Bls == "NC" | GTMother == "NC" | GTFather == "NC" | GTRef == "NC" | GTMother == "AA" | GTMother == "BB")

	B1 <- matrix(0, nrow(Bls), 1)
	B2 <- matrix(0, nrow(Bls), 1)

	if(is.na(ScSexes[1, 2])){
		print(paste("The Chr.X haplotype of Sample", EmbryoID, "could not be reconstructed as the origin of this chromosome could not be determined..."))
	} else if (ScSexes[1, 2] == "male") {
		Bls[Bls == "AB"] <- "NC"
		BlHaps <- do.call("rbind", strsplit(Bls, split = ""))

		B1[M1 == BlHaps[, 1]] <- 1
		B1[M2 == BlHaps[, 1]] <- 2
		B1[UnInf] <- 0
		print(paste("The", EmbryoID, "blastomere is from a male Embryo"))
	} else if (ScSexes[1, 2] == "female") {
		Bls[PatHap == "B" & Bls == "AB"] = "BA"#Phasing of blastomere based on the father's genotype
		BlHaps <- do.call("rbind", strsplit(Bls, split = ""))
		B1[M1 == BlHaps[, 2]] <- 1
		B1[M2 == BlHaps[, 2]] <- 2
		B2[, 1] <- 1
		B1[UnInf] <- 0
		print(paste("The", EmbryoID, "blastomere is from a female Embryo"))
	}
	
	# Combine haplotypes for maternal and paternal chromosome
	PhBl <- cbind(B1, B2)
	PhBls <- PhBl

	matHaps <- PhBls[, seq(1, ncol(PhBls), 2)]
	patHaps <- PhBls[, seq(2, ncol(PhBls), 2)]

	matHaps <- as.matrix(matHaps)
	patHaps <- as.matrix(patHaps)
	colnames(matHaps) <- paste(EmbryoID, "_Mat", sep = "")
	colnames(patHaps) <- paste(EmbryoID, "_Pat", sep = "")

	# Combine haplotypes with other relevant information like names, chromsome and position
	HapsChrX <- cbind(Gtypes[Gtypes$Chr == "X", c("Names", "Chr", "Position")], patHaps, matHaps)
	
	# Return the Haplotype from chromosome X
	return(HapsChrX)
}#end function

