#' Haplarithmisis_qcbyparents - Function to perform quality control on genotype data utilziing parents genotype data
#' Allele drop in (ADI) and allele drop out (ADO) are calculated

Haplarithmisis_qcbyparents <- function(Father, Mother, Gtypes, EmbryoID, Chroms, outPath){
	print("QC by parents analysis...")
  
  # Initialize a list to store quality control results for ADO and ADI
	QC_par <- vector("list", 2)
	names(QC_par)<-c("ADO","ADI")
	
	# Genome-wide quality control by parents
	# Filter genotype data for autosomes
	GTFather <- Gtypes[, Father][Gtypes$Chr != "X" | Gtypes$Chr == "XY" | Gtypes$Chr == "Y"]
	GTMother <- Gtypes[, Mother][Gtypes$Chr != "X" | Gtypes$Chr == "XY" | Gtypes$Chr == "Y"]
	GTSib <- Gtypes[Gtypes$Chr != "X" | Gtypes$Chr == "XY" | Gtypes$Chr == "Y", EmbryoID]
	
	# Compute ADO and ADI for autosomes
	ADO_AutGenome <- (sum((GTFather=="AB" & GTMother=="BB" & GTSib == "AA") |
						   				   (GTFather == "AB" & GTMother == "AA" & GTSib == "BB") |
						    			   (GTFather == "BB" & GTMother == "AB" & GTSib == "AA") |
						    			   (GTFather == "AA" & GTMother == "AB" & GTSib == "BB") |
						    			   (GTFather == "AA" & GTMother == "BB" & (GTSib == "AA" |  GTSib == "BB")) |
						    			   (GTFather == "BB" & GTMother == "AA" & (GTSib == "AA" |  GTSib == "BB"))) /
						   				     sum((GTSib == "AA" | GTSib == "BB") & (GTFather != "NC" | GTMother != "NC"))) * 100
	
	ADI_AutGenome <- (sum((GTFather == "AA" & GTMother == "AA" & GTSib == "AB") |
   						 				  (GTFather == "BB" & GTMother == "BB" & GTSib == "AB")) / 
   						 				    sum(GTSib == "AB" & (GTFather != "NC" | GTMother != "NC"))) * 100

	
	# Chromosome-specific QC by parents
	for(chr in Chroms){
	
		FatherChr <- Gtypes[, Father][Gtypes$Chr == chr]
		MotherChr <- Gtypes[, Mother][Gtypes$Chr == chr]
		SibChr <- Gtypes[Gtypes$Chr == chr, EmbryoID]
		
		ADO <- (sum((FatherChr == "AB" & MotherChr == "BB" & SibChr == "AA") |
						    (FatherChr == "AB" & MotherChr == "AA" & SibChr == "BB") | 
						    (FatherChr == "BB" & MotherChr == "AB" & SibChr == "AA") |
						    (FatherChr == "AA" & MotherChr == "AB" & SibChr == "BB") |
						    (FatherChr == "AA" & MotherChr == "BB" & (SibChr == "AA" |  SibChr == "BB")) |
						    (FatherChr == "BB" & MotherChr == "AA" & (SibChr == "AA" |  SibChr == "BB")))/sum((GTSib == "AA" | GTSib == "BB") & (GTFather != "NC" | GTMother != "NC")))*100
	
		ADI <- (sum((FatherChr == "AA" & MotherChr == "AA" & SibChr == "AB") | 
		              (FatherChr == "BB" & MotherChr == "BB" & SibChr == "AB")) /
		          sum(GTSib == "AB" & (GTFather != "NC" | GTMother != "NC"))) * 100
		
	# If it is the first chromosome, initialize matrices
	if (chr == Chroms[1]){ 
	  ADOGenome <- ADO 
	  ADIGenome <- ADI 
	  } else {
	  # If not, concatenate to existing matrix
	  ADOGenome <- cbind(ADOGenome,ADO)
	  ADIGenome <- cbind(ADIGenome,ADI)
	  }

	}# End chromosome loop
	
	# Combine genome-wide ADO and ADI with chromosome-specific results
	ADOGenome <- cbind(ADOGenome,ADO_AutGenome)
	ADIGenome <- cbind(ADIGenome,ADI_AutGenome)

	ADOs <- ADOGenome
	ADIs <- ADIGenome

	# Set column and rownames
	colnames(ADOs) <- c(Chroms,"GenomeAut")
	colnames(ADIs) <- c(Chroms,"GenomeAut")
	rownames(ADOs) <- EmbryoID
	rownames(ADIs) <- EmbryoID
	
	# Write ADO matrix to file
	write.table(ADOs, paste(outPath, paste(EmbryoID, "ChrSpec_ADO.txt", sep = "_"), sep = "/"), 
	            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
	# Write ADI matrix to file
	write.table(ADIs, paste(outPath, paste(EmbryoID, "ChrSpec_ADI.txt", sep = "_"), sep = "/"), 
	            sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
	
	# Store results in the list and return
	QC_par[["ADO"]] <- ADOs
	QC_par[["ADI"]] <- ADIs
	
	
	return(QC_par)
} # End Haplarithmisis_qcbyparents function 
