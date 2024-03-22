# Haplarithmisis_htypingAutOpt_SiblingMale - Function to haplotyping the autosomes when male sibling is referent

Haplarithmisis_htypingAutOpt_SiblingMale <- function(script, Gtypes, Father, Mother, RefSampleID, GTEmbryo, EmbryoID){
  # Source script to phase parents
  source(paste(script, "analyses", "Haplarithmisis", "functions","Haplarithmisis_Phasing_parents.R", sep = "/"))
  
  # Print messages to indicate the start of the function and its processes
  print("*******************************************")
	print("****     Option2 Htype analysis...     ****")
	print("-------------------------------------------")
	print("#1Phasing parental genotypes...")
	
	# Phase parental genotypes
	Parents <- Haplarithmisis_Phasing_parents(Father, Mother, Gtypes[, Father], Gtypes[, Mother], Gtypes[, RefSampleID])
	GTFather <- Parents[,Father]
	GTMother <- Parents[,Mother]
	
	# Print message indicating the phasing of autosomes
	print("-------------------------------------------")
	print("#2 Phasing autosomes of the siblings...")
	InfSNPs <- rep(0, length(GTFather))
	
	# Determination of parental informative SNPs, showing patern of homozygousity in one parent and heterozygousity in the other parent
	InfSNPs[GTFather == "AB" & GTMother == "AA"] <- 1
	InfSNPs[GTFather == "AB" & GTMother == "BB"] <- 2
	InfSNPs[GTFather == "BA" & GTMother == "AA"] <- 3
	InfSNPs[GTFather == "BA" & GTMother == "BB"] <- 4
	InfSNPs[GTFather == "AA" & GTMother == "AB"] <- 5
	InfSNPs[GTFather == "BB" & GTMother == "AB"] <- 6
	InfSNPs[GTFather == "AA" & GTMother == "BA"] <- 7
	InfSNPs[GTFather == "BB" & GTMother == "BA"] <- 8

	# Initialize matrices for paternal and maternal haplotypes
	PatHap <- rep(0, length(GTFather))
	MatHap <- rep(0, length(GTMother))
	
	cond1 <- 1
	cond2 <- 2

	# Determine paternal and maternal haplotypes based on informative SNPs and embryo genotype
	PatHap[(InfSNPs == 2 & GTEmbryo == "AB") | (InfSNPs == 2 & GTEmbryo == "AA") | (InfSNPs == 1 & GTEmbryo == "AA") |
	       (InfSNPs == 3 & GTEmbryo == "AB") | (InfSNPs == 3 & GTEmbryo == "BB") | (InfSNPs == 4 & GTEmbryo == "BB") ] <- cond1
				
	PatHap[(InfSNPs == 1 & GTEmbryo == "AB") | (InfSNPs == 1 & GTEmbryo == "BB") | (InfSNPs == 2 & GTEmbryo == "BB") |
	       (InfSNPs == 4 & GTEmbryo == "AB") | (InfSNPs == 4 & GTEmbryo == "AA") | (InfSNPs == 3 & GTEmbryo == "AA") ] <- cond2
	
		
	MatHap[(InfSNPs == 6 & GTEmbryo == "AB") | (InfSNPs == 6 & GTEmbryo == "AA") | (InfSNPs == 5 & GTEmbryo == "AA") |
	       (InfSNPs == 7 & GTEmbryo == "AB") | (InfSNPs == 7 & GTEmbryo == "BB") | (InfSNPs == 8 & GTEmbryo == "BB") ] <- cond1
				
	MatHap[(InfSNPs == 5 & GTEmbryo == "AB") | (InfSNPs == 5 & GTEmbryo == "BB") | (InfSNPs == 6 & GTEmbryo == "BB") |
	       (InfSNPs == 8 & GTEmbryo == "AB") | (InfSNPs == 8 & GTEmbryo == "AA") | (InfSNPs == 7 & GTEmbryo == "AA") ] <- cond2
	
	# Print message indicating successful phasing of the sibling
	print(paste("Sibling", EmbryoID, "is phased"))
	
	# Combine haplotypes for paternal and maternal chromosomes
	MatHaps <- MatHap;
	PatHaps <- PatHap

	PatHaps <- as.matrix(PatHaps)
	MatHaps <- as.matrix(MatHaps)

	colnames(PatHaps)<-paste(EmbryoID, "_Pat", sep = "")
	colnames(MatHaps)<-paste(EmbryoID, "_Mat", sep = "")

	# Combine haplotypes with other relevant information
	HapsAut <- cbind(Gtypes[Gtypes$Chr != "X", c("Names", "Chr", "Position")], PatHaps[Gtypes$Chr != "X", ], MatHaps[Gtypes$Chr != "X", ])
	colnames(HapsAut) <- c("Names", "Chr", "Position", colnames(PatHaps), colnames(MatHaps))
	
	# Print end-of-function messages
	print("-------------------------------------------")
	print("*******************************************")

return(HapsAut)
}# end function
