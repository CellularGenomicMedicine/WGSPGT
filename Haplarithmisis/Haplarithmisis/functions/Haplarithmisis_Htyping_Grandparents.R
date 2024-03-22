# Haplarithmisis_Htyping_Grandparents - Function for haplotyping using grandparents as reference

Haplarithmisis_Htyping_Grandparents <- function(script,Father, Mother, Grandfather, Grandmother, Gtypes, ParScore, parent, Parent1, EmbryoID, Chroms){

  # Source the functions used in this function
  # Source the function used for sex determination of the embryo
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_sexdetermination.R",sep="/"))
  # Source the function for handling AB SNPs on chromosome X of Father
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_ChrXFatherAB.R",sep="/"))
  # Source the function for haplotyping chromosome X using grandparents
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_chrxhtyping_Grandparents.R",sep="/"))
  # Source the function to segment the haplotypes
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_inthapnew1.R",sep="/"))
  
  # Print a message indicating the start of haplotyping
	print("*******************************************")
	print("****     Option1Htype analysis...     ****")
	
	# Matrix to store sex determination results
	ScSexes <- matrix(NA, length(ParScore), 2)
	rownames(ScSexes) <- names(ParScore)
	
	# Perform sex determination using the ParScore
	ScSexes <- Haplarithmisis_sexdetermination(ScSexes, ParScore)

	# Change "AB" SNPs on chromosome X of father to NoCalls "NC"
	Gtypes[Gtypes[,"Chr"] == "X", Father] <- Haplarithmisis_ChrXFatherAB(Gtypes, Father)

	# Perform haplotyping for chromosome X using grandparents as reference
	HapXMat <- Haplarithmisis_chrxhtyping_Grandparents(Father, Mother, Grandfather, Grandmother, Gtypes, parent, ScSexes, EmbryoID)

	# Phasing for parental haplotypes of the parent of indication
	PhasedPar <- Gtypes[, Parent1]
	
	if (parent == "Mother") {
	  Parent2 = Father 
	} else { 
	  Parent2 = Mother
	  }

	# Phase parental genotypes based genotype inherited from grandfather as first and from Grandmother as second
	PhasedPar[(Gtypes[, Grandfather] == "AB" & Gtypes[, Grandmother] == "AA" & Gtypes[, Parent1] == "AB") |
				    (Gtypes[, Grandfather] == "BB" & Gtypes[, Grandmother] == "AB" & Gtypes[, Parent1] == "AB") |
				    (Gtypes[, Grandfather] == "BB" & Gtypes[, Grandmother] == "AA" & Gtypes[, Parent1] == "AB") ] <- "BA"

	# Phase parental genotypes based on errorenous conditions to NC (NoCall)
	PhasedPar[(Gtypes[, Grandfather] == "AA" & Gtypes[, Grandmother] == "AA" & Gtypes[, Parent1] == "AB") |
	          (Gtypes[, Grandfather] == "BB" & Gtypes[, Grandmother] == "BB" & Gtypes[, Parent1] == "AB") |
	          (Gtypes[, Grandfather] == "AA" & Gtypes[, Grandmother] == "BB" & Gtypes[, Parent1] == "BB") |
	          (Gtypes[, Grandfather] == "AA" & Gtypes[, Grandmother] == "BB" & Gtypes[, Parent1] == "AA") |
	          ((Gtypes[, Grandfather] == "AA" | Gtypes[, Grandmother] == "AA") & Gtypes[, Parent1] == "BB") |
	          ((Gtypes[, Grandfather] == "AA" | Gtypes[, Grandmother] == "AB") & Gtypes[, Parent1] == "BB") |
	          ((Gtypes[, Grandfather] == "BB" | Gtypes[, Grandmother] == "AA") & Gtypes[, Parent1] == "BB") |
	          ((Gtypes[, Grandfather] == "BB" | Gtypes[, Grandmother] == "AA") & Gtypes[, Parent1] == "AA") |
	          ((Gtypes[, Grandfather] == "BB" | Gtypes[, Grandmother] == "AB") & Gtypes[, Parent1] == "AA") |
	          ((Gtypes[, Grandfather] == "BB" | Gtypes[, Grandmother] == "BB") & Gtypes[, Parent1] == "AA") |
	          ((Gtypes[, Grandfather] == "NC" | Gtypes[, Grandmother] == "NC") & Gtypes[, Parent1] == "AB")] <- "NC"

	
	Gtypes[, Parent1] <- PhasedPar

	# Initialize variable for haplotyping
	Hap1 <- rep(0, nrow(Gtypes))
	Hap2 <- rep(0, nrow(Gtypes))
	print(EmbryoID)
	GtypeChild <- Gtypes[, EmbryoID]

	cond1 <- 1
	cond2 <- 2

	# Assign haplotypes based on parental genotypes and specific conditions
	Hap1[(Gtypes[, Parent1] == "BA" & Gtypes[, Parent2] == "AA" & Gtypes[, EmbryoID] == "AB") |
	       (Gtypes[, Parent1] == "BA" & Gtypes[, Parent2] == "AA" & Gtypes[, EmbryoID] == "BB") |
	       (Gtypes[, Parent1] == "BA" & Gtypes[, Parent2] == "BB" & Gtypes[, EmbryoID] == "BB") |
	       (Gtypes[, Parent1] == "AB" & Gtypes[, Parent2] == "BB" & Gtypes[, EmbryoID] == "AB") |
	       (Gtypes[, Parent1] == "AB" & Gtypes[, Parent2] == "BB" & Gtypes[, EmbryoID] == "AA") |
	       (Gtypes[, Parent1] == "AB" & Gtypes[, Parent2] == "AA" & Gtypes[, EmbryoID] == "AA")] <- cond1

	Hap1[(Gtypes[, Parent1] == "BA" & Gtypes[, Parent2] == "AA" & Gtypes[, EmbryoID] == "AA") |
	       (Gtypes[, Parent1] == "BA" & Gtypes[, Parent2] == "BB" & Gtypes[, EmbryoID] == "AB") |
	       (Gtypes[, Parent1] == "BA" & Gtypes[, Parent2] == "BB" & Gtypes[, EmbryoID] == "AA") |
	       (Gtypes[, Parent1] == "AB" & Gtypes[, Parent2] == "BB" & Gtypes[, EmbryoID] == "BB") |
	       (Gtypes[, Parent1] == "AB" & Gtypes[, Parent2] == "AA" & Gtypes[, EmbryoID] == "AB") |
	       (Gtypes[, Parent1] == "AB" & Gtypes[, Parent2] == "AA" & Gtypes[, EmbryoID] == "BB")] <- cond2
	
	# Assign NoCall to specific genotypes		 
	Hap1[Gtypes[, Grandfather] == "AB" & Gtypes[, Grandmother] == "AB"] <- 0
	ParHaps <- Hap1
	names(ParHaps) <- EmbryoID

	# Initialize a list to store haplotype matrices
	Haps <- vector("list", 2)
	names(Haps) <- c("dataHapRaw", "dataHap")

	# Replace haplotypes for chromosome X in ParHaps by the HapXMat that was generated
	ParHaps[Gtypes$Chr == "X"] <- HapXMat

	ParHapsAll <- data.frame(Gtypes[,c("Names", "Chr", "Position")], ParHaps, stringsAsFactors = FALSE, check.names = FALSE)
	names(ParHapsAll) <- c("Names","Chr","Position", EmbryoID)
	dataHapRaw <- ParHapsAll

	# Perform haplotype interpretation 
	dataHap1 <- Haplarithmisis_inthapnew1(script, dataHapRaw, Window, Int, Chroms)

	if (parent == "Father") { 
	  par1 = "_Pat"; 
	  par2 = "_Mat" 
	  } else { 
	  par1 = "_Mat"; 
	  par2 = "_Pat" 
	  }

	Haps[["dataHapRaw"]] <- ParHapsAll
	Haps[["dataHapRaw"]]$mat <- 0
	colnames(Haps[["dataHapRaw"]]) <- c("Names", "Chr", "Position", paste(EmbryoID, par1, sep = ""), paste(EmbryoID, par2, sep = ""))

	Haps[["dataHap"]] <- dataHap1
	Haps[["dataHap"]]$mat <- 0
	colnames(Haps[["dataHap"]])	<- colnames(Haps[["dataHapRaw"]])

	# Return the haplotypes
	return(Haps)
}#end function
