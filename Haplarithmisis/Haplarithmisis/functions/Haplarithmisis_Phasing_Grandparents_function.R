# Haplarithmisis_phgrandparents - Function to phase parental genotypes using the grandparents

Haplarithmisis_phgrandparents <- function(Gtypes, parent, Parent1, Father, Mother, Grandfather, Grandmother){
  # Subset genetic data for chromosome X
	MotherChrX <- Gtypes[Gtypes$Chr == "X", Mother]
	GtypesChrX <- Gtypes[Gtypes$Chr == "X", ]
	
	# Initialize PhasedPar with Parent1 (indication parent) genotypes 
	PhasedPar <- Gtypes[, Parent1]

	# Update the PhasedPar based on phased genotypes
	PhasedPar[(Gtypes[, Grandfather] == "AB" & Gtypes[, Grandmother] == "AA" & Gtypes[, Parent1] == "AB") |
				   (Gtypes[, Grandfather] == "BB" & Gtypes[, Grandmother] == "AB" & Gtypes[, Parent1] == "AB") |
				   (Gtypes[, Grandfather] == "BB" & Gtypes[, Grandmother] == "AA" & Gtypes[, Parent1] == "AB")] <- "BA"

	# Update the PhasedPar based on erroneous genotypes
	PhasedPar[(Gtypes[, Grandfather] == "AA" & Gtypes[, Grandmother] == "AA" & Gtypes[, Parent1] == "AB") |
				   (Gtypes[, Grandfather] == "BB" & Gtypes[, Grandmother] == "BB" & Gtypes[, Parent1] == "AB") |
				   (Gtypes[, Grandfather] == "AA" & Gtypes[, Grandmother] == "BB" & Gtypes[, Parent1] == "BB") |
				   (Gtypes[, Grandfather] == "AA" & Gtypes[, Grandmother] == "BB" & Gtypes[, Parent1] == "AA") |
				   ((Gtypes[, Grandfather] == "NC" | Gtypes[, Grandmother] == "NC") & Gtypes[, Parent1] == "AB") |
				   (Gtypes[, Grandfather] == "AB" & Gtypes[, Grandmother] == "AB")] <- "NC"

	# Update Parent1 genotypes in Gtypes with phased genotypes
	Gtypes[, Parent1] <- PhasedPar

	# Handle phasing based on the indication parent 
	if(parent == "Mother"){
	  # Update MotherChrX based on phased genotypes
		MotherChrX[(GtypesChrX[, Grandfather] == "BB" & GtypesChrX[, Mother] == "AB")] <- "BA"

		MotherChrX[(GtypesChrX[, Grandfather] == "NC" & GtypesChrX[, Mother] == "AB") |
					   (GtypesChrX[, Grandfather] == "BB" & GtypesChrX[, Grandmother] == "BB" & GtypesChrX[, Mother] == "AB") |

					   (GtypesChrX[, Grandfather] == "AA" & GtypesChrX[, Grandmother] == "AA" & GtypesChrX[, Mother] == "AB") ] <- "NC"
		
		# update Parent1 genotypes in chromsome X with phased MotherChrX
		Gtypes[Gtypes$Chr == "X", Parent1] <- MotherChrX
		# Assign genotypes directly for Father and Mother based on Parent1
		GTFather = Gtypes[, Father]
		GTMother = Gtypes[, Parent1]
	}	   

	if(parent == "Father"){
	  GTFather = Gtypes[, Parent1]
	  GTMother = Gtypes[, Mother]
	}
	
	# Create a data frame to store phased parental genotypes
	Parents <- data.frame(GTFather, GTMother, stringsAsFactors = FALSE, check.names = FALSE)
	names(Parents) <- c(Father, Mother)
	
	# Return phased parental genotypes
	return(Parents)

}#end function
