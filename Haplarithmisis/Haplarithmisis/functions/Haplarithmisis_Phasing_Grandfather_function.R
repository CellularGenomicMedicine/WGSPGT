# Haplarithmisis_phgrandfather - Function to phase parental genotypes using the grandfather

Haplarithmisis_phgrandfather <- function(Gtypes, parent, Parent1, Father, Mother, RefSampleID){
  # Subset genetic data for chromosome X
  MotherChrX <- Gtypes[Gtypes$Chr == "X", Mother]
  GtypesChrX <- Gtypes[Gtypes$Chr == "X", ]
  
  # Initialize PhasedPar with Parent1 (indication parent) genotypes 
  PhasedPar <- Gtypes[, Parent1]
  
  # Update the PhasedPar based on phased genotypes
  PhasedPar[(Gtypes[, RefSampleID] == "BB" & Gtypes[,Parent1] == "AB")] <- "BA"
  PhasedPar[(Gtypes[, RefSampleID] == "AA" & Gtypes[,Parent1] == "AB")] <- "AB"

  # Update the PhasedPar based on erroneous genotypes
  PhasedPar[(Gtypes[, RefSampleID] == "NC" & Gtypes[,Parent1] == "AB") |
            (Gtypes[, RefSampleID] == "AB" & Gtypes[,Parent1] == "AB") ] <- "NC"
 
  # Update Parent1 genotypes in Gtypes with phased genotypes
  Gtypes[, Parent1] <- PhasedPar

  # Handle phasing based on the indication parent 
  if (parent == "Mother") {
    # Update MotherChrX based on phased genotypes
    MotherChrX[(GtypesChrX[, RefSampleID] == "BB" & GtypesChrX[,Mother] == "AB")] <- "BA"
    
    # update Parent1 genotypes in chromsome X with phased MotherChrX
    Gtypes[Gtypes$Chr == "X", Parent1] <- MotherChrX
    # Assign genotypes directly for Father and Mother based on Parent1
    GTFather <- Gtypes[, Father]
    GTMother <- Gtypes[, Parent1]
  }
  
  if (parent == "Father" ) {
    GTFather = Gtypes[, Parent1]
    GTMother = Gtypes[, Mother]
    }

  # Create a data frame to store phased parental genotypes
  Parents <- data.frame(GTFather, GTMother, stringsAsFactors = FALSE, check.names = FALSE)
  names(Parents) <- c(Father, Mother)
  
  # Return phased parental genotypes
  return(Parents)
  
}#end function
