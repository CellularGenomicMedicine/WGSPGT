# Haplarithmisis_Phasing_parents - Function to phase parental genotypes

Haplarithmisis_Phasing_parents <- function(Father, Mother, GTFather, GTMother, GTRef){
  # Print messages to indicate the start of the function and its processes
  print("Phasing parental genotypes...")
  print("SNP-calls of the reference sibling that are in violation with Mendelian inheritance, are corrected and will be treated as NoCalls...")
  
  # Initialize a vector to store informative SNPs
  InfSNPs <- rep(0,length(GTFather))

  # Determine  parental informative SNPs, showing paterns of homozygosity in one parent and heterozygosity in the other parent
  InfSNPs[GTFather == "AA" & GTMother == "AB"] <- 1
  InfSNPs[GTFather == "BB" & GTMother == "AB"] <- 2
  InfSNPs[GTFather == "AB" & GTMother == "AA"] <- 3
  InfSNPs[GTFather == "AB" & GTMother == "BB"] <- 4
  InfSNPs[GTFather == "AB" & GTMother == "AB"] <- 5

  # Assuming that Pat1 and Mat1 homologous chromosomes were transmitted to the reference sibling
  
  # Corrections for the mother's genotype
  GTMother[GTRef == "AB" & InfSNPs == 1] <- "BA"
  GTMother[GTRef == "BB" & InfSNPs == 1] <- "NC" # Not possible
  GTMother[GTRef == "BB" & InfSNPs == 2] <- "BA"
  GTMother[GTRef == "AA" & InfSNPs == 2] <- "NC" # Not possible
  GTMother[GTRef == "BB" & InfSNPs == 5] <- "BA"
  GTMother[GTRef == "AB" & InfSNPs == 5] <- "NC"

  # Corrections for the father's genotype
  GTFather[GTRef == "AB" & InfSNPs == 3 ] <- "BA"
  GTFather[GTRef == "BB" & InfSNPs == 3 ] <- "NC" # Not possible
  GTFather[GTRef == "BB" & InfSNPs == 4 ] <- "BA"
  GTFather[GTRef == "AA" & InfSNPs == 4 ] <- "NC" # Not possible
  GTFather[GTRef == "BB" & InfSNPs == 5 ] <- "BA"
  GTFather[GTRef == "AB" & InfSNPs == 5 ] <- "NC"

  # Replace "NC" (NoCall) in the genotype vectors
  GTMother[GTRef == "NC"] <- "NC"
  GTFather[GTRef == "NC"] <- "NC"

  # Create a data frame to store phased parental genotypes
  Parents <- data.frame(GTFather,GTMother,stringsAsFactors=F,check.names=F)
  names(Parents) <- c(Father,Mother)
  
  # Return phased parental genotypes
  return(Parents)
  
} # End of function

