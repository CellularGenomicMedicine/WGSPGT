# Haplarithmisis_Phasing_SiblingMale - Function to phase when male sibling is referent

Haplarithmisis_Phasing_SiblingMale <- function(script, Father, Mother, RefSampleID, Gtypes, BAFs, EmbryoID, flip){
  # Source the function for handling AB SNPs on chromosome X of Father
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_ChrXFatherAB.R", sep = "/"))
  # Source the function to phase chromosome X using parents
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_parentsX.R", sep = "/"))
  # Source the function to phase using parents
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_parents.R", sep = "/"))
  
  # Initialize a list to store phased BAFs
  PhBAF <- vector("list", 8)
  names(PhBAF) <- c("P1", "P2", "M1", "M2", "P1Seg", "P2Seg", "M1Seg", "M2Seg")
  
  # print message indicating the start of phased BAF computation
  print("Computing phased BAF...")
  
  # Change AB SNPs on chromosome X of father to NoCalls "NC"
  GTFatherX <- Haplarithmisis_ChrXFatherAB(Gtypes, Father)

  GTFather <- Gtypes[, Father]
  GTMother <- Gtypes[, Mother]
  GTRef <- Gtypes[, RefSampleID]

  GTMotherX <- GTMother[Gtypes$Chr == "X"]
  GTRefX    <- GTRef[Gtypes$Chr == "X"]

  RefSex<- "male"
  
  # Phase the chromosome X using the parents as referent
  ParentsX <- Haplarithmisis_Phasing_parentsX (Father, Mother, GTFatherX, GTMotherX, GTRefX, RefSex)
  # Phase using the parents as referent
  Parents <- Haplarithmisis_Phasing_parents(Father, Mother, GTFather, GTMother, GTRef)
  
  GTFather <- Parents[,Father]
  GTMother <- Parents[,Mother]
  GTFather[Gtypes$Chr == "X"] <- ParentsX[,Father]
  GTMother[Gtypes$Chr == "X"] <- ParentsX[,Mother]

  # Flip BAFs based on the flip parameter
  if (flip == 0) {
    BAFs[GTFather == "BA" | GTMother== "BA",!names(BAFs) %in% c("Names","Chr","Position")] <- 1 - BAFs[GTFather == "BA" | GTMother== "BA", !names(BAFs) %in% c("Names", "Chr", "Position")]
    }
  if (flip == 1) {
    BAFs[GTFather == "AB" | GTMother== "AB",!names(BAFs) %in% c("Names","Chr","Position")] <- 1 - BAFs[GTFather == "AB" | GTMother== "AB", !names(BAFs) %in% c("Names", "Chr", "Position")]
    }

  # Filter BAFs and genotypes based on NoCalls
  BAFs <- BAFs[GTRef != "NC",]
  GTFather <- GTFather[GTRef != "NC"]
  GTMother <- GTMother[GTRef != "NC"]

  # Prepare BAFs for output
  BAFs2 <- BAFs[,c("Names", "Chr", "Position", EmbryoID)]
  
  # Assign phased BAFs into categories based on parental genotypes 
  PhBAF[["P1"]] <- BAFs2[((GTFather == "AB" & GTMother == "AA") | (GTFather == "BA" & GTMother == "BB")), ]
  PhBAF[["P2"]] <- BAFs2[((GTFather == "AB" & GTMother == "BB") | (GTFather == "BA" & GTMother == "AA")), ]
  PhBAF[["M1"]] <- BAFs2[((GTMother == "AB" & GTFather == "AA") | (GTMother == "BA" & GTFather == "BB")), ]
  PhBAF[["M2"]] <- BAFs2[((GTMother == "AB" & GTFather == "BB") | (GTMother == "BA" & GTFather == "AA")), ]

  # Return phased BAFs
  return(PhBAF)
}#end function

