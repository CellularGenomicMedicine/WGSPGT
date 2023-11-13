Haplarithmisis_Phasing_SiblingFemale <- function(script,Father,Mother,RefSampleID,Gtypes,BAFs,EmbryoID,flip){
  PhBAF <- vector("list",8)
  names(PhBAF) <- c("P1","P2","M1","M2","P1Seg","P2Seg","M1Seg","M2Seg")
  print("Computing phased BAF...")
  GTFather <- Gtypes[,Father]
  GTMother <- Gtypes[,Mother]
  GTRef    <- Gtypes[,RefSampleID]

  GTMotherX <- GTMother[Gtypes$Chr=="X"]
  GTRefX    <- GTRef[Gtypes$Chr=="X"]

  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_ChrXFatherAB.R",sep="/"))
  GTFatherX <- Haplarithmisis_ChrXFatherAB(Gtypes,Father)

  RefSex <- "female"
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_parentsX.R",sep="/"))
  ParentsX <- Haplarithmisis_Phasing_parentsX (Father,Mother,GTFatherX,GTMotherX,GTRefX,RefSex)
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_parents.R",sep="/"))
  Parents <- Haplarithmisis_Phasing_parents(Father,Mother,GTFather,GTMother,GTRef)
  
  GTFather <- Parents[,Father]
  GTMother <- Parents[,Mother]
  GTFather[Gtypes$Chr=="X"] <- ParentsX[,Father]
  GTMother[Gtypes$Chr=="X"] <- ParentsX[,Mother]

  if (flip == 0) {BAFs[GTFather=="BA" | GTMother== "BA",!names(BAFs)%in%c("Names","Chr","Position")] <- 1 - BAFs[GTFather=="BA" | GTMother== "BA",!names(BAFs)%in%c("Names","Chr","Position")]}
  if (flip == 1) {BAFs[GTFather=="AB" | GTMother== "AB",!names(BAFs)%in%c("Names","Chr","Position")] <- 1 - BAFs[GTFather=="AB" | GTMother== "AB",!names(BAFs)%in%c("Names","Chr","Position")]}

  BAFs     <- BAFs[GTRef!="NC",]
  GTFather <- GTFather[GTRef!="NC"]
  GTMother <- GTMother[GTRef!="NC"]

  BAFs2       <- BAFs[,c("Names","Chr","Position",EmbryoID)]
  PhBAF[["P1"]] <- BAFs2[((GTFather=="AB"& GTMother=="AA") | (GTFather=="BA"& GTMother=="BB")),]
  PhBAF[["P2"]] <- BAFs2[((GTFather=="AB"& GTMother=="BB") | (GTFather=="BA"& GTMother=="AA")),]
  PhBAF[["M1"]] <- BAFs2[((GTMother=="AB"& GTFather=="AA") | (GTMother=="BA"& GTFather=="BB")),]
  PhBAF[["M2"]] <- BAFs2[((GTMother=="AB"& GTFather=="BB") | (GTMother=="BA"& GTFather=="AA")),]

  return(PhBAF)
}#end function

