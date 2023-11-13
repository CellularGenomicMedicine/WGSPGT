Haplarithmisis_Phasing_Grandfather <- function(script,Father,Mother,RefSampleID,Gtypes,BAFs,parent,Parent1,EmbryoID,flip){
  PhBAF <- vector("list",8)
  names(PhBAF) <- c("P1","P2","M1","M2","P1Seg","P2Seg","M1Seg","M2Seg")
  print("Computing phased BAF...")
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_ChrXFatherAB.R",sep="/"))
  Gtypes[Gtypes[,"Chr"]=="X",Father] <- Haplarithmisis_ChrXFatherAB(Gtypes,Father)

  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_Grandfather_function.R",sep="/"))
  Parents <- Haplarithmisis_phgrandfather(Gtypes,parent,Parent1,Father,Mother,RefSampleID)
  GTFather <- Parents[,Father]
  GTMother <- Parents[,Mother]
  
  if (flip == 0) {BAFs[GTFather=="BA" | GTMother== "BA",!names(BAFs)%in%c("Names","Chr","Position")] <- 1 - BAFs[GTFather=="BA" | GTMother== "BA",!names(BAFs)%in%c("Names","Chr","Position")]}
  if (flip == 1) {BAFs[GTFather=="AB" | GTMother== "AB",!names(BAFs)%in%c("Names","Chr","Position")] <- 1 - BAFs[GTFather=="AB" | GTMother== "AB",!names(BAFs)%in%c("Names","Chr","Position")]}

  BAFs      <- BAFs[Gtypes[,RefSampleID]!="NC",]
  GTFather  <- GTFather[Gtypes[,RefSampleID]!="NC"]
  GTMother  <- GTMother[Gtypes[,RefSampleID]!="NC"]
  
  BAFs2 <- BAFs[,c("Names","Chr","Position",EmbryoID)]

  PhBAF[["P1"]] <-  BAFs2[((GTFather=="AB"& GTMother=="AA") | (GTFather=="BA"& GTMother=="BB")),]
  PhBAF[["P2"]] <-  BAFs2[((GTFather=="AB"& GTMother=="BB") | (GTFather=="BA"& GTMother=="AA")),]
  
  PhBAF[["M1"]] <- BAFs2[((GTMother=="AB"& GTFather=="AA") | (GTMother=="BA"& GTFather=="BB")),]
  PhBAF[["M2"]] <- BAFs2[((GTMother=="AB"& GTFather=="BB") | (GTMother=="BA"& GTFather=="AA")),]

  ToRep <- BAFs2
  ToRep[,EmbryoID] <- -1
  
  if(parent=="Father"){
    PhBAF[["M1"]]<-ToRep;PhBAF[["M2"]]<-ToRep;print("Parent1 is the Father")
  }
  if(parent=="Mother"){
    PhBAF[["P1"]]<-ToRep;PhBAF[["P2"]]<-ToRep;print("Parent1 is the Mother")
  }

  return(PhBAF)
}#end function

