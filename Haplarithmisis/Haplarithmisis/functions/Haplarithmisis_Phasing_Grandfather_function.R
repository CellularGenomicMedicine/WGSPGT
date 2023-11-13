Haplarithmisis_phgrandfather <- function(Gtypes,parent,Parent1,Father,Mother,RefSampleID){
  MotherChrX <- Gtypes[Gtypes$Chr=="X",Mother]
  GtypesChrX <- Gtypes[Gtypes$Chr=="X",]
  
  PhasedPar <- Gtypes[,Parent1]
  
  PhasedPar[(Gtypes[,RefSampleID]=="BB" & Gtypes[,Parent1]=="AB")]<-"BA"
  PhasedPar[(Gtypes[,RefSampleID]=="AA" & Gtypes[,Parent1]=="AB")]<-"AB"

  PhasedPar[(Gtypes[,RefSampleID]=="NC" & Gtypes[,Parent1]=="AB") |
            (Gtypes[,RefSampleID]=="AB" & Gtypes[,Parent1]=="AB") ] <- "NC"
 
  Gtypes[,Parent1] <- PhasedPar

  if(parent=="Mother"){
    
    MotherChrX[(GtypesChrX[,RefSampleID]=="BB" & GtypesChrX[,Mother]=="AB")]<-"BA"
    
    Gtypes[Gtypes$Chr=="X",Parent1] <- MotherChrX
    GTFather <- Gtypes[,Father]
    GTMother <- Gtypes[,Parent1]
  }
  
  if(parent=="Father" ){GTFather=Gtypes[,Parent1];GTMother=Gtypes[,Mother]}


  Parents <- data.frame(GTFather,GTMother,stringsAsFactors=F,check.names=F)
  names(Parents)  <- c(Father,Mother)
  return(Parents)
  
}#end function
