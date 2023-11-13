Haplarithmisis_Htyping <- function(script,Father,Mother,REF, RefSampleID, Gtypes, ParScore,Int,Window,parent, Parent1, EmbryoID,Chroms){
  if(REF=="Grandfather") {
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Htyping_Grandfather.R",sep="/"))
    Haps <- Haplarithmisis_Htyping_Grandfather(script,Father,Mother,RefSampleID, Gtypes, ParScore, Window,Int,Chroms,parent,Parent1, EmbryoID)
  }
  if (REF=="Grandmother") {
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Htyping_Grandmother.R",sep="/"))
    Haps <- Haplarithmisis_Htyping_Grandmother(script,Father,Mother,RefSampleID, Gtypes, ParScore, Window,Int,Chroms,parent,Parent1, EmbryoID)
  }
  if (REF=="SiblingM") {
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Htyping_SiblingMale.R",sep="/"))
    Haps <- Haplarithmisis_Htyping_SiblingMale(script,Father,Mother,RefSampleID, ParScore,Gtypes, Window, Int,Chroms,EmbryoID)
  }
  if (REF=="SiblingF") {
    source(paste(script, "analyses", "Haplarithmisis", "functions","Haplarithmisis_Htyping_SiblingFemale.R",sep="/"))
    Haps <- Haplarithmisis_Htyping_SiblingFemale(script,Father,Mother,RefSampleID, ParScore,Gtypes, Window, Int,Chroms,EmbryoID)
  }
  return(Haps)
}