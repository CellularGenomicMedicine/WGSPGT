Haplarithmisis_Phasing <- function(script,REF,Father,Mother,RefSampleID,Gtypes,BAFs,parent,Parent1,EmbryoID,flip){
  if (REF == "Grandfather") {
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_Grandfather.R",sep="/"))
    PhBAF <- Haplarithmisis_Phasing_Grandfather(script,Father,Mother,RefSampleID,Gtypes,BAFs,parent,Parent1,EmbryoID,flip)
  }
  if (REF == "Grandmother") {
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_Grandmother.R",sep="/"))
    PhBAF <- Haplarithmisis_Phasing_Grandmother(script,Father,Mother,RefSampleID,Gtypes,BAFs,parent,Parent1,EmbryoID,flip)
  }
  if (REF == "SiblingM") {
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_SiblingMale.R",sep="/"))
    PhBAF <-Haplarithmisis_Phasing_SiblingMale(script,Father,Mother,RefSampleID,Gtypes,BAFs,EmbryoID,flip)
  }
  if (REF == "SiblingF") {
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_SiblingFemale.R",sep="/"))
    PhBAF <- Haplarithmisis_Phasing_SiblingFemale(script,Father,Mother,RefSampleID,Gtypes,BAFs,EmbryoID,flip)
  }
  return(PhBAF)
}