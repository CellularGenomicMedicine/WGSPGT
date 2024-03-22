# Haplarithmisis_Htyping - function to perform haplotyping and select which function to use based on the referent sample

Haplarithmisis_Htyping <- function(script,Father,Mother,REF, RefSampleID, Gtypes, ParScore, Int, Window, parent, Parent1, EmbryoID, Chroms){
  
  # Source the functions used in this function
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Htyping_Grandfather.R", sep = "/"))
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Htyping_Grandmother.R", sep = "/"))
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Htyping_SiblingMale.R", sep = "/"))
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Htyping_SiblingFemale.R", sep = "/"))
  
  # Condition to select the appropriate function based on the reference sample type
  if(REF == "Grandfather") {
    Haps <- Haplarithmisis_Htyping_Grandfather(script, Father, Mother, RefSampleID, Gtypes, ParScore, Window, Int, Chroms, parent, Parent1, EmbryoID)
  }
  if (REF == "Grandmother") {
    Haps <- Haplarithmisis_Htyping_Grandmother(script, Father, Mother, RefSampleID, Gtypes, ParScore, Window, Int,Chroms, parent, Parent1, EmbryoID)
  }
  if (REF == "SiblingM") {
    Haps <- Haplarithmisis_Htyping_SiblingMale(script, Father, Mother, RefSampleID, ParScore, Gtypes, Window, Int, Chroms, EmbryoID)
  }
  if (REF == "SiblingF") {
    Haps <- Haplarithmisis_Htyping_SiblingFemale(script, Father, Mother, RefSampleID, ParScore, Gtypes, Window, Int, Chroms, EmbryoID)
  }
  
  # Return the haplotypes
  return(Haps)
}