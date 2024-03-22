# Haplarithmisis_Phasing - function to perform phasing of B-allele frequencies and select which function to use based on the referent sample

Haplarithmisis_Phasing <- function(script, REF, Father, Mother, RefSampleID, Gtypes, BAFs, parent, Parent1, EmbryoID, flip){
  
  # Source the functions used in this function
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_Grandfather.R", sep = "/"))
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_Grandmother.R", sep = "/"))
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_SiblingMale.R", sep = "/"))
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_SiblingFemale.R", sep = "/"))
  
  # Condition to select the appropriate function based on the reference sample type
  if (REF == "Grandfather") {
    PhBAF <- Haplarithmisis_Phasing_Grandfather(script, Father, Mother, RefSampleID, Gtypes, BAFs, parent, Parent1, EmbryoID, flip)
  }
  if (REF == "Grandmother") {
    PhBAF <- Haplarithmisis_Phasing_Grandmother(script, Father, Mother, RefSampleID, Gtypes, BAFs, parent, Parent1, EmbryoID, flip)
  }
  if (REF == "SiblingM") {
    PhBAF <- Haplarithmisis_Phasing_SiblingMale(script, Father, Mother, RefSampleID, Gtypes, BAFs, EmbryoID, flip)
  }
  if (REF == "SiblingF") {
    PhBAF <- Haplarithmisis_Phasing_SiblingFemale(script, Father, Mother, RefSampleID, Gtypes, BAFs, EmbryoID, flip)
  }
  
  # return the phased B-allele frequencies
  return(PhBAF)
}