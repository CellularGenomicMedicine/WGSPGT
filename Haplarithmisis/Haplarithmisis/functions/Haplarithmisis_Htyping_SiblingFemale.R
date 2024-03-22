# Haplarithmisis_Htyping_SiblingFemale - Function for haplotyping when female sibling is referent

Haplarithmisis_Htyping_SiblingFemale <- function(script,Father,Mother,RefSampleID, ParScore,Gtypes, Window, Int,Chroms,EmbryoID){
  
  # Source function for haplotyping autosomes for sibling female
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_htypingAutOpt_SiblingFemale.R",sep="/"))
  # Source function for haplotyping chromosome X for sibling female
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_chrxhtypingOpt_SiblingFemale.R",sep="/"))
  # Source function for segmenting haplotypes
  source(paste(script, "analyses", "Haplarithmisis", "functions","Haplarithmisis_inthapnew1.R",sep="/"))
  
  # Initialize a list to store the haplotypes
  Haps <- vector("list", 2)
  names(Haps) <- c("dataHapRaw", "dataHap")
  
  # Extract genotype data for the specific embryo
  GTEmbryo <- Gtypes[,EmbryoID]
  GTEmbryo <- as.matrix(GTEmbryo)
  colnames(GTEmbryo) <- EmbryoID

  # Haplotyping autosomes for sibling female
  HapsAut <- Haplarithmisis_htypingAutOpt_SiblingFemale(script, Gtypes, Father, Mother, RefSampleID, GTEmbryo, EmbryoID)
  print("Autosomes are haplotyped")

  # Haplotyping of chromosome X for sibling female
  HapsChrX <- Haplarithmisis_chrxhtypingOpt_SiblingFemale(script,Father,Mother,RefSampleID, Gtypes,ParScore,EmbryoID)
  
  # Combine autosomes and chromsome X haplotyeps into dataHapRaw list element
  Haps[["dataHapRaw"]] <- rbind(HapsAut,HapsChrX)

  # Segment haplotypes for the combined data and store in dataHap
  Haps[["dataHap"]] <- Haplarithmisis_inthapnew1(script,Haps[["dataHapRaw"]],Window,Int,Chroms)
 
  # Return the haplotype data in Haps
  return(Haps)
}#end function

