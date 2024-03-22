# Haplarithmisis_Htyping_SiblingMale - Function for haplotyping when male sibling is referent

Haplarithmisis_Htyping_SiblingMale <- function(script, Father, Mother, RefSampleID, ParScore, Gtypes, Window, Int, Chroms, EmbryoID){
  
  # Source function for haplotyping autosomes for sbiling male
  source(paste(script, "analyses", "Haplarithmisis", "functions","Haplarithmisis_htypingAutOpt_SiblingMale.R",sep="/"))
  
  # Source function for haplotyping chromosome X for sibling male
  source(paste(script, "analyses", "Haplarithmisis", "functions","Haplarithmisis_chrxhtypingOpt_SiblingMale.R",sep="/"))
  
  # Source function for segmenting haplotypes
  source(paste(script, "analyses", "Haplarithmisis", "functions","Haplarithmisis_inthapnew1.R",sep="/"))
  
  # Initialize a list to store the haplotypes
  Haps <- vector("list",2)
  names(Haps) <- c("dataHapRaw","dataHap")
  
  # Extract genotype data for the specific embryo
  GTEmbryo <- Gtypes[,EmbryoID]
  GTEmbryo <- as.matrix(GTEmbryo)
  colnames(GTEmbryo)<-EmbryoID

  # Haplotyping autosomes for sibling male
  HapsAut <- Haplarithmisis_htypingAutOpt_SiblingMale(script,Gtypes,Father,Mother,RefSampleID,GTEmbryo,EmbryoID)
  print("Autosomes are haplotyped")
  
  # Haplotyping of chromosome X for sibling male
  HapsChrX <- Haplarithmisis_chrxhtypingOpt_SiblingMale(Father,Mother,RefSampleID,Gtypes,ParScore,EmbryoID)
  
  # Combine autosomes and chromsome X haplotyeps into dataHapRaw list element
  Haps[["dataHapRaw"]] <- rbind(HapsAut,HapsChrX)
  
  # Segment haplotypes for the combined data and store in dataHap
  Haps[["dataHap"]]    <- Haplarithmisis_inthapnew1(script,Haps[["dataHapRaw"]],Window,Int,Chroms)
  
  # Return the haplotype data in Haps
  return(Haps)
}#end function

