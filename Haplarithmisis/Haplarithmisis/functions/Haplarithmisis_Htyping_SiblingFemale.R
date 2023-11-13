Haplarithmisis_Htyping_SiblingFemale <- function(script,Father,Mother,RefSampleID, ParScore,Gtypes, Window, Int,Chroms,EmbryoID){
  Haps <- vector("list",2)
  names(Haps) <- c("dataHapRaw","dataHap")
  GTEmbryo <- Gtypes[,EmbryoID]
  GTEmbryo <- as.matrix(GTEmbryo)
  colnames(GTEmbryo) <- EmbryoID

  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_htypingAutOpt_SiblingFemale.R",sep="/"))
  HapsAut <- Haplarithmisis_htypingAutOpt_SiblingFemale(script,Gtypes,Father,Mother,RefSampleID,GTEmbryo,EmbryoID)
  print("Autosmes are haplotyped")

  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_chrxhtypingOpt_SiblingFemale.R",sep="/"))
  HapsChrX <- Haplarithmisis_chrxhtypingOpt_SiblingFemale(script,Father,Mother,RefSampleID, Gtypes,ParScore,EmbryoID)
  Haps[["dataHapRaw"]] <- rbind(HapsAut,HapsChrX)

  source(paste(script, "analyses", "Haplarithmisis", "functions","Haplarithmisis_inthapnew1.R",sep="/"))
  Haps[["dataHap"]]    <- Haplarithmisis_inthapnew1(script,Haps[["dataHapRaw"]],Window,Int,Chroms)
  return(Haps)

}#end function

