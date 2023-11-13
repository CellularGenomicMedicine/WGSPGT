# Haplarithmisis by Masoud
# Adapted by Jeroen
# Adapted by Kasper
# Version 4.0 (last edit: 26-04-2022)

# *** Setting parameters
args <- commandArgs(TRUE)
config_file_fam <- args[1]
EmbryoID <- args[2]
Gamma_value <- as.numeric(args[3])
errorFilePath <- args[4]

print(paste("Potential errors written to:", errorFilePath, sep = " "))
# clear any existing error file (in case of rerun)
if (file.exists(errorFilePath)) {
  warning("Error file already exists. Removing file and resuming program...")
  file.remove(errorFilePath)
}

# *** Load needed Libraries
library(limma)
library(signal)
library(plotrix)
library(MASS)
library(data.table)
library(QDNAseq)
library(Biobase)
library(plotrix)

tryCatch({
  source(as.character(config_file_fam))
  # *** Loading reference values and sources
  # *** Set default parameter settings
  gtype_window <- as.numeric(gtypemodulator_window)
  load(paste(script, "analyses", "Rda", ideogram, sep = "/"))
  load(paste(script, "analyses", "Rda", "REF_24h_QC_illuminaCytoSNP12.rda", sep = "/"))
  options(scipen = 999)
  Func <- "mean"

  samplesheet <- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
  fam_members <- samplesheet[, c("SampleID", "Sample_MetaInfo", "Sample_Status")]
  fam_members_parents <- fam_members[fam_members[, "Sample_MetaInfo"] %in% c("Father", "Mother"),]
  fam_members_nonEmbryo <- fam_members[!fam_members[, "Sample_Status"] %in% "E",]
  fam_members_Embryo <- fam_members[fam_members[, "SampleID"] %in% EmbryoID,]
  fam_members <- rbind(fam_members_nonEmbryo, fam_members_Embryo)
  colID <- fam_members[, "SampleID"]
  source(paste(script, "analyses", "functions", "checkDirExistsAndCreate.R", sep = "/"))
  if(exists("HaplarithmisisData")==FALSE){
  outPathData <- checkDirExistsAndCreate(paste(inputPath, "Haplarithmisis", sep = "/"))
  write(paste0("HaplarithmisisData=", '"', outPathData, '"'), config_file_fam, append = TRUE)
  } else { outPathData <- HaplarithmisisData }
  SharePath <- checkDirExistsAndCreate(paste(inputPath, EmbryoID, sep = "/"))
  outPathPlots <- checkDirExistsAndCreate(paste(SharePath, Gamma_value, sep = "/"))

  source(paste(script, "analyses", "functions", "checkFileExistsAndFread.R", sep = "/"))
  GC <- checkFileExistsAndFread(paste(dataPath, paste(gtype_window, "GCcontent.txt", sep = "_"), sep = "/"))

  for (parent in parents) {
    # *** Import data
    print("Loading parameters and intervals...")
    Params <- read.table(paste(inputPath, paste(family, parent, "parameters.txt", sep = "_"), sep = "/"), sep = " ", header = T, stringsAsFactors = F)
    Int <- read.table(paste(inputPath, paste(family, parent, "intervals.txt", sep = "_"), sep = "/"), sep = "\t", header = T, stringsAsFactors = F)
    GT <- checkFileExistsAndFread(paste(dataPath, paste(parent, EmbryoID, family, "GT.txt", sep = "_"), sep = "/"))
    BAF <- checkFileExistsAndFread(paste(dataPath, paste(parent, EmbryoID, family, "BAF.txt", sep = "_"), sep = "/"))

    if (parent == "Mother") { Parent1 <- Mother }
    if (parent == "Father") { Parent1 <- Father }

    # *** Data filtering
    Gtypes <- na.omit(GT)
    BAFs <- na.omit(BAF)
    Gtypes <- Gtypes[Gtypes[, "Chr"] != "Y",]
    Gtypes <- Gtypes[Gtypes[, "Chr"] != "XY",]
    #GC      <- GC[order(GC[,4]),]
    GC <- GC[GC[, 4] %in% Gtypes$Names,]
    Gtypes <- Gtypes[Gtypes$Names %in% GC[, 4],]
    #Filter Gtypes on BAFs without NA
    Gtypes <- Gtypes[Gtypes$Names %in% BAFs$Names,]
    Gtypes <- Gtypes[order(Gtypes$Chr, Gtypes$Position),]
    #Filter BAFs to match Gtypes
    BAFs <- BAFs[BAFs$Names %in% Gtypes$Names,]
    BAFs <- BAFs[order(BAFs$Chr, BAFs$Position),]
    write.table(Gtypes, paste(outPathData, paste0(parent, "_", EmbryoID, "_", family,"_",Gamma_value, "_Gtypes", ".txt"), sep = "/"), row.names = F, quote = F, sep = "\t")

    source(paste(script, "analyses", "functions", "validateGtype.R", sep = "/"))
    validateGtype(Gtypes, Chroms, outPathPlots)
    #filterout chromosomes that lack data, check for chromosome of interest is done in ConvertGenotype
    Chroms <- unique(Gtypes$Chr)
    # *** QC
    ChrPos <- Gtypes[, c("Chr", "Position")]
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_qcgtype.R", sep = "/"))
    QC <- Haplarithmisis_qcgtype(script, colID, Father, Mother, Gtypes, ChrPos, EmbryoID, outPathPlots, Chroms)
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_qcbyparents.R", sep = "/"))
    QCbyParents <- Haplarithmisis_qcbyparents(Father, Mother, Gtypes, EmbryoID, Chroms)
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_po2.R", sep = "/"))
    dataPo <- Haplarithmisis_po2(Father, Mother, Gtypes, family, EmbryoID, outPathData)
    write.table(dataPo, paste(outPathData, paste0(parent, "_", EmbryoID, "_", family,"_",Gamma_value, "_dataPo", ".txt"), sep = "/"), row.names = F, quote = F, sep = "\t")
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_patscore2.R", sep = "/"))
    ParScore <- Haplarithmisis_patscore2(Father, Mother, dataPo, QC, Chroms, Gtypes, EmbryoID)

    if (REF == "Grandparents") {
      source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Htyping_Grandparents.R", sep = "/"))
      Haps <- Haplarithmisis_Htyping_Grandparents(script, Father, Mother, Grandfather, Grandmother, Gtypes, ParScore, parent, Parent1, EmbryoID, Chroms)
      source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_Grandparents.R", sep = "/"))
      PhBAF <- Haplarithmisis_Phasing_Grandparents(script, Father, Mother, Grandfather, Grandmother, Gtypes, BAFs, parent, Parent1, EmbryoID, flip)
    } else {
      source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Htyping.R", sep = "/"))
      Haps <- Haplarithmisis_Htyping(script, Father, Mother, REF, RefSampleID, Gtypes, ParScore, Int, Window, parent, Parent1, EmbryoID, Chroms)
      source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing.R", sep = "/"))
      PhBAF <- Haplarithmisis_Phasing(script, REF, Father, Mother, RefSampleID, Gtypes, BAFs, parent, Parent1, EmbryoID, flip)
    }
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_callpcfBAF.R", sep = "/"))
    print("Started calculating SegBAF")
    PhBAF[["P1Seg"]] <- Haplarithmisis_callpcfBAF(script, PhBAF[["P1"]], Gamma_value, plateau, EmbryoID)
    PhBAF[["P2Seg"]] <- Haplarithmisis_callpcfBAF(script, PhBAF[["P2"]], Gamma_value, plateau, EmbryoID)
    PhBAF[["M1Seg"]] <- Haplarithmisis_callpcfBAF(script, PhBAF[["M1"]], Gamma_value, plateau, EmbryoID)
    PhBAF[["M2Seg"]] <- Haplarithmisis_callpcfBAF(script, PhBAF[["M2"]], Gamma_value, plateau, EmbryoID)

    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_writePhBAF.R", sep = "/"))
    Haplarithmisis_writePhBAF(PhBAF, outPathData, parent, EmbryoID, family,Gamma_value)
    write.table(Haps[["dataHap"]], paste(outPathData, paste0(parent, "_", EmbryoID, "_", family,"_",Gamma_value, "_dataHap", ".txt"), sep = "/"), row.names = F, quote = F, sep = "\t")
    write.table(Haps[["dataHapRaw"]], paste(outPathData, paste0(parent, "_", EmbryoID, "_", family,"_",Gamma_value,"_dataHapRaw", ".txt"), sep = "/"), row.names = F, quote = F, sep = "\t")
    source(paste(script, "analyses", "functions", "writeData.R", sep = "/"))
    writeData(BAFs, family, fam_members, parents, paste(Gamma_value,"BAF",sep="_"), outPathData)

    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_intphappropser3.R", sep = "/"))
    Intp <- Haplarithmisis_intphappropser3(Haps, ParScore, EmbryoID, Chroms)
    source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_disthap.R", sep = "/"))
    Haplarithmisis_disthap(Intp, Int, outPathData, EmbryoID, Chrom)
  }
}, error = function(e) {
  write(paste0("Haplarithmisis gave the following error: ", e), errorFilePath, sep = "")
  stop(e)
})