library(VariantAnnotation)
library(vcfR)
library(snpStats)
library(data.table)
library(QDNAseq)
library(Biobase)
library(plotrix)
options(scipen = 999)

args <- commandArgs(TRUE)
config_file_fam <- args[1]
dbsnp_path <- args[2]
EmbryoID <- args[3]
Gamma_value <- as.numeric(args[4])
errorFilePath <- args[5]

if (file.exists(errorFilePath)) {
  warning("Error file already exists. Removing file and resuming program...")
  file.remove(errorFilePath)
}

tryCatch({
  source(as.character(config_file_fam))
  load(paste(script, "analyses", "Rda", ideogram, sep = "/"))
  load(paste(script, "analyses", "Rda", "REF_24h_QC_illuminaCytoSNP12.rda", sep = "/"))

  Window <- as.numeric(gtypemodulator_window)

  source(paste(script, "analyses", "functions", "checkDirExistsAndCreate.R", sep = "/"))
  SharePath <- checkDirExistsAndCreate(paste(inputPath, EmbryoID, sep = "/"))
  outPathPlots <- checkDirExistsAndCreate(paste(SharePath, Gamma_value, "", sep = "/"))

  samplesheet <- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
  fam_members <- samplesheet[, c("SampleID", "Sample_MetaInfo", "Sample_Status")]
  fam_members_parents <- fam_members[fam_members[, "Sample_MetaInfo"] %in% c("Father", "Mother"),]
  fam_members_nonEmbryo <- fam_members[!fam_members[, "Sample_Status"] %in% "E",]
  fam_members_Embryo <- fam_members[fam_members[, "SampleID"] %in% EmbryoID,]
  fam_members <- rbind(fam_members_nonEmbryo, fam_members_Embryo)
  colID <- fam_members[, "SampleID"]

  for (parent in parents) {
    if (parent == "Mother") {
      Parent1 <- Mother
      segs <- c("M1", "M2")
    }
    if (parent == "Father") {
      Parent1 <- Father
      segs <- c("P1", "P2")
    }
    # *** Configuration  validation
    Int <- read.table(paste(inputPath, paste(family, parent, "intervals.txt", sep = "_"), sep = "/"), sep = "\t", header = T, stringsAsFactors = F)
    Int <- Int[complete.cases(Int),]
    source(paste(script, "analyses", "functions", "checkFileExistsAndFread.R", sep = "/"))
    Gtypes <- checkFileExistsAndFread(paste(HaplarithmisisData, paste(parent, EmbryoID, family,Gamma_value, "Gtypes.txt", sep = "_"), sep = "/"))
    dataHapRaw <- checkFileExistsAndFread(paste(HaplarithmisisData, paste(parent, EmbryoID, family,Gamma_value, "dataHapRaw.txt", sep = "_"), sep = "/"))
    PhBAF <- vector("list", 8)
    for (seg in c("P1", "P2", "M1", "M2", "P1Seg", "P2Seg", "M1Seg", "M2Seg")) {
      PhBAF[[seg]] <- checkFileExistsAndFread(paste(HaplarithmisisData, paste(parent, EmbryoID, family,Gamma_value, paste0(seg, ".txt"), sep = "_"), sep = "/"))
    }
    source(paste(script, "analyses", "EmbryoTestReportData", "functions", "EmbryoTestReportData_inf_snps_Rs_gtype_inheretance_Chrom_int_2mb.R", sep = "/"))
    EmbryoTestReportData_inf_snps_Rs_gtype_inheretance_Chrom_int_2mb(dbsnp_path, outPathPlots, Gtypes, PhBAF, dataHapRaw, segs, Int, Chrom, flip, parent, HelixoutPath, Gamma_value, EmbryoID)
  }
}, error = function(e) {
  write(paste0("EmbryoTestReportData gave the following error: ", e), errorFilePath, sep = "")
  stop(e)
})