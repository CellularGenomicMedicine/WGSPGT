#' EmbryoTestReportPlot - Read in all the output data from haplarithmisis and plot it

# Retrieve command line arguments
args <- commandArgs(TRUE)
config_file_fam <- args[1]
dbsnp_path <- args[2]
EmbryoID <- args[3]
Gamma_value <- as.numeric(args[4])
errorFilePath <- args[5]

# Clear any existing error file (in case of rerun)
if (file.exists(errorFilePath)) {
  warning("Error file already exists. Removing file and resuming program...")
  file.remove(errorFilePath)
}

# Load libraries
library(VariantAnnotation)
library(vcfR)
library(snpStats)
library(data.table)
library(QDNAseq)
library(Biobase)
library(plotrix)

# Set options
options(scipen = 999)

# Execute the main script with error handling
tryCatch({
  # Source configuration file
  source(as.character(config_file_fam))
  load(paste(script, "analyses", "Rda", ideogram, sep = "/"))
  load(paste(script, "analyses", "Rda", "REF_24h_QC_illuminaCytoSNP12.rda", sep = "/"))

  # Load all the scripts needed
  source(paste(script, "analyses", "functions", "checkDirExistsAndCreate.R", sep = "/"))
  source(paste(script, "analyses", "functions", "checkFileExistsAndFread.R", sep = "/"))
  source(paste(script, "analyses", "EmbryoTestReportData", "functions", "EmbryoTestReportData_inf_snps_Rs_gtype_inheritance_Chrom_int_2mb.R", sep = "/"))
  
  # Check file paths exist and create
  SharePath <- checkDirExistsAndCreate(paste(inputPath, EmbryoID, sep = "/"))
  outPathPlots <- checkDirExistsAndCreate(paste(SharePath, Gamma_value, "", sep = "/"))

  # Read Metainfo from samplesheet
  samplesheet <- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
  # Process family members and parents
  fam_members <- samplesheet[, c("SampleID", "Sample_MetaInfo", "Sample_Status")]
  fam_members_parents <- fam_members[fam_members[, "Sample_MetaInfo"] %in% c("Father", "Mother"), ]
  fam_members_nonEmbryo <- fam_members[!fam_members[, "Sample_Status"] %in% "E", ]
  fam_members_Embryo <- fam_members[fam_members[, "SampleID"] %in% EmbryoID, ]
  fam_members <- rbind(fam_members_nonEmbryo, fam_members_Embryo)
  # Retrieve the sampleIDs
  colID <- fam_members[, "SampleID"]

  # Iterate over parents for analysis
  for (parent in parents) {
    if (parent == "Mother") {
      Parent1 <- Mother
      segs <- c("M1", "M2")
    }
    if (parent == "Father") {
      Parent1 <- Father
      segs <- c("P1", "P2")
    }

    # Load interval data and genetic data for the parent
    Int <- read.table(paste(inputPath, paste(family, parent, "intervals.txt", sep = "_"), sep = "/"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    Int <- Int[complete.cases(Int), ]
    Gtypes <- checkFileExistsAndFread(paste(HaplarithmisisData, paste(parent, EmbryoID, family, Gamma_value, "Gtypes.txt", sep = "_"), sep = "/"))
    dataHapRaw <- checkFileExistsAndFread(paste(HaplarithmisisData, paste(parent, EmbryoID, family,Gamma_value, "dataHapRaw.txt", sep = "_"), sep = "/"))
    PhBAF <- vector("list", 8)
    for (seg in c("P1", "P2", "M1", "M2", "P1Seg", "P2Seg", "M1Seg", "M2Seg")) {
      PhBAF[[seg]] <- checkFileExistsAndFread(paste(HaplarithmisisData, paste(parent, EmbryoID, family,Gamma_value, paste0(seg, ".txt"), sep = "_"), sep = "/"))
    }
    # Perform analysis and generate table with number of informative SNP information
    EmbryoTestReportData_inf_snps_Rs_gtype_inheritance_Chrom_int_2mb(dbsnp_path, outPathPlots, Gtypes, PhBAF, dataHapRaw, segs, Int, Chrom, flip, parent, HelixoutPath, Gamma_value, EmbryoID)
  }
}, error = function(e) {
  write(paste0("EmbryoTestReportData gave the following error: ", e), errorFilePath, sep = "")
  stop(e)
})