#' Haplarithmisis

#' Input:
  #' - Config file from family
  #' - Samplesheet family
  #' - Genotype data 
  #' - B-allele frequency data
  #' - Genetic interval 

#' Output:
  #' - QC data
  #'   - Call Rates
  #'   - Allele Drop-in, Allele Drop-outs
  #'   - Mendelian inconsistency
  #'   - Parent of origin information
  #'   - Parental scores
  #'   - (segmented) Haplotypes
  #'   - (segmented) Phased B-allele frequencies

# Retrieve command line arguments
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

# Load libraries
library(limma)
library(signal)
library(plotrix)
library(MASS)
library(data.table)
library(QDNAseq)
library(Biobase)
library(plotrix)

# Set options
options(scipen = 999)

# Execute the main script with error handling
tryCatch({
  # Source the configuration file
  source(as.character(config_file_fam))
  
  # Load all the scripts
  # Source utility script to check directory existence and creation
  source(paste(script, "analyses", "functions", "checkDirExistsAndCreate.R", sep = "/"))
  # Source utility script to check if file exists and read file
  source(paste(script, "analyses", "functions", "checkFileExistsAndFread.R", sep = "/"))
  # Source script to count the number of genotypes per chromosome
  source(paste(script, "analyses", "functions", "validateGtype.R", sep = "/"))
  # Source script to write data
  source(paste(script, "analyses", "functions", "writeData.R", sep = "/"))
  
  # Source script to calculate the genotype call rates and genotype error rates (Mendelian inconsistency rates)
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_qcgtype.R", sep = "/"))
  # Source script to calculate the allele drop in (ADI) and allele drop out (ADO) in the embryo determined using genotype data from parents
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_qcbyparents.R", sep = "/"))
  # Source script to calculate the parent of origin
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_parentoforigin.R", sep = "/"))
  # Source script to calculate the parental score
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_ParScore.R", sep = "/"))
  
  # Source script to perform haplotyping when reference is grandparents
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Htyping_Grandparents.R", sep = "/"))
  # Source script to Phase the parental genotypes when reference is grandparents
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing_Grandparents.R", sep = "/"))
  # Source script to perform haplotyping when reference is not grandparents
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Htyping.R", sep = "/"))
  # Source script to phase the B-allele frequencies
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_Phasing.R", sep = "/"))
  # Source script to segment the phased B-allele frequencies 
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_callpcfBAF.R", sep = "/"))
  # Source script to write the segmented phased B-allele frequencies
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_writePhBAF.R", sep = "/"))
  # Source script to interpret the hapl
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_intphappropser3.R", sep = "/"))
  source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_disthap.R", sep = "/"))
  
  
  # Loading reference values and sources
  # Set default parameter settings
  load(paste(script, "analyses", "Rda", ideogram, sep = "/"))
  # Load reference .. data
  load(paste(script, "analyses", "Rda", "REF_24h_QC_illuminaCytoSNP12.rda", sep = "/"))

  # Read Metainfo from samplesheet path
  samplesheet <- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
  
  # Process family members and parents
  fam_members <- samplesheet[, c("SampleID", "Sample_MetaInfo", "Sample_Status")]
  fam_members_parents <- fam_members[fam_members[, "Sample_MetaInfo"] %in% c("Father", "Mother"),]
  fam_members_nonEmbryo <- fam_members[!fam_members[, "Sample_Status"] %in% "E",]
  fam_members_Embryo <- fam_members[fam_members[, "SampleID"] %in% EmbryoID,]
  fam_members <- rbind(fam_members_nonEmbryo, fam_members_Embryo)
  
  # Retrieve the sampleIDs
  colID <- fam_members[, "SampleID"]
  
  # Set paths for EmbryoID that is analysed
  SharePath <- checkDirExistsAndCreate(paste(inputPath, EmbryoID, sep = "/"))
  outPathPlots <- checkDirExistsAndCreate(paste(SharePath, Gamma_value, sep = "/"))
  
  # Loop through each parent of indication
  for (parent in parents) {
    # Import data
    print("Loading parameters and intervals...")
    
    # Read interval file
    Int <- read.table(paste(inputPath, paste(family, parent, "intervals.txt", sep = "_"), sep = "/"), sep = "\t", header = T, stringsAsFactors = F) # Read intervals file
    
    # Read GT file
    GT <- checkFileExistsAndFread(paste(dataPath, paste(parent, EmbryoID, family, "GT.txt", sep = "_"), sep = "/")) # Read GT file
    
    # Read BAF file
    BAF <- checkFileExistsAndFread(paste(dataPath, paste(parent, EmbryoID, family, "BAF.txt", sep = "_"), sep = "/")) # Read BAF file

    # Set Parent1 variable based on the parent of the indication
    if (parent == "Mother") { Parent1 <- Mother }
    if (parent == "Father") { Parent1 <- Father }

    # Data filtering
    Gtypes <- na.omit(GT) # Remove rows with NA values
    BAFs <- na.omit(BAF) # Remove rows with NA values
    Gtypes <- Gtypes[Gtypes[, "Chr"] != "Y", ] # Filter out Y chromosome data
    Gtypes <- Gtypes[Gtypes[, "Chr"] != "XY", ] # Filter out XY chromosome data
    
    # Filter Gtypes on BAFs without NA
    Gtypes <- Gtypes[Gtypes$Names %in% BAFs$Names, ] # Filter Gtypes based on BAFs
    # Order Gtypes data
    Gtypes <- Gtypes[order(Gtypes$Chr, Gtypes$Position), ]
    
    # Filter BAFs to match Gtypes
    BAFs <- BAFs[BAFs$Names %in% Gtypes$Names,]
    # Order BAfs data
    BAFs <- BAFs[order(BAFs$Chr, BAFs$Position), ]
    # Write Gtypes data to file
    write.table(Gtypes, paste(HaplarithmisisData, paste0(parent, "_", EmbryoID, "_", family,"_",Gamma_value, "_Gtypes", ".txt"), sep = "/"), row.names = FALSE, quote = FALSE, sep = "\t")

    # Count the number of Gtype data points for each of the chromosomes and throw error if the number is lower than 25.
    validateGtype(Gtypes, Chroms, outPathPlots, Chrom)
    
    Chroms <- unique(Gtypes$Chr)
    
    # Quality Control
    
    # Retrieve the chromosomes and positions
    ChrPos <- Gtypes[, c("Chr", "Position")]
    
    # Perform genotype quality control on genotype call rates and genotype error rates (Mendelian inconsistency rates) and write file
    QC <- Haplarithmisis_qcgtype(script, colID, Father, Mother, Gtypes, ChrPos, EmbryoID, outPathPlots, Chroms)

    # Perform genotype quality control on genotype data using genotype data from parents to calculate ADI and ADO and write file
    QCbyParents <- Haplarithmisis_qcbyparents(Father, Mother, Gtypes, EmbryoID, Chroms, outPathPlots)
    
    # Determine the parent of origin and write table 
    dataPo <- Haplarithmisis_parentoforigin(Father, Mother, Gtypes, family, EmbryoID, HaplarithmisisData, Gamma_value, parent)
    
    # Calculate the parental score 
    ParScore <- Haplarithmisis_ParScore(Father, Mother, dataPo, QC, Chroms, Gtypes, EmbryoID)

    # Haplotyping and Phasing
    if (REF == "Grandparents") { # Check if reference is grandparents
      # Perform haplotyping for grandparents as reference
      Haps <- Haplarithmisis_Htyping_Grandparents(script, Father, Mother, Grandfather, Grandmother, Gtypes, ParScore, parent, Parent1, EmbryoID, Chroms)
      # Perform phasing of parental B-allele frequencies for grandparents as reference
      PhBAF <- Haplarithmisis_Phasing_Grandparents(script, Father, Mother, Grandfather, Grandmother, Gtypes, BAFs, parent, Parent1, EmbryoID, flip)
    } else { # If referent is not grandparents
      # Perform haplotyping
      Haps <- Haplarithmisis_Htyping(script, Father, Mother, REF, RefSampleID, Gtypes, ParScore, Int, Window, parent, Parent1, EmbryoID, Chroms)
      # Perform phasing of parental B-allele frequencies 
      PhBAF <- Haplarithmisis_Phasing(script, REF, Father, Mother, RefSampleID, Gtypes, BAFs, parent, Parent1, EmbryoID, flip)
    }
    
    # Calculate the segmented phased B-allele frequencies 
    print("Started calculating SegBAF")
    PhBAF[["P1Seg"]] <- Haplarithmisis_callpcfBAF(script, PhBAF[["P1"]], Gamma_value, plateau, EmbryoID)
    PhBAF[["P2Seg"]] <- Haplarithmisis_callpcfBAF(script, PhBAF[["P2"]], Gamma_value, plateau, EmbryoID)
    PhBAF[["M1Seg"]] <- Haplarithmisis_callpcfBAF(script, PhBAF[["M1"]], Gamma_value, plateau, EmbryoID)
    PhBAF[["M2Seg"]] <- Haplarithmisis_callpcfBAF(script, PhBAF[["M2"]], Gamma_value, plateau, EmbryoID)

    # Write table of the (segmented) phased B-allele frequencies into separate files
    Haplarithmisis_writePhBAF(PhBAF, HaplarithmisisData, parent, EmbryoID, family, Gamma_value)
    # Write table of the interpreted / segmented Haplotypes
    write.table(Haps[["dataHap"]], paste(HaplarithmisisData, paste0(parent, "_", EmbryoID, "_", family,"_",Gamma_value, "_dataHap", ".txt"), sep = "/"), row.names = FALSE, quote = FALSE, sep = "\t")
    # Write table of the interpreted / segmented Haplotypes
    write.table(Haps[["dataHapRaw"]], paste(HaplarithmisisData, paste0(parent, "_", EmbryoID, "_", family,"_",Gamma_value,"_dataHapRaw", ".txt"), sep = "/"), row.names = FALSE, quote = FALSE, sep = "\t")
    # Write table of the B-allele frequencies
    writeData(BAFs, family, fam_members, parents, paste(Gamma_value, "BAF", sep = "_"), HaplarithmisisData)

    # Process and filter computed data
    Intp <- Haplarithmisis_intphappropser3(Haps, ParScore, EmbryoID, Chroms)
    
    # Determine the distance to the genetic interval of interest
    Haplarithmisis_disthap(Intp, Int, HaplarithmisisData, EmbryoID, Chrom)
  }
}, error = function(e) {
  write(paste0("Haplarithmisis gave the following error: ", e), errorFilePath, sep = "")
  stop(e)
})
