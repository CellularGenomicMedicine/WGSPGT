#' QDNASeq - Calculate the copy number variation for each sample
#' 
#' Input:
#'  - config_file_fam.txt as PGT_config_{family}
#'  - errorFilePath: provide error file path
  #' The config_file_fam links to filepath where 
#' Output:
#' - dataframe with LogR values for each sample per column
#' - dataframe with Segmented LogR values for each sample per column

# Retrieve command line arguments
args <- commandArgs(TRUE)
config_file_fam <- args[1]
errorFilePath <- args[2]

# Clear any existing error file (in case of rerun)
print(paste("Potential errors written to:", errorFilePath, sep = " "))
# clear any existing error file (in case of rerun)
if (file.exists(errorFilePath)) {
  warning("Error file already exists. Removing file and resuming program...")
  file.remove(errorFilePath)
}

# Load libraries
library(data.table)
library(QDNAseq)
library(Biobase)

# Set options
options(scipen=999)

# Execute the main script with error handling
tryCatch({
  # Source configuration file
	source(as.character(config_file_fam))

  # Load all the scripts
  # Source utility script to check if dir exists and create
  source(paste(script, "analyses", "functions", "checkDirExistsAndCreate.R", sep = "/"))
  # source script to filter the bam files from the family members in this PGT family
  source(paste(script, "analyses", "QDNASeq", "functions", "filterBams.R", sep = "/"))
  # Source script to calculate the logR values
  source(paste(script, "analyses", "QDNASeq", "functions", "calculateLogR.R", sep = "/"))
  # Source utility script to write data
  source(paste(script, "analyses", "functions","writeData.R", sep = "/"))
  # Source script to segmented the logR values to Segmented logR values (SegLogRs)
  source(paste(script, "analyses", "QDNASeq", "functions", "SegLogRs.R", sep = "/"))
  
  # Load the reference dataframes
	load(paste(script, "analyses", "Rda", ideogram, sep = "/"))
	load(paste(script, "analyses", "Rda", QDNASeqBins, sep = "/"))

	# Append "Y" to the Chroms vector
	Chroms <- c(Chroms, "Y")
	
	outPath <- checkDirExistsAndCreate(paste(inputPath, "QDNASeq", sep = "/"))
	write(paste0("logRData=", '"', outPath, '"'), config_file_fam, append = TRUE)
	
	# Read Metainfo from samplesheet
	samplesheet	<- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
	
	# Process family members and parents
	fam_members	<- samplesheet[, c("SampleID", "Sample_MetaInfo", "Sample_Status")]
	fam_members_nonEmbryo <- fam_members[!fam_members[, "Sample_Status"] %in% "E",]
	fam_members_Embryo <- fam_members[fam_members[, "Sample_Status"] %in% "E",]
	
	# Retrieve the sampleIDs
	colID <- fam_members[, "SampleID"] 
	
	# Read the list of bamfiles that were used to generate the GATK file
	bamFiles <- read.table(paste(GATK, paste(family, "list", sep = "."), sep = "/"), sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

	# Filter BAM files based on family members
	filterbamFiles <- filterBams(bamFiles, fam_members)

	print("Started calculating LogR")
	# Calculate LogR values
	logRs <- calculateLogR(bins, colID, filterbamFiles)
	# Write LogR data
	writeData(logRs, family, fam_members, parents, "LogR", outPath)

	print(paste("Started calculating SegLogR for gamma", gammaMC, sep = " "))
	# Calculate segmemted LogR values
	segLogRs <- SegLogRs(script, logRs, gammaMC, gammaSC, plateau, inputPath, fam_members_Embryo)
	
	# Write segmented LogR data
	writeData(segLogRs, family, fam_members, parents, "SegLogR", outPath)
}, error = function(e) {
	write(paste0("QDNASeq gave the following error: ", e), errorFilePath, sep="")
	stop(paste0("QDNASeq gave the following error: ", e))
})
