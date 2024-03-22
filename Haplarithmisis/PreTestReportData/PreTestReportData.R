#' PreTestReportData - Read in all the genetic data and generate mendelian inconsistency dataframe and number of informative SNPs of the region of interest

# Retrieve command line arguments
args <- commandArgs(TRUE)
config_file_fam <- args[1]
errorFilePath <- args[2]
dbsnp_path <- args[3]

# Clear any existing error file (in case of rerun)
print(paste("Potential errors written to:", errorFilePath, sep = " "))
# clear any existing error file (in case of rerun)
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
options(scipen=999)

# Execute the main script with error handling
tryCatch({
  # Source configuration file
  source(as.character(config_file_fam))

	load(paste(script,"analyses","Rda",ideogram,sep="/"))

	# Load all the scripts needed
	source(paste(script, "analyses", "functions","checkDirExistsAndCreate.R", sep="/"))
	source(paste(script, "analyses", "functions","checkFileExistsAndFread.R", sep="/"))
	source(paste(script, "analyses", "PreTestReportData","functions","PreTestReportData_mendInc.R",sep="/"))
	source(paste(script, "analyses", "PreTestReportData","functions","PreTestReportData_phasing.R",sep="/"))
	
	# Check file paths exist and create
	outPath  <- checkDirExistsAndCreate(paste(inputPath,"PreTestReport",sep="/"))

	# Load genetic data for the parents
	GT <- checkFileExistsAndFread(paste(dataPath, paste(family, "_GT.txt", sep = ""), sep = "/"), inputPath)
	GT <- na.omit(GT)
	Chroms 	<- unique(GT[, "Chr"])
	
	# Read metainfo from samplesheet
	samplesheet <- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
	# Process family members
	fam_members	<- samplesheet[,c("SampleID","Sample_MetaInfo","Sample_Status")]
	# Retrieve the sampleIDs
	colID <- fam_members[, "SampleID"]

	# Create data list for informative and genome-wide SNPs and mendelian inconsistency
	NGS_PGD_GENOMEWIDE_list <- list()
	for (parent in parents){
		# Load interval data
	  Interval <- paste(inputPath, paste(family, parent, "intervals.txt",sep = "_"), sep = "/")
		Int <- read.table(Interval, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
		Parent <- fam_members[fam_members[, "Sample_MetaInfo"] %in% parent, "SampleID"]

		# Calculate the Mendelian inconsistency
		PreTestReportData_mendinc(script, config_file_fam, Parent, REF, Father, Mother, GT, outPath, Chroms)
		# Phasing and extracting informative SNPs
		NGS_PGD_GENOMEWIDE_list <- PreTestReportdata_phasing(script, config_file_fam, GT, Parent, Father, Mother, REF, Chroms, Int, fam_members, dbsnp_path, parent, HelixoutPath, Chrom, NGS_PGD_GENOMEWIDE_list)
	}
	
	NGS_PGD_GENOMEWIDE_df <- rbind(NGS_PGD_GENOMEWIDE_list[["TOT_AANT"]], NGS_PGD_GENOMEWIDE_list[["NGS_PGD_GENOMEWIDE_M_df"]], NGS_PGD_GENOMEWIDE_list[["NGS_PGD_GENOMEWIDE_P_df"]])
	write.table(NGS_PGD_GENOMEWIDE_df, paste(paste(HelixoutPath, "", sep = "/"), Mother, "-NGS_PGD_GENOMEWIDE", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")

}, error = function(e) {
	write(paste0("PreTestReport gave the following error: ", e), errorFilePath, sep = "")
	stop(paste0("PreTestReport gave the following error: ", e))
})