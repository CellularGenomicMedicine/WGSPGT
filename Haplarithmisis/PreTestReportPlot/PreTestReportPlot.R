#' PreTestReportPlot - Read in all the genetic data and generate mendelian inconsistency dataframe and number of informative SNPs of the region of interest

# Retrieve command line arguments
args	<- commandArgs(TRUE)
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
	source(paste(script, "analyses", "PreTestReportPlot", "functions","PreTestReportPlot_bafplots.R", sep="/"))
	source(paste(script, "analyses", "functions", "XYplot.R", sep = "/"))
	
	Chroms	<- c(Chroms,"Y")

	outPath  <- checkDirExistsAndCreate(paste(inputPath, "PreTestReport", sep = "/"))

	#Check and load data files
	logRs <- checkFileExistsAndFread(paste(logRData, paste(family, "_LogR.txt", sep = ""), sep="/"), inputPath)
	SeglogRs <- checkFileExistsAndFread(paste(logRData, paste(family, "_SegLogR.txt", sep = ""), sep = "/"), inputPath)
	baf	<- checkFileExistsAndFread(paste(dataPath, paste(family, "_BAF.txt", sep=""), sep = "/"), inputPath)

	# Read metainfo from samplesheet
	samplesheet <- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
	# Process family members
	fam_members	<- samplesheet[,c("SampleID","Sample_MetaInfo","Sample_Status")]
	# Retrieve the sampleIDs
	colID <- fam_members[,"SampleID"]

	# Process data for each parent
	for(parent in parents){
	  # Load interval data
		Interval <- paste(inputPath,paste(family, parent, "intervals.txt", sep = "_"), sep = "/")
		Int	<- read.table(Interval, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
		upstream_boundary <- as.numeric(Int[Int[, 1] == Chrom, 3])
		downstream_boundary	<- as.numeric(Int[Int[, 1] == Chrom, 2])

		Father <- fam_members[grepl("Father", fam_members[, "Sample_MetaInfo"]), 1]
		Mother <- fam_members[grepl("Mother", fam_members[, "Sample_MetaInfo"]), 1]
		Refs <- fam_members[!grepl("Mother", fam_members[, "Sample_MetaInfo"]) & !grepl("Father", fam_members[,"Sample_MetaInfo"]), 1]
		Parent <- fam_members[grepl(parent, fam_members[, "Sample_MetaInfo"]), 1]
		baf_int <- subset(baf, Chr == Chrom)
		
		# Generate BAF plots
		PreTestReport_bafplots(script, Chrom, upstream_boundary, downstream_boundary, ChrsLength, ideogram, outPath, colID, fam_members, baf_int, parent)

		# Plot XY logRs
		XYplot(script, logRs, SeglogRs, ChrsLength, ideogram, outPath, colID, fam_members)
	}
}, error = function(e) {
  # Write error information to file and stop execution
	write(paste0("PreTestReport gave the following error: ", e), errorFilePath, sep = "")
	stop(paste0("PreTestReport gave the following error: ", e))
})