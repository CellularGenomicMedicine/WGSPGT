#' MetaInfo - Extract MetaInfo from PGT-samplesheet
#' 
#' Input: 
#' - config_file_fam.txt as PGT_config_{family}
  #' Example provided: ExampleConfigFile.txt and includes the following 
  #' - family = {family number}
  #' - script = "/path/to/scripts"
  #' - indication = "GENE_PGD" (Affected gene of interest)
  #' - onderzoeksnr = (Research number used to identify the indication including year of analysis)
  #' - analyses = "/path/to/analysisfolder" (Folder with output for each PGD family is stored)
  #' - ideogram = "Ideogram_GRCh38.rda" (Ideogram for GRCh38)
  #' - Window = 22 (Window for haplotype interpretation)
  #' - gammaSC = 300 (gamma for segmentation of logR for Embryo)
  #' - gammaMC = 50 (gamma for segmentation of logR for parents)
  #' - plateau = 100 (value for window determination of standard error)
#' - errorFilePath: provide error file path

#' Output: 
  #' - Genomic coordinates (Crd)
  #' - Genotypes (GT)
  #' - B-Allele Frequencies (BAF)
  #' - Read Depth (DP)

# Retrieve command line arguments
args <- commandArgs(TRUE)
config_file_fam <- args[1]
errorFilePath <- args[2]

print(paste("Potential errors written to:", errorFilePath, sep = " "))

# Clear any existing error file (in case of rerun)
if (file.exists(errorFilePath)) {
	warning("Error file already exists. Removing file and resuming program...")
	file.remove(errorFilePath)
}

# Execute the main script with error handling
tryCatch({
  # Source the configuration file
	source(as.character(config_file_fam))

  # Load all the scripts
  # Source utility script to check if file exists
  source(paste(script, "analyses", "functions","checkFileExists.R", sep = "/"))
  # Source utility script to check directory existence and creation
  source(paste(script, "analyses", "functions","checkDirExistsAndCreate.R", sep="/"))
  # Source script to extract information which parent is the carrier parent and who is the seed for phasing
  source(paste(script, "analyses", "MetaInfo","functions","Parameters_extraction.R", sep = "/"))
  # Source script to extract information on the genomic interval of interest
  source(paste(script, "analyses", "MetaInfo","functions","Interval_extraction.R", sep = "/"))
  # Source script to check inheritance
  source(paste(script, "analyses", "MetaInfo", "functions", "Inheritance_check.R",sep = "/"))
  
  # Replace the '/' with '_' in onderzoeksnr
	onderzoeksnummer <- gsub("/", "_", onderzoeksnr)
	write(paste0("onderzoeksnummer=", '"', onderzoeksnummer, '"'), config_file_fam, append = TRUE)

	# Construct input path
	inputPath <- paste(analyses, family, paste(onderzoeksnummer, indication, sep = "-"), sep = "/")
	
	# Add input path to config file
	write(paste0("inputPath=", '"', inputPath, '"'), config_file_fam, append = TRUE)
	
	# Construct samplesheet path
	samplesheet_path <- paste(inputPath, paste(paste(family, paste(onderzoeksnummer, indication, sep = "-"), sep = "_"), "csv", sep = "."), sep = "/")
	# Check if samplesheet exists
	checkFileExists(samplesheet_path,inputPath)
	# Add samplesheet path to config file
	write(paste0("samplesheet_path=", '"', samplesheet_path, '"'), config_file_fam, append = TRUE)
	
	# Read samplesheet data
	samplesheet <- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
	
	# Define chromosomes
	Chroms <- c(1:22, "X")
	write(paste0("Chroms=", 'c(1:22,"X")'), config_file_fam, append = TRUE)
	
	# Create output directory for converted genotype data
	outPath <- checkDirExistsAndCreate(paste(inputPath, "ConvertGenotype", sep = "/"))
	# Add output path to config file
	write(paste0("dataPath=", '"', outPath, '"'), config_file_fam, append = TRUE)
	
	# Create output directory for patient portal Helix
	HelixoutPath <- checkDirExistsAndCreate(paste(inputPath, "HelixData", sep = "/"))
	write(paste0("HelixoutPath=", '"', HelixoutPath, '"'), config_file_fam, append = TRUE)

	# Define path where GATK output is stored
	GATK <- paste(inputPath, "HaplotypeCaller", sep = "/")
	write(paste0("GATK=", '"', GATK, '"'), config_file_fam, append = TRUE)

	# Process family members and parents
	fam_members	<- samplesheet[,c("SampleID", "Sample_MetaInfo", "Sample_Status")]
	fam_members_parents	<- fam_members[fam_members[, "Sample_MetaInfo"] %in% c("Father", "Mother"),]
	if (nrow(fam_members_parents) < 2){
		write(paste0("MetaInfo gave the following error: single parent in samplesheet", e), errorFilePath, sep="")
		stop(paste0("MetaInfo gave the following error: single parent in samplesheet", e))
	}
	
	fam_members_nonEmbryo <- fam_members[!fam_members[, "Sample_Status"] %in% "E",]
	fam_members_Embryo <- fam_members[fam_members[, "Sample_Status"] %in% "E",]
	
	# Retrieve the sampleIDs
	colID <- fam_members[,"SampleID"]

	# Determine the parent of indication based on the Sample Status in the samplesheet
	if (length(fam_members[grepl("AF", fam_members[, "Sample_Status"]) | grepl("AM", fam_members[, "Sample_Status"]), "Sample_Status"]) == 2) {
		parents <- c("Father", "Mother");
		write('parents=c("Father","Mother")', config_file_fam, append = TRUE)
	} else {
		if (fam_members[fam_members[, "Sample_MetaInfo"] %in% "Father", "Sample_Status"] == "AF") {
			parents <- "Father";
			write('parents="Father"', config_file_fam, append = TRUE)
		}
		if (fam_members[fam_members[, "Sample_MetaInfo"] %in% "Mother", "Sample_Status"] == "AM") {
			parents <- "Mother";
			write('parents="Mother"', config_file_fam, append = TRUE)
		}
	}

	# Write parent IDs to config file
	write(paste0("Father=", '"', fam_members[fam_members[,"Sample_MetaInfo"] %in% "Father", "SampleID"], '"'), config_file_fam, append = TRUE)
	write(paste0("Mother=", '"', fam_members[fam_members[,"Sample_MetaInfo"] %in% "Mother", "SampleID"], '"'), config_file_fam, append = TRUE)

	# Assign the referent from the samplesheet to REF
	REF <- Parameters_extraction(inputPath, fam_members)
	write(paste0("REF=", '"', REF, '"'), config_file_fam, append = TRUE)
	if (REF == "Grandparents") {
		RefID <- fam_members[grepl("AGF", fam_members[, "Sample_Status"]) | grepl("AGM", fam_members[, "Sample_Status"]), "Sample_MetaInfo"]
		write(paste0("RefID=", '"', RefID, '"'), config_file_fam, append = TRUE)
		GrandFather <- fam_members[grepl("GF", fam_members[, "Sample_Status"]), "SampleID"]
		write(paste0("Grandfather=", '"', GrandFather, '"'), config_file_fam, append = TRUE)
		GrandMother <- fam_members[grepl("GM", fam_members[, "Sample_Status"]), "SampleID"]
		write(paste0("Grandmother=", '"', GrandMother,'"'), config_file_fam, append = TRUE)
	} else {
		RefID <- REF
		write(paste0("RefID=", '"', REF, '"'), config_file_fam, append = TRUE)
		write(paste0("RefSampleID=", '"', fam_members[fam_members[, "Sample_MetaInfo"] %in% REF, "SampleID"], '"'), config_file_fam, append = TRUE)
	}
	
	# Determine flipping based on the carrier status of the referent
	flip <- 0
	if (fam_members[grepl(RefID, fam_members[, "Sample_MetaInfo"]), "Sample_Status"] %in% "AGM" & REF != "Grandmother") { flip <- 1 }
	if (fam_members[grepl(RefID, fam_members[, "Sample_MetaInfo"]), "Sample_Status"] %in% "US") { flip <- 1 }
	if (fam_members[grepl(RefID, fam_members[, "Sample_MetaInfo"]), "Sample_Status"] %in% "UGM") { flip <- 1 }
	if (fam_members[grepl(RefID, fam_members[, "Sample_MetaInfo"]), "Sample_Status"] %in% "UGF") { flip <- 1 }
	write(paste0("flip=", flip), config_file_fam, append = TRUE)

	# Extract genomic interval for indication
	Chrom <- Intervals_extraction(inputPath, samplesheet)
	write(paste0("Chrom=", '"', Chrom, '"'), config_file_fam, append = TRUE)

	# Perform inheritance check
	Inheritance_check(fam_members, REF, Chrom)
	
	# Set paths
	if (exists("HaplarithmisisData") == FALSE) {
	  outPathData <- checkDirExistsAndCreate(paste(inputPath, "Haplarithmisis", sep = "/"))
	  write(paste0("HaplarithmisisData=", '"', outPathData, '"'), config_file_fam, append = TRUE)
	} else { outPathData <- HaplarithmisisData }
	
	
}, error = function(e) {
	write(paste0("MetaInfo gave the following error: ", e), errorFilePath, sep = "")
	stop(paste0("MetaInfo gave the following error: ", e))
})
