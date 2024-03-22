# ConvertGenotype - read in the family-wise vcf file generated via GATK Joint-genotyping and convert it into a dataframe containing genotypes

# Input: family-wise vcf file
# Output: 
#   Genomic coordinates (Crd)
#   Genotypes (GT)
#   B-Allele Frequencies (BAF)
#   Read Depth (DP)

# Retrieve command line arguments
args <- commandArgs(TRUE)
config_file_fam <- args[1]
errorFilePath <- args[2]

# Print information about error file
print(paste("Potential errors written to:", errorFilePath, sep = " "))

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
library(Biobase)

# Set options
options(scipen=999)

# Execute the main script with error handling
tryCatch({
  # Source configuration file
	source(as.character(config_file_fam))
  
  # Load all the scripts
  # Source utility script to check if file exists
  source(paste(script, "analyses", "functions","checkFileExists.R", sep="/"))
  # Source the script to calculate the BAF.
  source(paste(script, "analyses", "ConvertGenotype", "functions","calculateBAF.R", sep = "/"))
  # Source script to count the number of genotypes per chromosome
  source(paste(script, "analyses", "functions", "validateGtype.R", sep = "/"))
  # Source utility script to write data
  source(paste(script, "analyses", "functions", "writeData.R", sep = "/"))
  
  # Append "Y" to the Chroms vector
	Chroms <- c(Chroms, "Y")
	
	# Construct the path to the vcf file
	gvcf_file <- paste(GATK, "/", family, ".g.vcf", sep = "")
	
	# Check if vcf file exists
	checkFileExists(gvcf_file,inputPath)

	# Read metainfo from samplesheet
	samplesheet <- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
	
	# Process family members and parents
	fam_members	<- samplesheet[,c("SampleID", "Sample_MetaInfo", "Sample_Status")]
	fam_members_nonEmbryo	<- fam_members[!fam_members[, "Sample_Status"] %in% "E", ]
	fam_members_Embryo <- fam_members[fam_members[, "Sample_Status"] %in% "E", ]
	
	# Retrieve the sampleIDs
	colID <- fam_members[,"SampleID"] 
	
	# Read VCF file using read.vcfR for Genomic Coordinates (Crd) and Alleles
	print("Started reading vcfR for Crd and Alleles")
	vcf <- read.vcfR (gvcf_file, verbose = FALSE) # Read vcf file
	Crd <- data.frame(Chr = gsub("chr","", getCHROM(vcf)), Position = as.numeric(getPOS(vcf)), stringsAsFactors = FALSE) # Extract genomic coordinates
	Names <- paste0("chr", paste(Crd[,1], Crd[,2], sep=":")) # Generate names from genomic coordinates
	Crd	<- data.frame(Crd, Names = Names,stringsAsFactors = FALSE) # Create dataframe for genomic coordinates
	Alleles <- data.frame(as.character(getREF(vcf)), getALT(vcf), stringsAsFactors = FALSE) # Extract reference and alternative alleles
	Alleles	<- data.frame(Alleles, Names = Names, stringsAsFactors = FALSE) # Create dataframe for alleles

	# Read VCF file using readVcf to retrieve Read Depth information (DP)
	print("Determine SNP coverage")
	DP <- extract.gt(vcf, element='DP') # Extract read depth information
	DP <- data.frame(Names = Names, DP, stringsAsFactors = FALSE, check.names = FALSE) # Create dataframe for read depth
	CrdDP	<- merge(Crd, DP, by = "Names") # Merge genomic coordinates with read depth information
	CrdDP	<- CrdDP[order(CrdDP$Chr,CrdDP$Position),] # Order the dataframe by chromosome and position

	# Read VCF file using readVcf to retrieve genotypes (GT)
	print("Started reading Vcf for GT")
	vcf1 <- readVcf(gvcf_file, verbose = FALSE) # Read vcf file
	res <- genotypeToSnpMatrix(vcf1, uncertain = FALSE) # Convert genotypes to SNP matrix
	GT <- as.data.frame(t(as(res$genotype, "character")), stringsAsFactors = FALSE, check.names = FALSE) # Convert SNP matrix to dataframe
	
	# Clean and process GT data
	for(col in 1:ncol(GT)){
		GT[,col] <- gsub("/", "", GT[,col]) # Remove slashes
		GT[,col] <- gsub("NA", "NC", GT[,col]) # Replace NA values with NC
	}
	GT$ID <- row.names(GT) # Add ID column as GT rownames
	GT$Names <- sapply(GT$ID,function(w) { unlist(strsplit(w, "_"))[1] }) # Extract names from ID column
	GTcombine	<- GT # Store GT dataframe
	CrdGT	<- merge(Crd, GT, by = "Names") # Merge genomic coordinates with genotypes
	CrdGT	<- CrdGT[order(CrdGT$Chr,CrdGT$Position),] # Order the dataframe by chromosome and position
	colnames(GT) <- gsub("$", ".GType", colnames(GT)) # Append ".GType" to column names of GT dataframe
	
	# Calculate BAF
	print("Started calculating BAF")
	AD <- extract.gt(vcf, element = 'AD') # Extract allelic depths
	baf <- calculateBAF(AD, colID, fam_members) # Calculate BAF
	bafcombine <- baf # Store BAF dataframe
	baf <- data.frame(Names = Names, baf, stringsAsFactors = FALSE, check.names = FALSE) # Create dataframe for BAF
	Crdbaf <- merge(Crd, baf, by = "Names")
	Crdbaf <- Crdbaf[order(Crdbaf$Chr,Crdbaf$Position),]#filter on GT
	colnames(baf) <- paste0(colnames(baf), ".B Allele Freq") # Append ".B Allele Freq" to column names of BAF dataframe
	
	# Combine data files
	Gtype	<- data.frame(Crd, GTcombine, bafcombine, stringsAsFactors = FALSE, check.names = FALSE)
	
	# Perform filtering
	Gtype <- Gtype[nchar(Alleles[,1]) == 1 & nchar(Alleles[,2]) == 1,] 	# Only select genotypes that have a reference and alternative allele that is one nucleotide in length
	
	# Perform Data validation
	# If the number of genotypes for a chromosome is lower than 25, a txt file will be generated.
	# If the number of genotypes for the chromosome of interest is lower than 25, the script will stop.
	validateGtype(Gtype,Chroms,inputPath,Chrom)
	
	
	Crd <- Crd[Crd$Names %in% Gtype$Names, ] # Filter Genomic coordinates based on filtered genotypes
	Crdbaf <- Crdbaf[Crdbaf$Names %in% Gtype$Names, ] # Filter BAF based on filtered genotypes
	CrdGT <- CrdGT[CrdGT$Names %in% Gtype$Names, ] # Filter genotypes based on filtered genotypes
	CrdDP <- CrdDP[CrdDP$Names %in% Gtype$Names, ] # Filter read depth based on filtered genotypes
	
	# Write Genotype data
	write.table(Crd, paste(dataPath, paste(family, "_", "Crd", ".txt", sep = ""), sep = "/"), row.names = FALSE, quote = FALSE, sep = "\t") # Write genomic coordinates to file
	writeData(CrdGT, family, fam_members, parents,"GT", dataPath) # Write genotypes to file
	writeData(Crdbaf, family, fam_members, parents,"BAF", dataPath) # Write BAF to file
	writeData(CrdDP, family, fam_members, parents,"DP", dataPath) # Write read depth to file

}, error = function(e) {
  # Handle errors
	write(paste0("ConvertGenotype gave the following error: ", e), errorFilePath, sep = "")
	stop(paste0("ConvertGenotype gave the following error: ", e))
})