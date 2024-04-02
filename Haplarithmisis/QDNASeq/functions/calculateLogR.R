# calculateLogR - function to calculate logR from .bam file using functions from the QDNAseq package
# bins: data defining the genomic regions. 
# colID: column identifier: SampleID. 
# filterbamFiles: matrix containing the bam file paths

calculateLogR <- function(bins, colID, filterbamFiles){
  
  # call the function binReadCounts to calculate read counts for each bin
	readCounts <- binReadCounts(bins, bamfiles = filterbamFiles[, 1], chunkSize = 1e8)
	
	# Apply filters to the read count data, filter residuals of 4.0 s.d. and do not overlap with blacklisted regions
	readCounts <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)
	
	# Estimate correction for the read counts data using GC content and mappability as predictors in a loess regression model 
	readCounts <- estimateCorrection(readCounts)
	
	# Apply filters and ignore specific chromosomes
	readCounts <- applyFilters(readCounts, chromosomes = NA)
	
	# Correct the read counts for GC content and mappability to obtain copy numbers
	copyNumbers <- correctBins(readCounts)
	
	# Normalize the copy numbers using median normalization
	copyNumbersNormalized <- normalizeBins(copyNumbers)
	
	# Smooth out any outlier bins in the normalized copy numbers
	copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
	
	# Extract data assay element from the smoothed copy numbers
	dat <- assayDataElement(copyNumbersSmooth, "copynumber")
	
	# Set the column names of data to the corresponding SampleID
	colnames(dat) <- colID
	
	# Create a data frame logRs containing the log ratios
	logRs <- data.frame(
	  Names = featureNames(copyNumbersSmooth), # feature names
	  Chr = fData(copyNumbersSmooth)$chromosome, # chromosome
	  Position = as.integer(fData(copyNumbersSmooth)$start), # Start position
	  as.data.frame(dat, stringsAsFactors = FALSE), # assay data
	  check.names = FALSE, 
	  stringsAsFactors = FALSE
	  )
	
	# Remove rows with missing values in the SampleID column
	logRs	<- logRs[rowSums(!is.na(logRs[, colID])) == ncol(dat), ]
	
	# Round the values in column 'colID' to three decimal places after applying log2 transformation to the copy number.
	# To avoid undefined logR in case of a copy number of 0, 0.0005 is added to all copy numbers before log transformation.
	logRs[, colID] <- round(log2(logRs[, colID] + 0.0005), digits = 3)
	
	# Return the resulting data frame 'logRs'
	return(logRs)
}
