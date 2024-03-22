# filterBams - Function to filter BAM files based on family members

filterBams <- function(bamFiles,fam_members){
  # Create a list of filter BAM file names based on family members
  filter_bams <- paste0("/", fam_members[, "SampleID"], ".bam")
	filterbamFiles <- data.frame(stringsAsFactors = FALSE, row.names = FALSE)
	
	# Iterate over each BAM file
	for(filter_bam in filter_bams){
		filterbamFile <- data.frame(bamfiles = bamFiles[grep(filter_bam,bamFiles[, 1]), ], stringsAsFactors = FALSE, row.names = grep(filter_bam, filter_bams))
		# Check if it is the first filtered BAM file
		if(filter_bam == filter_bams[1]) {
			filterbamFiles <- filterbamFile
		} else {
		  # Append the filtered BAM file to the existing data frame
			filterbamFiles <- rbind(filterbamFiles,filterbamFile)
		}
	}
	# Return the filtered BAM files data frame
	return(filterbamFiles)
}
