# Haplarithmisis_disthap - Function to calculate the distance to the genetic interval of interest

Haplarithmisis_disthap <- function(Intp, Int, outPathData, EmbryoID, Chrom){
  # Extract the haplotype data from the Intp object
	HapsIntp <- Intp[["dataHap"]]
	
	# Initialize an empty data frame to store the distance information
	DistFam <- NULL
	
	# Select column in HapsInt corresponding to maternal and paternal haplotypes based on EmbryoID
	Scs <- colnames(HapsIntp)[colnames(HapsIntp) %in% paste(EmbryoID, "Mat", sep = "_") | colnames(HapsIntp) %in% paste(EmbryoID, "Pat", sep = "_")]
	
	# loop through each selected column (haplotype)
	for (ind in Scs) {
	  # Subset the Int data frame for the specified chromosome
		IntChr <- Int[Int[, 1] == Chrom, ]
		# Subset the haplotype data for the specified chromosome and selected haplotypes
		HapChr <- HapsIntp[HapsIntp[, "Chr"] == Chrom, c("Position", ind)]
		# Subset haplotype data for position up- and downstream of the interval
		UP <- HapChr[HapChr$Position <= as.numeric(IntChr[2]), ]
		Down <- HapChr[HapChr$Position >= as.numeric(IntChr[3]), ]

		# Compute run-length encoding (RLE) for haplotype values up- and downstream of the interval
		UPRle <- rle(UP[, ind])
		DownRle <- rle(Down[, ind])
		
		# Combine distance information into a matrix
		Dist <- cbind(ind, Chrom, UPRle$lengths[length(UPRle$values)], DownRle$lengths[1], UPRle$values[length(UPRle$values)], DownRle$values[1])
		Dists <- Dist

		# Append Dist to DistFam data frame
		if (ind == Scs[1]) { 
		  DistFam <- Dists 
		} else { 
		  DistFam <- rbind(DistFam, Dists)
		}
		# Print the current haplotype being processed
		print(ind)
	}# End ind loop
	
	# Set column names for the Distance data frame
	colnames(DistFam) <- c("Haplotype", "Chr", "LengthUp", "LengthDown", "ValueUp", "ValueDown")
	# Write table with distance data to file
	write.table(DistFam,paste(outPathData, "DistanceToHR.txt", sep = "/"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}
