# validateGtype - function to validate genotype data from each chromosome
validateGtype <- function(Gtype,Chroms,outPath,Chrom){
  # Create an empty dataframe to store chromosomes with less than 25 genotype values
	GtypeLess <- data.frame("chr" = 0, length = 0, stringsAsFactors = FALSE)
	# Loop through each chromosome 
	for (chr in Chroms) {
		Gtypes <- Gtype[Gtype$Chr == chr,]
		
		# Check if the number of genotype values for the current chromosome is less than 25
		if(nrow(Gtypes) < 25) {
		  
		  # Record chromosome and number of genotype values less than 25
			GtypeLess <- rbind(GtypeLess, c(chr, length(Gtypes)));
			
			# Write chromosome information to a file
			write.table(chr, paste(outPath, paste("Less_then_25_Gtype_values_chr", chr, ".txt", sep = ""),sep = "/"), sep = "\t")
			
			# Check if the current chromosome is the chromosome of interest and has less than 25 genotype values 
			if (nrow(Gtypes) < 25 & chr==Chrom) {
			  # print message indicating chromosome of interest has less than 25 genotype values
				cat("Chromosome of interest ", chr, " has less then 25 Gtype values")
				stop(paste0("Chromosome of interest ",chr," has less then 25 Gtype"))
			}
		}
	}
}
