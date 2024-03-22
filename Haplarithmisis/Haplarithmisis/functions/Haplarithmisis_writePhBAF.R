# Haplarithmisis_writePhBAF - Function to write phased BAF data to text files

Haplarithmisis_writePhBAF <- function(PhBAF, outPathData, parent, EmbryoID, family, Gamma_value){
  # Iterate through each segment in PhBAF
  for(PhBAF_segs in names(PhBAF)){
		print(PhBAF_segs)
		write.table(PhBAF[[PhBAF_segs]], # Data frame to write
		            paste(outPathData, paste0(parent, "_", EmbryoID, "_", family, "_", Gamma_value, "_", PhBAF_segs, ".txt"), sep = "/"), # File path
		            row.names = FALSE, 
		            quote = FALSE, 
		            sep = "\t"
		            )
	}
}
