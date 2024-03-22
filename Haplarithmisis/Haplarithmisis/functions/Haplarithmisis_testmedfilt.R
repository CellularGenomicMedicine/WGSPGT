# Haplarithmisis_testmedfilt - Function to perform median filtering for haplotypes

Haplarithmisis_testmedfilt <- function(Hap, Window, Chroms){
	# Loop through each chromosome
  for(chr in Chroms){
    # Extract haplotype data for the current chromosome
		HapChr <- Hap[Hap$Chr == chr, 4]
		
		# Define chromosome-specific window size for median filtering
		WindowChr <- round((Window * length(Hap$Chr == chr)) / length(Hap$Chr == 1))
		if((WindowChr %% 2) == 0) { 
		  WindowChr <- WindowChr + 1
		}
		
		# Perform median filtering using runmed function
    HapIntChr <- runmed(HapChr, WindowChr, "median")
    
    # Append the filtered haplotype data to the genome-wide list
		if(chr == Chroms[1]) {
		  HapIntGenome <- HapIntChr
		} else { 
		  HapIntGenome <- c(HapIntGenome, HapIntChr)}
	   
   }# End chromosome loop
  
  # Return the segmented haplotypes for the entire genome
   return(HapIntGenome)
} # End function
