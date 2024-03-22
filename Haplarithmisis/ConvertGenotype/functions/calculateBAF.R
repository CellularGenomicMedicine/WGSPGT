# Function to calculate B-allele frequencies (BAF) from allele depth
calculateBAF <- function(AD, colID, fam_members){
	# Convert input matrix to a dataframe
  AD2 <- data.frame(AD,stringsAsFactors=FALSE,check.names=F)
  
  # Loop through each sample in the family
	for (fam_member in colID) {
	  # Calculate BAF for the current sample
		baf2 <- sapply(AD2[, fam_member], function(w) {  
		  Y <- as.numeric(unlist(strsplit(w,",")))
		  Y2 <- Y[2]/(Y[1]+Y[2])
		  Y2
		  })
		
		# Assign BAF values to the corresponding column in the output dataframe
		if (fam_member == fam_members[1,1]) {
		  baf <- data.frame(baf2, stringsAsFactors = FALSE) 
		  names(baf) <- fam_member
		  } else {
			colnmsbaf3 <- names(baf)
			baf <- data.frame(baf, baf2, stringsAsFactors = FALSE)
			names(baf) <- c(colnmsbaf3, fam_member)
		}
	}
  
  # Return the dataframe containing the BAF values
	return(baf)
}
