# Haplarithmisis_ChrXFatherAB - Function to process chromosome X genetic data from the father

Haplarithmisis_ChrXFatherAB <- function(Gtypes, Father){
  # Subset the genetic data for the father on chromosome X
	GTFatherX <- Gtypes[Gtypes[, "Chr"] == "X", Father]
	
	# Check if there are any heterozygous SNP calls ("AB") in the father's genetic data on chromosome X
	if(sum(GTFatherX == "AB") >= 1){
	  # Issue a warning about the presence of heterozygous SNP calls on chromosome X
		warning(paste("On the chr. X  of the father", sum(GTFatherX == "AB"),"(",sum(GTFatherX == "AB")/sum(Gtypes$Chr == "X"), "%) of the SNP-calls are heterozygous!!"))
		
	  # Inform about treating heterozygous calls as NoCalls
	  print("These will be treated as NoCalls")
	  
	  # Change heterozygous calls ("AB") to "NC" (NoCalls) to reflect uncertainty / missing data
		GTFatherX[GTFatherX == "AB"] <- "NC"# There could not be heterozygous SNP calls on chromosome X of the father
	}
	
	# Return the modified father's genetic data on chromosome X
	return(GTFatherX)
}
