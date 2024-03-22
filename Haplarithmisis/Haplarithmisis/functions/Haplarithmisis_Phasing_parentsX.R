# Haplarithmisis_Phasing_parentsX - Function to phase parental genotypes on chromosome X

Haplarithmisis_Phasing_parentsX <- function(Father, Mother, GTFatherX, GTMotherX, GTRefX, RefSex){
  # Check the sex of the reference individual
  if(RefSex == "male"){
    # Update GTMotherX based on specific conditions for male sibling reference
		GTMotherX[GTRefX == "BB" & GTMotherX == "AB"] <- "BA"
	} else if (RefSex == "female"){
	  # Update GTMotherX based on specific conditions for female sibling reference
		GTMotherX[GTRefX == "BB" & GTFatherX == "BB" & GTMotherX == "AB"] <- "BA"
		GTMotherX[GTRefX == "AB" & GTFatherX == "AA" & GTMotherX == "AB"] <- "BA"
	} else {
	  # Print a message if the sex of the reference individual is not determined
		print("Sex of ref. individual is not detemined and chr. X htype is not reliable")
	}
  
  # Update GTMotherX for other conditions
	GTMotherX[GTRefX == "NC" & GTMotherX == "AB"] <- "NC"

	# Create a data frame with phased parental genotypes
	Parents <- data.frame(GTFatherX, GTMotherX, stringsAsFactors = FALSE, check.names = FALSE)
	names(Parents) <- c(Father, Mother)
	
	return(Parents)
}#end function
