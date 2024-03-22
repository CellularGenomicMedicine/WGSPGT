# Haplarithmisis_parentoforigin - Function to perform parent of origin analysis

Haplarithmisis_parentoforigin <- function(Father, Mother, Gtypes, Family, EmbryoID, outPath, Gamma_value, parent){
	
	print("Parent of origin analysis...")
  
  # Extract genotypes of parents and embryo
	PatGtype <- Gtypes[,Father]
	MatGtype <- Gtypes[,Mother]
	ChildGtype <- Gtypes[,EmbryoID]
	
	# Initialize matrix to store parent of origin information
	P <- matrix(0,nrow = nrow(Gtypes), ncol = 1)
	
	# Assign parent of origin labels on genotype combinations
	P[MatGtype == "AA" & PatGtype == "BB" & ChildGtype == "AA"] <- 1
	P[MatGtype == "AA" & PatGtype == "BB" & ChildGtype == "BB"] <- -1
	P[MatGtype == "BB" & PatGtype == "AA" & ChildGtype == "AA"] <- -1
	P[MatGtype == "BB" & PatGtype == "AA" & ChildGtype == "BB"] <- 1
	P[MatGtype == "AB" & PatGtype == "AA" & ChildGtype == "BB"] <- 0.5
	P[MatGtype == "AA" & PatGtype == "AB" & ChildGtype == "BB"] <- -0.5
	P[MatGtype == "BB" & PatGtype == "AB" & ChildGtype == "AA"] <- -0.5
	P[MatGtype == "AB" & PatGtype == "BB" & ChildGtype == "AA"] <- 0.5

	# Set column names of the matrix to the IDs of the embryos
  colnames(P) <- EmbryoID
  
  # Store matrix P in variable Ps
  Ps <- P
  
  # Print the IDs of the embryos
  print(EmbryoID)

  # Create a dataframe combining genotype data with parent of origin information
	dataPo <- data.frame(Gtypes[, c("Names","Chr","Position")], Ps, stringsAsFactors = FALSE, check.names = FALSE)
	
	# Write a data frame of dataPo
	write.table(dataPo, paste(outPath, paste0(parent, "_", EmbryoID, "_", Family,"_", Gamma_value, "_dataPo", ".txt"), sep = "/"), 
	            row.names = FALSE, quote = FALSE, sep = "\t")
	
	# Return the data frame containing parent of origin information
	return(dataPo)
} #end Haplarithmisis_po2 function