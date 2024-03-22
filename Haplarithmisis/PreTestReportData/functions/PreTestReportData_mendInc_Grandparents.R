# PreTestReportData_mendInc_Grandparents - Function to calculate Mendelian Inconsistency Percentage for when grandparents are reference

PreTestReportData_mendInc_Grandparents <- function(Grandmother, Grandfather, Parent, GT, Chroms){
  # Iterate over each chromosome
  for(Chr in Chroms){
    # Subset genetic data for the current chromosome
    GT_chrom <- GT[GT$Chr %in% Chr, ]
    
    # Remove rows with missing values for relevant columns
    GT_chrom <- GT_chrom[GT_chrom[, Parent] != "NC" & GT_chrom[, Grandmother] != "NC" & GT_chrom[, Grandfather] != "NC", ]
    
    # Calculate Mendelian inconsistencies based on specific genotype patterns
    GT_mendInc <- GT_chrom[(GT_chrom[,Parent] == "AB" & GT_chrom[, Grandmother] == "AA" & GT_chrom[, Grandfather] == "AA") |
                           (GT_chrom[,Parent] == "AB" & GT_chrom[, Grandmother] == "BB" & GT_chrom[, Grandfather] == "BB") |
                           (GT_chrom[,Parent] == "AA" & GT_chrom[, Grandmother] == "BB" & GT_chrom[, Grandfather] == "BB") |
                           (GT_chrom[,Parent] == "AA" & GT_chrom[, Grandmother] == "AA" & GT_chrom[, Grandfather] == "BB") |
                           (GT_chrom[,Parent] == "AA" & GT_chrom[, Grandmother] == "BB" & GT_chrom[, Grandfather] == "AA") |
                           (GT_chrom[,Parent] == "AA" & GT_chrom[, Grandmother] == "AB" & GT_chrom[, Grandfather] == "BB") |
                           (GT_chrom[,Parent] == "AA" & GT_chrom[, Grandmother] == "BB" & GT_chrom[, Grandfather] == "AB") |
                           (GT_chrom[,Parent] == "BB" & GT_chrom[, Grandmother] == "AA" & GT_chrom[, Grandfather] == "AA") |
                           (GT_chrom[,Parent] == "BB" & GT_chrom[, Grandmother] == "AA" & GT_chrom[, Grandfather] == "BB") |
                           (GT_chrom[,Parent] == "BB" & GT_chrom[, Grandmother] == "BB" & GT_chrom[, Grandfather] == "AA") |
                           (GT_chrom[,Parent] == "BB" & GT_chrom[, Grandmother] == "AB" & GT_chrom[, Grandfather] == "AA") |
                           (GT_chrom[,Parent] == "BB" & GT_chrom[, Grandmother] == "AA" & GT_chrom[, Grandfather] == "AB"), ]
    # Calculate the percentage of Mendelian inconsistencies
    Percent_MendInc <- (nrow(GT_mendInc) / nrow(GT_chrom) * 100)
    
    # Store the percentage for the current chromosome
    if(Chr == Chroms[1]){
      tot_Percent_MendInc <- Percent_MendInc
    } else {
        tot_Percent_MendInc <- rbind(tot_Percent_MendInc, Percent_MendInc)
    }
  }
  
  # Combine chromosome names with Mendelian inconsistency percentage
  tot_Percent_MendInc <- cbind(Chroms, tot_Percent_MendInc)
  # Return the data frame containing Mendelian inconsistency percentages
  return(tot_Percent_MendInc)
}
