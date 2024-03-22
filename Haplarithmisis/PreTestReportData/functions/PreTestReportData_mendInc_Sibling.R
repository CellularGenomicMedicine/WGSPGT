# PreTestReportData_mendInc_Sibling - Function to calculate Mendelian Inconsistency Percentage when sibling is Reference

PreTestReportData_mendInc_Sibling <- function(Father, Mother, Ref, GT, Chroms){
  for(Chr in Chroms){
    # Subset genetic data for the current chromosome
    GT_chrom <- GT[GT$Chr %in% Chr, ]
    
    # Remove rows with missing values for relevant columns
    GT_chrom <- GT_chrom[GT_chrom[, Ref] != "NC" & GT_chrom[, Mother] != "NC" & GT_chrom[, Father] != "NC", ]
    
    # Calculate Mendelian inconsistencies based on specific genotype patterns
    GT_mendInc <- GT_chrom[(GT_chrom[, Ref] == "AB" & GT_chrom[, Mother] == "AA" & GT_chrom[, Father] == "AA") |
                           (GT_chrom[, Ref] == "AB" & GT_chrom[, Mother] == "BB" & GT_chrom[, Father] == "BB") |
                           (GT_chrom[, Ref] == "AA" & GT_chrom[, Mother] == "BB" & GT_chrom[, Father] == "BB") |
                           (GT_chrom[, Ref] == "AA" & GT_chrom[, Mother] == "AA" & GT_chrom[, Father] == "BB") |
                           (GT_chrom[, Ref] == "AA" & GT_chrom[, Mother] == "BB" & GT_chrom[, Father] == "AA") |
                           (GT_chrom[, Ref] == "AA" & GT_chrom[, Mother] == "AB" & GT_chrom[, Father] == "BB") |
                           (GT_chrom[, Ref] == "AA" & GT_chrom[, Mother] == "BB" & GT_chrom[, Father] == "AB") |
                           (GT_chrom[, Ref] == "BB" & GT_chrom[, Mother] == "AA" & GT_chrom[, Father] == "AA") |
                           (GT_chrom[, Ref] == "BB" & GT_chrom[, Mother] == "AA" & GT_chrom[, Father] == "BB") |
                           (GT_chrom[, Ref] == "BB" & GT_chrom[, Mother] == "BB" & GT_chrom[, Father] == "AA") |
                           (GT_chrom[, Ref] == "BB" & GT_chrom[, Mother] == "AB" & GT_chrom[, Father] == "AA") |
                           (GT_chrom[, Ref] == "BB" & GT_chrom[, Mother] == "AA" & GT_chrom[, Father] == "AB"), ]
    # Calculate the percentage of Mendelian inconsistencies
    Percent_MendInc <- (nrow(GT_mendInc) / nrow(GT_chrom) * 100)
    
    # Store the percentage for the current chromosome
    if(Chr == Chroms[1]) {
      tot_Percent_MendInc <- Percent_MendInc
      } else {
        tot_Percent_MendInc <- rbind(tot_Percent_MendInc, Percent_MendInc)
      }
  }
  
  # Combine chromosome names with Mendelian inconsistency percentages
  tot_Percent_MendInc <- cbind(Chroms, tot_Percent_MendInc)
  
  # Return the data frame containing Mendelian inconsistency percentages
  return(tot_Percent_MendInc)
}
