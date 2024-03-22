# PreTestReportData_mendInc_Grandparent - Function to calculate Mendelian Inconsistency Percentage for Grandparent Reference

PreTestReportData_mendInc_Grandparent <- function(Ref, Parent, GT, Chroms){
  # Iterate over each chromosome
  for(Chr in Chroms){
    # Subset genetic data for the current chromosome
    GT_chrom <- GT[GT$Chr %in% Chr, ]
    
    # Filter out SNPs where parent and reference information is available
    GT_chrom <- GT_chrom[GT_chrom[, Parent] != "NC" & GT_chrom[, Ref] != "NC", ]
    
    # Filter Mendelian inconsistencies
    GT_mendInc <- GT_chrom[(GT_chrom[, Parent] == "AA" & GT_chrom[, Ref] == "BB") |
                           (GT_chrom[, Parent] == "BB" & GT_chrom[, Ref] == "AA"), ]
    
    # Calculate percentage of Mendelian inconsistencies
    Percent_MendInc <- (nrow(GT_mendInc) / nrow(GT_chrom) * 100)
    
    # Store percentage of Mendelian inconsistencies for this chromosome
    if(Chr == Chroms[1]) {
      tot_Percent_MendInc <- Percent_MendInc 
    } else {
      tot_Percent_MendInc <- rbind(tot_Percent_MendInc, Percent_MendInc)}
  }
  
  # Combine chromosome names with corresponding Mendelian inconsistency percentages
  tot_Percent_MendInc <- cbind(Chroms, tot_Percent_MendInc)
  return(tot_Percent_MendInc)
}