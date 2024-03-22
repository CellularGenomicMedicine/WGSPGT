# helixData_ngs_pgd_chrom - Function to create a data frame for the informative SNPs

helixData_ngs_pgd_chrom <- function(inf_snps_pretest, Chrom){
  # Create data frames for different SNP counts
    SNP_AANT <- data.frame(ID = "Tot. aant. SNPs", Value = inf_snps_pretest[, paste("chr", Chrom, "total", sep = "_")], stringsAsFactors = FALSE)
    INF_AANT <- data.frame(ID = "Aant. inf. SNPs", Value = inf_snps_pretest[, paste("chr", Chrom, "inf", sep = "_")], stringsAsFactors = FALSE)
    CHROM <- data.frame(ID = "Chromosoom", Value = Chrom, stringsAsFactors = FALSE)
    
    # Combine the data frames to create the final data frame
    NGS_PGD_CHROM_df <- rbind(CHROM,SNP_AANT, INF_AANT)
    
    # Return the resulting data frame
    return(NGS_PGD_CHROM_df)
}
