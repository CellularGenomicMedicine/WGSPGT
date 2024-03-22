# NGS_PGD_INFSNP - Function to generate and save INFSNP data for Helix

NGS_PGD_INFSNP <- function(GT_chr_inf_pos_start_merge,GT_chr_inf_pos_stop_merge,parent,Mother,HelixoutPath){
  # Calculate counts for SNP types within the 2Mb interval
  AFm2MB <- data.frame(ID = "Af SNP -2Mb", Value = length(which(GT_chr_inf_pos_start_merge[, "inheritance"] == "Af")), stringsAsFactors = FALSE)
  WTm2MB <- data.frame(ID = "WT SNP -2Mb", Value = length(which(GT_chr_inf_pos_start_merge[, "inheritance"] == "WT")), stringsAsFactors = FALSE)
  TOTm2MB <- data.frame(ID = "Tot SNP -2Mb", Value = nrow(GT_chr_inf_pos_start_merge), stringsAsFactors = FALSE)
  AFp2MB <- data.frame(ID = "Af SNP +2Mb", Value = length(which(GT_chr_inf_pos_stop_merge[, "inheritance"] == "Af")), stringsAsFactors = FALSE)
  WTp2MB <- data.frame(ID = "WT SNP +2Mb", Value = length(which(GT_chr_inf_pos_stop_merge[, "inheritance"] == "WT")), stringsAsFactors = FALSE)
  TOTp2MB <- data.frame(ID = "Tot SNP +2Mb", Value = nrow(GT_chr_inf_pos_stop_merge), stringsAsFactors = FALSE)
  
  # Combine SNP counts into a data frame
  NGS_PGD_INFSNP_df <- rbind(AFm2MB, WTm2MB, TOTm2MB, AFp2MB, WTp2MB, TOTp2MB)
  
  # Save data frame to a CSV file based on the parent type (Mother/Father)
  if (parent == "Mother"){
		write.table(NGS_PGD_INFSNP_df, paste(paste(HelixoutPath, "", sep = "/"), Mother, "-NGS_PGD_INFSNP_MAT",sep=""),col.names = TRUE, row.names = FALSE, quote = FALSE, sep=",")
	}
  if (parent == "Father"){
    write.table(NGS_PGD_INFSNP_df, paste(paste(HelixoutPath, "", sep = "/"), Mother, "-NGS_PGD_INFSNP_PAT",sep=""),col.names = TRUE, row.names = FALSE, quote = FALSE, sep=",")
	}
}