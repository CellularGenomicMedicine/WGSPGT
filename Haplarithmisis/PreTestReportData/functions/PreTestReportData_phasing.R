# PreTestReportdata_phasing - Function to perform phasing and extract informative SNPs

PreTestReportdata_phasing <- function(script, config_file_fam, GT, Parent, Father, Mother, REF, Chroms, Int, fam_members, dbsnp_path, parent, HelixoutPath, Chrom, NGS_PGD_GENOMEWIDE_list){
  # Load all the scripts needed
  source(paste(script, "analyses", "PreTestReportData", "functions", "PreTestReportData_Grandparents_inf_snps.R", sep = "/"))
  source(paste(script, "analyses", "PreTestReportData", "functions", "PreTestReportData_recessive_inf_snps.R", sep = "/"))
  source(paste(script, "analyses", "PreTestReportData", "functions", "PreTestReportData_Sibling_inf_snps.R", sep = "/"))
  source(paste(script, "analyses", "PreTestReportData", "functions", "PreTestReportData_Grandparent_inf_snps.R", sep = "/"))
  source(paste(script, "analyses", "PreTestReportData", "functions", "PreTestReportData_helixData_ngs_pgd_chrom.R", sep = "/"))
  
  # Check the reference sample type and perform necessary operations
  if (REF == "Grandparents"){
			source(as.character(config_file_fam))
			inf_snps_pretest_ref_list <- PreTestReportData_Grandparents_inf_snps(script, GT, Parent, Father, Mother, Grandfather, Grandmother, Chroms, Int, fam_members, dbsnp_path, parent, HelixoutPath, Chrom)
		}

	if (REF == "SiblingF" | REF == "SiblingM"){
	  if (length(fam_members[grepl("AF", fam_members[, "Sample_Status"]) | grepl("AM", fam_members[, "Sample_Status"]), "Sample_Status"]) == 2) {
	    source(as.character(config_file_fam))
			inf_snps_pretest_ref_list <- PreTestReportData_recessive_inf_snps(script, GT, Parent, Father, Mother, REF, RefSampleID, Chroms, Int, fam_members, dbsnp_path, parent,HelixoutPath,Chrom)
			} else {
			  source(as.character(config_file_fam))
			  inf_snps_pretest_ref_list <- PreTestReportData_Sibling_inf_snps(script, GT, Parent, Father, Mother, REF, RefSampleID, Chroms, Int, fam_members, dbsnp_path, parent, HelixoutPath, Chrom)
			}
	}
  
  if (REF == "Grandmother" | REF == "Grandfather"){
    source(as.character(config_file_fam))
		inf_snps_pretest_ref_list <- PreTestReportData_Grandparent_inf_snps(script, GT, Parent, Father, Mother, REF, RefSampleID, Chroms, Int, fam_members, dbsnp_path, parent, HelixoutPath, Chrom)
  }
  
  
  # Extract necessary information from the results
  pretest_Ref <- grep(paste("pretest", REF, sep = "_"), names(inf_snps_pretest_ref_list))
	inf_snps_pretest <- inf_snps_pretest_ref_list[[pretest_Ref]]
	Chrom_gtype <- grep(paste("chr", Chrom, sep = "_"), names(inf_snps_pretest_ref_list))
	
	# Write results to file
	write.table(inf_snps_pretest_ref_list[[Chrom_gtype]], paste(paste(outPath, "", sep = "/"), parent, "_Informative_Gtypes_RS_nrs", REF, "_2MB_chroms_of_interest_chr_", Chrom, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

	# Create dataframe for informative SNPs per chromosome and link to Helix
	NGS_PGD_CHROM_df <- helixData_ngs_pgd_chrom(inf_snps_pretest, Chrom)
	if (parent == "Mother"){
	  write.table(NGS_PGD_CHROM_df, paste(paste(HelixoutPath, "", sep = "/"), Mother, "-NGS_PGD_CHROM_MAT", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")
	}
	if (parent == "Father"){
	  write.table(NGS_PGD_CHROM_df, paste(paste(HelixoutPath, "", sep = "/"), Mother, "-NGS_PGD_CHROM_PAT", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ",")
	}

	# Create dataframe for genome-wide number of informative SNPs and link to Helix
	TOT_AANT <- data.frame(ID = "Tot. aantal SNPs", Value = inf_snps_pretest[, "total_genome"], stringsAsFactors = FALSE)
	if (parent == "Mother") {
	  INF_AT_M <- data.frame(ID = "Inform. SNPs M",Value = inf_snps_pretest[,"total_inf"], stringsAsFactors = FALSE)
		INF_PC_M <- data.frame(ID = "% Inf. SNPs M",Value = round(inf_snps_pretest[,"ratio_inf"] * 100, 2), stringsAsFactors = FALSE)
		NGS_PGD_GENOMEWIDE_list[["NGS_PGD_GENOMEWIDE_M_df"]] <- rbind(INF_AT_M, INF_PC_M)
	}
	if (parent == "Father") {
	  INF_AT_P <- data.frame(ID = "Inform. SNPs P", Value = inf_snps_pretest[, "total_inf"], stringsAsFactors = FALSE)
		INF_PC_P <- data.frame(ID = "% Inf. SNPs P", Value = round(inf_snps_pretest[, "ratio_inf"] * 100, 2), stringsAsFactors = FALSE)
		NGS_PGD_GENOMEWIDE_list[["NGS_PGD_GENOMEWIDE_P_df"]] <- rbind(INF_AT_P,INF_PC_P)
	}
	
	# Store total SNP count data in the genome-wide list
	NGS_PGD_GENOMEWIDE_list[["TOT_AANT"]] <- TOT_AANT
	# Return the updated NGS_PGD_GENOMEWIDE_list
	return(NGS_PGD_GENOMEWIDE_list)
}