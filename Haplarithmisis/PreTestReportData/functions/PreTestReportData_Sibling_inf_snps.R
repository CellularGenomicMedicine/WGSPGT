# PreTestReportData_Sibling_inf_snps - Function to calulcate informative SNPs when sibling is reference

PreTestReportData_Sibling_inf_snps <- function(script, GT, Parent, Father, Mother, REF, Ref, Chroms, Int, fam_members, dbsnp_path, parent, HelixoutPath, Chrom){
	
  # Load all the scripts needed
  source(paste(script, "analyses", "PreTestReportData", "functions", "PreTestReportData_snp_origin_call_sibling.R", sep = "/"))
  source(paste(script, "analyses", "PreTestReportData", "functions", "chrom_int_2mb.R", sep = "/"))
  
  # Initialize lists and data frames
  pretest_list <- list()
	inf_snps_pretest_ref <- data.frame(Ref = Ref, stringsAsFactors = FALSE)
	
	# Iterate over each chromosome
	for (Chr in Chroms){
	  
	  # Subset genetic data for the current chromosome
	  GT_chr <- GT[GT$Chr %in% Chr, ]
		
	  # Filter out SNPs where parent information is not available
		GT_chr 	<- GT_chr[GT_chr[, Father] != "NC" & GT_chr[, Mother] != "NC", ]
		
		# Calculate total SNPs for the current chromosome
		inf_snps_pretest_ref_total <- data.frame(SNPS = nrow(GT_chr), stringsAsFactors = FALSE)
		colnames(inf_snps_pretest_ref_total)	<- paste("chr", Chr, "total", sep = "_")
		inf_snps_pretest_ref <- data.frame(inf_snps_pretest_ref, inf_snps_pretest_ref_total)

		# Calculate semi-informative SNPs for chromosome X
		if( Chr == "X") {
			GT_chr_semi_inf <- GT_chr[(GT_chr[, Father] == "AB" & 
			                           GT_chr[, Mother] == "AA") | (GT_chr[, Father] == "AB" & 
			                           GT_chr[, Mother] == "BB") | (GT_chr[, Father] == "AA" & 
			                           GT_chr[, Mother] == "AB") | (GT_chr[, Father] == "BB" & 
			                           GT_chr[, Mother] == "AB") & (GT_chr[, Ref] != "NC"), ]
		} else {
		  # Calculate semi-informative SNPs for all autosomal chromosomes 
			GT_chr_semi_inf <- GT_chr[(GT_chr[, Father] == "AB" & 
			                           GT_chr[, Mother] == "AA" & 
			                           GT_chr[, Ref]!="BB" & 
			                           GT_chr[, Ref]!="NC") | (GT_chr[, Father] == "AB" & 
			                           GT_chr[, Mother] == "BB" & 
			                           GT_chr[, Ref] != "AA" & 
			                           GT_chr[, Ref] != "NC") | (GT_chr[, Father] == "AA" & 
			                           GT_chr[, Mother] == "AB" & 
			                           GT_chr[, Ref] != "BB" & 
			                           GT_chr[, Ref] != "NC") | (GT_chr[, Father] == "BB" & 
			                           GT_chr[, Mother] == "AB" & 
			                           GT_chr[, Ref] != "AA" & 
			                           GT_chr[, Ref] != "NC"), ]
		}
		
		# Add semi-informative SNP count to the reference data frame
		inf_snps_pretest_chr_semi_inf <- data.frame(SNPS = nrow(GT_chr_semi_inf), stringsAsFactors = FALSE)
		colnames(inf_snps_pretest_chr_semi_inf) <- paste("chr", Chr, "semi_inf", sep = "_")
		inf_snps_pretest_ref <- data.frame(inf_snps_pretest_ref, inf_snps_pretest_chr_semi_inf)

		# Filter informative SNPs based on parental genotypes
		if (Parent == Father) {
		  Parent2 = Mother
		  }
		if (Parent == Mother) {
		  Parent2 = Father
		  }
		GT_chr_inf <- GT_chr_semi_inf[(GT_chr_semi_inf[, Parent] == "AB" ) & (GT_chr_semi_inf[, Parent2] != "AB") & (GT_chr_semi_inf[, Ref] != "NC"), ]
		
		# Assign inheritance status for the informative SNPs
		GT_chr_inf$inheritance <- "WT"
		
		Indication_Ref <- fam_members[fam_members[,"SampleID"]%in%Ref,"Sample_Status"]
		if (grepl("AS",Indication_Ref)){ GT_chr_inf <- snp_origin_call_Affected_sibling(GT_chr_inf, Parent, Parent2, Ref) }
		if (grepl("US",Indication_Ref)){ GT_chr_inf <- snp_origin_call_Unaffected_sibling(GT_chr_inf, Parent, Parent2, Ref) } 
 
		inf_snps_pretest_chr_inf <- data.frame(SNPS = nrow(GT_chr_inf),stringsAsFactors = FALSE)
		colnames(inf_snps_pretest_chr_inf) <- paste("chr", Chr, "inf", sep = "_")
		inf_snps_pretest_ref <- data.frame(inf_snps_pretest_ref, inf_snps_pretest_chr_inf)
		
		# If the chromosome is of interest, perform additional processing
		if( any(Chr == Chrom)){
		  Chr_id <- paste("Gtype_inf_snps_chr", Chr, "of", REF, sep = "_")
		  pretest_list[[ Chr_id]] <- chrom_int_2mb(script, dbsnp_path, inf_snps_pretest_ref, Chrom, Int, GT_chr_inf, parent, Mother, HelixoutPath)
		  }
	}
	
	# Calculate total genome-wide SNP count
	inf_snps_pretest_ref$total_genome <- nrow(GT)
	# Calculate total semi-informative SNPs and their ratio
	total_semi_inf_ref <- rowSums(inf_snps_pretest_ref[, grepl("semi_inf", names(inf_snps_pretest_ref))])
	inf_snps_pretest_ref <- data.frame(inf_snps_pretest_ref, total_semi_inf = total_semi_inf_ref, ratio_semi_inf = total_semi_inf_ref / nrow(GT), stringsAsFactors = FALSE)
	
	# Calculate total informative SNPs and their ratio
	total_inf_ref <- inf_snps_pretest_ref[, grepl("inf", colnames(inf_snps_pretest_ref))]
	total_inf_ref <- rowSums(total_inf_ref[, !grepl("semi", names(total_inf_ref))])
	inf_snps_pretest_ref <- data.frame(inf_snps_pretest_ref, total_inf = total_inf_ref, ratio_inf = total_inf_ref / nrow(GT), stringsAsFactors = FALSE)
	
	# Define the reference ID
	ref_id <- paste("Inf_snps_pretest", REF, sep = "_")
	
	# Store the results in the pretest_list
	pretest_list[[ref_id]] <- inf_snps_pretest_ref
	return(pretest_list)
}
