# PreTestReportData_Grandparents_inf_snps - Function to calculate informative SNPs when grandparents are referent

PreTestReportData_Grandparents_inf_snps <- function(script, GT, Parent, Father, Mother, Grandfather, Grandmother, Chroms, Int, fam_members, dbsnp_path, parent, HelixoutPath, Chrom){
  # Load all the scripts needed
  source(paste(script, "analyses", "PreTestReportData", "functions", "PreTestReportData_snp_origin_call_grandparents.R", sep = "/"))
  source(paste(script, "analyses", "PreTestReportData", "functions", "chrom_int_2mb.R", sep = "/"))
  
  # Initialize an empty list to store results
  pretest_list <- list()
  # Create a data frame to store reference information
	inf_snps_pretest_ref <- data.frame(Ref = REF, stringsAsFactors = FALSE)
	
	# Iterate over each chromosome
	for (Chr in Chroms){
	  # Subset genetic data for the current chromosome
	  GT_chr <- GT[GT$Chr %in% Chr, ]
		# Filter out SNPs where parent information is not available
		GT_chr <- GT_chr[GT_chr[, Father] != "NC" & GT_chr[, Mother] != "NC", ]
		
		# Calculate total SNps for the current chromosome
		inf_snps_pretest_ref_total <- data.frame(SNPS = nrow(GT_chr), stringsAsFactors = FALSE)
		colnames(inf_snps_pretest_ref_total) <- paste("chr", Chr, "total", sep = "_")
		# Add total SNP count to the reference data frame
		inf_snps_pretest_ref <- data.frame(inf_snps_pretest_ref,inf_snps_pretest_ref_total)

		# Filter semi-informative SNPs based on parental genotype
		GT_chr_semi_inf <- GT_chr[(GT_chr[, Father] == "AB" & GT_chr[, Mother] == "AA") | 
		                          (GT_chr[, Father] == "AB" & GT_chr[, Mother] == "BB") | 
		                          (GT_chr[, Father] == "AA" & GT_chr[, Mother] == "AB") | 
		                          (GT_chr[, Father] == "BB" & GT_chr[, Mother] == "AB"), ]
		# Filter out SNPs where grandparent information is not available
		GT_chr_semi_inf <-  GT_chr_semi_inf[GT_chr_semi_inf[, Grandfather] != "NC" & GT_chr_semi_inf[, Grandmother] != "NC", ]
		inf_snps_pretest_chr_semi_inf <- data.frame(SNPS = nrow(GT_chr_semi_inf), stringsAsFactors = FALSE)
		colnames(inf_snps_pretest_chr_semi_inf) <- paste("chr", Chr, "semi_inf", sep = "_")
		# Add semi-informative SNP count to the reference data frame
		inf_snps_pretest_ref <- data.frame(inf_snps_pretest_ref, inf_snps_pretest_chr_semi_inf)

		# informative snps
		Gtypes <- data.frame(matrix(ncol = 0, nrow = nrow(GT_chr_semi_inf)), row.names = row.names(GT_chr_semi_inf), stringsAsFactors = FALSE)
		Gtypes$Grandfather <- GT_chr_semi_inf[,Grandfather]
		Gtypes$Grandmother <- GT_chr_semi_inf[,Grandmother]
		Gtypes$AfParent <- GT_chr_semi_inf[,Parent]
		GT_chr_inf_pre <- GT_chr_semi_inf[!((Gtypes$Grandfather=="AA" & Gtypes$Grandmother=="AA" & Gtypes$AfParent=="AB") |
		                                    (Gtypes$Grandfather=="BB" & Gtypes$Grandmother=="BB" & Gtypes$AfParent=="AB") | 
		                                    (Gtypes$Grandfather=="AB" & Gtypes$Grandmother=="AB" & Gtypes$AfParent=="AB") ),]

		# Filter informative SNPs based on grandparental genotypes
		if (Parent == Father) {
		  Parent2 = Mother
		  }
		if (Parent == Mother) {
		  Parent2 = Father
		  }
		GT_chr_inf <- GT_chr_inf_pre[(GT_chr_inf_pre[, Parent] == "AB" ) & (GT_chr_inf_pre[, Parent2] != "AB"), ]
		GT_chr_inf$inheritance <- "WT"
		
		# Define the informative and uninformative SNPs based on the sample status
		Ref_inf <- fam_members[grep("AG", fam_members[, "Sample_Status"]), "SampleID"]
		Ref_Uninf <- fam_members[grep("UG", fam_members[, "Sample_Status"]), "SampleID"]
		GT_chr_inf <- snp_origin_call_grandparents(GT_chr_inf, Parent, Parent2, Ref_inf, Ref_Uninf)

		# Calculate informative SNP count for the current chromosome
		inf_snps_pretest_chr_inf  <- data.frame(SNPS = nrow(GT_chr_inf), stringsAsFactors = FALSE)
		colnames(inf_snps_pretest_chr_inf) <- paste("chr", Chr, "inf", sep = "_")
		# Add informative SNP count ot the reference data frame
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
	inf_snps_pretest_ref <- data.frame(inf_snps_pretest_ref, total_semi_inf = total_semi_inf_ref,ratio_semi_inf = total_semi_inf_ref / nrow(GT), stringsAsFactors = FALSE)
	
	# Calculate total informative SNPs and their ratio
	total_inf_ref <- inf_snps_pretest_ref[, grepl("inf", colnames(inf_snps_pretest_ref))]
	total_inf_ref <- rowSums(total_inf_ref[, !grepl("semi", names(total_inf_ref))])
	inf_snps_pretest_ref <- data.frame(inf_snps_pretest_ref,total_inf = total_inf_ref, ratio_inf = total_inf_ref / nrow(GT), stringsAsFactors = FALSE)
	# Define the reference ID
	ref_id <- "Inf_snps_pretest_Grandparents"
	# Store the results in the pretest_list
	pretest_list[[ref_id]] <- inf_snps_pretest_ref
	return(pretest_list)
}
