# chrom_int_2mb - Function to process the informative SNP information within a 2 megabase (Mb) interval on a chromosome

chrom_int_2mb <- function(script, dbsnp_path, inf_snps_pretest_ref, Chrom, Int, GT_chr_inf, parent, Mother, HelixoutPath){
	
  # Load all the scripts needed
  source(paste(script, "analyses", "PreTestReportData", "functions", "NGS_PGD_INFSNP.R", sep = "/"))
  
  # Read dbSNP data for the chromsome
  dbsnp_chr <- fread(paste(dbsnp_path, paste("chr", Chrom, sep = ""), ".bed", sep = ""), data.table = FALSE, stringsAsFactors = FALSE, colClasses = c("character", "numeric", "numeric", "character"), sep = "\t")
  names(dbsnp_chr) <- c("chr", "start", "Position", "RS_nr")
  
  # Extract start and stop positions for the 2-megabase interval
  Pos_start_int <- Int[Int[, "Chr"] %in% Chrom, "Start"]
	Pos_stop_int <- Int[Int[, "Chr"] %in% Chrom, "Stop"]
	
	# if the start position of the 2 mb downstream would be before the start position of the chromosome, the start position becomes 0 (end of the chromosome)
	if ((Pos_start_int - 2000000) < 0) { 
	  Pos2mb_start_int <- 0
	} else {
	  Pos2mb_start_int <- Pos_start_int - 2000000
	}
	
	Pos2mb_stop_int <- Pos_stop_int + 2000000
	
	# Subset genetic data within the 2-megabase interval
	GT_chr_inf_pos_start <- GT_chr_inf[GT_chr_inf[, "Position"] > Pos2mb_start_int, ]
	GT_chr_inf_pos_start <- GT_chr_inf_pos_start[GT_chr_inf_pos_start[, "Position"] < Pos_start_int, ]
	GT_chr_inf_pos_stop <- GT_chr_inf[GT_chr_inf[, "Position"] < Pos2mb_stop_int, ]
	GT_chr_inf_pos_stop <- GT_chr_inf_pos_stop[GT_chr_inf_pos_stop[, "Position"]> Pos_stop_int,]
	
	# Merge dbSNP data with genetic data within the interval
	GT_chr_inf_pos_start_merge <- merge(GT_chr_inf_pos_start,dbsnp_chr, "Position")
	GT_chr_inf_pos_stop_merge <- merge(GT_chr_inf_pos_stop,dbsnp_chr, "Position")
	
	# Update the informative SNPs data frame with counts within the interval
	inf_snps_pretest_chr_int_down <- data.frame(SNPS = nrow(GT_chr_inf_pos_start_merge), stringsAsFactors = FALSE)
	colnames(inf_snps_pretest_chr_int_down) <- paste("chr", Chrom, "2mb", "down", sep = "_")
	inf_snps_pretest_ref <- data.frame(inf_snps_pretest_ref, inf_snps_pretest_chr_int_down)
	inf_snps_pretest_chr_int_up <- data.frame(SNPS = nrow(GT_chr_inf_pos_stop_merge), stringsAsFactors = FALSE)
	colnames(inf_snps_pretest_chr_int_up) <- paste("chr", Chrom, "2mb", "up", sep = "_")
	inf_snps_pretest_ref <- data.frame(inf_snps_pretest_ref,inf_snps_pretest_chr_int_up)
	
	# Generate information about SNP loci within the interval
	loci_int <- GT_chr_inf_pos_start_merge[1:2, ]
	rownames(loci_int) <- c("loci_inf_region_start", "loci_inf_region_stop")
	loci_int[1, c("Position")] <- Pos_start_int
	loci_int[2, c("Position")] <- Pos_stop_int
	loci_int[, !names(loci_int) %in% c("Position", "inheritance")] <- "ND"
	loci_int[1, "inheritance"] <- paste(length(which(GT_chr_inf_pos_start_merge[, "inheritance"] == "WT")), "WT_SNPS_found_", length(which(GT_chr_inf_pos_start_merge[, "inheritance"] == "Af")), "Af_SNPS_found_of", nrow(GT_chr_inf_pos_start_merge), "total_up", sep = "_")
	loci_int[2, "inheritance"] <- paste(length(which(GT_chr_inf_pos_stop_merge[, "inheritance"] == "WT")), "WT_SNPS_found_", length(which(GT_chr_inf_pos_stop_merge[, "inheritance"] == "Af")), "Af_SNPS_found_of", nrow(GT_chr_inf_pos_stop_merge), "total_down", sep = "_")
	loci_int <- rbind(GT_chr_inf_pos_start_merge,loci_int,GT_chr_inf_pos_stop_merge)
	
	# Generate and save INFSNP data for Helix
	NGS_PGD_INFSNP(GT_chr_inf_pos_start_merge, GT_chr_inf_pos_stop_merge, parent, Mother, HelixoutPath)

  return(loci_int) # Return the modified genetic data frame within the interval
}