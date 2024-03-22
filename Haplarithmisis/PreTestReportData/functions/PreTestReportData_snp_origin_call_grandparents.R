# snp_origin_call_grandparents - Function to assign inheritance status for informative SNPs based on grandparental genotypes 

snp_origin_call_grandparents <- function(GT_chr_inf, Parent, Parent2, Ref_inf, Ref_Uninf){
  # Assign affected ("Af") inheritance for specific parental genotypes and SNP alleles
	GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "AA" & GT_chr_inf[, Ref_inf] == "BB" & GT_chr_inf[, Ref_Uninf] == "AA" ), "inheritance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "AA" & GT_chr_inf[, Ref_inf] == "AB" & GT_chr_inf[, Ref_Uninf] == "AA" ), "inheritance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "BB" & GT_chr_inf[, Ref_inf] == "AA" & GT_chr_inf[, Ref_Uninf] == "BB" ), "inheritance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "BB" & GT_chr_inf[, Ref_inf] == "AB" & GT_chr_inf[, Ref_Uninf] == "BB" ), "inheritance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "BB" & GT_chr_inf[, Ref_inf] == "AA" & GT_chr_inf[, Ref_Uninf] == "AB" ), "inheritance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "AA" & GT_chr_inf[, Ref_inf] == "BB" & GT_chr_inf[, Ref_Uninf] == "AB" ), "inheritance"] <- "Af"
	return(GT_chr_inf) # Return the modified DataFrame
}