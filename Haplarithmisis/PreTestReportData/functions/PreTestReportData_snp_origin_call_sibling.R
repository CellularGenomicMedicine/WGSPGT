# snp_origin_call_Affected_sibling - Function to assign inheritance status for informative SNPs in affected sibling as referent

snp_origin_call_Affected_sibling <- function(GT_chr_inf, Parent, Parent2, Ref){
  # Assign affected ("Af") inheritance for specific parental genotypes and SNP alleles
  GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "BB" & GT_chr_inf[, Ref] == "AA"), "inheritance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "AA" & GT_chr_inf[, Ref] == "BB"), "inheritance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "AA" & GT_chr_inf[, Ref] == "AB"), "inheritance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "BB" & GT_chr_inf[, Ref] == "AB"), "inheritance"] <- "Af"
	return(GT_chr_inf) # Return the modified DataFrame
}


# snp_origin_call_Unaffected_sibling - Function to assign inheritance status for informative SNPs in unaffected sibling as referent

snp_origin_call_Unaffected_sibling <- function(GT_chr_inf, Parent, Parent2, Ref){
  # Assign affected ("Af") inheritance for specific parental genotypes and SNP alleles
  GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "BB" & GT_chr_inf[, Ref] == "BB"), "inheritance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[, Parent] == "AB" & GT_chr_inf[, Parent2] == "AA" & GT_chr_inf[, Ref] == "AA"), "inheritance"] <- "Af"
	return(GT_chr_inf) # Return the modified DataFrame
}
