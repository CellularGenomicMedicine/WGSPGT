snp_origin_call_Affected_sibling <- function(GT_chr_inf,Parent,Parent2,Ref){
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="BB" & GT_chr_inf[,Ref]=="AA"),"inheretance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="AA" & GT_chr_inf[,Ref]=="BB"),"inheretance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="AA" & GT_chr_inf[,Ref]=="AB"),"inheretance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="BB" & GT_chr_inf[,Ref]=="AB"),"inheretance"] <- "Af"
	return(GT_chr_inf)
}

snp_origin_call_Unaffected_sibling <- function(GT_chr_inf,Parent,Parent2,Ref){
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="BB" & GT_chr_inf[,Ref]=="BB"),"inheretance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="AA" & GT_chr_inf[,Ref]=="AA"),"inheretance"] <- "Af"
	return(GT_chr_inf)
}
