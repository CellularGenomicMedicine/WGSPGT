snp_origin_call_grandparents <- function(GT_chr_inf,Parent,Parent2,Ref_inf,Ref_Uninf){
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="AA" & GT_chr_inf[,Ref_inf]=="BB" & GT_chr_inf[,Ref_Uninf]=="AA" ),"inheretance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="AA" & GT_chr_inf[,Ref_inf]=="AB" & GT_chr_inf[,Ref_Uninf]=="AA" ),"inheretance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="BB" & GT_chr_inf[,Ref_inf]=="AA" & GT_chr_inf[,Ref_Uninf]=="BB" ),"inheretance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="BB" & GT_chr_inf[,Ref_inf]=="AB" & GT_chr_inf[,Ref_Uninf]=="BB" ),"inheretance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="BB" & GT_chr_inf[,Ref_inf]=="AA" & GT_chr_inf[,Ref_Uninf]=="AB" ),"inheretance"] <- "Af"
	GT_chr_inf[(GT_chr_inf[,Parent]=="AB" & GT_chr_inf[,Parent2]=="AA" & GT_chr_inf[,Ref_inf]=="BB" & GT_chr_inf[,Ref_Uninf]=="AB" ),"inheretance"] <- "Af"
	return(GT_chr_inf)
}