PreTestReportData_Grandparent_inf_snps <- function(script,GT,Parent,Father,Mother,REF,Ref,Chroms,Int,fam_members,dbsnp_path,parent,HelixoutPath,Chrom){
	pretest_list <- list()
	inf_snps_pretest_ref <- data.frame(Ref=Ref,stringsAsFactors=FALSE)
	for (Chr in Chroms){
		GT_chr <- GT[GT$Chr%in%Chr,]
		#total snps
		GT_chr 						<- GT_chr[GT_chr[,Father]!="NC" & GT_chr[,Mother]!="NC",]
		inf_snps_pretest_ref_total	<- data.frame(SNPS=nrow(GT_chr),stringsAsFactors=FALSE)
		colnames(inf_snps_pretest_ref_total) <- paste("chr",Chr,"total",sep="_")
		inf_snps_pretest_ref    <- data.frame(inf_snps_pretest_ref,inf_snps_pretest_ref_total)

		# semi informative snps
		GT_chr_semi_inf <- GT_chr[(GT_chr[,Father]=="AB" & GT_chr[,Mother]=="AA") | (GT_chr[,Father]=="AB" & GT_chr[,Mother]=="BB") | (GT_chr[,Father]=="AA" & GT_chr[,Mother]=="AB") | (GT_chr[,Father]=="BB" & GT_chr[,Mother]=="AB") & (GT_chr[,Ref]!="NC"),]
		inf_snps_pretest_chr_semi_inf  <- data.frame(SNPS=nrow(GT_chr_semi_inf),stringsAsFactors=FALSE)
		colnames(inf_snps_pretest_chr_semi_inf) <- paste("chr",Chr,"semi_inf",sep="_")
		inf_snps_pretest_ref    <- data.frame(inf_snps_pretest_ref,inf_snps_pretest_chr_semi_inf)

		# informative snps
		if (Parent==Father) {Parent2=Mother}
		if (Parent==Mother) {Parent2=Father}
		if (any(grep("G",fam_members[,3]))){
			GT_chr_inf <- GT_chr_semi_inf[(GT_chr_semi_inf[,Parent]== "AB" ) & (GT_chr_semi_inf[,Ref]!= "AB") & (GT_chr_semi_inf[,Parent2]!= "AB") & (GT_chr_semi_inf[,Ref]!="NC"),]
		}
		GT_chr_inf$inheretance <- "WT"
		Indication_Ref    <- fam_members[fam_members[,"SampleID"]%in%Ref,"Sample_Status"]
		source(paste(script, "analyses","PreTestReportData","functions","PreTestReportData_snp_origin_call_grandparent.R", sep="/"))
		if (grepl("AG",Indication_Ref)){ GT_chr_inf <- snp_origin_call_Affected_grandparent(GT_chr_inf,Parent,Parent2,Ref) }
		if (grepl("UG",Indication_Ref)){ GT_chr_inf <- snp_origin_call_Unaffected_grandparent(GT_chr_inf,Parent,Parent2,Ref) } 
		inf_snps_pretest_chr_inf  			<- data.frame(SNPS=nrow(GT_chr_inf),stringsAsFactors=FALSE)
		colnames(inf_snps_pretest_chr_inf) 	<- paste("chr",Chr,"inf",sep="_")
		inf_snps_pretest_ref    			<- data.frame(inf_snps_pretest_ref,inf_snps_pretest_chr_inf)
		if( any(Chr==Chrom)){
		    Chr_id <- paste("Gtype_inf_snps_chr",Chr,"of",REF,sep="_")
			source(paste(script, "analyses", "functions","chrom_int_2mb.R", sep="/"))
			pretest_list[[ Chr_id]] <- chrom_int_2mb(script,dbsnp_path,inf_snps_pretest_ref,Chrom,Int,GT_chr_inf,parent,Mother,HelixoutPath)
		}
 	}
	inf_snps_pretest_ref$total_genome <- nrow(GT)
	total_semi_inf_ref                <- rowSums(inf_snps_pretest_ref[,grepl("semi_inf",names(inf_snps_pretest_ref))])
	inf_snps_pretest_ref              <- data.frame(inf_snps_pretest_ref,total_semi_inf=total_semi_inf_ref,ratio_semi_inf=total_semi_inf_ref/nrow(GT),stringsAsFactors=FALSE)
	total_inf_ref                     <- inf_snps_pretest_ref[,grepl("inf",colnames(inf_snps_pretest_ref))]
	total_inf_ref                     <- rowSums(total_inf_ref[,!grepl("semi",names(total_inf_ref))])
	inf_snps_pretest_ref              <- data.frame(inf_snps_pretest_ref,total_inf=total_inf_ref, ratio_inf=total_inf_ref/nrow(GT),stringsAsFactors=FALSE)
	ref_id <- paste("Inf_snps_pretest",REF,sep="_")
	pretest_list[[ref_id]] <- inf_snps_pretest_ref
	return(pretest_list)
}