PreTestReportData_mendinc <- function(script,config_file_fam,Parent,REF,Father,Mother,GT,outPath,Chroms){
		#calculate informative SNPS and generate stoftests
		#mendelian inconsistency, phasing and extracting informative SNPs
		if (REF=="Grandparents"){
			source(as.character(config_file_fam))
			source(paste(script, "analyses", "PreTestReportData", "functions","PreTestReportData_mendInc_Grandparents.R", sep="/"))
			tot_Percent_MendInc       <- PreTestReportData_mendInc_Grandparents(Grandmother,Grandfather,Parent,GT,Chroms)
		}

		if (REF=="SiblingF" | REF=="SiblingM"){
			source(as.character(config_file_fam))
			source(paste(script, "analyses", "PreTestReportData","functions","PreTestReportData_mendInc_Sibling.R", sep="/"))
			tot_Percent_MendInc       <- PreTestReportData_mendInc_Sibling(Father,Mother,RefSampleID,GT,Chroms)
		}

		if (REF=="Grandmother" | REF=="Grandfather"){
			source(as.character(config_file_fam))
			source(paste(script, "analyses", "PreTestReportData","functions","PreTestReportData_mendInc_Grandparent.R", sep="/"))
			tot_Percent_MendInc       <- PreTestReportData_mendInc_Grandparent(RefSampleID,Parent,GT,Chroms)
		}
		write.table(tot_Percent_MendInc,paste(paste(outPath,"",sep="/"),parent,"_MendInc_",REF,".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")

}