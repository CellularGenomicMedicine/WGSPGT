NGS_PGD_INFSNP <- function(GT_chr_inf_pos_start_merge,GT_chr_inf_pos_stop_merge,parent,Mother,HelixoutPath){
	#INFSNP Helix
    AFm2MB                  <- data.frame(ID="Af SNP -2Mb",Value=length(which(GT_chr_inf_pos_start_merge[,"inheretance"]=="Af")),stringsAsFactors=FALSE)
    WTm2MB                  <- data.frame(ID="WT SNP -2Mb",Value=length(which(GT_chr_inf_pos_start_merge[,"inheretance"]=="WT")),stringsAsFactors=FALSE)
    TOTm2MB                 <- data.frame(ID="Tot SNP -2Mb",Value=nrow(GT_chr_inf_pos_start_merge),stringsAsFactors=FALSE)
    AFp2MB                  <- data.frame(ID="Af SNP +2Mb",Value=length(which(GT_chr_inf_pos_stop_merge[,"inheretance"]=="Af")),stringsAsFactors=FALSE)
    WTp2MB                  <- data.frame(ID="WT SNP +2Mb",Value=length(which(GT_chr_inf_pos_stop_merge[,"inheretance"]=="WT")),stringsAsFactors=FALSE)
    TOTp2MB                 <- data.frame(ID="Tot SNP +2Mb",Value=nrow(GT_chr_inf_pos_stop_merge),stringsAsFactors=FALSE)
    NGS_PGD_INFSNP_df       <- rbind(AFm2MB,WTm2MB,TOTm2MB ,AFp2MB ,WTp2MB,TOTp2MB)
    if (parent=="Mother"){
		write.table(NGS_PGD_INFSNP_df,paste(paste(HelixoutPath,"",sep="/"),Mother,"-NGS_PGD_INFSNP_MAT",sep=""),col.names=T,row.names=F,quote=F,sep=",")
	}
    if (parent=="Father"){
		write.table(NGS_PGD_INFSNP_df,paste(paste(HelixoutPath,"",sep="/"),Mother,"-NGS_PGD_INFSNP_PAT",sep=""),col.names=T,row.names=F,quote=F,sep=",")
	}
}