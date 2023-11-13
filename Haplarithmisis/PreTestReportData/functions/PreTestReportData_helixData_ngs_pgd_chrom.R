helixData_ngs_pgd_chrom <- function(inf_snps_pretest,Chrom){
    #CHROM SNP Helix
    SNP_AANT                <- data.frame(ID="Tot. aant. SNPs",Value=inf_snps_pretest[,paste("chr",Chrom,"total",sep="_")],stringsAsFactors=FALSE)
    INF_AANT                <- data.frame(ID="Aant. inf. SNPs",Value=inf_snps_pretest[,paste("chr",Chrom,"inf",sep="_")],stringsAsFactors=FALSE)
    CHROM                   <- data.frame(ID="Chromosoom",Value=Chrom,stringsAsFactors=FALSE)
    NGS_PGD_CHROM_df    <- rbind(CHROM,SNP_AANT,INF_AANT)
    NGS_PGD_CHROM_df
}
