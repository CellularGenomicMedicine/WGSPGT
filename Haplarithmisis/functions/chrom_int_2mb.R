chrom_int_2mb <- function(script,dbsnp_path,inf_snps_pretest_ref,Chrom,Int,GT_chr_inf,parent,Mother,HelixoutPath){
	dbsnp_chr           <- fread(paste(dbsnp_path,paste("chr",Chrom,sep=""),".bed",sep=""),data.table=FALSE,stringsAsFactors=FALSE,colClasses=c("character","numeric","numeric","character"),sep="\t")
    names(dbsnp_chr)    <- c("chr","start","Position","RS_nr")
    Pos_start_int       <- Int[Int[,"Chr"]%in%Chrom,"Start"]
	Pos_stop_int        <- Int[Int[,"Chr"]%in%Chrom,"Stop"]
	if ( (Pos_start_int-2000000)<0) { Pos2mb_start_int     <- 0} else {Pos2mb_start_int <- Pos_start_int-2000000}
	Pos2mb_stop_int         <- Pos_stop_int+2000000
	GT_chr_inf_pos_start    <- GT_chr_inf[GT_chr_inf[,"Position"]>Pos2mb_start_int,]
	GT_chr_inf_pos_start    <- GT_chr_inf_pos_start[GT_chr_inf_pos_start[,"Position"]< Pos_start_int,]
	GT_chr_inf_pos_stop     <- GT_chr_inf[GT_chr_inf[,"Position"]< Pos2mb_stop_int,]
	GT_chr_inf_pos_stop     <- GT_chr_inf_pos_stop[GT_chr_inf_pos_stop[,"Position"]> Pos_stop_int,]
	#merge dbsnp
	GT_chr_inf_pos_start_merge 		<- merge(GT_chr_inf_pos_start,dbsnp_chr,"Position")
	GT_chr_inf_pos_stop_merge  		<- merge(GT_chr_inf_pos_stop,dbsnp_chr,"Position")
	inf_snps_pretest_chr_int_down  	<- data.frame(SNPS=nrow(GT_chr_inf_pos_start_merge),stringsAsFactors=FALSE)
	colnames(inf_snps_pretest_chr_int_down) <- paste("chr",Chrom,"2mb","down",sep="_")
	inf_snps_pretest_ref    <- data.frame(inf_snps_pretest_ref,inf_snps_pretest_chr_int_down)
	inf_snps_pretest_chr_int_up  <- data.frame(SNPS=nrow(GT_chr_inf_pos_stop_merge),stringsAsFactors=FALSE)
	colnames(inf_snps_pretest_chr_int_up) <- paste("chr",Chrom,"2mb","up",sep="_")
	inf_snps_pretest_ref    <- data.frame(inf_snps_pretest_ref,inf_snps_pretest_chr_int_up)
	loci_int     <- GT_chr_inf_pos_start_merge[1:2,]
	rownames(loci_int) <- c("loci_inf_region_start","loci_inf_region_stop")
	loci_int[1,c("Position")] <- Pos_start_int
	loci_int[2,c("Position")] <- Pos_stop_int
	loci_int[,!names(loci_int)%in%c("Position","inheretance")] <- "ND"
	loci_int[1,"inheretance"] <- paste(length(which(GT_chr_inf_pos_start_merge[,"inheretance"]=="WT")),"WT_SNPS_found_",length(which(GT_chr_inf_pos_start_merge[,"inheretance"]=="Af")),"Af_SNPS_found_of",nrow(GT_chr_inf_pos_start_merge),"total_up",sep="_")
	loci_int[2,"inheretance"] <- paste(length(which(GT_chr_inf_pos_stop_merge[,"inheretance"]=="WT")),"WT_SNPS_found_",length(which(GT_chr_inf_pos_stop_merge[,"inheretance"]=="Af")),"Af_SNPS_found_of",nrow(GT_chr_inf_pos_stop_merge),"total_down",sep="_")
	loci_int        <- rbind(GT_chr_inf_pos_start_merge,loci_int,GT_chr_inf_pos_stop_merge)
	#generate INFSNP stoftest Helix
	source(paste(script, "analyses", "functions","NGS_PGD_INFSNP.R", sep="/"))
	NGS_PGD_INFSNP(GT_chr_inf_pos_start_merge,GT_chr_inf_pos_stop_merge,parent,Mother,HelixoutPath)

    return(loci_int)
}