EmbryoTestReportData_inf_snps_Rs_gtype_inheretance_Chrom_int_2mb <- function(dbsnp_path,outPathPlots,Gtypes,PhBAF,dataHapRaw,segs,Int,Chrom,flip,parent,HelixoutPath,Gamma_value,EmbryoID){
  inf_snps_seg_df <- data.frame(segs=character(),chr=character(),up=numeric(),down=numeric(),stringsAsFactors=FALSE)
  seg_df_chr_pos <- data.frame(segs=character(),chr=numeric(),down=numeric(),up=numeric(),origin=character(),stringsAsFactors=FALSE)
  dbsnp_chr <- fread(paste(dbsnp_path,paste("chr",Chrom,sep=""),".bed",sep=""),data.table=FALSE,stringsAsFactors=FALSE,colClasses=c("character","numeric","numeric","character"),sep="\t")
  names(dbsnp_chr) <- c("chr","start","end","RS_nr")
  for (seg in segs){
    seg_df <- PhBAF[[seg]]
    dataHap_df_chr <- dataHapRaw[dataHapRaw[,"Chr"]%in%Chrom,]
    if(any(grepl("M1",segs))){
      parent_dataHap_remove=names(dataHap_df_chr)[grepl("Pat",names(dataHap_df_chr))];
      parent_dataHap=names(dataHap_df_chr)[grepl("Mat",names(dataHap_df_chr))]
    }
    if(any(grepl("P1",segs))){
      parent_dataHap_remove=names(dataHap_df_chr)[grepl("Mat",names(dataHap_df_chr))];
      parent_dataHap=names(dataHap_df_chr)[grepl("Pat",names(dataHap_df_chr))]
    }
    dataHap_df_chr <- dataHap_df_chr[dataHap_df_chr[,parent_dataHap]!=0,]
    seg_df_chr           <- seg_df[seg_df[,"Chr"]%in%Chrom,]
    seg_df_chr$inheretance <- "WT"
    Pos_start_int        <- Int[Int[,1]%in%Chrom,2]
    Pos_stop_int         <- Int[Int[,1]%in%Chrom,3]
    if ( (Pos_start_int-2000000)<0) {
      Pos2mb_start_int     <- 0
    } else {
      Pos2mb_start_int <- Pos_start_int-2000000
    }
    Pos2mb_stop_int      <- Pos_stop_int+2000000
    seg_df_chr_pos_start <- seg_df_chr[seg_df_chr[,"Position"]>Pos2mb_start_int,]
    seg_df_chr_pos_start <- seg_df_chr_pos_start[seg_df_chr_pos_start[,"Position"]< Pos_start_int,]
    seg_df_chr_pos_stop  <- seg_df_chr[seg_df_chr[,"Position"]< Pos2mb_stop_int,]
    seg_df_chr_pos_stop  <- seg_df_chr_pos_stop[seg_df_chr_pos_stop[,"Position"]> Pos_stop_int,]
    inf_snps_seg_df_chr  <- data.frame(segs=seg,chr=Chrom,down=nrow(seg_df_chr_pos_start),up=nrow(seg_df_chr_pos_stop),stringsAsFactors=FALSE)
    inf_snps_seg_df      <- rbind(inf_snps_seg_df,inf_snps_seg_df_chr)
    seg_df_chr_pos_start_merge <- merge(seg_df_chr_pos_start,dbsnp_chr,3)
    seg_df_chr_pos_stop_merge  <- merge(seg_df_chr_pos_stop,dbsnp_chr,3)
    if (nrow(seg_df_chr_pos_start_merge)>0){
      seg_df_chr_pos_start_merge$origin <- seg
      seg_df_chr_pos       <- rbind(seg_df_chr_pos,seg_df_chr_pos_start_merge)
    }
    if (nrow(seg_df_chr_pos_stop_merge)>0){
     seg_df_chr_pos_stop_merge$origin <- seg
     seg_df_chr_pos       <- rbind(seg_df_chr_pos,seg_df_chr_pos_stop_merge)
    }
  }
  if (sum(inf_snps_seg_df[,"down"]>0 | inf_snps_seg_df[,"up"]>0)>0) {
    names(seg_df_chr_pos) <- c("Position","Names","Chr","BAF","inheretance","chr","start","RS_nr","origin")
    seg_df_chr_pos        <- seg_df_chr_pos[,c("Names","BAF","start","RS_nr","origin","inheretance")]
    write.table(inf_snps_seg_df,paste(outPathPlots,paste("Informative_snps_2MB_chr",Chrom,"phasing",parent,"for",EmbryoID,sep="_"),".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    write.table(seg_df_chr_pos,paste(outPathPlots,paste("Informative_snps_2MB_with_RS_chr",Chrom,"phasing",parent,"for",EmbryoID,sep="_"),".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    seg_df_chr_pos_gtype <- merge(seg_df_chr_pos,Gtypes,by="Names")
    seg_df_chr_pos_gtype <- seg_df_chr_pos_gtype[!duplicated(seg_df_chr_pos_gtype[,"Names"]),]
    seg_df_chr_pos_gtype <- merge(seg_df_chr_pos_gtype,dataHap_df_chr[,c("Names",parent_dataHap)],by="Names")
    cond1=1
    if (flip==1) {cond1=2}
    if (any(grepl("M1",segs))) {
      seg_df_chr_pos_gtype[seg_df_chr_pos_gtype[,grepl("Mat",names(seg_df_chr_pos_gtype))]==cond1,"inheretance"] <- "Af"
    }
    if (any(grepl("P1",segs))) {
      seg_df_chr_pos_gtype[seg_df_chr_pos_gtype[,grepl("Pat",names(seg_df_chr_pos_gtype))]==cond1,"inheretance"] <- "Af"
    }
    seg_df_chr_pos_gtype <- seg_df_chr_pos_gtype[!duplicated(seg_df_chr_pos_gtype[,"Names"]),]
    loci_int     <- seg_df_chr_pos_gtype[1:2,]
    loci_int[1,c("Names","BAF","start","RS_nr","origin")] <-  data.frame(Names="ROI",BAF=0.5,start=Pos_start_int,RS_nr="loci_int",origin=parent,stringsAsFactors=FALSE)
    loci_int[2,c("Names","BAF","start","RS_nr","origin")] <-  data.frame(Names="ROI",BAF=0.5,start=Pos_stop_int,RS_nr="loci_int",origin=parent,stringsAsFactors=FALSE)
    loci_int[,!names(loci_int)%in%c("Name","BAF","start","RS_nr","origin","Chr","Position")] <- "ND"
    seg_df_chr_pos_gtype_start <- seg_df_chr_pos_gtype[seg_df_chr_pos_gtype[,"Position"]< Pos_start_int,]
    seg_df_chr_pos_gtype_stop  <- seg_df_chr_pos_gtype[seg_df_chr_pos_gtype[,"Position"]> Pos_stop_int,]
    loci_int[1,"inheretance"] <- paste(length(which(seg_df_chr_pos_gtype_start[,"inheretance"]=="WT")),"WT_SNPS_found_",length(which(seg_df_chr_pos_gtype_start[,"inheretance"]=="Af")),"Af_SNPS_found(",round((length(which(seg_df_chr_pos_gtype_start[,"inheretance"]=="Af"))/nrow(seg_df_chr_pos_gtype_start))*100,1),"pct)_up",sep="_")
    loci_int[2,"inheretance"] <- paste(length(which(seg_df_chr_pos_gtype_stop[,"inheretance"]=="WT")),"WT_SNPS_found_",length(which(seg_df_chr_pos_gtype_stop[,"inheretance"]=="Af")),"Af_SNPS_found(",round((length(which(seg_df_chr_pos_gtype_stop[,"inheretance"]=="Af"))/nrow(seg_df_chr_pos_gtype_stop))*100,1),"pct)_up",sep="_")
    seg_df_chr_pos_gtype        <- rbind(seg_df_chr_pos_gtype,loci_int)
    seg_df_chr_pos_gtype        <- seg_df_chr_pos_gtype[order(seg_df_chr_pos_gtype[,"start"]),]
    write.table(seg_df_chr_pos_gtype,paste(outPathPlots,paste("Informative_snps_2MB_with_RS_gtype_chr",Chrom,"phasing",parent,"for",EmbryoID,sep="_"),".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    #INFSNP Helix
    AFm2MB                  <- data.frame(ID="Af SNP -2Mb",Value=length(which(seg_df_chr_pos_gtype_start[,"inheretance"]=="Af")),stringsAsFactors=FALSE)
    WTm2MB                  <- data.frame(ID="WT SNP -2Mb",Value=length(which(seg_df_chr_pos_gtype_start[,"inheretance"]=="WT")),stringsAsFactors=FALSE)
    TOTm2MB                 <- data.frame(ID="Tot SNP -2Mb",Value=nrow(seg_df_chr_pos_gtype_start),stringsAsFactors=FALSE)
    AFp2MB                  <- data.frame(ID="Af SNP +2Mb",Value=length(which(seg_df_chr_pos_gtype_stop[,"inheretance"]=="Af")),stringsAsFactors=FALSE)
    WTp2MB                  <- data.frame(ID="WT SNP +2Mb",Value=length(which(seg_df_chr_pos_gtype_stop[,"inheretance"]=="WT")),stringsAsFactors=FALSE)
    TOTp2MB                 <- data.frame(ID="Tot SNP +2Mb",Value=nrow(seg_df_chr_pos_gtype_stop),stringsAsFactors=FALSE)
    CCAFm2MB                <- data.frame(ID="Conc Aff -2Mb",Value=round(length(which(seg_df_chr_pos_gtype_start[,"inheretance"]=="Af"))/nrow(seg_df_chr_pos_gtype_start)*100,2),stringsAsFactors=FALSE)
    CCAFp2MB                <- data.frame(ID="Conc Aff +2Mb",Value=round(length(which(seg_df_chr_pos_gtype_stop[,"inheretance"]=="Af"))/nrow(seg_df_chr_pos_gtype_stop)*100,2),stringsAsFactors=FALSE)
    CCWTm2MB                <- data.frame(ID="Conc WT -2Mb",Value=round(length(which(seg_df_chr_pos_gtype_start[,"inheretance"]=="WT"))/nrow(seg_df_chr_pos_gtype_start)*100,2),stringsAsFactors=FALSE)
    CCWTp2MB                <- data.frame(ID="Conc WT +2Mb",Value=round(length(which(seg_df_chr_pos_gtype_stop[,"inheretance"]=="WT"))/nrow(seg_df_chr_pos_gtype_stop)*100,2),stringsAsFactors=FALSE)
    NGS_PGD_INFSNP_df       <- rbind(AFm2MB,WTm2MB,TOTm2MB ,AFp2MB ,WTp2MB,TOTp2MB,CCAFm2MB,CCWTm2MB,CCAFp2MB,CCWTp2MB)
    if ( Gamma_value == 50 ){
      if ( parent=="Father"){
        write.table(NGS_PGD_INFSNP_df,paste(HelixoutPath,paste0(EmbryoID,"-NGS_PGD_INFSNP_PAT_EMBRYO"),sep="/"),col.names=T,row.names=F,quote=F,sep=",")
      }
      if ( parent=="Mother"){
        write.table(NGS_PGD_INFSNP_df,paste(HelixoutPath,paste0(EmbryoID,"-NGS_PGD_INFSNP_MAT_EMBRYO"),sep="/"),col.names=T,row.names=F,quote=F,sep=",")
      }
    }
  }else {
    write.table(inf_snps_seg_df,paste(outPathPlots,paste("NO informative snps 2MB chr",Chrom,"phasing",parent,"for",EmbryoID,sep=" "),".txt",sep=""),col.names=T,row.names=F,quote=F,sep="\t")
    stop(paste("NO informative snps up and or down 2MB chr",Chrom,"phasing",parent,"for",EmbryoID,sep=" "))
  }
}

