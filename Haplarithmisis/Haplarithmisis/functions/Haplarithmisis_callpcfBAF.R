Haplarithmisis_callpcfBAF <- function(script,PhBAFseg,Gamma_value,plateau,EmbryoID){
  source(paste(script, "analyses", "functions", "getMad.R",sep="/"))
  source(paste(script, "analyses",  "functions","selectFastPcf.R",sep="/"))
  PhBAFseg <- PhBAFseg[PhBAFseg[,"Chr"]!="Y",]
  print("PCF segmentation is applying...")
  SegGenomes <- NULL
  SegGenome <- NULL
  sdevGenome <- getMad(na.omit(PhBAFseg[,colnames(PhBAFseg)%in%EmbryoID]),k=plateau)
  if (sdevGenome==0) {
    sdevCat <- NULL;
    for(chr in unique(PhBAFseg[,"Chr"])){
      logRChr <- PhBAFseg[as.character(PhBAFseg$Chr)==chr,EmbryoID];
      while(sum(is.na(logRChr))>=1){logRChr[which(is.na(logRChr))]<-logRChr[which(is.na(logRChr))-1]}
      sdev <- getMad(logRChr,k=plateau);
      sdevCat <- c(sdevCat,sdev)
    }
    sdevGenome <- mean(sdevCat)
  }
  print(paste("sdevGenome",sdevGenome,sep="_"))
  print(paste("gamma",Gamma_value,sep="_"))
  for(chr in unique(PhBAFseg[,"Chr"])){
    logRChr <- PhBAFseg[as.character(PhBAFseg$Chr)==chr,EmbryoID]
    while(sum(is.na(logRChr))>=1){logRChr[which(is.na(logRChr))]<-logRChr[which(is.na(logRChr))-1]}
    sdev <- getMad(logRChr,k=plateau)
    if ( is.na(sdev) | sdev == 0  ) {
      print(paste(chr,"sdev",sdev,sep="_")); sdev <- sdevGenome
    }
    if ( length(logRChr)<100 ){
      SegChr <- logRChr; print(paste(chr,"logRChr",length(logRChr),sep="_"))
    } else { res <- selectFastPcf(logRChr,3,Gamma_value*sdev,T)
      SegChr <- res$yhat
    }
    SegGenome <- c(SegGenome,SegChr)
  }#end chr loop

  SegGenomes<-cbind(SegGenomes,SegGenome)
  print(paste(EmbryoID,"==> gamma",Gamma_value, "is applied"))

  PhBAFSeg <- data.frame(PhBAFseg[,c("Names","Chr","Position")],SegGenomes,stringsAsFactors=FALSE)
  colnames(PhBAFSeg)<- colnames(PhBAFseg)

  return(PhBAFSeg)
}#end function
