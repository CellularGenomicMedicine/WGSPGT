SegLogRs <- function(script,logRs,gammaMC,gammaSC,plateau,outPath,fam_members_Embryo){

print("PCF segmentation is applying...")

	SegGenomes <- NULL
	segLogRs <- NULL
	source(paste(script, "analyses", "functions", "getMad.R",sep="/"))
	source(paste(script, "analyses",  "functions","selectFastPcf.R",sep="/"))
	for(i in 4:ncol(logRs)){
		if(any(colnames(logRs)[i])%in%fam_members_Embryo[,"SampleID"]) {Gamma=gammaSC}else{Gamma=gammaMC}
        SegGenome <- NULL
		for(chr in unique(logRs[,"Chr"])){
			logRChr <- logRs[as.character(logRs$Chr)==chr,i]
			while(sum(is.na(logRChr))>=1){logRChr[which(is.na(logRChr))]<-logRChr[which(is.na(logRChr))-1]}
			sdev <- getMad(logRChr,k=plateau)
			if ( is.na(sdev) | sdev == 0 | length(logRChr)<100 ) {
				SegChr <- logRChr; write.table(logRChr,paste(paste(outPath,"",sep="/"),"Chr",chr,"_for_",colnames(logRs)[i],"_has_less_then_100_logRseg_values",".txt",sep=""),sep="\t",quote=F,col.names=T,row.names=T)
			} else {
                 	res <- selectFastPcf(logRChr,3,Gamma*sdev,T)
                 	SegChr <- res$yhat
			}
			SegGenome <- c(SegGenome,SegChr)
		
		}#end chr loop
		SegGenomes<-cbind(SegGenomes,SegGenome)
		print(paste(colnames(logRs)[i],"==> gamma",Gamma, "is applied"))
	}#end file loop

	segLogRs <- data.frame(logRs[,c("Names","Chr","Position")],SegGenomes,stringsAsFactors=FALSE,check.names=F)
	colnames(segLogRs)<- colnames(logRs)

	return(segLogRs)
}#end function

