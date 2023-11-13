Haplarithmisis_intphappropser3 <- function(Haps,ParScore,EmbryoID,Chroms){
	dataHap    <- Haps[["dataHap"]]
	dataHapRaw <- Haps[["dataHapRaw"]]
	MetaInfoCols <- c("Names","Chr","Position")
	MinLength = 300
	MinFrac = 0.6
	dataHap1 <- dataHap
	#Chroms <- c(1:22,"X")
	Samps 			<- vector("list",ncol(dataHap1)-length(MetaInfoCols))
	names(Samps)	<- colnames(dataHap1[,!names(dataHap1)%in%MetaInfoCols])

	ParScore2 				<- ParScore
	ParScore[[EmbryoID]] 	<- ParScore[[EmbryoID]][!is.na(ParScore[[EmbryoID]][,"Par"]),]
	if(sum(nrow(ParScore[[EmbryoID]]))>1){
		PatChr <- rownames(ParScore[[EmbryoID]])[ParScore[[EmbryoID]][,"Par"] ==11]
		MatChr <- rownames(ParScore[[EmbryoID]])[ParScore[[EmbryoID]][,"Par"] ==22]
		NullChr <- rownames(ParScore[[EmbryoID]])[ParScore[[EmbryoID]][,"Par"] ==0 ]
		if(length(PatChr)!=0){for(chr in PatChr){dataHap1[dataHap1[,"Chr"]==chr,paste(EmbryoID,"_Mat",sep="")]<- 0}}
		if(length(MatChr)!=0){for(chr in MatChr){dataHap1[dataHap1[,"Chr"]==chr,paste(EmbryoID,"_Pat",sep="")]<- 0}}
		if(length(NullChr)!=0){for(chr in NullChr){dataHap1[dataHap1[,"Chr"]==chr,paste(EmbryoID,"_Pat",sep="")]<- 0
			dataHap1[dataHap1[,"Chr"]==chr,paste(EmbryoID,"_Mat",sep="")]<- 0}}
		print(EmbryoID)
	}else{
		print(paste(EmbryoID,"sample sbould be excluded for diagnosis, as the origin of the chromosomes could not be determined"))
	}


	ParScore <- ParScore2
	print("--------------------------------------------------------------")
	print("Probability computation...")

	Blasts <- colnames(dataHap1)[grepl(EmbryoID,colnames(dataHap1))]
	for(samp in Blasts){
		Indiv <- EmbryoID
		for(chr in Chroms){
			if(!is.na(ParScore[[EmbryoID]][chr,"Par"])){
				dataHapChr <- dataHap1[dataHap1$Chr==chr,]

				IntBlock <- rle(dataHap1[dataHap1$Chr==chr,samp])
				IntBlock <- cbind(IntBlock$values,IntBlock$lengths)
				colnames(IntBlock)<- c("values","lengths")
				IntBlock <- data.frame(rbind(IntBlock,c(0,0)))

				CumLength <- rbind(cumsum(IntBlock$lengths),IntBlock$values)

				if(nrow(IntBlock)==0){
					CumLength<- cbind(0,0,0,0)
				}else{
					InfSNPs <- rep(NA,ncol(CumLength))
					InfSNPsDisc <- rep(NA,ncol(CumLength))
					Frac <- rep(NA,ncol(CumLength))
					Annot <- matrix(NA,ncol(CumLength),3)

					for(i in 1:ncol(CumLength)){
						if(i==1){start=1}else{start=CumLength[1,i-1]+1}
						if(start>CumLength[1,i]){start=CumLength[1,i]}
						InfSNPs[i] <- sum(dataHapRaw[dataHap1$Chr==chr,samp][start:CumLength[1,i]]==CumLength[2,i])
						InfSNPsDisc[i] <- sum(dataHapRaw[dataHap1$Chr==chr,samp][start:CumLength[1,i]]!=CumLength[2,i] & dataHapRaw[dataHap1$Chr==chr,samp][start:CumLength[1,i]]!=0)
						Frac[i] <- round(InfSNPs[i]/(InfSNPs[i]+InfSNPsDisc[i]),digits=3)
						Annot[i,] <- cbind(chr,dataHapChr$Position[start],dataHapChr$Position[CumLength[1,i]])
						if(Frac[i]=="NaN"){Frac[i] <- 0}#when no informative SNP is present for the smoothed haplotype block
						if( Frac[i] < MinFrac){dataHap1[dataHap1$Chr==chr,samp][start:CumLength[1,i]] <- 0}
					}#end i loop
					CumLength<- cbind.data.frame(Annot,IntBlock$lengths,t(CumLength),InfSNPs,InfSNPsDisc,Frac)
				}
				Chroms2 <- rownames(ParScore[[EmbryoID]])
				if(chr==Chroms2[which(ParScore[[EmbryoID]][,"Par"]!="NA")[1]]){CumLengths<-CumLength}else{CumLengths<-rbind(CumLengths,CumLength)}
			}else{print(paste("The origin of chr.",chr, "of", EmbryoID,"could not be determined"))}
		}#end chr loop

		if(exists("CumLengths")==FALSE){CumLengths <- matrix(0,2,9)}
		colnames(CumLengths) <- c("Chr","Start","Stop","Length","CumLength","Value","#InfSNPs","#InfSNPsDisc","Frac")
		CumLengths[CumLengths[,"Frac"]<MinFrac,"Value"]<-0
		CumLengths <- CumLengths[CumLengths[,"Value"]!=0,]
		Samps[[samp]]<-CumLengths
		print(samp)
		rm(list="CumLengths")
	}#end samp loop

	MDAs <- do.call(rbind,Samps[grepl(EmbryoID,names(Samps))])

	X <-MDAs[which(MDAs$Chr=="X"),]
	MDAs <- rbind(MDAs[MDAs$Chr!="X",],X[grep("Mat",rownames(X)),])
	MDAs <- MDAs[MDAs[,"Value"]!=0,]
	MDAs <- MDAs[MDAs[,"Length"]>=MinLength,]

	MDAs[,"Chr"] <- as.character(MDAs[,"Chr"])
	MDAs[,"Start"] <- as.numeric(as.character(MDAs[,"Start"]))
	MDAs[,"Stop"] <- as.numeric(as.character(MDAs[,"Stop"]))


	LmFit <- rlm(MDAs[,"#InfSNPs"]~MDAs[,"Length"], data=MDAs)
	DataCompl<-cbind(MDAs,LmFit$fitted.values,LmFit$resid,rep(NA,nrow(MDAs)))
	colnames(DataCompl)[(ncol(DataCompl)-2):ncol(DataCompl)]<-cbind("Fitted","Resid","Prob")
	NegRes <- subset(DataCompl, Resid < 0)


	zscores<-scale(DataCompl$Resid,center=TRUE,scale=TRUE)
	datawithzscores<-cbind(DataCompl,zscores)
	DataCompl[,"Prob"] <- round(pnorm(-abs(datawithzscores$zscores))+0.5,digits=3)
	DataCompl[DataCompl[,"Resid"]>0,"Prob"] <- 1

	Intp <- vector("list",2)
	names(Intp) <- c("dataHap","DataCompl")
	Intp[["dataHap"]] <- dataHap1
	Intp[["DataCompl"]] <- DataCompl
	return(Intp)

}

