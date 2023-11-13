Haplarithmisis_patscore2 <- function(Father,Mother,dataPo,QC,Chroms,Gtypes,EmbryoID){
	#change to different values with WGS remove before putting in production
	AvgMDA		<-2.2862648
	SdMDA		<-0.4721314
	MinCallRate <-16.38

	print("Computing parental scores...")
	dataCallRate <- QC[["CallRateChrsInd"]]
			
	CallRates <- dataCallRate[,Chroms]

	Pos <- dataPo[,EmbryoID]
	Pos <- as.matrix(Pos)
	colnames(Pos)<-EmbryoID
	SamplesPo <- vector("list",ncol(Pos))
	names(SamplesPo)<- EmbryoID
		
	PoScore <- matrix(NA,length(Chroms),6)
	rownames(PoScore) <- Chroms
	for(chr in Chroms){
		
		PatGtype <- Gtypes[Gtypes[,"Chr"]==chr,Father]
		MatGtype <- Gtypes[Gtypes[,"Chr"]==chr,Mother]
		InfPo <- sum((MatGtype== "AA" & PatGtype == "BB") |
							    (MatGtype== "BB" & PatGtype == "AA") | 
							    (MatGtype== "AB" & PatGtype == "AA") | 
							    (MatGtype== "AA" & PatGtype == "AB") | 
							    (MatGtype== "BB" & PatGtype == "AB") | 
							    (MatGtype== "AB" & PatGtype == "BB"))

		dataPoChr <- Pos[dataPo[,"Chr"]==chr ,EmbryoID]
		Pat <- (abs(sum(dataPoChr[dataPoChr<0]))/InfPo)*100
		Mat <- (sum(dataPoChr[dataPoChr>0])/InfPo)*100
		if(Pat==0){Pat=0.0001} ; if(Mat==0){Mat=0.0001}
	
		PatMat <- Pat/(Mat+Pat)
		MatPat <- Mat/(Mat+Pat)
		
		if(Pat<=(AvgMDA-(3*SdMDA)) & Mat<=(AvgMDA-(3*SdMDA))){PatMat<-0;MatPat<-0}

		PoScore[chr,]<-cbind(Pat,Mat,PatMat,MatPat,CallRates[EmbryoID,chr],NA)
	}
	
	colnames(PoScore) <- c("Pat","Mat","PatMat","MatPat","CallRate","Par")
	SamplesPo[[EmbryoID]]<-PoScore
	print(EmbryoID)


	SamplesPo2 <- vector("list",ncol(Pos))
	names(SamplesPo2)<- colnames(Pos)

	for(ind in names(SamplesPo)){

		Sample <- SamplesPo[[ind]]
		Sample[Sample[,"CallRate"]>MinCallRate & ((Sample[,"Pat"]< (AvgMDA+SdMDA)& Sample[,"Mat"]< (AvgMDA+SdMDA)) | (Sample[,"PatMat"]<0.9 & Sample[,"MatPat"]<0.9)),"Par"]<- 12
		Sample[Sample[,"Pat"]> (AvgMDA+SdMDA) & Sample[,"CallRate"]>MinCallRate & Sample[,"PatMat"]>=0.9 ,"Par"]<- 11
		Sample[Sample[,"Mat"]> (AvgMDA+SdMDA) & Sample[,"CallRate"]>MinCallRate & Sample[,"MatPat"]>=0.9  ,"Par"]<- 22
		Sample[Sample[,"CallRate"]<MinCallRate & (Sample[,"PatMat"]<0.5 & Sample[,"MatPat"]<0.5) ,"Par"]<- 00
 		SamplesPo2[[ind]] <- Sample

	}# end ind loop

return(SamplesPo2)

}#end function