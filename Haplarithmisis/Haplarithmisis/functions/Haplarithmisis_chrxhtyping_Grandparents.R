Haplarithmisis_chrxhtyping_Grandparents <- function(Father,Mother,Grandfather,Grandmother,Gtypes,parent,ScSexes,EmbryoID){
	GtypesChrX <- Gtypes[Gtypes$Chr=="X",]
	MotherChrX <- GtypesChrX[,Mother]
	FatherChrX <- GtypesChrX[,Father]
	GrandfatherChrX <- GtypesChrX[,Grandfather]
	GrandmotherChrX <- GtypesChrX[,Grandmother]

	HapXMat <- matrix(0,nrow(GtypesChrX),ncol(as.matrix(GtypesChrX[,which(colnames(Gtypes)%in%EmbryoID)])))
	colnames(HapXMat)<-EmbryoID

	if(parent=="Mother"){

		MotherChrX[(GtypesChrX[,Grandfather]=="BB" & GtypesChrX[,Mother]=="AB")]<-"BA"
	
		MotherChrX[(GtypesChrX[,Grandfather]=="NC" & GtypesChrX[,Mother]=="AB") |
					   (GtypesChrX[,Grandfather]=="BB" & GtypesChrX[,Grandmother]=="BB" & GtypesChrX[,Mother]=="AB") |
					   (GtypesChrX[,Grandfather]=="AA" & GtypesChrX[,Grandmother]=="AA" & GtypesChrX[,Mother]=="AB") ]<-"NC"
	}	   

	na.omit(ScSexes[ScSexes[,2]=="male",1])

	if(ScSexes[,2]=="male"){
		print(paste(EmbryoID,"is a male sample but has",sum(GtypesChrX[,EmbryoID]=="AB"),"heterozygous SNP-calls; these are treated as NoCalls"))
		GtypesChrX[GtypesChrX[,EmbryoID]=="AB",EmbryoID]<-"NC"
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,EmbryoID] == "BB") |
				   (MotherChrX=="BA" & GtypesChrX[,EmbryoID] == "AA"),EmbryoID] <- 2
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,EmbryoID] == "AA") |
					(MotherChrX=="BA" & GtypesChrX[,EmbryoID] == "BB"),EmbryoID] <- 1
		HapXMat[GrandfatherChrX=="NC" | GrandmotherChrX=="NC" ,EmbryoID] <- 0
		print(EmbryoID)
	}

	if(ScSexes[,2]=="female"){
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,EmbryoID] == "BB") |
				   (MotherChrX=="BA" & GtypesChrX[,EmbryoID] == "AA") |
				   (MotherChrX=="AB" & FatherChrX == "AA" & GtypesChrX[,EmbryoID] == "AB") |
				   (MotherChrX=="BA" & FatherChrX == "BB" & GtypesChrX[,EmbryoID] == "AB") ,EmbryoID] <- 2
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,EmbryoID] == "AA") |
				    (MotherChrX=="BA" & GtypesChrX[,EmbryoID] == "BB") |
					(MotherChrX=="BA" & FatherChrX == "AA" & GtypesChrX[,EmbryoID] == "AB")|
					(MotherChrX=="AB" & FatherChrX == "BB" & GtypesChrX[,EmbryoID] == "AB"),EmbryoID] <- 1
		HapXMat[GrandfatherChrX=="NC" | GrandmotherChrX=="NC" | FatherChrX=="NC" ,EmbryoID] <- 0
		print(EmbryoID)
 	}
	return(HapXMat)
}#end chrxhtypingOpt1 function