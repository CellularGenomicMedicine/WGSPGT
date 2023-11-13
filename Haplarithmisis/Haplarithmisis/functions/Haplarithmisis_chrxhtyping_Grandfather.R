Haplarithmisis_chrxhtyping_Grandfather <- function(Father,Mother,RefSampleID,Gtypes,parent,ScSexes,EmbryoID){

	GtypesChrX <- Gtypes[Gtypes$Chr=="X",]
	MotherChrX <- GtypesChrX[,Mother]
	FatherChrX <- GtypesChrX[,Father]
	GrandfatherChrX <- GtypesChrX[,RefSampleID]

	HapXMat           <- matrix(0,nrow(GtypesChrX),ncol(as.matrix(GtypesChrX[,colnames(GtypesChrX)%in%EmbryoID])))
	colnames(HapXMat) <- colnames(Gtypes)[colnames(Gtypes)%in%EmbryoID]

	if(parent=="Mother"){
		MotherChrX[(GtypesChrX[,RefSampleID]=="BB" & GtypesChrX[,Mother]=="AB")]<-"BA"
	}

	if(ScSexes[,2]=="male"){
		print(paste(EmbryoID,"is a male sample but has",sum(GtypesChrX[,EmbryoID]=="AB"),"heterozygous SNP-calls; these are treated as NoCalls"))
		GtypesChrX[GtypesChrX[,EmbryoID]=="AB",EmbryoID]<-"NC"
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,EmbryoID] == "BB") |
				   (MotherChrX=="BA" & GtypesChrX[,EmbryoID] == "AA"),EmbryoID] <- 2
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,EmbryoID] == "AA") |
					(MotherChrX=="BA" & GtypesChrX[,EmbryoID] == "BB"),EmbryoID] <- 1
		HapXMat[GrandfatherChrX=="NC" ,EmbryoID] <- 0
		print(EmbryoID)
 	}#end s loop

	if(ScSexes[,2]=="female"){
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,EmbryoID] == "BB") |
				   (MotherChrX=="BA" & GtypesChrX[,EmbryoID] == "AA") |
				   (MotherChrX=="AB" & FatherChrX == "AA" & GtypesChrX[,EmbryoID] == "AB") |
				   (MotherChrX=="BA" & FatherChrX == "BB" & GtypesChrX[,EmbryoID] == "AB") ,EmbryoID] <- 2
		HapXMat[(MotherChrX=="AB" & GtypesChrX[,EmbryoID] == "AA") |
				    (MotherChrX=="BA" & GtypesChrX[,EmbryoID] == "BB") |
					(MotherChrX=="BA" & FatherChrX == "AA" & GtypesChrX[,EmbryoID] == "AB")|
					(MotherChrX=="AB" & FatherChrX == "BB" & GtypesChrX[,EmbryoID] == "AB"),EmbryoID] <- 1
		HapXMat[GrandfatherChrX=="NC" | FatherChrX=="NC" ,EmbryoID] <- 0
		print(EmbryoID)
 	}#end s loop
	
	return(HapXMat)

}#end chrxhtypingOpt1 function