Haplarithmisis_chrxhtypingOpt_SiblingFemale <- function(script,Father,Mother,RefSampleID, Gtypes,ParScore,EmbryoID){
	print("Haplotype reconstrucion of the sex chromosome...")
	print("Ammended chrxhtyping")

	source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_ChrXFatherAB.R",sep="/"))
	Gtypes[Gtypes[,"Chr"]=="X",Father] <- Haplarithmisis_ChrXFatherAB(Gtypes,Father)

	GTFather <- Gtypes[,Father][Gtypes$Chr=="X"]
	GTMother <- Gtypes[,Mother][Gtypes$Chr=="X"]

	GTRef <- Gtypes[Gtypes$Chr=="X",RefSampleID]
	ScSexes <- matrix(NA,length(ParScore),2)
	rownames(ScSexes)<-names(ParScore)

	source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_sexdetermination.R",sep="/"))
	ScSexes <- Haplarithmisis_sexdetermination(ScSexes,ParScore)

	GTRef[GTRef=="AB" & GTFather=="BB"]="BA"#There could not be htz SNPs on chromosome X of male affected
	GTMother[GTRef=="AB" & GTFather== "AA" & GTMother=="AB"] = "BA" #M1 is the carrier chromosome
	GTMother[GTRef=="BB" & GTMother=="AB"] = "BA" #M1 is the carrier chromosome
	GTMother[GTRef=="BB" & GTFather== "AA" & GTMother=="AB"] = "NC"
	GTMother[GTRef=="AA" & GTFather== "BB" & GTMother=="AB"] = "NC"

	MatHaps <- do.call("rbind",strsplit(GTMother,split=""))
	PatHap <-	as.matrix(do.call("rbind",strsplit(GTFather,split=""))[,1])

	M1<-as.matrix(MatHaps[,1])
	M2<-as.matrix(MatHaps[,2])

	Bls<- Gtypes[Gtypes$Chr=="X",EmbryoID]#Blastomeres
	Bls <- as.matrix(Bls)

	UnInf<-which(Bls=="NC" | GTMother=="NC"| GTFather=="NC"| GTRef=="NC" | GTMother=="AA" | GTMother=="BB" )

	B1<-matrix(0,nrow(Bls),1)
	B2<-matrix(0,nrow(Bls),1)

	if(is.na(ScSexes[1,2])){
		print(paste("The Chr.X haplotype of Sample",EmbryoID,"could not be reconstructed as the origin of this chromosome could not be determined..."))
	}else if(ScSexes[1,2]=="male"){
		Bls[Bls=="AB"]<-"NC"
		BlHaps<- do.call("rbind",strsplit(Bls,split=""))

		B1[M1==BlHaps[,1]]<-1
		B1[M2==BlHaps[,1]]<-2
		B1[UnInf]<-0
		print(paste("The",EmbryoID,"blastomere is from a male Embryo"))
	}else if (ScSexes[1,2]=="female"){
		Bls[PatHap=="B" & Bls=="AB"]="BA"#Phasing of blastomere based on the father's genotype
		BlHaps<- do.call("rbind",strsplit(Bls,split=""))
		B1[M1==BlHaps[,2]]<-1
		B1[M2==BlHaps[,2]]<-2
		B2[,1]<-1
		B1[UnInf]<-0
		print(paste("The",EmbryoID,"blastomere is from a female Embryo"))
	}
	PhBl<-cbind(B1,B2)
	PhBls<-PhBl

	matHaps <-  PhBls[,seq(1,ncol(PhBls),2)]
	patHaps <-  PhBls[,seq(2,ncol(PhBls),2)]

	matHaps <- as.matrix(matHaps)
	patHaps <- as.matrix(patHaps)
	colnames(matHaps)<-paste(EmbryoID,"_Mat",sep="")
	colnames(patHaps)<-paste(EmbryoID,"_Pat",sep="")

	HapsChrX <- cbind(Gtypes[Gtypes$Chr=="X",c("Names","Chr","Position")],patHaps,matHaps)
	return(HapsChrX)
}#end function

