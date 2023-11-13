Haplarithmisis_Htyping_Grandparents <- function(script,Father,Mother,Grandfather,Grandmother, Gtypes, ParScore, parent,Parent1, EmbryoID,Chroms){

	print("*******************************************")
	print("****     Option1Htype analysis...     ****")
	
	ScSexes <- matrix(NA,length(ParScore),2)
	rownames(ScSexes)<-names(ParScore)
	source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_sexdetermination.R",sep="/"))
	ScSexes <- Haplarithmisis_sexdetermination(ScSexes,ParScore)

	source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_ChrXFatherAB.R",sep="/"))
	Gtypes[Gtypes[,"Chr"]=="X",Father] <- Haplarithmisis_ChrXFatherAB(Gtypes,Father)

	#------------------------------------------------------------
	#------------			Phasing chr X			-------------
	source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_chrxhtyping_Grandparents.R",sep="/"))
	HapXMat <- Haplarithmisis_chrxhtyping_Grandparents(Father,Mother,Grandfather,Grandmother,Gtypes,parent,ScSexes,EmbryoID)

	PhasedPar <- Gtypes[,Parent1]
	if(parent=="Mother"){Parent2 = Father}else{Parent2 = Mother}

	PhasedPar[ (Gtypes[,Grandfather]=="AB" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="AB")   |
				   (Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="AB" & Gtypes[,Parent1]=="AB") |
				   (Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="AB") ] <-"BA"

	PhasedPar[(Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="AB") |
				   (Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="AB") |
				   (Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="BB") |
				   (Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="AA") |
				   ((Gtypes[,Grandfather]=="AA" | Gtypes[,Grandmother]=="AA") & Gtypes[,Parent1]=="BB") |
				   ((Gtypes[,Grandfather]=="AA" | Gtypes[,Grandmother]=="AB") & Gtypes[,Parent1]=="BB") |
                                   ((Gtypes[,Grandfather]=="BB" | Gtypes[,Grandmother]=="AA") & Gtypes[,Parent1]=="BB") |
                                   ((Gtypes[,Grandfather]=="BB" | Gtypes[,Grandmother]=="AA") & Gtypes[,Parent1]=="AA") |
                                   ((Gtypes[,Grandfather]=="BB" | Gtypes[,Grandmother]=="AB") & Gtypes[,Parent1]=="AA") |
				   ((Gtypes[,Grandfather]=="BB" | Gtypes[,Grandmother]=="BB") & Gtypes[,Parent1]=="AA") |
				   ((Gtypes[,Grandfather]=="NC" | Gtypes[,Grandmother]=="NC") & Gtypes[,Parent1]=="AB") 	] <- "NC"

	Gtypes[,Parent1] <- PhasedPar

	Hap1 <- rep(0,nrow(Gtypes))
	Hap2 <- rep(0,nrow(Gtypes))
	print(EmbryoID)
	GtypeChild <- Gtypes[,EmbryoID]

	cond1 <- 1
	cond2 <- 2

	Hap1[(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,EmbryoID]=="AB") |
		 (Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,EmbryoID]=="BB") |
		 (Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,EmbryoID]=="BB") |
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,EmbryoID]=="AB") |
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,EmbryoID]=="AA") |
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,EmbryoID]=="AA")] <- cond1

	Hap1[(Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="AA" & Gtypes[,EmbryoID]=="AA") |
		 (Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,EmbryoID]=="AB") |
		 (Gtypes[,Parent1]=="BA" & Gtypes[,Parent2]=="BB" & Gtypes[,EmbryoID]=="AA") |
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="BB" & Gtypes[,EmbryoID]=="BB") |
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,EmbryoID]=="AB") |
		 (Gtypes[,Parent1]=="AB" & Gtypes[,Parent2]=="AA" & Gtypes[,EmbryoID]=="BB")] <- cond2
			 
	Hap1[Gtypes[,Grandfather]=="AB" & Gtypes[,Grandmother]=="AB"] <-0
	ParHaps <- Hap1
	names(ParHaps)<- EmbryoID

	Haps <- vector("list",2)
	names(Haps) <- c("dataHapRaw","dataHap")

	ParHaps[Gtypes$Chr=="X"] <- HapXMat

	ParHapsAll        	<- data.frame(Gtypes[,c("Names","Chr","Position")],ParHaps,stringsAsFactors=FALSE,check.names=F)
	names(ParHapsAll) 	<- c("Names","Chr","Position",EmbryoID)
	dataHapRaw 			<- ParHapsAll

	source(paste(script, "analyses", "Haplarithmisis", "functions", "Haplarithmisis_inthapnew1.R",sep="/"))
	dataHap1			<- Haplarithmisis_inthapnew1(script,dataHapRaw,Window,Int,Chroms)

	if(parent=="Father"){ par1="_Pat"; par2="_Mat" }else{ par1="_Mat"; par2="_Pat" }

	Haps[["dataHapRaw"]] 			<- ParHapsAll
	Haps[["dataHapRaw"]]$mat 		<- 0
	colnames(Haps[["dataHapRaw"]])	<- c("Names","Chr","Position",paste(EmbryoID,par1,sep=""),paste(EmbryoID,par2,sep=""))

	Haps[["dataHap"]] 			<- dataHap1
	Haps[["dataHap"]]$mat 		<- 0
	colnames(Haps[["dataHap"]])	<- colnames(Haps[["dataHapRaw"]])

	return(Haps)
}#end function
