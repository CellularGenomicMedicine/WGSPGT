Haplarithmisis_phgrandparents <- function(Gtypes,parent,Parent1,Father,Mother,Grandfather,Grandmother){

	MotherChrX <- Gtypes[Gtypes$Chr=="X",Mother]
	GtypesChrX <- Gtypes[Gtypes$Chr=="X",]
	PhasedPar <- Gtypes[,Parent1]

	PhasedPar[(Gtypes[,Grandfather]=="AB" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="AB") |
				   (Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="AB" & Gtypes[,Parent1]=="AB") |
				   (Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="AB")]<-"BA"

	PhasedPar[(Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="AA" & Gtypes[,Parent1]=="AB") |
				   (Gtypes[,Grandfather]=="BB" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="AB") |
				   (Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="BB") |
				   (Gtypes[,Grandfather]=="AA" & Gtypes[,Grandmother]=="BB" & Gtypes[,Parent1]=="AA") |
				   ((Gtypes[,Grandfather]=="NC" | Gtypes[,Grandmother]=="NC") & Gtypes[,Parent1]=="AB") |
				   (Gtypes[,Grandfather]=="AB" & Gtypes[,Grandmother]=="AB")] <- "NC"

	Gtypes[,Parent1] <- PhasedPar

	if(parent=="Mother"){

		MotherChrX[(GtypesChrX[,Grandfather]=="BB" & GtypesChrX[,Mother]=="AB")]<-"BA"

		MotherChrX[(GtypesChrX[,Grandfather]=="NC" & GtypesChrX[,Mother]=="AB") |
					   (GtypesChrX[,Grandfather]=="BB" & GtypesChrX[,Grandmother]=="BB" & GtypesChrX[,Mother]=="AB") |

					   (GtypesChrX[,Grandfather]=="AA" & GtypesChrX[,Grandmother]=="AA" & GtypesChrX[,Mother]=="AB") ]<-"NC"
		Gtypes[Gtypes$Chr=="X",Parent1] <- MotherChrX
		GTFather=Gtypes[,Father]
		GTMother=Gtypes[,Parent1]
	}	   

	if(parent=="Father" ){GTFather=Gtypes[,Parent1];GTMother=Gtypes[,Mother]}

	Parents <- data.frame(GTFather,GTMother,stringsAsFactors=F,check.names=F)
	names(Parents)  <- c(Father,Mother)
	return(Parents)

}#end function
