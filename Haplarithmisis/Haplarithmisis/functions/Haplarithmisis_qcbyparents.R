Haplarithmisis_qcbyparents <- function(Father,Mother,Gtypes,EmbryoID,Chroms){
	print("QC by parents analysis...")
	QC_par <- vector("list",2)
	names(QC_par)<-c("ADO","ADI")
	#Genome-wide QC by parents
	GTFather <- Gtypes[,Father][Gtypes$Chr!="X" | Gtypes$Chr=="XY" | Gtypes$Chr=="Y"]
	GTMother <- Gtypes[,Mother][Gtypes$Chr!="X" | Gtypes$Chr=="XY" | Gtypes$Chr=="Y"]
	GTSib <- Gtypes[Gtypes$Chr!="X" | Gtypes$Chr=="XY" | Gtypes$Chr=="Y",EmbryoID]
	
	ADO_AutGenome <- (sum((GTFather=="AB" & GTMother=="BB" & GTSib=="AA") |
						   				   (GTFather=="AB" & GTMother=="AA" & GTSib=="BB") |
						    			   (GTFather=="BB" & GTMother=="AB" & GTSib=="AA") |
						    			   (GTFather=="AA" & GTMother=="AB" & GTSib=="BB") |
						    			   (GTFather=="AA" & GTMother=="BB" & (GTSib=="AA" |  GTSib=="BB")) |
						    			   (GTFather=="BB" & GTMother=="AA" & (GTSib=="AA" |  GTSib=="BB")))/sum((GTSib=="AA" | GTSib=="BB") & (GTFather!="NC" | GTMother!="NC")))*100
	
	ADI_AutGenome <- (sum((GTFather=="AA" & GTMother=="AA" & GTSib=="AB") |
   						 				  (GTFather=="BB" & GTMother=="BB" & GTSib=="AB") )/sum(GTSib=="AB" & (GTFather!="NC" | GTMother!="NC")))*100

	
	#Chr-specific QC by parents
	for(chr in Chroms){
	
		FatherChr <- Gtypes[,Father][Gtypes$Chr==chr]
		MotherChr <- Gtypes[,Mother][Gtypes$Chr==chr]
		SibChr <- Gtypes[Gtypes$Chr==chr,EmbryoID]
		
		ADO <- (sum((FatherChr=="AB" & MotherChr=="BB" & SibChr=="AA") |
						    (FatherChr=="AB" & MotherChr=="AA" & SibChr=="BB") | 
						    (FatherChr=="BB" & MotherChr=="AB" & SibChr=="AA") |
						    (FatherChr=="AA" & MotherChr=="AB" & SibChr=="BB") |
						    (FatherChr=="AA" & MotherChr=="BB" & (SibChr=="AA" |  SibChr=="BB")) |
						    (FatherChr=="BB" & MotherChr=="AA" & (SibChr=="AA" |  SibChr=="BB")))/sum((GTSib=="AA" | GTSib=="BB") & (GTFather!="NC" | GTMother!="NC")))*100
	
		ADI <- (sum((FatherChr=="AA" & MotherChr=="AA" & SibChr=="AB") |
   						   (FatherChr=="BB" & MotherChr=="BB" & SibChr=="AB") )/sum(GTSib=="AB" & (GTFather!="NC" | GTMother!="NC")))*100
						  
	if(chr==Chroms[1]){ADOGenome <- ADO; ADIGenome <- ADI }else{ADOGenome <- cbind(ADOGenome,ADO); ADIGenome <- cbind(ADIGenome,ADI) }

	}#end chr loop
	
	ADOGenome <- cbind(ADOGenome,ADO_AutGenome)
	ADIGenome <- cbind(ADIGenome,ADI_AutGenome)

	ADOs <- ADOGenome; ADIs <- ADIGenome

	colnames(ADOs) <- c(Chroms,"GenomeAut")
	colnames(ADIs) <- c(Chroms,"GenomeAut")
	rownames(ADOs) <- EmbryoID
	rownames(ADIs) <- EmbryoID
	
	QC_par[["ADO"]] <- ADOs
	QC_par[["ADI"]] <- ADIs
	return(QC_par)
}#end function 
