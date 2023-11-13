Haplarithmisis_Phasing_parentsX <- function(Father,Mother,GTFatherX,GTMotherX,GTRefX,RefSex){
	if(RefSex=="male"){
		GTMotherX[GTRefX=="BB" & GTMotherX=="AB"] <- "BA"
	}else if(RefSex=="female"){
		GTMotherX[GTRefX=="BB" & GTFatherX=="BB" & GTMotherX=="AB"] <- "BA"
		GTMotherX[GTRefX=="AB" & GTFatherX=="AA" & GTMotherX=="AB"] <- "BA"
	}else{
		print("Sex of ref. individual is not detemined and chr. X htype is not reliable")
	}
	GTMotherX[GTRefX=="NC" & GTMotherX=="AB"] <- "NC"

	Parents <- data.frame(GTFatherX,GTMotherX,stringsAsFactors=F,check.names=F)
	names(Parents) <- c(Father,Mother)
	return(Parents)
}#end function
