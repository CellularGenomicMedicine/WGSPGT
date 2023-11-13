calculateBAF <-function(A,colID,fam_members){
	A2 <- data.frame(A,stringsAsFactors=FALSE,check.names=F)
	for(fam_member in colID){
		baf2 <- sapply(A2[,fam_member],function(w) {  Y <- as.numeric(unlist(strsplit(w,","))); Y2 <- Y[2]/(Y[1]+Y[2]);	Y2 })
		if( fam_member == fam_members[1,1]) {baf <- data.frame(baf2,stringsAsFactors=FALSE); names(baf) <- fam_member} else {
			colnmsbaf3 <- names(baf)
			baf <- data.frame(baf,baf2,stringsAsFactors=FALSE)
			names(baf) <- c(colnmsbaf3,fam_member)
		}
	}
	return(baf)
}
