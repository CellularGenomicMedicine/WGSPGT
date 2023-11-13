calculateLogR <-function(bins,colID,filterbamFiles){
	readCounts				<- binReadCounts(bins,bamfiles=filterbamFiles[,1],chunkSize=100000000)
	readCounts 				<- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
	readCounts 				<- estimateCorrection(readCounts)
	readCounts 				<- applyFilters(readCounts, chromosomes=NA)
	copyNumbers				<- correctBins(readCounts)
	copyNumbersNormalized 	<- normalizeBins(copyNumbers)
	copyNumbersSmooth 		<- smoothOutlierBins(copyNumbersNormalized)
	dat   					<- assayDataElement(copyNumbersSmooth, "copynumber")
	colnames(dat) 			<- colID
	logRs 					<- data.frame(Names=featureNames(copyNumbersSmooth) , Chr=fData(copyNumbersSmooth)$chromosome, Position=as.integer(fData(copyNumbersSmooth)$start) , as.data.frame(dat,stringsAsFactors=FALSE), check.names=FALSE, stringsAsFactors=FALSE)
	logRs					<- logRs[rowSums(!is.na(logRs[,colID]))==ncol(dat),]
	logRs[,colID]  			 <- round(log2(logRs[,colID]+0.0005), digits=3)
	return(logRs)
}
