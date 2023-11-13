Intervals_extraction <- function(outPath,samplesheet){
	status_comp     <- samplesheet[,"Sample_Status"]
	Value1_stat  	<- status_comp[grepl("AF",status_comp)|grepl("AM",status_comp)]
	Family_Interval <- unique(samplesheet[,"Family_Interval"])
	if(length(Family_Interval)>1) {	write.table(Family_Interval,paste(outPath,paste("Multiple intervals.txt",sep="_"),sep="/"),row.names=FALSE,quote=FALSE,na=''); stop("Multiple intervals")}
	interval_1 <- data.frame(t(unlist(strsplit(Family_Interval,"_"))),stringsAsFactors=FALSE)
	Family_second_Interval <- unique(samplesheet[,"Family_second_Interval"])
	if(length(Family_second_Interval)>1) {	write.table(Family_Interval,paste(outPath,paste("Multiple intervals.txt",sep="_"),sep="/"),row.names=FALSE,quote=FALSE,na=''); stop("Multiple intervals")}
	interval_2 <- data.frame(t(unlist(strsplit(Family_second_Interval,"_"))),stringsAsFactors=FALSE)
	if(ncol(interval_2)>1) { intervals <- rbind(interval_1,interval_2)}
	if(ncol(interval_2)==1) { intervals <- interval_1}
	for(interval in 1:nrow(intervals)){
		for(v1 in Value1_stat) {
			if (intervals[interval,4]=="Pat") {parent="Father"}
			if (intervals[interval,4]=="Mat") {parent="Mother"}
			Interval <- data.frame(Chr=c(0,gsub("chr","",intervals[interval,1])),Start=c(0,intervals[interval,2]),Stop=c(0,intervals[interval,3]),Origin=intervals[interval,4],stringsAsFactors=FALSE)
			write.table(Interval,paste(outPath,paste(family,parent,"intervals.txt",sep="_"),sep="/"),sep="\t",row.names=FALSE,quote=FALSE,na='')
		}
	}
}
