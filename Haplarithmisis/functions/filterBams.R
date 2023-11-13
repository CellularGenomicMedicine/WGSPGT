filterBams <-function(bamFiles,fam_members){
	filter_bams <- paste0("/",fam_members[,1],".bam")
	filterbamFiles <- data.frame(stringsAsFactors=F,row.names=F)
	for(filter_bam in filter_bams){
		filterbamFile <- data.frame(bamfiles=bamFiles[grep(filter_bam,bamFiles[,1]),],stringsAsFactors=F,row.names=grep(filter_bam,filter_bams))
		if(filter_bam==filter_bams[1]) {
			filterbamFiles <- filterbamFile
		} else {
			filterbamFiles <- rbind(filterbamFiles,filterbamFile)
		}
	}
	filterbamFiles
}
