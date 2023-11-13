Haplarithmisis_writePhBAF <-function(PhBAF,outPathData,parent,EmbryoID,family,Gamma_value){
	for(PhBAF_segs in names(PhBAF)){
		print(PhBAF_segs)
		write.table(PhBAF[[PhBAF_segs]],paste(outPathData,paste0(parent,"_",EmbryoID,"_",family,"_",Gamma_value,"_",PhBAF_segs,".txt"),sep="/"),row.names=F, quote=F, sep="\t")
	}
}
