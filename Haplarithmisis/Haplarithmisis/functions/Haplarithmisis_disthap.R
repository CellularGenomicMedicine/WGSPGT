Haplarithmisis_disthap <-function(Intp,Int,outPathData,EmbryoID,Chrom){
	HapsIntp <- Intp[["dataHap"]]
	DistFam <- NULL
	Scs <- colnames(HapsIntp)[colnames(HapsIntp)%in%paste(EmbryoID,"Mat",sep="_")|colnames(HapsIntp)%in%paste(EmbryoID,"Pat",sep="_")]
	for(ind in Scs){
		IntChr <- Int[Int[,1]==Chrom,]
		HapChr <- HapsIntp[HapsIntp[,"Chr"]==Chrom,c("Position",ind)]
		UP <- HapChr[HapChr$Position<=as.numeric(IntChr[2]),]
		Down <- HapChr[HapChr$Position>=as.numeric(IntChr[3]),]

		UPRle <- rle(UP[,ind])
		DownRle <- rle(Down[,ind])
		Dist <- cbind(ind,Chrom,UPRle$lengths[length(UPRle$values)],DownRle$lengths[1],UPRle$values[length(UPRle$values)], DownRle$values[1])#}
		Dists <- Dist

		if(ind==Scs[1]){DistFam <- Dists }else{DistFam <-rbind(DistFam,Dists)}
		print(ind)
	}#end ind loop
	colnames(DistFam) <- c("Haplotype","Chr","LengthUp","LengthDown","ValueUp","ValueDown")
	write.table(DistFam,paste(outPathData,"DistanceToHR.txt",sep="/"),quote=F,sep="\t",col.names=T,row.names=F)
}
