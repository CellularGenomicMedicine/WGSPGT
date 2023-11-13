Haplarithmisis_inthapnew1 <- function(script,dataHapRaw,Window,Int,Chroms){
	source(paste(script, "analyses", "Haplarithmisis", "functions","Haplarithmisis_testmedfilt.R",sep="/"))
	HapsInt <- vector("list",(ncol(dataHapRaw)-length(c("Names","Chr","Position"))))
	names(HapsInt) <- colnames(dataHapRaw)[!colnames(dataHapRaw)%in%c("Names","Chr","Position")]

	#Naively exclude the putative abbarant regions
	for(i in 1:nrow(Int)){

		dataHapRaw[dataHapRaw$Chr==Int[i,1] & as.numeric(dataHapRaw$Position) >= as.numeric(Int[i,2]) & as.numeric(dataHapRaw$Position) <= as.numeric(Int[i,3]) ,grep(Int[i,4],colnames(dataHapRaw))] <- 0

	}#end i loop

	for(h in colnames(dataHapRaw)[!colnames(dataHapRaw)%in%c("Names","Chr","Position")]){

		Hap <- dataHapRaw[dataHapRaw[,h]!=0,c("Names","Chr","Position",h)]
		rownames(Hap)<-Hap$Names
		HapMed <- Haplarithmisis_testmedfilt(Hap,Window,Chroms)
		HapsInt[[h]]<- cbind(Hap,HapMed)
	}#end h loop
	HapsInt2<- vector("list",(ncol(dataHapRaw)-3))
	names(HapsInt2) <- colnames(dataHapRaw)[!colnames(dataHapRaw)%in%c("Names","Chr","Position")]

	for(h in names(HapsInt)){
		#Mapping back the smoothed haplotypes
		HapMed <- matrix(0,nrow(dataHapRaw),1)
		rownames(HapMed)<- dataHapRaw$Names
		HapMed[rownames(HapsInt[[h]]),]<-HapsInt[[h]][,"HapMed"]
		HapMed2 <- cbind(dataHapRaw[,c("Names","Chr","Position")],dataHapRaw[,h],HapMed)

		for(chr in Chroms){
			#Defining chr-specific window size
			HapBlocks <- rle(as.numeric( HapMed2[HapMed2$Chr==chr,"HapMed"]))
			if(length(HapBlocks$values)>2){
				for(k in 2:(length(HapBlocks$values)-1)){
					if(HapBlocks$values[k]==0 & HapBlocks$values[k+1] == HapBlocks$values[k-1]){
						HapBlocks$values[k] <- HapBlocks$values[k-1]}
				}#end k loop
			}#end if statement
			HapIntChr <- inverse.rle(HapBlocks)
			if(chr==Chroms[1]){HapIntGenome <- HapIntChr}else{HapIntGenome <- c(HapIntGenome,HapIntChr)}
		}#end chr loop

		if(h==names(HapsInt2)[1]){HapInds <- HapIntGenome}else{HapInds <- cbind(HapInds,HapIntGenome)}
		HapsInt2[[h]] <- cbind(HapMed2,HapIntGenome)
		colnames(HapsInt2[[h]])[4] <- "Raw"
		print(h)

	}#end h loop
	dataHap <- cbind(dataHapRaw[,c("Names","Chr","Position")],HapInds)
	colnames(dataHap) <- c("Names","Chr","Position",names(HapsInt2))
	return(dataHap)

}#end function