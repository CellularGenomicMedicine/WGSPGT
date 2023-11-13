library(data.table)
options(scipen=999)

args			<- commandArgs(TRUE)
config_file_fam <- args[1]
errorFilePath 	<- args[2]

if (file.exists(errorFilePath)) {
	warning("Error file already exists. Removing file and resuming program...")
	file.remove(errorFilePath)
}
tryCatch({
	source(as.character(config_file_fam))
	load(paste(script,"analyses","Rda",ideogram,sep="/"))

	Chroms 				<- c(Chroms,"Y")

	source(paste(script, "analyses", "functions","checkFileExistsAndFread.R", sep="/"))
	Crd		<- checkFileExistsAndFread(paste(dataPath,paste(family,"_Crd.txt",sep=""),sep="/"),inputPath)

	for (parent in parents){
		#remove Mitochondrial Chromosome, can cause pileup of strange haplotypes and is a contamination
		Crd           	<- Crd[Crd[,"Chr"]!="M",]
		SnpAnotFile		<- Crd[,c("Names","Chr","Position")]
		Interval		<- data.frame(start=SnpAnotFile$Position-round(Window/2),end=SnpAnotFile$Position+round(Window/2),stringsAsFactors=FALSE)
		ChrsLength		<- data.frame(ChrsLength)
		ChrsLength[,1]	<- as.character(ChrsLength[,1])
		ChrsLength[,2]	<- as.numeric(as.character(ChrsLength[,2]))
		Bed				<- data.frame(Chr=SnpAnotFile$Chr,Interval,Name=SnpAnotFile$Name,stringsAsFactors=FALSE)

		for(chr in Chroms){
			BedChr    <- Bed[Bed[,1]==chr,]
			if(nrow(BedChr)>1){
				ChrLength <- ChrsLength[ChrsLength[,1]==paste("chr",chr,sep=""),2]
				if(sum(BedChr[,2]<0)>=1){
					BedChr[BedChr[,2]<0,2]<-0
				}
				if(BedChr[nrow(BedChr),3]>ChrLength){
					BedChr[nrow(BedChr),3]<-ChrLength
				}
				if(BedChr[nrow(BedChr),2]>ChrLength){
					BedChr[nrow(BedChr),2]<-ChrLength
				}
				if(chr==Chroms[1]){
					BedGenome<-BedChr
				}else{
					BedGenome<-rbind(BedGenome,BedChr)
				}
			}
		}#end chr loop

		BedGenome[,1] <- paste0("chr",BedGenome[,1],sep="")
		write.table(BedGenome,paste(dataPath,paste("Window",as.numeric(gtypemodulator_window),".bed",sep=""),sep="/"),sep="\t",quote=F,col.names=F,row.names=F)
	}
}, error = function(e) {
	write(paste0("NucBedPrep gave the following error: ", e), errorFilePath, sep="")
	stop(paste0("NucBedPrep gave the following error: ", e))
})