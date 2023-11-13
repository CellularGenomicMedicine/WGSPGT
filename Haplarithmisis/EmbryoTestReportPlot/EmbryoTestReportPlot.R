library(VariantAnnotation)
library(vcfR)
library(snpStats)
library(data.table)
library(QDNAseq)
library(Biobase)
library(plotrix)
options(scipen=999)

args		<- commandArgs(TRUE)
config_file_fam <- args[1]
EmbryoID        <- args[2]
Gamma_value     <- as.numeric(args[3])
errorFilePath   <- args[4]

if (file.exists(errorFilePath)) {
	warning("Error file already exists. Removing file and resuming program...")
	file.remove(errorFilePath)
}
tryCatch({
	source(as.character(config_file_fam))
	load(paste(script,"analyses","Rda",ideogram,sep="/"))
	Window   <- as.numeric(gtypemodulator_window)

	source(paste(script, "analyses", "functions","checkDirExistsAndCreate.R", sep="/"))
	SharePath		<- checkDirExistsAndCreate(paste(inputPath,EmbryoID,sep="/"))
	outPathPlots	<- checkDirExistsAndCreate(paste(SharePath,Gamma_value,"",sep="/"))

	samplesheet           <- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
	fam_members           <- samplesheet[, c("SampleID", "Sample_MetaInfo", "Sample_Status")]
	fam_members_parents   <- fam_members[fam_members[, "Sample_MetaInfo"] %in% c("Father", "Mother"),]
	fam_members_nonEmbryo <- fam_members[!fam_members[, "Sample_Status"] %in% "E",]
	fam_members_Embryo    <- fam_members[fam_members[, "SampleID"] %in% EmbryoID,]
	fam_members           <- rbind(fam_members_nonEmbryo, fam_members_Embryo)
	colID                 <- fam_members[,"SampleID"]

	for (parent in parents){
		source(paste(script, "analyses", "functions","checkFileExistsAndFread.R", sep="/"))
		dataHap		<- checkFileExistsAndFread(paste(HaplarithmisisData,paste(parent,EmbryoID,family,Gamma_value,"dataHap.txt",sep="_"), sep = "/"),inputPath)
		dataHapRaw	<- checkFileExistsAndFread(paste(HaplarithmisisData,paste(parent,EmbryoID,family,Gamma_value,"dataHapRaw.txt",sep="_"),sep="/"),inputPath)
		dataPo		<- checkFileExistsAndFread(paste(HaplarithmisisData,paste(parent,EmbryoID,family,Gamma_value,"dataPo.txt",sep="_"),sep="/"),inputPath)
		logRs		<- checkFileExistsAndFread(paste(logRData,paste(parent,EmbryoID,family,"LogR.txt",sep="_"),sep="/"),inputPath)
		SegLogRs	<- checkFileExistsAndFread(paste(logRData,paste(parent,EmbryoID,family,"SegLogR.txt",sep="_"),sep="/"),inputPath)
		BAF			<- checkFileExistsAndFread(paste(HaplarithmisisData,paste(parent,EmbryoID,family,Gamma_value,"BAF.txt",sep="_"),sep="/"),inputPath)
		Gtypes		<- checkFileExistsAndFread(paste(HaplarithmisisData,paste(parent,EmbryoID,family,Gamma_value,"Gtypes.txt",sep="_"),sep="/"),inputPath)
		PhBAF <- vector("list",8)
		for(seg in c("P1","P2","M1","M2","P1Seg","P2Seg","M1Seg","M2Seg")){
			PhBAF[[seg]] <- checkFileExistsAndFread(paste(HaplarithmisisData,paste(parent,EmbryoID,family,Gamma_value,paste0(seg,".txt"),sep="_"), sep = "/"),inputPath)
		}
		Int	    <- read.table(paste(inputPath, paste(family, parent, "intervals.txt", sep = "_"), sep = "/"),sep="\t",header=T,stringsAsFactors=F)

		#filter BAFs for Gtypes overlap
		BAFs    <- BAF[BAF$Names%in%Gtypes$Names,]
		if( exists("Chrom")){
			if (Chrom=="X" & Gamma_value == 50){
				source(paste(script, "analyses", "functions","XYplot.R",sep="/"))
				XYplot(script,logRs,SegLogRs,ChrsLength,ideogram,outPathPlots,fam_members_Embryo[,"SampleID"],fam_members)
			}
			up  <- if(any(Int[Int[,1]!=0,1]==Chrom)){Int[Int[,1]==Chrom,3]}
			down <- if(any(Int[Int[,1]!=0,1]==Chrom)){Int[Int[,1]==Chrom,2]}
			upcalc <- up+2000000
			downcalc <- down-2000000
			if (downcalc < 0) {downcalc=0}
				source(paste(script,"analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_combinedplot.R",sep="/"))
				EmbryoTestReportPlot_combinedplot(script,dataHap,dataHapRaw,dataPo,BAFs,logRs,SegLogRs,PhBAF,ChrsLength,ideogram,parent,EmbryoID,outPathPlots,Chrom,Int,flip)
		} else {
			source(paste(script,"analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_genomebafplot.R",sep="/"))
			EmbryoTestReportPlot_genomebafplot(script,BAFs,logRs,SegLogRs,PhBAF,ideogram,parent,Father,Mother,EmbryoID,outPathPlots)
		}
	}
}, error = function(e) {
	write(paste0("EmbryoTestReportPlot gave the following error: ", e), errorFilePath, sep="")
	stop(paste0("EmbryoTestReportPlot gave the following error: ", e))
})