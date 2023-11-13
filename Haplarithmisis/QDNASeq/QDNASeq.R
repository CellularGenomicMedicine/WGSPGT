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
errorFilePath 	<- args[2]

print(paste("Potential errors written to:", errorFilePath, sep = " "))
# clear any existing error file (in case of rerun)
if (file.exists(errorFilePath)) {
	warning("Error file already exists. Removing file and resuming program...")
	file.remove(errorFilePath)
}
tryCatch({
	source(as.character(config_file_fam))

	load(paste(script,"analyses","Rda",ideogram,sep="/"))
	load(paste(script,"analyses","Rda",QDNASeqBins,sep="/"))

	Chroms 				<- c(Chroms,"Y")

	source(paste(script, "analyses", "functions","checkDirExistsAndCreate.R", sep="/"))
	outPath  <- checkDirExistsAndCreate(paste(inputPath, "QDNASeq", sep = "/"))
	write(paste0("logRData=",'"',outPath,'"'),config_file_fam,append=TRUE)
	#read metainfo from samplesheet
	samplesheet           	<- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
	fam_members				<- samplesheet[,c("SampleID","Sample_MetaInfo","Sample_Status")]
	fam_members_nonEmbryo 	<- fam_members[!fam_members[,"Sample_Status"]%in%"E",]
	fam_members_Embryo 		<- fam_members[fam_members[,"Sample_Status"]%in%"E",]
	colID <- fam_members[,1]
	bamFiles 				<- read.table(paste(GATK,paste(family,"list",sep="."),sep="/"),sep=",",check.names=FALSE,stringsAsFactors=FALSE)

	source(paste(script, "analyses", "functions","filterBams.R", sep="/"))
	filterbamFiles <- filterBams(bamFiles,fam_members)

	print("Started calculating LogR")
	source(paste(script, "analyses", "QDNASeq", "functions","calculateLogR.R", sep="/"))
	logRs <- calculateLogR(bins,colID,filterbamFiles)
	source(paste(script, "analyses", "functions","writeData.R", sep="/"))
	writeData(logRs,family,fam_members,parents,"LogR",outPath)

	print(paste("Started calculating SegLogR for gamma",gammaMC,sep=" "))
	source(paste(script, "analyses", "functions","SegLogRs.R", sep="/"))
	segLogRs 		<- SegLogRs(script,logRs,gammaMC,gammaSC,plateau,inputPath,fam_members_Embryo)
	writeData(segLogRs,family,fam_members,parents,"SegLogR",outPath)
}, error = function(e) {
	write(paste0("QDNASeq gave the following error: ", e), errorFilePath, sep="")
	stop(paste0("QDNASeq gave the following error: ", e))
})