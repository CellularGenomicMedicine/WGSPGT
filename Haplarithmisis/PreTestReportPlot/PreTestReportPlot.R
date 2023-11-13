library(VariantAnnotation)
library(vcfR)
library(snpStats)
library(data.table)
library(QDNAseq)
library(Biobase)
library(plotrix)
options(scipen=999)

args			<- commandArgs(TRUE)
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
	Chroms 				<- c(Chroms,"Y")

	source(paste(script, "analyses", "functions","checkDirExistsAndCreate.R", sep="/"))
	outPath  <- checkDirExistsAndCreate(paste(inputPath,"PreTestReport",sep="/"))

	#Check and load data files
	source(paste(script, "analyses", "functions","checkFileExistsAndFread.R", sep="/"))
	logRs 		<- checkFileExistsAndFread(paste(logRData,paste(family,"_LogR.txt",sep=""),sep="/"),inputPath)
	SeglogRs 	<- checkFileExistsAndFread(paste(logRData,paste(family,"_SegLogR.txt",sep=""),sep="/"),inputPath)
	baf			 <- checkFileExistsAndFread(paste(dataPath,paste(family,"_BAF.txt",sep=""),sep="/"),inputPath)

	#read metainfo from samplesheet
	samplesheet           <- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
	fam_members				<- samplesheet[,c("SampleID","Sample_MetaInfo","Sample_Status")]
	colID <- fam_members[,1]

	#Create helix data informative and genomewide SNPS etc
	for(parent in parents){
		Interval			<- paste(inputPath,paste(family,parent,"intervals.txt",sep="_"),sep="/")
		Int					<- read.table(Interval,sep="\t",header=T,stringsAsFactors=F)
		upstream_boundary 	<- as.numeric(Int[Int[,1]==Chrom,3])
		downstream_boundary	<- as.numeric(Int[Int[,1]==Chrom,2])

		Father <- fam_members[grepl("Father",fam_members[,2]),1]
		Mother <- fam_members[grepl("Mother",fam_members[,2]),1]
		Refs   <- fam_members[!grepl("Mother",fam_members[,2])&!grepl("Father",fam_members[,2]),1]
		Parent <- fam_members[grepl(parent,fam_members[,2]),1]
		baf_int <- subset(baf, Chr==Chrom)
		source(paste(script, "analyses", "PreTestReportPlot", "functions","PreTestReportPlot_bafplots.R", sep="/"))
		PreTestReport_bafplots(script,Chrom,upstream_boundary,downstream_boundary,ChrsLength,ideogram,outPath,colID,fam_members,baf_int,parent)

		#plot XY logRs
		source(paste(script, "analyses","functions","XYplot.R",sep="/"))
		XYplot(script,logRs,SeglogRs,ChrsLength,ideogram,outPath,colID,fam_members)
	}
}, error = function(e) {

	write(paste0("PreTestReport gave the following error: ", e), errorFilePath, sep="")
	stop(paste0("PreTestReport gave the following error: ", e))
})