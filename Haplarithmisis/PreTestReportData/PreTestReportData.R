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
dbsnp_path      <- args[3]

print(paste("Potential errors written to:", errorFilePath, sep = " "))
# clear any existing error file (in case of rerun)
if (file.exists(errorFilePath)) {
	warning("Error file already exists. Removing file and resuming program...")
	file.remove(errorFilePath)
}

tryCatch({
	source(as.character(config_file_fam))

	load(paste(script,"analyses","Rda",ideogram,sep="/"))
	#Chroms 				<- c(Chroms,"Y")

	source(paste(script, "analyses", "functions","checkDirExistsAndCreate.R", sep="/"))
	outPath  <- checkDirExistsAndCreate(paste(inputPath,"PreTestReport",sep="/"))

	source(paste(script, "analyses", "functions","checkFileExistsAndFread.R", sep="/"))
	GT<- checkFileExistsAndFread(paste(dataPath,paste(family,"_GT.txt",sep=""),sep="/"),inputPath)
	GT<- na.omit(GT)
	Chroms 	<- unique(GT[,"Chr"])
	#read metainfo from samplesheet
	samplesheet          	<- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
	fam_members				<- samplesheet[,c("SampleID","Sample_MetaInfo","Sample_Status")]
	colID 					<- fam_members[,1]

	NGS_PGD_GENOMEWIDE_list <- list()
	#Create helix data informative and genomewide SNPS etc
	for(parent in parents){
		Interval			<- paste(inputPath,paste(family,parent,"intervals.txt",sep="_"),sep="/")
		Int					<- read.table(Interval,sep="\t",header=T,stringsAsFactors=F)
		Parent <- fam_members[fam_members[,"Sample_MetaInfo"]%in%parent,"SampleID"]

		#calculate informative SNPS and generate stoftests
		#mendelian inconsistency
		source(paste(script, "analyses", "PreTestReportData","functions","PreTestReportData_mendInc.R",sep="/"))
		PreTestReportData_mendinc(script,config_file_fam,Parent,REF,Father,Mother,GT,outPath,Chroms)
		#phasing and extracting informative SNPs
		source(paste(script, "analyses", "PreTestReportData","functions","PreTestReportData_phasing.R",sep="/"))
		NGS_PGD_GENOMEWIDE_list <- PreTestReportdata_phasing(script,config_file_fam,GT,Parent,Father,Mother,REF,Chroms,Int,fam_members,dbsnp_path,parent,HelixoutPath,Chrom,NGS_PGD_GENOMEWIDE_list)
	}
	NGS_PGD_GENOMEWIDE_df <- rbind(NGS_PGD_GENOMEWIDE_list[["TOT_AANT"]],NGS_PGD_GENOMEWIDE_list[["NGS_PGD_GENOMEWIDE_M_df"]],NGS_PGD_GENOMEWIDE_list[["NGS_PGD_GENOMEWIDE_P_df"]])
	write.table(NGS_PGD_GENOMEWIDE_df,paste(paste(HelixoutPath,"",sep="/"),Mother,"-NGS_PGD_GENOMEWIDE",sep=""),col.names=T,row.names=F,quote=F,sep=",")

}, error = function(e) {
	write(paste0("PreTestReport gave the following error: ", e), errorFilePath, sep="")
	stop(paste0("PreTestReport gave the following error: ", e))
})