library(VariantAnnotation)
library(vcfR)
library(snpStats)
library(data.table)
library(Biobase)
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
	Chroms 		<- c(Chroms,"Y")
	gvcf_file 	<- paste(GATK,"/",family,".g.vcf",sep="")
	source(paste(script, "analyses", "functions","checkFileExists.R", sep="/"))
	checkFileExists(gvcf_file,inputPath)

	#read metainfo from samplesheet
	samplesheet           <- read.table(samplesheet_path, sep = ",", check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
	fam_members				<- samplesheet[,c("SampleID","Sample_MetaInfo","Sample_Status")]
	fam_members_nonEmbryo 	<- fam_members[!fam_members[,"Sample_Status"]%in%"E",]
	fam_members_Embryo 		<- fam_members[fam_members[,"Sample_Status"]%in%"E",]
	colID 					<- fam_members[,1]

	print("Started reading vcfR for Crd and DP")
	vcf           	<- read.vcfR (gvcf_file, verbose=FALSE)
	Crd           	<- data.frame(Chr=gsub("chr","",getCHROM(vcf)),Position=as.numeric(getPOS(vcf)),stringsAsFactors=FALSE)
	Names         	<- paste0("chr",paste(Crd[,1],Crd[,2],sep=":"))
	Crd				<- data.frame(Crd,Names=Names,stringsAsFactors=F)
	Alleles       	<- data.frame(as.character(getREF(vcf)),getALT(vcf),stringsAsFactors=F)
	Alleles		 	<- data.frame(Alleles,Names=Names,stringsAsFactors=F)

	print("Determine SNP coverage")
	DP            	<- extract.gt(vcf,element='DP')
	DP				<- data.frame(Names=Names,DP,stringsAsFactors=F,check.names=F)
	CrdDP			<- merge(Crd,DP,by="Names")
	CrdDP			<- CrdDP[order(CrdDP$Chr,CrdDP$Position),]

	print("Started reading Vcf for GT")
	vcf1          <- readVcf(gvcf_file, verbose=FALSE)
	res           <- genotypeToSnpMatrix(vcf1,uncertain=FALSE)
	GT            <- as.data.frame(t(as(res$genotype, "character")) ,stringsAsFactors=F,check.names=F)

	for(col in 1:ncol(GT)){
		GT[,col] <- gsub("/","",GT[,col])
		GT[,col] <- gsub("NA","NC",GT[,col])
	}

	GT$ID 		<- row.names(GT)
	GT$Names 	<- sapply(GT$ID,function(w) {unlist(strsplit(w,"_"))[1]})
	GTcombine	<- GT
	CrdGT		<- merge(Crd,GT,by="Names")
	CrdGT		<- CrdGT[order(CrdGT$Chr,CrdGT$Position),]
	#CALCULATE BAF
	print("Started calculating BAF")
	source(paste(script, "analyses", "ConvertGenotype", "functions","calculateBAF.R", sep="/"))
	A           <- extract.gt(vcf,element='AD')
	baf 		<- calculateBAF(A,colID,fam_members)
	bafcombine	<- baf
	baf 		<- data.frame(Names=Names,baf,stringsAsFactors=FALSE,check.names=F)
	Crdbaf		<- merge(Crd,baf,by="Names")
	Crdbaf		<- Crdbaf[order(Crdbaf$Chr,Crdbaf$Position),]#filter on GT
	#combine data files
	colnames(GT)  	<- gsub("$",".GType",colnames(GT))
	colnames(baf) 	<- paste0(colnames(baf),".B Allele Freq")
	Gtype			<- data.frame(Crd,GTcombine,bafcombine,stringsAsFactors=F,check.names=F)
	# *** Data validation
	Gtype 		<- Gtype[nchar(Alleles[,1])==1 & nchar(Alleles[,2])==1,]
	source(paste(script, "analyses", "functions","validateGtype.R", sep="/"))
	validateGtype(Gtype,Chroms,inputPath,Chrom)
	#Write data
	Crd 		<- Crd[Crd$Names%in%Gtype$Names,] #filter on GT
	Crdbaf 		<- Crdbaf[Crdbaf$Names%in%Gtype$Names,]
	CrdGT		<- CrdGT[CrdGT$Names%in%Gtype$Names,]
	CrdDP		<- CrdDP[CrdDP$Names%in%Gtype$Names,]
	source(paste(script, "analyses", "functions","writeData.R", sep="/"))
	write.table(Crd,paste(dataPath,paste(family,"_","Crd",".txt",sep=""),sep="/"),row.names=F, quote=F, sep="\t")
	writeData(CrdGT,family,fam_members,parents,"GT",dataPath)
	writeData(Crdbaf,family,fam_members,parents,"BAF",dataPath)
	writeData(CrdDP,family,fam_members,parents,"DP",dataPath)

}, error = function(e) {
	write(paste0("ConvertGenotype gave the following error: ", e), errorFilePath, sep="")
	stop(paste0("ConvertGenotype gave the following error: ", e))
})