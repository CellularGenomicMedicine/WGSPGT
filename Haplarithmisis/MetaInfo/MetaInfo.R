#MetaInfo
#Extract MetaInfo from PGT-samplesheet
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

	onderzoeksnummer 	<- gsub("/","_",onderzoeksnr)
	write(paste0("onderzoeksnummer=",'"',onderzoeksnummer,'"'),config_file_fam,append=TRUE)

	inputPath 			<- paste(analyses,family,paste(onderzoeksnummer,indication,sep="-"),sep="/")
	#add to path to config_file
	write(paste0("inputPath=",'"',inputPath,'"'),config_file_fam,append=TRUE)

	source(paste(script, "analyses", "functions","checkFileExists.R", sep="/"))
	samplesheet_path <- paste(inputPath,paste(paste(family,paste(onderzoeksnummer,indication,sep="-"),sep="_"),"csv",sep="."),sep="/")
	checkFileExists(samplesheet_path,inputPath)
	write(paste0("samplesheet_path=",'"',samplesheet_path,'"'),config_file_fam,append=TRUE)
	samplesheet			<- read.table(samplesheet_path,sep=",",check.names=FALSE,stringsAsFactors=FALSE,header=TRUE)

	Chroms 				<- c(1:22,"X")
	write(paste0("Chroms=",'c(1:22,"X")'),config_file_fam,append=TRUE)

	#change to different values with WGS in Haplarithmissis_patscore2
	AvgMDA		<-2.2862648
	write("AvgMDA=2.2862648",config_file_fam,append=TRUE)
	SdMDA		<-0.4721314
	write("SdMDA=0.4721314",config_file_fam,append=TRUE)
	MinCallRate <-16.38
	write("MinCallRate=16.38",config_file_fam,append=TRUE)

	source(paste(script, "analyses", "functions","checkDirExistsAndCreate.R", sep="/"))
	outPath  <- checkDirExistsAndCreate(paste(inputPath,"ConvertGenotype",sep="/"))
	#add to paths to config_file
	write(paste0("dataPath=",'"',outPath,'"'),config_file_fam,append=TRUE)

	HelixoutPath  <- checkDirExistsAndCreate(paste(inputPath,"HelixData",sep="/"))
	write(paste0("HelixoutPath=",'"',HelixoutPath,'"'),config_file_fam,append=TRUE)

	GATK     <- paste(inputPath,"HaplotypeCaller",sep="/")
	write(paste0("GATK=",'"',GATK,'"'),config_file_fam,append=TRUE)

	fam_members				<- samplesheet[,c("SampleID","Sample_MetaInfo","Sample_Status")]
	fam_members_parents 	<- fam_members[fam_members[,"Sample_MetaInfo"]%in%c("Father","Mother"),]
	if (nrow(fam_members_parents)<2){
		write(paste0("MetaInfo gave the following error: single parent in samplesheet", e), errorFilePath, sep="")
		stop(paste0("MetaInfo gave the following error: single parent in samplesheet", e))
	}
	fam_members_nonEmbryo 	<- fam_members[!fam_members[,"Sample_Status"]%in%"E",]
	fam_members_Embryo 		<- fam_members[fam_members[,"Sample_Status"]%in%"E",]
	colID <- fam_members[,1]

	if (length(fam_members[grepl("AF",fam_members[,"Sample_Status"])|grepl("AM",fam_members[,"Sample_Status"]),"Sample_Status"])==2) {
		parents <- c("Father","Mother");
		write('parents=c("Father","Mother")',config_file_fam,append=TRUE)
	} else {
		if (fam_members[fam_members[,"Sample_MetaInfo"]%in%"Father","Sample_Status"]=="AF") {
			parents <- "Father";
			write('parents="Father"',config_file_fam,append=TRUE)
		}
		if (fam_members[fam_members[,"Sample_MetaInfo"]%in%"Mother","Sample_Status"]=="AM") {
			parents <- "Mother";
			write('parents="Mother"',config_file_fam,append=TRUE)
		}
	}

	write(paste0("Father=",'"',fam_members[fam_members[,"Sample_MetaInfo"]%in%"Father","SampleID"],'"'),config_file_fam,append=TRUE)
	write(paste0("Mother=",'"',fam_members[fam_members[,"Sample_MetaInfo"]%in%"Mother","SampleID"],'"'),config_file_fam,append=TRUE)
	source(paste(script, "analyses", "MetaInfo","functions","Parameters_extraction.R", sep="/"))
	REF <- Parameters_extraction(inputPath,fam_members)
	write(paste0("REF=",'"',REF,'"'),config_file_fam,append=TRUE)
	if (REF == "Grandparents") {
		RefID       <- fam_members[grepl("AGF", fam_members[, "Sample_Status"]) | grepl("AGM", fam_members[, "Sample_Status"]), "Sample_MetaInfo"]
		write(paste0("RefID=",'"',RefID,'"'),config_file_fam,append=TRUE)
		GrandFather <- fam_members[grepl("GF", fam_members[, "Sample_Status"]), "SampleID"]
		write(paste0("Grandfather=",'"',GrandFather,'"'),config_file_fam,append=TRUE)
		GrandMother <- fam_members[grepl("GM", fam_members[, "Sample_Status"]), "SampleID"]
		write(paste0("Grandmother=",'"',GrandMother,'"'),config_file_fam,append=TRUE)
	} else{
		RefID <- REF
		write(paste0("RefID=",'"',REF,'"'),config_file_fam,append=TRUE)
		write(paste0("RefSampleID=",'"',fam_members[fam_members[, "Sample_MetaInfo"]%in%REF, "SampleID"],'"'),config_file_fam,append=TRUE)
	}
	flip <- 0
	if (fam_members[grepl(RefID, fam_members[, "Sample_MetaInfo"]), "Sample_Status"] %in% "AGM" & REF != "Grandmother") { flip <- 1 }
	if (fam_members[grepl(RefID, fam_members[, "Sample_MetaInfo"]), "Sample_Status"] %in% "US") { flip <- 1 }
	if (fam_members[grepl(RefID, fam_members[, "Sample_MetaInfo"]), "Sample_Status"] %in% "UGM") { flip <- 1 }
	if (fam_members[grepl(RefID, fam_members[, "Sample_MetaInfo"]), "Sample_Status"] %in% "UGF") { flip <- 1 }
	write(paste0("flip=",flip),config_file_fam,append=TRUE)

	source(paste(script, "analyses", "MetaInfo","functions","Interval_extraction.R",sep="/"))
	Chrom <- Intervals_extraction(inputPath,samplesheet)
	write(paste0("Chrom=",'"',Chrom,'"'),config_file_fam,append=TRUE)

	source(paste(script, "analyses", "MetaInfo", "functions","inheritance_check.R",sep="/"))
	inheritance_check(fam_members,REF,Chrom)
}, error = function(e) {
	write(paste0("MetaInfo gave the following error: ", e), errorFilePath, sep="")
	stop(paste0("MetaInfo gave the following error: ", e))
})
