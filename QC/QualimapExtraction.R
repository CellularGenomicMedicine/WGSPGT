###########################################################################################################################
# Author: Anouk Janssen
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: Extract the QC metrics from Qualimap 

# load libraries
library(dplyr)

# collect the data from all embryos and combine
samplesheet <- read.table("sampleSheet.txt", header = T, sep = "\t", stringsAsFactors=T) 

results <- data.frame()

for (i in 1:nrow(samplesheet)){
Folder = paste0(samplesheet$Folder[i])
Name = samplesheet$Name[i]
Files = list.files(Folder, recursive = T, full.names=T)
SampleFiles = Files[grep(pattern = Name, Files)]

# filter out ETS or MT DNA qualimap output and retrieve the output from the "genome_results" txt file
qualimapfile = SampleFiles[!grepl(pattern = "(ets|MT)", SampleFiles, ignore.case = TRUE) & grepl(pattern = "genome_results.txt", SampleFiles, ignore.case = TRUE)]

# Read the contents of the file
file_contents = readLines(qualimapfile)  

for (line in file_contents) {
  if (grepl ("mean coverageData", line)) {
    meanDepth <- as.numeric(gsub("mean coverageData =", "", gsub("X", "", line)))
  } 
  
  if (grepl ("number of reads", line)) {
  numReads <- as.numeric(gsub(",", "", gsub("number of reads =", "", line)))
  } 
  
  if (grepl ("number of mapped bases", line)) {
   mappedBases <- as.numeric(gsub(",", "", gsub("number of mapped bases =", "", gsub(" bp", "",line))))
  } 
  
  if (grepl ("coverageData >= 1X", line)){
    meanBreadth <- as.numeric(gsub("There is a ","", gsub("% of reference with a coverageData >= 1X","",line)))
  } 
  }


#load the other characteristics for this specific sample
PGT_number <- samplesheet_long$PGT_number[i]
Cycle <- samplesheet_long$Cycle[i]
Sample <- samplesheet_long$Sample[i]
SampleType <- samplesheet_long$SampleType[i]
StudyPart <- samplesheet_long$StudyPart[i]
Family <- samplesheet_long$Family[i]
EmbryoStatus <- samplesheet_long$EmbryoStatus[i]
PaperFamName <- samplesheet_long$PaperFamName[i]
Method <- samplesheet_long$Method[i]
HaplaOutput <- samplesheet_long$HaplaOutput[i]

#append the results to the dataframe
result <- data.frame(PGT_number = PGT_number,
                     Cycle = Cycle,
                     Sample = Sample,
                     SampleType = SampleType,
                     StudyPart = StudyPart,
                     Family = Family,
                     EmbryoStatus = EmbryoStatus,
                     PaperFamName = PaperFamName,
                     Method = Method,
                     Name = Name, 
                     Folder = Folder,
                     HaplaOutput = HaplaOutput,
                     meanDepth = meanDepth, 
                     numReads = numReads, 
                     mappedBases = mappedBases, 
                     meanBreadth = meanBreadth, 
                     stringsAsFactors = F)

results <- rbind(results, result)
}

write.table(results, "outDir/QualimapExtractionAllSamples.txt", sep = "\t")
  