###########################################################################################################################
# Author: Anouk Janssen
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: liftover the hg19 coordinates to hg38 coordinates from onePGT (GBS) data

#delete all objects from memory
rm(list=ls(all=TRUE))

# set parameters
input_dir <- "path/to/input"
output_dir <- "path/to/output"
family <- "" # select the family folder to be liftover

# load packages
library("rtracklayer")
library(devtools)
library(GenomicRanges)
library(IRanges)

#Download the hg19 to hg38 liftover chain file from: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

# load the chain file for the hg19 to hg38 liftover
hg19ToHg38_chain <- import.chain("hg19ToHg38.over.chain")

# Family folder
filedir <- paste0(input_dir,"/",family,"/Haplarithmisis")

output <- paste0(output_dir,"/", family, "/hg38_Haplarithmisis")
if (!file.exists(output)){
  dir.create(output)
}

# iterate over all the files in family analyses folder
for (i in 1:length(list.files(filedir))){
inputfile <- list.files(filedir)[i]

data_hg19 <- read.table(inputfile, header = T, stringsAsFactors = F)

# Extract the chromosome, start position and end position from the position column
data_hg19$NewChrColumn <- paste0("chr", data_hg19$Chr)
data_hg19$old_hg19_bp <- data_hg19$Position
data_hg19$start <- as.numeric(sub("\\..*$", "", data_hg19$Position))
data_hg19$end <- data_hg19$start

#Convert the input data with hg19 coordinates into a GRanges object
input_GRanges <- makeGRangesFromDataFrame(
  data_hg19,
  seqnames.field = "NewChrColumn", # name of the column with chromosome information including "chr" prefix
  start.field = "start", # name of the column with start.field (start of position)
  end.field = "end", # name of the column with end.field (end of position)
  keep.extra.columns = TRUE) #strand information is unknown

# Call the liftOver function from the rtracklayer package to convert coordinates form hg19 to hg38.
output_GRanges <- as.data.frame(liftOver(input_GRanges, hg19ToHg38_chain))

write.table(output_GRanges, paste0(output,"/hg38_",inputfile), quote=F, sep = "\t", row.names=F)
}#end list.files loop
