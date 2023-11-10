###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: calculate the number of informative SNPs per sample ber 1000bp bin

# input: parentally phased informative SNP positions (from haplarithmisis) 

# output:  number of informative SNPs per 1000bp bin

###########################################################################################################################

#load libraries
library(data.table)
library(dplyr)

# load the bins object
bin_object <- data.table::fread("bins100000.txt", stringsAsFactors = F) %>%
  as.data.frame()

# load the sample sheet
samplesheet <- read.table("sampleSheet.csv", sep = ";", header=T,
                          stringsAsFactors = F)

for (s in 1:nrow(samplesheet)){

  # set the sample specific parameters needed for file loading
  EmbryoID <- samplesheet$sampleName[s]
  FamilyID <- samplesheet$Family[s])
  
  WGSorOnePGT <- "WGS"

  Gamma_value = 50

  setwd(samplesheet$dataDir)
  files = list.files()
  
  setwd(dataDir)

  M1file = grep(paste0(mBAFs, "_", EmbryoID, "_", FamilyID, "_", Gamma_value, "_M1.txt"), files, value=T)

   M2file = grep(paste0(mBAFs, "_", EmbryoID, "_", FamilyID, "_", Gamma_value, "_M2.txt"), files, value=T)

   P1file = grep(paste0(pBAFs, "_", EmbryoID, "_", FamilyID, "_", Gamma_value, "_P1.txt"), files, value=T)

   P2file = grep(paste0(pBAFs, "_", EmbryoID, "_", FamilyID, "_", Gamma_value, "_P2.txt"), files, value=T)

  # read all the informative SNP files per embryo
  M1 <- fread(M1file, header=T, sep="\t",data.table=F, check.names=F) %>%
        filter(.[ , 4] >= 0)

  M2 <- fread(M2file, header=T, sep="\t",data.table=F, check.names=F) %>%
        filter(.[ , 4] >= 0)

  P1 <- fread(P1file, header=T, sep="\t",data.table=F, check.names=F) %>%
        filter(.[ , 4] >= 0)

  P2 <- fread(P2file, header=T, sep="\t",data.table=F, check.names=F) %>%
        filter(.[ , 4] >= 0)

  # combine the objects
  InfSnPs <- rbind(M1,M2,P1,P2)

  #bin_object for loop
  for (i in 1:nrow(bin_object)){

    # select bin locations
    chr = bin_object$chr[i]
    start = bin_object$start[i]
    end = bin_object$end[i]

    # calculate number of informative SNPs in the bin
    count = InfSnPs %>% filter(Chr == chr, Position >= start, Position <= end) %>% nrow(.)

    # add the value to the data frame
    bin_object[i, EmbryoID] = count

      }

}

# write the full object
write.table(bin_object,
            file = "InfSNPsPerBin100000Bin_WGS_phasedOnly.txt",
            row.names = F, quote = F, sep = "\t")

## this script can easily be adjusted to calculate the informative SNPs per phased parent if desired.


