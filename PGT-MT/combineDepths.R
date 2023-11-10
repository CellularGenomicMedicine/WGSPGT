###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: combine individual sample depth per position files

# input: depth per position (1x .txt file per sample - output of samtoolsDepth.sh), sample sheet (for file IDs)

# output:  multi-sapmple depth per mitochondrial position file (.csv)

###########################################################################################################################

# loadPackages
library(dplyr)

# collect the data from all embryos & combine

toRun = data.table::fread("sampleSheet.csv", stringsAsFactors = F) 

for(i in 1:nrow(toRun)){
  
  sample = toRun[i, "sampleID"]
  
  setwd(paste0(toRun[i, "dataDir"], "/", sample))
  
  setwd(list.files()[1])
  
  data = data.table::fread(paste0("MT_depth/", sample, "_MT.txt"))
  
  colnames(data) = gsub(".bam", "", colnames(data))
  
  if(i == 1){
    mtCoverage = data
  }else{
    
    mtCoverage = full_join(mtCoverage, data, by = c("#CHROM", "POS"))
  }
  
}

# write out the data
write.csv(mtCoverage, "outDir/depthPerPos.csv", row.names = F)
