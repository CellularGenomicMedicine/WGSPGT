###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: heteroplasmy level calculation

# input: mitochondrial variant information (.vcf) - output from gatk.sh

# output:  heteroplasmy level for 1 sample (.csv)

###########################################################################################################################

## load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(VariantAnnotation))

## set the indication site
indication = 3243

## load the sample data
sampleID = "nameOfSample"

vcf = readVcf(paste0(dataDir, "/", sampleID, ".g.vcf"))

## extract the annotation information from the vcf file
ranges = as.data.frame(rowRanges(vcf)) %>%
    filter(start == indication)

# create the results object
heteroplasmy = data.frame(ref = ranges$REF,
                         alt1 = unlist(ranges$ALT)[1],
                         alt2 = unlist(ranges$ALT)[2],
                         alt3 = unlist(ranges$ALT)[3])

# extract the genotype
GT = geno(vcf)$GT %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column("position") %>%
          filter(grepl(indication, position))

heteroplasmy = heteroplasmy %>% mutate(gt = GT[ , sampleID])

# extract the depth per allele
AD = geno(vcf)$AD %>% 
          as.data.frame() %>% 
          tibble::rownames_to_column("position") %>%
          filter(grepl(indication, position))

heteroplasmy = heteroplasmy %>% 
          mutate(allele0 = unlist(AD[ , sampleID])[1],
                 allele1 = unlist(AD[ , sampleID])[2],
                 allele3 = unlist(AD[ , sampleID])[3],
                 allele4 = unlist(AD[ , sampleID])[4])

# extract the total depth at the position
DP = geno(vcf)$DP %>%
                as.data.frame() %>% 
                tibble::rownames_to_column("position") %>%
                filter(grepl(indication, position))
  
heteroplasmy = heteroplasmy %>% 
                mutate(total = DP[ , sampleID)])

# calculate the heteroplasmy percentage
heteroplasmy = heteroplasmy %>% 
                mutate(heteroplasmy = allele1 / (allele0 + allele1) * 100)

# write the output
write.csv(heteroplasmy, 
          file = paste0("outDir/", sampleID, "_heteroplasmy.csv", 
          row.names = F)

                       

