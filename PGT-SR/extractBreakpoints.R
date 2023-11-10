###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: extract relevant breakpoint information from the manta output

# input: structural variant information (.vcf) - output from manta.sh

# output: filtered structural variants per sample (.csv)

###########################################################################################################################

# load packages
library(dplyr)
library(VariantAnnotation)
library(tibble)

# load the data
sampleID = "sampleName"

vcf = readVcf("sampleDataDir/diploidSV.vcf")

# set the chromosomes of interest
firstChr = "chr2"
secondChr = "chr11"

# extract and combine the relevant data from the vcf file
pos = rowRanges(vcf) %>% as.data.frame() %>% dplyr::select(-REF, -ALT)
var = info(vcf) %>% as.data.frame() %>% dplyr::select(-HOMSEQ, - LEFT_SVINSSEQ, -RIGHT_SVINSSEQ)

data = cbind(pos, var)
data = data %>% rownames_to_column("ID")

# filter based on the chromosomes expected to be involved
# (retain only breakpoints with matching MATEIDs)

data_firtChr = data %>% filter(seqnames == firstChr, SVTYPE == "BND")

data_secondChr = data %>% filter(seqnames == secondChr, MATEID %in% data_firstChr$ID)

data_firstChr = data_firstChr %>% filter(MATEID %in% data_secondChr$ID)

results = rbind(data_firstChr, data_secondChr) %>% mutate(sample = sampleID)

# select relevant columns
results = results %>% dplyr::select(ID, seqnames, start, end, QUAL, FILTER, IMPRECISE, SVTYPE, MATEID, JUNCTION_QUAL, sample)

# substitute problematic characters
results$FILTER = gsub(";", ".", results$FILTER)

# save the output
write.csv(results, "outDIR/sampleID_breakpoints.csv", row.names = F)

## Quality filters could also be applied at theis point using the filter, qual and imprecise columns. These were not applied in this
## project to demonstrate the unfiltered returns and the effect of possible filtering conditions.









