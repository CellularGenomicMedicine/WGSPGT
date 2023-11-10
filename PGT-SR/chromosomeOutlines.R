###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: create basic chromosome coordinates for plotting

# input: chromosome lengths

# output: basic chromosome coordinates for plotting

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))

# define the positions for the chromosome outlines

load("Ideogram_GRCh38.rda")

centros = ideogram %>% filter(Color == "acen")

centromeres = data.frame(chromosome = NA, xmin = NA, xmax = NA,
                         stringsAsFactors = F)

for(chr in unique(centros$Chromosome)){
  
  dat = centros %>% filter(Chromosome == chr)
  
  xmin = min(dat$Start)
  
  xmax = max(dat$Stop)
  
  chromosome = gsub("chr", "", chr)
  
  row = c(chromosome, xmin, xmax)
  
  centromeres = rbind(centromeres, row)
  
}

centromeres = na.omit(centromeres)

centromeres = centromeres %>% filter(chromosome != "Y")
centromeres$chromosome = factor(centromeres$chromosome, levels = c("X", seq(1, by = 1, length.out = 22)))
centromeres = centromeres %>% arrange(chromosome)

centromeres_gbs = centromeres %>% mutate(ymin = seq(0.9, by = 3.5, length.out = 23),
                                     ymax = seq(1.1, by = 3.5, length.out = 23),
                                     method = "gbs")

centromeres_wgs = centromeres %>% mutate(ymin = seq(2.4, by = 3.5, length.out = 23),
                                         ymax = seq(2.6, by = 3.5, length.out = 23),
                                         method = "wgs")


parms = data.frame(chromosome = NA, xmin = NA, xmax = NA,
                   stringsAsFactors = F)

for(chr in unique(centromeres$chromosome)){
  
  chromosome = chr
  
  xmin = 0
  
  xmax = centromeres %>% filter(chromosome == chr) %>% .$xmin
  
  row = c(chromosome, xmin, xmax)
  
  parms = rbind(parms, row)
}

parms = na.omit(parms)

parms = parms %>% filter(chromosome != "Y")
parms$chromosome = factor(parms$chromosome, levels = c("X", seq(1, by = 1, length.out = 22)))
parms = parms %>% arrange(chromosome)

parms_gbs = parms %>% mutate(ymin = seq(0.5, by = 3.5, length.out = 23),
                         ymax = seq(1.5, by = 3.5, length.out = 23),
                         method = "gbs")

parms_wgs = parms %>% mutate(ymin = seq(2, by = 3.5, length.out = 23),
                             ymax = seq(3, by = 3.5, length.out = 23),
                             method = "wgs")



qarms = data.frame(chromosome = NA, xmin = NA, xmax = NA,
                   stringsAsFactors = F)

for(chr in unique(centromeres$chromosome)){
  
  chromosome = chr
  
  xmin = centromeres %>% filter(chromosome == chr) %>% .$xmax
  
  xmax = ChrsLength %>% filter(Chromosome == paste0("chr", chr)) %>% .$Length
  
  row = c(chromosome, xmin, xmax)
  
  qarms = rbind(qarms, row)
}

qarms = na.omit(qarms)

qarms = qarms %>% filter(chromosome != "Y")
qarms$chromosome = factor(qarms$chromosome, levels = c("X", seq(1, by = 1, length.out = 22)))
qarms = qarms %>% arrange(chromosome)

qarms_gbs = qarms %>% mutate(ymin = seq(0.5, by = 3.5, length.out = 23),
                         ymax = seq(1.5, by = 3.5, length.out = 23),
                         method = "gbs")
qarms_wgs = qarms %>% mutate(ymin = seq(2, by = 3.5, length.out = 23),
                         ymax = seq(3, by = 3.5, length.out = 23),
                         method = "wgs")

outlines = rbind(centromeres_gbs, centromeres_wgs, parms_gbs, parms_wgs, qarms_gbs, qarms_wgs)

outlines = outlines %>% mutate(xmin = as.numeric(xmin), xmax = as.numeric(xmax))

write.csv(outlines, "chromosomeOutlines.csv", row.names = F)
