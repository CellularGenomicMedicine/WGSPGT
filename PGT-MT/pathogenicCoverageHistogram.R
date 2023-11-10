###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: combine individual sample depth per position files

# input: multi-sapmple depth per mitochondrial position file (.csv), sample sheet

# output: pathogenic position coverage histogram (.jpg)
###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

### load the data

depth = data.table::fread("dataDir/depthPerPos.csv.csv",
                          header = T, stringsAsFactors = F) %>%
  as.data.frame() %>%
  select(chr, POS)

### load the sample sheet

sampleSheet = data.table::fread("sampleSheet.csv", stringsAsFactors = F) 

### load the pathological positions

pathMITO = read.csv("pathMITO.csv")

### calculate the mean depth per condition
# conditions: validation (PGT-M / PGT-SR embryos), sub10x (deep-sequencing and subsampling to 10x), seq10x (direct 10x sequencing)
# the code code chunk below can be updated with the rowMins function from the MatrixStats package to calculate the minimum
# coverage of any sample at a given site.

validation = sampleSheet %>% filter(condition == "validation") %>% select(sampleID)
sub40x = sampleSheet %>% filter(condition == "sub40x") %>% select(sampleID)
sub30x = sampleSheet %>% filter(condition == "sub30x") %>% select(sampleID)
sub20x = sampleSheet %>% filter(condition == "sub20x") %>% select(sampleID)
sub10x = sampleSheet %>% filter(condition == "sub10x") %>% select(sampleID)
sub5x = sampleSheet %>% filter(condition == "sub5x") %>% select(sampleID)
seq10x = sampleSheet %>% filter(condition == "seq10x") %>% select(sampleID)

depth$mean_validation = rowMeans(depth %>% select(all_of(validation)))
depth$mean_sub40x = rowMeans(depth %>% select(all_of(sub40x)))
depth$mean_sub30x = rowMeans(depth %>% select(all_of(sub30x)))
depth$mean_sub20x = rowMeans(depth %>% select(all_of(sub20x)))
depth$mean_sub10x = rowMeans(depth %>% select(all_of(sub10x)))
depth$mean_sub5x = rowMeans(depth %>% select(all_of(sub5x)))
depth$mean_seq10x = rowMeans(depth %>% select(all_of(seq10x)))

### filter the positions to only contain the pathological positions
indication = depth %>% filter(POS %in% pathMITO$Position)

### select only the relevant columns
indication = indication %>% select(POS, all_of((contains("mean")))

### Re-format the data for plotting
depth.longer = indication %>% pivot_longer(cols = c(2:8),
                                    names_to = "conditionID",
                                    values_to = "depth") %>%
  as.data.frame()

# order the conditions for plotting
depth.longer$conditionID = gsub("mean_", "", depth.longer$conditionID)      
                                   
depth.longer$conditionID = factor(depth.longer$conditionID, 
                               levels = c("sub40x", "sub30x", "sub20x", "sub10x", "sub5x", "seq10x", "validation"))
                                   
# calculate log10 of the depth                                  
depth.longer = depth.longer %>% mutate(depthLog10 = log10(depth)) 

# set the plot labels                                   
labels.coverage = c("full PGT-MT", "30X PGT-MT", "20X PGT-MT", "10X PGT-MT", "5X PGT-MT", "10X family 21", "10X PGT-M/PGT-SR")
names(labels.coverage) = c("sub40x", "sub30x", "sub20x", "sub10x", "sub5x", "seq10x", "validation")

# make the histogram
path = depth.longer %>% ggplot(aes(x = depthLog10)) +
  geom_histogram(bins = 50) +
  geom_hline(yintercept = 0, colour = "black", size = 0.5) +
  theme_classic() +
  facet_wrap(~conditionID, ncol = 1, strip.position = "right", labeller = labeller(conditionID = labels.coverage)) +
  scale_x_continuous(breaks = log10(c(1500, 3000, 5000, 10000)), labels = c(1500, 3000, 5000, 10000)) +
  scale_y_continuous(expand = c(0, 0.5), breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme(axis.line.x = element_blank(),
        axis.text = element_text(size = 5),
        strip.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        panel.grid.major.y  = element_line(colour = "grey80", size = 0.2, linetype = "dashed")) +
  xlab("Reads per position (log10)")

# save the plot
ggsave(path, file = "outDIR/covgHistogram_pathMITO105.jpg",
       width = 100, height = 150, device = "jpg", units = "mm", dpi = 900)


