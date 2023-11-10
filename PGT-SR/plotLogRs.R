###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualise (segmented) logRs from subsampled sequencing data of PGT-SR samples

# input: (segmented) logR values

# output: logR plot (.jpeg)

###########################################################################################################################

# load packages & set options
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot))
suppressPackageStartupMessages(library(tibble))
options(scipen=999)

setwd("dataDir")

chrs = c(8, 16)

sampleID = "sampleName"

# load, filter and combine data

logR = data.table::fread("full/sampleID.txt") %>%
  filter(Chr %in% chrs)

full = data.table:: fread(paste0("full/", "sampleID_gammaSC50_SegLogR.txt")) %>%
  filter(Chr %in% chrs) %>% 
  select(Names, SampleID) %>% 
  rename("sampleID.full" = sampleID)

X10 = data.table:: fread(paste0("10X/", "sampleID_gammaSC50_SegLogR.txt")) %>%
  filter(Chr %in% chrs) %>% 
  select(Names, SampleID) %>% 
  rename("sampleID.10" = sampleID)

X20 = data.table:: fread(paste0("20X/", "sampleID_gammaSC50_SegLogR.txt")) %>%
  filter(Chr %in% chrs) %>% 
  select(Names, SampleID) %>% 
  rename("sampleID.20" = sampleID)

X30 = data.table:: fread(paste0("30X/", "sampleID_gammaSC50_SegLogR.txt")) %>%
  filter(Chr %in% chrs) %>% 
  select(Names, SampleID) %>% 
  rename("sampleID.30" = sampleID)

X5 = data.table:: fread(paste0("5X/", "sampleID_gammaSC50_SegLogR.txt")) %>%
  filter(Chr %in% chrs) %>%
  select(Names, SampleID) %>% 
  rename("sampleID.5" = sampleID)

data = list(logR, full, X10, X20, X30, X5) %>% reduce(full_join, by = "Names")

## re-format the data for plotting
toplot = pivot_longer(data, cols = colnames(data)[grep("sampleID.", colnames(data))],
                            names_to = "depth",
                            values_to = "SegLogR") %>%
  as.data.frame() %>%
  select(Names, Position, Chr, sampleID, depth, SegLogR) %>%
  rename("logR" = sampleID) %>%
  mutate(embryo = "sampleID", identifier = paste0(embryo, "_", Chr))

## adjust the variables for plotting

toPlot$depth = gsub("sampleID.", "", toPlot$depth)
toPlot$depth = factor(toPlot$depth, levels = c ("full", 30, 20, 10, 5))
toPlot$identifier = factor(toPlot$identifier, levels = c("sampleID_8", "sampleID_16"))


### plot 
labels.embryo = c("Embryo 1: chr8", "Embryo1: chr16")
names(labels.embryo) = c("sampleID_8", "sampleID_16")
labels.coverage = c("Full", "30X", "20X", "10X", "5X")
names(labels.coverage) = c("full", 30, 20, 10, 5)

adjust_x = function(x){x/1000000}

logRs =  toPlot %>% ggplot() +
  geom_hline(yintercept = c(-1.5, -1, -0.5, 0.5, 1, 1.5), linetype = "dashed", colour = "gray30", size = 0.1) +
  geom_hline(yintercept = 0, colour = "gray30") +
  geom_point(aes(x = Position, y = logR), size = 0.01, alpha = 0.5) +
  geom_point(aes(x = Position, y = SegLogR), colour = "red", shape = "square", size = 0.1) +
  facet_grid(rows = vars(depth), cols = vars(identifier), scales = "free_x", space = "free_x",
             labeller = labeller(depth = labels.coverage, identifier = labels.embryo)) +
  scale_y_continuous(limits = c(-1.6, 1.6), breaks = c(-1, -0.5, 0, 0.5, 1), expand = c(0, 0)) +
  scale_x_continuous(labels = adjust_x) +
  xlab("Position (Mb)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5))



ggsave(logRs, filename = "logRsPlot.jpg", device = "jpg",
       width = 180, height = 180, units = "mm", dpi = 900)
