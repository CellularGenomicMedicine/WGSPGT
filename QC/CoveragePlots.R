###########################################################################################################################
# author: Anouk Janssen
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: plotting the coverage QC metrics

## plotting the breadth of coverage of subsampled and validation samples: Fig. 1d 
## plotting the depth of coverage from validation samples: Fig. 1e
## plotting the depth of coverage from subsampled samples: Supp Fig. 1

# delete all objects from memory
rm(list=ls(all=TRUE))

options(scipen = 0, digits = 5)

# load packages
library(ggplot2)
library(dplyr)
library(ggsignif)
library(tidyr)
library(cowplot)

# load the geom_flat_violin function
source("https://raw.githubusercontent.com/rkoeck/saliva_methylome/main/scripts/rainCloudPlot.R")

# Set parameters reference genome length
OnePGTRef <- 3095694854
WGSPGTRef <- 3215267892

# Set the directories
file.dir <- ""

out.dir <- ""

# load the QC table
QualimapDataTable <- read.table("QualimapExtractionAllSamples.txt", stringsAsFactors = F)

#change the character values to numeric
QualimapDataTable[c("meanDepth", "numReads", "mappedBases", "meanBreadth")] <- lapply(QualimapDataTable[c("meanDepth", "numReads", "mappedBases", "meanBreadth")], as.numeric)

# calculate the mean depth of coverage for the covered positions
QualimapDataTable <- QualimapDataTable %>%
  mutate(meanDepthCovPositions = ifelse(Method =="OnePGT", (mappedBases / ((OnePGTRef*meanBreadth)/100)), NA)) %>%
  mutate(meanDepthCovPositions = ifelse(Method == "WGS", (mappedBases / ((WGSPGTRef*meanBreadth)/100)), meanDepthCovPositions))

# Separate the Pilot and validation data GBS and WGS data for breadth of coverage 
CombinedPilotValidation <- QualimapDataTable %>%
  mutate(DownSampled = ifelse(endsWith(Folder, "X"), "Yes", "No")) %>%
  mutate(TargetCoverage = ifelse(DownSampled == "Yes", sub(".*/", "", Folder), "FULL")) %>%
  mutate(Method = ifelse(Method == "OnePGT", "GBS", Method)) %>%
  mutate(TargetCoverage = ifelse(Method == "GBS", "GBS", TargetCoverage)) %>%
  mutate(TargetCoverage = ifelse(StudyPart == "Validation" & Method == "GBS", "GBSVal", TargetCoverage)) %>%
  mutate(TargetCoverage = ifelse(StudyPart == "Validation" & Method == "WGS", "WGSVal", TargetCoverage)) %>%
  mutate(TargetCoverage = factor(TargetCoverage, levels = c("GBS", "5X", "10X", "20X", "30X", "FULL", "GBSVal", "WGSVal")))

###########################################################################################################################
################################### Breadth of Coverage Plot ##############################################################
###########################################################################################################################

BreadthofCoveragePlot <- CombinedPilotValidation %>%
  mutate(TargetCoverage = factor(TargetCoverage, levels = rev(c("GBS", "5X", "10X", "20X", "30X", "FULL", "GBSVal", "WGSVal")))) %>%
  filter(TargetCoverage != "FULL") %>%
  ggplot() +
  aes(x = TargetCoverage, y = meanBreadth, fill = TargetCoverage) +
  geom_boxplot(lwd=0.4,outlier.color = "white", width = 0.8) +
  geom_point(position = position_jitter(0.3), size = 0.1, col = "grey18") +
  labs(x = "Target Coverage", y = "breadth of coverage (%)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7)) +
  coord_flip() + #flip coordinates
  guides(fill=F) + #remove the legend
  scale_fill_manual(values = c("GBS" = "darkgoldenrod3",
                               "5X" = "skyblue",
                               "10X" = "deepskyblue",
                               "20X" = "deepskyblue3",
                               "30X" = "deepskyblue4",
                               "WGSVal" = "deepskyblue",
                               "GBSVal" = "darkgoldenrod3"
  )) +
  scale_y_continuous(breaks = seq(0,100,10), 
                     limits = c(0,100), 
                     expand = c(0,0)) +
  scale_x_discrete(expand = c(0.1,0.1), labels = c("GBS" = "GBS", 
                                                   "5X" = "5X", 
                                                   "10X" = "10X", 
                                                   "20X" = "20X", 
                                                   "30X" = "30X", 
                                                   "GBSVal" = "GBS", 
                                                   "WGSVal" = "WGS")) +
  geom_signif(comparisons = list(c("GBS", "5X"), c("5X", "10X"), c("10X", "20X"), c("20X", "30X"), c("10X", "30X")),
              map_signif_level = T, 
              vjust = 3,
              y_position = c(10,76,78,80,72),
              tip_length = -0.01,
              textsize = 1.5,
              size = 0.3)

# Save the plot
ggsave(filename = file.path(out.dir,"Fig1d.png"), 
       plot = BreadthofCoveragePlot, 
       width = 110, 
       height = 50, 
       dpi = 600, 
       units = "mm")

###########################################################################################################################
################################### Depth of Coverage Plot Validation #####################################################
###########################################################################################################################

Validation <- QualimapDataTable %>%
  filter(StudyPart == "Validation") %>% # Select data from the validation study
  mutate(Method = ifelse(Method == "OnePGT", "GBS", Method)) %>%
  mutate(MDAvsBULK = ifelse(SampleType == "Embryo", "MDA", "BULK")) # create a new column indicating MDA-based 

Validation$Method <- factor(Validation$Method, levels = c("GBS", "WGS"))
Validation$MDAvsBULK <- factor(Validation$MDAvsBULK, levels = c("MDA", "BULK"))

# Compare the GBS and WGS 
ggplot(Validation) +
  aes(x = Method, y = meanDepthCovPositions) +
  geom_boxplot() +
  geom_signif(test = "wilcox.test",
              comparisons = list(c("GBS", "WGS")),
              map_signif_level = T, 
              y_position = 17,
              tip_length = 0.01,
              textsize = 2.5,
              vjust = -2.5)

# mean depth at covered positions for the validation study separated by method (GBS or WGS) and input sample (MDA or bulk)
plot <- ggplot(Validation) +
  aes(x = MDAvsBULK, y = meanDepthCovPositions) +
  geom_flat_violin(data = filter(Validation, MDAvsBULK == "MDA"), aes(fill=Method, col = Method), width = 0.7) +
  geom_boxplot(data = filter(Validation, MDAvsBULK == "MDA"), lwd=0.4, fill = "white", width = 0.1, outlier.color = "NA") +
  geom_point(data = filter(Validation, MDAvsBULK == "MDA"), aes( x = as.numeric(MDAvsBULK) + 0.2, col = Method), position = position_jitter(width = 0.1), size = 0.1, alpha = 0.5) +
  geom_flat_violin(data = filter(Validation, MDAvsBULK == "BULK"), aes(fill = Method, col = Method), width = -0.7) +
  geom_boxplot(data = filter(Validation, MDAvsBULK == "BULK"), lwd=0.4, fill = "white", width = -0.1, outlier.color = "NA") +
  geom_point(data = filter(Validation, MDAvsBULK == "BULK"), aes( x = as.numeric(MDAvsBULK) - 0.2, col = Method), position = position_jitter(width = 0.1), size = 0.1, alpha = 0.5) +
  
  # Add geom_signif layer for significance testing between MDA and BULK
  geom_signif(test = "wilcox.test",
              comparisons = list(c("MDA", "BULK")),
              map_signif_level = T, 
              y_position = 4,
              tip_length = -0.01,
              textsize = 2,
              vjust = 2.5,
              size = 0.25) +
    
  labs(y = "depth of coverage (X)", fill = "Method", title = NULL) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 7),
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(0, "lines")) +
  guides(fill=F, col = F) + #remove the legend
  scale_fill_manual(values = c("darkgoldenrod3","deepskyblue")) +
  scale_color_manual(values = c("darkgoldenrod3","deepskyblue")) +
  scale_y_continuous(breaks = seq(0,20,5), limits = c(0,20), expand = c(0,0)) +
  coord_flip() +
  facet_grid(Method ~ ., scales = "free_y")
  
#Save the plot
ggsave(filename = file.path(out.dir,"Fig1e.png"), 
       plot = plot, 
       width = 110, 
       height = 40, 
       dpi = 600, 
       units = "mm")
# add error bar between WGS and GBS according to boxplot with comparison WGS and GBS

###########################################################################################################################
################################### Depth of Coverage Plot Pilot ##########################################################
###########################################################################################################################

# select the pilot data
Pilot <- QualimapDataTable %>%
  filter(StudyPart == "Pilot") %>%
  filter(PGT_number %in% c("PGD4621", "PGD4267")) %>%
  #Create a new column indicating if the file path ends with "X"
  mutate(DownSampled = ifelse(endsWith(Folder, "X"), "Yes", "No")) %>%
  #Extract the last part of the folder filepaths and assign it to the "TargetCoverage"
  mutate(TargetCoverage = ifelse(DownSampled == "Yes", sub(".*/", "", Folder), "FULL")) %>%
  mutate(Method = ifelse(Method == "OnePGT", "GBS", Method)) %>%
  mutate(TargetCoverage = ifelse(Method == "GBS", "GBS", TargetCoverage)) %>%
  mutate(TargetCoverage = factor(TargetCoverage, levels = c("GBS", "5X", "10X", "20X", "30X", "FULL")))

# plot the depthofcoverage
DepthPilot <- ggplot(Pilot) +
  aes(x = TargetCoverage, y = meanDepth, fill = TargetCoverage) +
  geom_boxplot(lwd=0.4,outlier.color = "white") +
  geom_point(position = position_jitter(0.3), size = 0.5, col = "grey18") +
  labs(x = "Target Coverage", y = "mean depth of coverage (X)", fill = "Target Coverage (X)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7)) +
  guides(fill=F) + #remove the legend
  scale_fill_manual(values = c("darkgoldenrod3","skyblue","deepskyblue","deepskyblue2","deepskyblue3","deepskyblue4")) +
  scale_y_continuous(breaks = seq(0,50,10), limits = c(0,56), expand = c(0,0)) +
  geom_signif(test = "wilcox.test",
              comparisons = list(c("GBS", "5X"), c("5X", "10X"), c("10X", "20X"), c("20X", "30X"), c("30X","FULL")),
              map_signif_level = T, 
              y_position = c(46.25,47,47.75,48.5,49.25,50),
              tip_length = 0.01,
              textsize = 2.5)

meanDepthCovPositionsPilot <- ggplot(Pilot) +
  aes(x = TargetCoverage, y = meanDepthCovPositions, fill = TargetCoverage) +
  geom_boxplot(lwd=0.3,outlier.color = "white") +
  geom_point(position = position_jitter(0.3), size = 0.3, col = "grey18") +
  labs(x = "Target Coverage", y = "mean depth of coverage (X) in covered positions", fill = "Target Coverage (X)") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7)) +
  guides(fill=F) + #remove the legend
  scale_fill_manual(values = c("darkgoldenrod3","skyblue","deepskyblue","deepskyblue2","deepskyblue3","deepskyblue4")) +
  scale_y_continuous(breaks = seq(0,50,10), limits = c(0,56), expand = c(0,0)) +
  geom_signif(test = "wilcox.test",
              comparisons = list(c("GBS", "5X"), c("5X", "10X"), c("10X", "20X"), c("20X", "30X"), c("30X","FULL")),
              map_signif_level = T, 
              y_position = c(46.25,47,47.75,48.5,49.25,50),
              tip_length = 0.01,
              textsize = 2.5)

combined_plot <- plot_grid(
  meanDepthCovPositionsPilot + ggtitle("a") + theme(plot.title = element_text(hjust = -0.05, size = 7)),
  DepthPilot + ggtitle("b") + theme(plot.title = element_text(hjust = -0.05, size = 7)),
  ncol = 2 # Display the plots in a single column
)

ggsave(filename = file.path(out.dir,"SuppFig1.png"), 
              plot = combined_plot, 
              width = 120, 
              height = 60, 
              dpi = 600, 
              units = "mm")
  
