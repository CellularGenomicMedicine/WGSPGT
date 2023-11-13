###########################################################################################################################
# Author: Anouk Janssen
# research group: Cellular Genomic Medicine, clinical genetics, Maastricht University Medical Centre + (MUMC+)

# purpose script: plotting the autosomal Mendelian inconsistency for the Validation study for GBS-PGT and WGS-PGT: Fig. 1f

#delete all objects from memory
rm(list=ls(all=TRUE))

#load packages
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(ggsignif)
library(cowplot)

# Set the directories
file.dir <- ""

out.dir <- ""

###########################################################################################################################
############################################## Prepare the data ###########################################################
###########################################################################################################################

# Set the path to the directory containing the Mendelian inconsistency files
data_directory <- ""

# Find all files with the pattern "MendInc.txt"
MendIncFiles <- list.files(data_directory, pattern = "MendInc.txt", recursive = T, full.names = T)

# Extract embryo names without path
embryo_names <- basename(MendIncFiles)

# Identify duplicated embryo names 
duplicated_names <- duplicated(embryo_names)

# Filter only the unique files
MendIncFiles <- MendIncFiles[!duplicated_names]

# create an empty data frame to store the combined data
MendIncPerChr <- data.frame()

# loop through each file and read its content, then bind it to 
for (i in MendIncFiles){
  # Read the file
  MendInc <- read.table(i, header = T, stringsAsFactors = F, check.names = F)
  
  # Bind the data to MendIncPerChr
  MendIncPerChr <- rbind(MendIncPerChr, MendInc)
}

###########################################################################################################################
############################################## if data is already combined ################################################


# load the tables with Mendelian inconsistency levels for WGS-PGT
MendIncPerChrWGSVal <- read.table("MendIncPerChrWGSVal.txt",sep = "\t", header=T, check.names=F)
MendIncPerChrWGSVal <- MendIncPerChrWGSVal %>%
  mutate(Method = "WGSVal") # add a column with the method WGS

# load the tables with Mendelian inconsistency levels for GBS-PGT
MendIncPerChrGBSVal <- read.table("MendIncPerChrGBSVal.txt",sep = "\t", header=T, check.names=F)
MendIncPerChrGBSVal <- MendIncPerChrGBSVal %>%
  mutate(Method = "GBSVal") # add a column with the method GBS

# Combine the data for GBS-PGT and WGS-PGT
MendIncPerChrVal <- rbind(MendIncPerChrWGSVal,MendIncPerChrGBSVal)        
MendIncPerChrVal$Method <- factor(MendIncPerChrVal$Method)

# Plot the autosomal Mendelian inconsistency levels for GBS-PGT and WGS-PGT
plot <- MendIncPerChrVal %>%
  ggplot() +
  aes(x = Method, y=chrAut) +
  geom_boxplot(data = filter(MendIncPerChrVal, Method == "WGSVal"), aes(fill = Method), lwd=0.4, width = 0.1, outlier.shape = NA) +
  geom_point(data = filter(MendIncPerChrVal, Method == "WGSVal"), aes(x = as.numeric(Method) + 0.2), size = 0.1, position = position_jitter(width = 0.1)) +
  geom_boxplot(data = filter(MendIncPerChrVal, Method == "GBSVal"), aes(fill = Method), lwd=0.4, width = 0.1, outlier.shape = NA) +
  geom_point(data = filter(MendIncPerChrVal, Method == "GBSVal"), aes(x = as.numeric(Method) - 0.2), size = 0.1, position = position_jitter(width = 0.1)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7),
    legend.position = "none") +
  labs(x = "Method", y = "Mendelian inconsistency (%)") +
  scale_fill_manual(values = c("GBSVal" = "darkgoldenrod3",
                               "WGSVal" = "deepskyblue")) +
  theme(legend.position = "none") +
  scale_y_continuous(expand=c(0,0), breaks = seq(0,20,5), limits = c(0,20)) +
  scale_x_discrete(labels = c("GBSVal" = "GBS", "WGSVal" = "WGS")) +
  geom_signif(comparisons = list(c("GBSVal", "WGSVal")),
              map_signif_level = T, 
              vjust = -1,
              y_position = 14,
              tip_length = 0.01,
              textsize = 1.5,
              size = 0.3)

ggsave(filename = file.path(out.dir, "Fig1f.png"),
       plot = plot, 
       width = 50, 
       height = 50, 
       dpi = 900, 
       units = "mm")

