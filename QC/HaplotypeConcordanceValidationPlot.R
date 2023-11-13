# WGS-PGT 
# author: Anouk Janssen
# research group: Cellular Genomic Medicine, clnical genetics, Maastricht University Medical Centre + (MUMC+)

# purpose script: plotting the haplotype concordance for maternal and paternal haplotypes 

# Fig. 1g

#delete all objects from memory
rm(list=ls(all=TRUE))

#load packages
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

# Set the directories
file.dir <- ""

out.dir <- ""

# load the QC table
ConcordanceTable <- file.path(file.dir, "ConcordanceTable.txt")

concordance <- read.table(ConcordanceTable, sep = "\t", stringsAsFactors = F)
concordance_longer <- concordance %>%
  pivot_longer(cols = ends_with("_concordance"),
               names_to = "Parent",
               values_to = "Concordance") %>%
  mutate(Parent = str_replace(Parent, "_concordance", ""))
concordance_longer$Parent <- factor(concordance_longer$Parent)

plot <- ggplot(concordance_longer, aes(x = Parent, y = as.numeric(Concordance), fill = Parent)) +
  geom_boxplot(data = filter(concordance_longer, Parent == "maternal"), lwd = 0.4, width = 0.1, outlier.colour = NA) +
  geom_point(data = filter(concordance_longer, Parent == "maternal"), aes(x = as.numeric(Parent) - 0.2, y =  as.numeric(Concordance)), size = 0.1, position = position_jitter(width = 0.1)) +
  geom_boxplot(data = filter(concordance_longer, Parent == "paternal"), lwd = 0.4, width = 0.1, outlier.colour = NA) +
  geom_point(data = filter(concordance_longer, Parent == "paternal"), aes(x = as.numeric(Parent) + 0.2, y =  as.numeric(Concordance)), size = 0.1, position = position_jitter(width = 0.1)) +
  scale_fill_manual(values = c(paternal = "cornflowerblue",
                               maternal = "#ffc0cb")) +
  scale_y_continuous(expand=c(0,0), breaks = seq(80,100,5), limits = c(80,100)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7),
    legend.position = "none") +
  labs(x = "Haplotype", y = "haplotype concordance (%)") 

ggsave(filename = file.path(out.dir, "Fig1g.png"),
       plot = plot, 
       width = 50, 
       height = 40, 
       dpi = 900, 
       units = "mm")

