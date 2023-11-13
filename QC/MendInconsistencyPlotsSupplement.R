# WGS-PGT 
# author: Anouk Janssen
# research group: Cellular Genomic Medicine, clinical genetics, Maastricht University Medical Centre + (MUMC+)

# purpose script: plotting the Mendelian inconsistency rates for the WGS and GBS methods for the pilot study (a) and the validation study per chromosome (b)

#delete all objects from memory
rm(list=ls(all=TRUE))

#load packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(ggsignif)
library(cowplot)

# set the directories
file.dir <- ""

out.dir <- ""

# load the data
MendIncPerChrPilot <- read.table("MendelianInPerChrPilot.txt", sep = "\t", header=T, check.names=F)
MendIncPerChrPilot[MendIncPerChrPilot[,"Target_coverage"]=="default","Target_coverage"] = "GBS"

# adjust the column names 
colnames(MendIncPerChrPilot) <- c("Family", 
                                  "Indication", 
                                  "Target_coverage", 
                                  "Embryo",
                                  "StudyPart",
                                  paste0("chr",rep(1:22)),
                                  "chrX",
                                  "chrAut")

#select the relevant columns
MendIncPerChrPilot <- MendIncPerChrPilot %>%
  select(c("Target_coverage","Embryo",paste0("chr",rep(1:22)),"chrX","chrAut"))

MendIncPerChrWGSVal <- read.table("MendIncPerChrWGSVal.txt",sep = "\t", header=T, check.names=F)
MendIncPerChrWGSVal <- MendIncPerChrWGSVal %>%
  mutate(Target_coverage = "WGSVal")

MendIncPerChrGBSVal <- read.table("MendIncPerChrGBSVal.txt",sep = "\t", header=T, check.names=F)
MendIncPerChrGBSVal <- MendIncPerChrGBSVal %>%
  mutate(Target_coverage = "GBSVal")

MendIncPerChrVal <- rbind(MendIncPerChrWGSVal,MendIncPerChrGBSVal)        
MendIncPerChrVal <- tibble::rownames_to_column(MendIncPerChrVal, var = "Embryo")         

MendIncPerChrAll <- rbind(MendIncPerChrPilot, MendIncPerChrVal)   

# To create the Violin Plot the data needs to be transformed to long format
MendIncPerChrAll <- MendIncPerChrAll %>%
  pivot_longer(cols = c( paste0("chr",rep(1:22)), "chrX", "chrAut"), names_to = "chr", values_to = "MendInc")

# Define the desired order of target_coverage levels
cov_order <- (c("GBS", "5", "10", "20", "30", "FULL", "GBSVal", "WGSVal"))
MendIncPerChrAll$Target_coverage <- factor(MendIncPerChrAll$Target_coverage, levels = cov_order)
# Define the desired order of chromosomes
chr_order <- paste0("chr", c(1:22, "X", "Aut"))
MendIncPerChrAll$chr <- factor(MendIncPerChrAll$chr, levels = chr_order)

# Supplemental Fig 2a: plot the mendelian inconsistency for the GBS data and corresponding pilot data per target coverage
MendIncPlotPilotVertical <- MendIncPerChrAll %>%
  filter(!Target_coverage %in% c("WGSVal", "GBSVal")) %>% # filter out the validation data
  filter(Target_coverage != "FULL",
         chr == "chrAut") %>% # look only at the GenomAut
  group_by(Target_coverage) %>%
  ggplot(aes(x = Target_coverage, y=MendInc)) +
  geom_boxplot(aes(fill = Target_coverage), lwd=0.4, width = 0.2, outlier.shape = NA) +
  geom_point(size = 0.1, position = position_jitter(width = 0.1)) +
  
  # Add geom_signif layer for significance testing between target coverages
  geom_signif(comparisons = list(c("GBS", "5"), c("5", "10"), c("10", "20"), c("20", "30"), c("10","30")),
                map_signif_level = T, 
                vjust = 0,
                y_position = c(13.5,13,6.5,4.5,8),
                tip_length = 0.01,
                textsize = 2,
                size = 0.25) +
  
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7)) +
  labs(x = "Target Coverage", y = "Mendelian inconsistency (%)") +
  scale_fill_manual(values = c("GBS" = "darkgoldenrod3",
                               "5" = "skyblue",
                               "10" = "deepskyblue",
                               "20" = "deepskyblue3",
                               "30" = "deepskyblue4")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 7)) +
  scale_y_continuous(expand=c(0,0), breaks = seq(0,20,5), limits = c(0,20)) +
  scale_x_discrete(labels = c("GBS" = "GBS", "5" = "5X", "10" = "10X", "20" = "20X", "30" = "30X"))

ggsave(filename = file.path(out.dir, "PublishPlots/SuppFig2a.png"),
       plot = MendIncPlotPilotVertical,
       width = 60,
       height = 60,
       dpi = 600,
       units = "mm")


# Fig 2b
#make a chr specific plot
MendIncPerChrPlotValidation <- MendIncPerChrAll %>%
  filter(Target_coverage%in% c("WGSVal", "GBSVal")) %>%
  filter(!(chr %in% c("chrX","chrAut"))) %>%
  group_by(chr) %>%
  ggplot() +
  aes(x = chr, y=MendInc) +
  #geom_flat_violin(aes(fill=Target_coverage,col=Target_coverage)) +
  geom_boxplot(lwd=0.4, width = 0.4, outlier.size = 0.1, aes(fill=Target_coverage, col =Target_coverage), position = position_dodge(width = 0.5)) +
  geom_boxplot(lwd=0.4, width = 0.4, outlier.size = 0.1, aes(fill=Target_coverage), position = position_dodge(width = 0.5), outlier.colour = NA) +
  theme_classic() +
  theme(legend.position="none") +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7)) +
  labs(x = "Chromosome", y = "Mendelian Inconsistency (%)") +
  scale_fill_manual(values = c("GBSVal" = "darkgoldenrod3",
                               "WGSVal" = "deepskyblue")) +
  scale_colour_manual(values = c("GBSVal" = "darkgoldenrod3",
                                 "WGSVal" = "deepskyblue")) +
  scale_x_discrete(labels = c(1:22))

ggsave(filename = file.path(out.dir, "PublishPlots/SuppFig2b.png"),
       plot = MendIncPerChrPlotValidation,
       width = 120,
       height = 60,
       dpi = 600,
       units = "mm")



# combine the plots
plot1 <- MendIncPlotPilotVertical +
  labs(tag = "a") +
  theme(plot.tag = element_text(face = "bold", size = 8))

plot2 <- MendIncPerChrPlotValidation +
  labs(tag = "b") +
  theme(plot.tag = element_text(face = "bold", size = 8))


combined_plot <- plot_grid(
  plot1,
  plot2,
  ncol = 2,
  rel_widths = c(1,2)
)

ggsave(filename = file.path(out.dir, "SuppFig2.png"),
plot = combined_plot,
width = 160,
height = 60,
dpi = 600,
units = "mm")

