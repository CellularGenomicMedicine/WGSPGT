###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: haplarithm plotting using ggplot

# input: output from haplarithmisis data processing

# output: haplarithm plot (.jpg)

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

# set the chromosome and position of interest (can include multiple positions for comlex heterozygous indications)
affectedChrs = data.frame(Chr = 8, pos = 71215510)

# set embryo parameters
EmbryoID = "sampleName"
flip = 0

opt1 <- 1
opt2 <- 2
if (flip == 1) {
  opt1 <- 2
  opt2 <- 1
}

# set plotting parameters
alphaBAF = 0.1
sizeBAF = 0.02
sizeAvg = 0.1
alphaParentalBAF = 0.2
linewidth = 0.2

margins = unit(c(0,1,0,0), units = "mm")

theme_xaxis = theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 5),
        axis.text.y = element_text(size = 5),
        text = element_text(size = 5),
        plot.margin = margins)

removey = theme(axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.line.y = element_blank(),
                text = element_text(size = 5),
                plot.margin = margins)

## load / set the chromsome variables
load("Ideogram_GRCh38.rda")
ChrsLength$Chromosome = gsub("chr", "", ChrsLength$Chromosome)

###############################################################################
###################### load the haplarithmisis data ###########################
###############################################################################

setwd("dataDir")

fileNameBase = "genericName"

BAFs = data.table::fread(paste0(fileNameBase, "BAF.txt")) %>% as.data.frame()

P1 = data.table::fread(paste0(fileNameBase, "P1.txt")) %>% as.data.frame()
P1Seg = data.table::fread(paste0(fileNameBase, "P1Seg.txt")) %>% as.data.frame()
P2 = data.table::fread(paste0(fileNameBase, "P2.txt")) %>% as.data.frame()
P2Seg = data.table::fread(paste0(fileNameBase, "P2Seg.txt")) %>% as.data.frame()

M1 = data.table::fread(paste0(fileNameBase, "M1.txt")) %>% as.data.frame()
M1Seg = data.table::fread(paste0(fileNameBase, "M1Seg.txt")) %>% as.data.frame()
M2 = data.table::fread(paste0(fileNameBase, "M2.txt")) %>% as.data.frame()
M2Seg = data.table::fread(paste0(fileNameBase, "M2Seg.txt")) %>% as.data.frame()

Hap = data.table::fread(paste0(fileNameBase, "dataHap.txt")) %>% as.data.frame()
HapRaw = data.table::fread(paste0(fileNameBase, "dataHapRaw.txt")) %>% as.data.frame()

logRs = data.table::fread(paste0(fileNameBase, "_LogR.txt") %>% as.data.frame()
segLogRs = data.table::fread(paste0(fileNameBase, "_gammaSC50_segLogRs.txt") %>% as.data.frame()

##################################################################################################################
################################### plot the relevant haplarithmisis segments ####################################
##################################################################################################################
                             
for(chr in unique(affectedChrs$Chr)){
  
  length = ChrsLength %>% filter(Chromosome == chr) %>% .$Length
  minx = ifelse(affectedChrs$pos[1] - 2000000 > 0, affectedChrs$pos[1] - 2000000, 0)
  maxx = affectedChrs$pos[1] + 2000000
  
  # calculate the axis labels
  ticks = c(minx, affectedChrs$pos[1], maxx)
  labels = c("", round(affectedChrs$pos[1] / 1000000, digits = 2), "")
  
  
  ## B allele frequency plot
  BAFs_chr = BAFs %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  
  BAFplot = ggplot(BAFs_chr) +
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
    geom_vline(xintercept = affectedChrs$pos, colour = "orange", size = linewidth * 2) +
    geom_point(aes(x = Position, y = BAFs_chr[ , EmbryoID]), size = sizeBAF, alpha = alphaBAF) +
    ylab("BAF") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
    ggtitle(paste0("Chr", chr)) +
    theme_xaxis 
  
  ## Paternal Haplotype
  dataHap_chr = Hap %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  
  pat_hap = ggplot() + 
    geom_linerange(data = dataHap_chr[dataHap_chr[, paste0(EmbryoID, "_Pat")] == opt1, ],
                   aes(x = Position, ymin = 0, ymax = 2), size = 0.1, colour = "blue") +
    geom_linerange(data = dataHap_chr[dataHap_chr[, paste0(EmbryoID, "_Pat")] == opt2, ],
                       aes(x = Position, ymin = 0, ymax = 2), size = 0.3, colour = "cornflowerblue") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
    theme_void()
  
  
  ## Paternal BAFs
  
  P1_chr = P1 %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  P1Seg_chr = P1Seg %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  P2_chr = P2 %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  P2Seg_chr = P2Seg %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  
  paternal_BAF <- ggplot() +
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
    geom_vline(xintercept = affectedChrs$pos, colour = "orange", size = linewidth * 2) +
    geom_point(aes(x = P1_chr$Position, y = P1_chr[ , EmbryoID]), colour = "blue1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = P2_chr$Position, y = P2_chr[ , EmbryoID]), colour = "red1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = P1Seg_chr$Position, y = P1Seg_chr[ , EmbryoID]), colour = "blue4", alpha = 1, shape = "square", size = sizeAvg) +
    geom_point(aes(x = P2Seg_chr$Position, y = P2Seg_chr[ , EmbryoID]), colour = "red4", alpha = 1, shape = "square", size = sizeAvg) +
    ylab("Pat-BAF") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.005,0.005)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
    theme_xaxis
  
  ## Maternal Haplotype
  
  mat_hap = ggplot() + 
     geom_linerange(data = dataHap_chr[dataHap_chr[, paste0(EmbryoID, "_Mat")] == opt1, ],
                     aes(x = Position, ymin = 0, ymax = 2), size = 0.1, colour = "red") +
    geom_linerange(data = dataHap_chr[dataHap_chr[, paste0(EmbryoID, "_Mat")] == opt2, ],
                   aes(x = Position, ymin = 0, ymax = 2), size = 0.1, colour = "pink") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
    theme_void()
  
  # Maternal BAFs
  
  M1_chr =M1 %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  M1Seg_chr =M1Seg %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  M2_chr =M2 %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  M2Seg_chr =M2Seg %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  
  maternal_BAF <- ggplot() +
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
    geom_vline(xintercept = affectedChrs$pos, colour = "orange", size = linewidth * 2) +
    geom_point(aes(x = M1_chr$Position, y = M1_chr[ , EmbryoID]), colour = "blue1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = M2_chr$Position, y = M2_chr[ , EmbryoID]), colour = "red1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = M1Seg_chr$Position, y = M1Seg_chr[ , EmbryoID]), colour = "blue4", alpha = 1, shape = "square", size = sizeAvg) +
    geom_point(aes(x = M2Seg_chr$Position, y = M2Seg_chr[ , EmbryoID]), colour = "red4", alpha = 1, shape = "square", size = sizeAvg) +
    ylab("Mat-BAF") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.005,0.005)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
    theme_xaxis
  
  ## Embryo logR plot (bottom plot)
  
  logRs_chr = logRs %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  segLogRs_chr = segLogRs %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  
  logRplots <- ggplot() + 
    geom_hline(yintercept = c(-1.5, -1, -0.5, 0.5, 1, 1.5), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
    geom_vline(xintercept = affectedChrs$pos, colour = "orange", size = linewidth * 2) +
    geom_point(aes(x = logRs_chr$Position, y = logRs_chr[ , EmbryoID]), alpha = alphaParentalBAF*20, size = sizeBAF) +
    geom_point(aes(x = segLogRs_chr$Position, y = segLogRs_chr[ , EmbryoID]), colour = "darkgoldenrod2", shape = "square", size = sizeAvg) +
    ylab("logR") +
    theme_classic() + 
    theme(axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          text = element_text(size = 5),
          plot.margin = margins) +
    scale_y_continuous(limits = c(-1.6, 1.6), breaks = c(-1, 0, 1), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx), breaks = ticks, labels = labels) #+
    #theme_xaxis
  
  # combine the plots and save
  
    plotwidth = 0.75 
  
  jpeg(paste0(getwd(), "/", EmbryoID, "_chr", chr, "_haplarithm.jpg"),
       width = plotwidth, height = 2.5, units = "in", res = 600)
  
  print(plot_grid(BAFplot + removey,
                  NULL,
                  pat_hap,
                  NULL,
                  paternal_BAF + removey,
                  NULL,
                 mat_hap,
                  NULL,
                  maternal_BAF + removey,
                  NULL,
                  logRplots + removey,
                  ncol = 1,
                  rel_heights = c(5, 0.3, 
                                 0.5, 0.1,
                                  4, 0.3,
                                 0.5, 0.1,
                                  4, 0.3, 
                                  5.5),
                  align = "v"))
  
  dev.off()
  
}

                             
# save the corresponding y axis plot as well
jpeg(paste0(getwd(), "/", EmbryoID, "_chr", chr, "_haplarithm.jpg"),
       width = plotwidth, height = 2.5, units = "in", res = 600)

print(plot_grid(BAFplot + removey,
                NULL,
                #pat_hap_raw,
                #NULL,
                pat_hap,
                NULL,
                paternal_BAF + removey,
                NULL,
                #mat_hap_raw,
                #NULL,
                mat_hap,
                NULL,
                maternal_BAF + removey,
                NULL,
                logRplots + removey,
                ncol = 1,
                rel_heights = c(5.5, 0.3, 
                                #0.5, 0.1,
                                0.5, 0.1,
                                4, 0.3,
                                #0.5, 0.1,
                                0.5, 0.1,
                                4, 0.3, 
                                5.5),
                align = "v"))

dev.off()
