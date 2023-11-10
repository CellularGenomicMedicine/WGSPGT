###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: haplarithm plotting with chromosome ideogram (for PGT-AO cases)

# input: haplarithm output and ideogram

# output: haplarithm plot

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

# chromosomes of intrest and break points
affectedChrs = c(2)

# Sample parameters
EmbryoID = "sampleName"

#### load the haplarithmisis data ###

setwd("dataDir")

fileNameBase = "genericName"

BAFs = data.table::fread(paste0(fileNameBase, "BAF.txt")) %>% as.data.frame()

P1 = data.table::fread(paste0(fileNameBase, "P1.txt")) %>% as.data.frame()
P1Seg = data.table::fread(paste0(fileNameBase, "P1Seg.txt")) %>% as.data.frame()
P2 = data.table::fread(paste0(fileNameBase, "P2.txt")) %>% as.data.frame()
P2Seg = data.table::fread(paste0(fileNameBase, "P2Seg.txt")) %>% as.data.frame()

P1_self = data.table::fread(paste0(fileNameBase, "P1.txt")) %>% as.data.frame()
P1Seg_self = data.table::fread(paste0(fileNameBase, "P1Seg.txt")) %>% as.data.frame()
P2_self = data.table::fread(paste0(fileNameBase, "P2.txt")) %>% as.data.frame()
P2Seg_self = data.table::fread(paste0( fileNameBase, "P2Seg.txt")) %>% as.data.frame()

M1_self = data.table::fread(paste0(fileNameBase, "M1.txt")) %>% as.data.frame()
M1Seg_self = data.table::fread(paste0(fileNameBase, "M1Seg.txt")) %>% as.data.frame()
M2_self = data.table::fread(paste0(fileNameBase, "M2.txt")) %>% as.data.frame()
M2Seg_self = data.table::fread(paste0(ileNameBase, "M2Seg.txt")) %>% as.data.frame()

logRs = data.table::fread(paste0(ileNameBase, "_LogR.txt") %>% as.data.frame()
segLogRs = data.table::fread(paste0(ileNameBase, "_gammaSC50_segLogRs.txt") %>% as.data.frame()

breaks = read.csv(paste0(ileNameBase, "_limitsP2Seg.csv")

################################### set the plotting options ####################################################

## load / set the chromsome variables
load("Ideogram_GRCh38.rda")
ChrsLength$Chromosome = gsub("chr", "", ChrsLength$Chromosome)
                  
alphaBAF = 0.005
sizeBAF = 0.02
sizeAvg = 0.1
alphaParentalBAF = 0.01
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

################################## make the chromosome ideograms ###############################################

source(ideogram.R)

################################### plot the relevant haplarithmisis segments ####################################

for(chr in unique(affectedChrs)){
  
  maxx = ChrsLength %>% filter(Chromosome == chr) %>% .$Length
  
  minx = 0
  
  # calculate the axis labels
  ticks = seq(0, maxx, maxx / 4)
  labels = round(seq(0, maxx / 1000000, maxx / 4000000), digits = 2)
  
  ## B allele frequency plot
  BAFs_chr = BAFs %>% filter(Chr == chr)
  
  BAFplot = ggplot(BAFs_chr) +
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
    geom_point(aes(x = Position, y = BAFs_chr[ , EmbryoID]), size = sizeBAF, alpha = alphaBAF) +
    ylab("BAF") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
 #   ggtitle(paste0("Chr", chr)) +
    theme_xaxis 
  
  ## Paternal BAFs
  
  P1_chr = P1 %>% filter(Chr == chr)
  P1Seg_chr = P1Seg %>% filter(Chr == chr)
  P2_chr = P2 %>% filter(Chr == chr)
  P2Seg_chr = P2Seg %>% filter(Chr == chr)
  
  paternal_BAF <- ggplot() +
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
    geom_point(aes(x = P1_chr$Position, y = P1_chr[ , EmbryoID]), colour = "blue1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = P2_chr$Position, y = P2_chr[ , EmbryoID]), colour = "red1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = P1Seg_chr$Position, y = P1Seg_chr[ , EmbryoID]), colour = "blue4", alpha = 1, shape = "square", size = sizeAvg) +
    geom_point(aes(x = P2Seg_chr$Position, y = P2Seg_chr[ , EmbryoID]), colour = "red4", alpha = 1, shape = "square", size = sizeAvg) +
    ylab("Pat-BAF") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.005,0.005)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
    theme_xaxis
  
  ## Paternal BAFs self
  
  P1_chr_self = P1_self %>% filter(Chr == chr)
  P1Seg_chr_self = P1Seg_self %>% filter(Chr == chr)
  P2_chr_self = P2_self %>% filter(Chr == chr)
  P2Seg_chr_self = P2Seg_self %>% filter(Chr == chr)
  
  paternal_BAF_self <- ggplot() +
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
    geom_point(aes(x = P1_chr_self$Position, y = P1_chr_self[ , EmbryoID]), colour = "blue1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = P2_chr_self$Position, y = P2_chr_self[ , EmbryoID]), colour = "red1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = P1Seg_chr_self$Position, y = P1Seg_chr_self[ , EmbryoID]), colour = "blue4", alpha = 1, shape = "square", size = sizeAvg) +
    geom_point(aes(x = P2Seg_chr_self$Position, y = P2Seg_chr_self[ , EmbryoID]), colour = "red4", alpha = 1, shape = "square", size = sizeAvg) +
    ylab("Pat-BAF-emb") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.005,0.005)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
    theme_xaxis
  
  # Maternal BAFs self
  
  M1_chr_self =M1_self %>% filter(Chr == chr)
  M1Seg_chr_self =M1Seg_self %>% filter(Chr == chr)
  M2_chr_self =M2_self %>% filter(Chr == chr)
  M2Seg_chr_self =M2Seg_self %>% filter(Chr == chr)
  
  maternal_BAF_self <- ggplot() +
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
    geom_point(aes(x = M1_chr_self$Position, y = M1_chr_self[ , EmbryoID]), colour = "blue1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = M2_chr_self$Position, y = M2_chr_self[ , EmbryoID]), colour = "red1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = M1Seg_chr_self$Position, y = M1Seg_chr_self[ , EmbryoID]), colour = "blue4", alpha = 1, shape = "square", size = sizeAvg) +
    geom_point(aes(x = M2Seg_chr_self$Position, y = M2Seg_chr_self[ , EmbryoID]), colour = "red4", alpha = 1, shape = "square", size = sizeAvg) +
    ylab("Mat-BAF-emb") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.005,0.005)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
    theme_xaxis
  
  ## Embryo logR plot (bottom plot)
  
  logRs_chr = logRs %>% filter(Chr == chr)
  segLogRs_chr = segLogRs %>% filter(Chr == chr)
  
  logRplots <- ggplot() + 
    geom_hline(yintercept = c(-1.5, -1, -0.5, 0.5, 1, 1.5), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
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
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx), breaks = ticks, labels = labels) +
    theme_xaxis
  
  # combine and save the plots
  plotwidth = 1.2

  jpeg(paste0("outDir", EmbryoID, "_chr", chr, ".jpg"),
       width = plotwidth, height = 3.8, units = "in", res = 600)
  
  print(plot_grid(maternal1 + removey,
                  NULL,
                  paternal1 + removey,
                  NULL,
                  paternal1 + removey,
                  NULL,
                  BAFplot + removey,
                  NULL,
                  paternal_BAF + removey,
                  NULL,
                  logRplots + removey,
                  NULL,
                  paternal_BAF_self + removey,
                  NULL,
                  maternal_BAF_self + removey,
                  ncol = 1,
                  rel_heights = c(0.5, 0.5, #ideogram1
                                  0.5, 0.1, #ideogram 2
                                  0.5, 1, #ideogram3
                                  5, 0.3, #BAFs
                                  4, 0.3, #matBAFs
                                  5, 1.5, # logRs
                                  4, 0.3, #matBAFs-self
                                  4), #pat-BAFs-self
                  align = "v"))
  
  dev.off()
  
}


# save the corresponding y axis plot 

jpeg(paste0("/outDir/", EmbryoID, "_yaxis.jpg"),
     width = 2, height = 3.8, units = "in", res = 600)

print(plot_grid(maternal1,
                NULL,
                paternal1,
                NULL,
                paternal1,
                NULL,
                BAFplot,
                NULL,
                paternal_BAF,
                NULL,
                logRplots,
                NULL,
                paternal_BAF_self,
                NULL,
                maternal_BAF_self,
                ncol = 1,
                rel_heights = c(0.5, 0.5, #ideogram1
                                0.5, 0.1, #ideogram 2
                                0.5, 1, #ideogram3
                                5, 0.3, #BAFs
                                4, 0.3, #matBAFs
                                5, 1.5, # logRs
                                4, 0.3, #matBAFs-self
                                4), #pat-BAFs-self
                align = "v"))

dev.off()
