###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualise haplarithms, including chromsome schematics and breakpoint information

# input: as per haplarithmisis plots, filtered break point information

# output: haplarithm plot per affected chromosome & corresponding y-axis (3x .jpg)

###########################################################################################################################

# load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

## set the relevant parameters

# chromosomes of intrest and break points (change as relevant)
affectedChrs = data.frame(Chr = c(8, 16), breakpoint = c(22092231, 89257105))

# sample parameters
EmbryoID = "sampleName"
flip = 0

opt1 <- 1
opt2 <- 2
if (flip == 1) {
  opt1 <- 2
  opt2 <- 1
}

# plotting parameters
linetype = "solid"
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

#####################################################################
################ load the haplarithmisis data ######################
#####################################################################

setwd("dataDir")

fileNameBase = "genericFileName"

BAFs = data.table::fread(paste0(fileNameBase, "BAF.txt")) %>% as.data.frame()

P1 = data.table::fread(paste0(fileNameBase, "P1.txt")) %>% as.data.frame()
P1Seg = data.table::fread(paste0(fileNameBase, "P1Seg.txt")) %>% as.data.frame()
P2 = data.table::fread(paste0(fileNameBase, "P2.txt")) %>% as.data.frame()
P2Seg = data.table::fread(paste0(fileNameBase, "P2Seg.txt")) %>% as.data.frame()

M1 = data.table::fread(paste0(fileNameBase, "M1.txt")) %>% as.data.frame()
M1Seg = data.table::fread(paste0(fileNameBase, "M1Seg.txt")) %>% as.data.frame()
M2 = data.table::fread(paste0(fileNameBase, "M2.txt")) %>% as.data.frame()
M2Seg = data.table::fread(paste0(fileNameBase, "M2Seg.txt")) %>% as.data.frame()

dataHapRaw = data.table::fread(paste0(fileNameBase, "dataHapRaw.txt")) %>% as.data.frame()
dataHap = data.table::fread(paste0(fileNameBase, "dataHap.txt")) %>% as.data.frame()

logRs = data.table::fread(paste0(fileNameBase, "_LogR.txt")) %>% as.data.frame()
segLogRs = data.table::fread(paste0(fileNameBase, "gammaSC50_SegLogR.txt")) %>% as.data.frame()

#####################################################################
################ load the plotting files ############################
#####################################################################

# the chromosome length information
load("Ideogram_GRCh38.rda")
ChrsLength$Chromosome = gsub("chr", "", ChrsLength$Chromosome)

# the chromosome schematics coordinates
outlines = read.csv("2mbTaperOutlines.csv") 
blanks = read.csv("2mbTaperBlanks.csv") 

outlines_fam = read.csv("extraChromosomeOutlines.csv")
blanks_fam = read.csv("extraChromosomeBlanks.csv")

#####################################################################
################ generate ideogram plots ############################
#####################################################################

ideograms = list()

# normal chr8

outlines8 = outlines %>% filter(chr == 8)
blanks8 = blanks %>% filter(chr == 8)

normal8 = ggplot() +
  geom_rect(aes(xmin = 0, xmax = max(outlines8$x), ymin = 0, ymax = 1), fill = "#FF9900") +
  geom_polygon(data = blanks8, aes(x = x, y = y, group = arm), fill = "white", color = NA) +
  geom_polygon(data = outlines8, aes(x = x, y= y, group = arm), fill = NA, colour = "black", size = 0.1) +
  theme_void() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(0, max(outlines8$x)))

# der8 #segment chr16 = 1081240
der8_outlines = outlines_fam %>% filter(chr == "der_8")
der8_outlines$x = gsub(124127645, 145138636, der8_outlines$x)
der8_outlines$x = gsub(122127645, 143138636, der8_outlines$x)
der8_outlines$x = as.numeric(der8_outlines$x)

der8_blanks = blanks_fam %>% filter(chr == "der_8")
der8_blanks$x = gsub(124127645, 145138636, der8_blanks$x)
der8_blanks$x = gsub(122127645, 143138636, der8_blanks$x)
der8_blanks$x = as.numeric(der8_blanks$x)

der8_fill = data.frame(xmin = c(21010991, 22092231), xmax = c(22092231, 145138636), origin = c("16", "8"))

der8 = ggplot() +
  geom_rect(data = der8_fill, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = origin, alpha = origin)) +
  scale_fill_manual(values = c("16" = "#99CCFF", "8" = "#FF9900")) +
  scale_alpha_manual(values = c("16" = 0.5, "8" = 1)) +
  geom_polygon(data = der8_blanks, aes(x = x, y = y, group = arm), fill = "white", color = NA) +
  geom_polygon(data = der8_outlines, aes(x = x, y= y, group = arm), fill = NA, colour = "black", size = 0.1) +
  geom_line(aes(x = 22092231, y = c(0.2,0.8)), size = 0.2, colour = "grey30", linetype = "dashed") +
  theme_void() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(0, max(der8_outlines$x)))

# chr16
outlines16 = outlines %>% filter(chr == 16)
blanks16 = blanks %>% filter(chr == 16)

normal16 = ggplot() +
  geom_rect(aes(xmin = 0, xmax = max(outlines16$x), ymin = 0, ymax = 1), fill = "#99CCFF") +
  geom_polygon(data = blanks16, aes(x = x, y = y, group = arm), fill = "white", color = NA) +
  geom_polygon(data = outlines16, aes(x = x, y= y, group = arm), fill = NA, colour = "black", size = 0.1) +
  theme_void() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(0, max(outlines16$x))) 

# extra ch16
extra16_outlines = outlines_fam %>% filter(chr == 16)
extra16_blanks = blanks_fam %>% filter(chr == 16) 

extra16 = ggplot() +
  geom_rect(aes(xmin = 89257105, xmax = max(extra16_outlines$x), ymin = 0, ymax = 1), fill = "#99CCFF") +
  geom_polygon(data = extra16_blanks, aes(x = x, y = y, group = arm), fill = "white", color = NA) +
  geom_polygon(data = extra16_outlines, aes(x = x, y= y, group = arm), fill = NA, colour = "black", size = 0.1) +
  theme_void() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(0, max(outlines16$x)))

# combine the ideogram parts
ideograms[["chr8"]] = plot_grid(NULL, NULL,
                                der8, NULL,
                                normal8, 
                                ncol = 1,
                                align = "v",
                                rel_heights = c(1, 0.05,
                                                1, 0.05,
                                                1))

ideograms[["chr16"]] = plot_grid(extra16, NULL,
                                 normal16, NULL,
                                 normal16, 
                                 ncol = 1,
                                 align = "v",
                                 rel_heights = c(1, 0.05,
                                                 1, 0.05,
                                                 1))

#####################################################################
########## plot the haplarithmisis segments #########################
#####################################################################

for(chr in unique(affectedChrs$Chr)){
  
  round = which(affectedChrs$Chr == chr)
  
  maxChr = ChrsLength %>% filter(Chromosome == chr) %>% .$Length
  
  breakpoint = affectedChrs %>% filter(Chr == chr) %>% .$breakpoint
  
  minx = 0
  maxx = maxChr
  
  breakpoint = affectedChrs %>% filter(Chr == chr) %>% .$breakpoint
  
  
  # calculate the axis labels
  ticks = c(minx, breakpoint, maxx)
  labels = c("", round(breakpoint / 1000000, digits = 1), "")
  
  ## B allele frequency plot
  BAFs_chr = BAFs %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  
  BAFplot = ggplot(BAFs_chr) +
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
    geom_vline(xintercept = breakpoint, colour = "chartreuse4", size = 0.3, linetype = linetype) +
    geom_point(aes(x = Position, y = BAFs_chr[ , EmbryoID]), size = sizeBAF, alpha = alphaBAF) +
    ylab("BAF") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
 #   ggtitle(paste0("Chr", chr)) +
    theme_xaxis 
  
  ## Paternal Haplotype
  
  dataHap_chr = dataHap %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  dataHapRaw_chr = dataHapRaw %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  
  pat_hap = ggplot() + 
    geom_linerange(data = dataHap_chr[dataHap_chr[, paste0(EmbryoID, "_Pat")] == opt1, ],
                   aes(x = Position, ymin = 0, ymax = 2), size = 0.01, colour = "blue") +
    geom_linerange(data = dataHap_chr[dataHap_chr[, paste0(EmbryoID, "_Pat")] == opt2, ],
                   aes(x = Position, ymin = 0, ymax = 2), size = 0.01, colour = "cornflowerblue") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
    theme_void()
  
  ## Paternal BAFs
  
  P1_chr = P1 %>% filter(Chr %in% chr, Position >= minx, Position <= maxx)
  P1Seg_chr = P1Seg %>% filter(Chr %in% chr, Position >= minx, Position <= maxx)
  P2_chr = P2 %>% filter(Chr %in% chr, Position >= minx, Position <= maxx)
  P2Seg_chr = P2Seg %>% filter(Chr %in% chr, Position >= minx, Position <= maxx)
  
  paternal_BAF <- ggplot() +
    geom_hline(yintercept = c(0.25, 0.5, 0.75, 1), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
    geom_vline(xintercept = breakpoint, colour = "chartreuse4", linetype = linetype, size = 0.3) +
    geom_point(aes(x = P1_chr$Position, y = P1_chr[ , EmbryoID]), colour = "blue1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = P2_chr$Position, y = P2_chr[ , EmbryoID]), colour = "red1", alpha = alphaParentalBAF, size = sizeBAF) +
    geom_point(aes(x = P1Seg_chr$Position, y = P1Seg_chr[ , EmbryoID]), colour = "blue4", alpha = 1, shape = "square", size = sizeAvg) +
    geom_point(aes(x = P2Seg_chr$Position, y = P2Seg_chr[ , EmbryoID]), colour = "red4", alpha = 1, shape = "square", size = sizeAvg) +
    ylab("Pat-BAF") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0.005,0.005)) +
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx)) +
    theme_xaxis
  
  ## Embryo logR plot (bottom plot)
  
  logRs_chr = logRs %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  segLogRs_chr = segLogRs %>% filter(Chr == chr, Position >= minx, Position <= maxx)
  
  logRplots <- ggplot() + 
    geom_hline(yintercept = c(-1.5, -1, -0.5, 0.5, 1, 1.5), linetype = "dashed", colour = "gray70", size = linewidth) +
    geom_hline(yintercept = 0, colour = "gray30") +
    geom_vline(xintercept = breakpoint, colour = "chartreuse4", linetype = linetype, size = 0.3) +
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
    scale_x_continuous(expand = c(0,0), limits = c(minx, maxx), breaks = ticks, labels = labels)
  
#####################################################################
########## combine the plots & save output #########################
#####################################################################
  
  plotwidth = (maxx / 235476981) * 1.5
  
  jpeg(paste0("outDir", EmbryoID, "_chr", chr, "breakpointHaplarithm.pdf"),
      width = plotwidth, height = 2.4, units = "in", res = 600) 
  
  print(plot_grid(ideograms[[paste0("chr", chr)]] + removey,
                  NULL,
                  BAFplot + removey,
                  NULL,
                  pat_hap + removey,
                  NULL,
                  paternal_BAF + removey,
                  NULL,
                  logRplots + removey,
                  ncol = 1,
                  rel_heights = c(1.5, 0.5,
                                  4, 0.3, 
                                  0.5, 0.2, 
                                  4, 0.3, 
                                  5.5),
                  align = "v"))
  
  dev.off()
  
}

# save the corresponding y axis plot separately

pdf(paste0("outDir", EmbryoID, "_yaxis", ".pdf"),
    width = 0.47, height = 2) 

print(plot_grid(BAFplot,
                NULL,
                pat_hap,
                NULL,
                paternal_BAF,
                NULL,
                logRplots,
                ncol = 1,
                rel_heights = c(5.5, 0.3, 0.5, 0.2, 4, 0.3, 5.5),
                align = "v"))

dev.off()

