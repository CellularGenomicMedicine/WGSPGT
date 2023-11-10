###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: generate an ideogram for triploid chromosomes depicting homolgous recombination sites as well

# input: homologous recombination site information from haplarithmisis

# output: chromosome ideogram

###########################################################################################################################

# load packages

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))


###################################### load / set the data ###############################################################

# chromosomes of intrest and break points
breaks = read.csv("segMentLimits.csv")

# annotate / adjust the breaks object
breaks$colour = c("cornflowerblue", "blue", "cornflowerblue")

# set chromosome
c = 2

################################### load the chromosome information ####################################################

outlines = read.csv("2mbTaperOutlines.csv") %>%
  filter(chr == c) 

blanks = read.csv("2mbTaperBlanks.csv") %>%
  filter(chr == c) 

################################# adjust the coordinate system according to the breaks ################################

## plot

paternal1 = ggplot() +
  geom_rect(data = breaks, aes(xmin = Start, xmax = end, ymin = 0, ymax = 1, fill = colour)) +
  geom_polygon(data = blanks, aes(x = x, y = y, group = arm), fill = "white") +
  geom_polygon(data = outlines, aes(x = x, y= y, group = arm), fill = NA, colour = "black") +
  scale_fill_manual(values = c("blue" = "blue", "cornflowerblue" = "lightblue")) +
  theme_void() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(0, max(outlines$x))) 

maternal1 = ggplot() +
  geom_rect(aes(xmin = 0, xmax = max(outlines$x), ymin = 0, ymax = 1), fill = "lightgrey") +
  geom_polygon(data = blanks, aes(x = x, y = y, group = arm), fill = "white") +
  geom_polygon(data = outlines, aes(x = x, y= y, group = arm), fill = NA, colour = "black") +
  theme_void() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0,0), limits = c(0, max(outlines$x))) 

## paste the plots together

plot_grid(paternal1,
          NULL,
          paternal1,
          NULL,
          maternal1,
          align = "v",
          ncol = 1,
          rel_heights = c(1, 0.05,
                          1, 0.3,
                          1))
