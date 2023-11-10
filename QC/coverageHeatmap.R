###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualise average Informative SNP distribution across the genome

# input: informative SNPs per 1000bp bins

# output: chromosome heatmap plot

###########################################################################################################################

## load packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(matrixStats))

## load the chromosome outlines & blanks

outlines = read.csv("chromosomeOutlines.csv",
                    stringsAsFactors = F) %>%
                    filter(chromosome != "Y")

blanks = read.csv("centromereBlanks.csv",
                  stringsAsFactors = F) %>%
                  filter(chromosome != "Y")

## load the data

#WGS
data1 = data.table::fread("InfSNPsPerBin100000Bin_WGS.txt", stringsAsFactors = F) 

data1 = data1 %>% mutate(mean = rowMeans(data1 %>% select(-chr, -start, -end), na.rm = T),
                         median = rowMedians(as.matrix(data1 %>% select(-chr, -start, -end)), na.rm = T))

data1 = data1 %>%
  select(chr, start, end, median, mean) %>%
  mutate(method = "WGS")

#GBS
data2 = data.table::fread("InfSNPsPerBin100000_GBS.txt", stringsAsFactors = F) 

data2 = data2 %>% mutate(mean = rowMeans(data2 %>% select(-chr, -start, -end), na.rm = T),
                         median = rowMedians(as.matrix(data2 %>% select(-chr, -start, -end)), na.rm = T))

data2 = data2 %>%
  select(chr, start, end, median, mean) %>%
  mutate(method = "GBS")

# join the data
data = rbind(data1, data2)

# set the new end position (sliding windows)
data = data %>% filter(chr != "Y")

# define the ymin and ymax coordinates for all chromosomes

coordinates = data.frame(chromosomes = c("X", 1:22),
                         ymin = seq(0.5, by = 3.5, length.out = 23),
                         ymax = seq(1.5, by = 3.5, length.out = 23),
                         stringsAsFactors = F)

# add the coordinates to the data

for(r in 1:nrow(data)){
  
  chromosome = data[r, "chr"]
  
  data[r, "ymin"] = coordinates %>% filter(chromosomes == as.character(chromosome)) %>% .$ymin
  
  data[r, "ymax"] = coordinates %>% filter(chromosomes == as.character(chromosome)) %>% .$ymax
}

# adjust the coordinates of the wgs data (move the columns over)

data$ymin = ifelse(data$method == "WGS", data$ymin + 1.5, data$ymin)
data$ymax = ifelse(data$method == "WGS", data$ymin + 1, data$ymax)

## plot the data side by side for gbs and wgs
ticks = seq(1.75, by = 3.5, length.out = 23)

 chrs = data %>% ggplot() + 
  geom_rect(aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = mean)) +
  geom_rect(data = blanks, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey98", colour = NA) +
  geom_rect(data = outlines, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = NA, colour = "black", size = 0.2) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 80.5), breaks = ticks, labels = c(c("X", 1:22))) +
  scale_x_continuous(limits = c(0, 248956422), expand = c(0, 0)) +
  scale_fill_viridis(discrete = F, direction = -1, na.value = "white", limits = c(NA, 300), trans = "log2", option = "viridis",
                     name = "log2 mean informative \n SNPs per 100kb window",
                     breaks = c(0.25, 4, 64)) +
  theme( axis.text.x = element_text(size = 6, face = "bold"),
         plot.title = element_text(size = 6),
         legend.text = element_text(size = 5),
         legend.title = element_text(size = 5),
         axis.line.y = element_blank(),
         axis.line.x = element_blank(),
         axis.ticks.y = element_blank(),
         axis.ticks.x = element_blank(),
         axis.text.y = element_blank(),
         legend.position = "bottom",
         legend.title.align = 0.5,
         legend.margin = margin(0,0,0,0),
         strip.text.x = element_text(size = 6),
         legend.key.height = unit(0.2, "cm"),
         legend.key.width = unit(0.5, "cm"),
         plot.background = element_rect(fill = "grey98"),
         panel.background = element_rect(fill = "grey98"),
         legend.background = element_rect(fill = "grey98")) +
  coord_flip() 

# save the plot
ggsave(chrs, filename = "coverageChrs.png",
       width = 14, height = 6, units = "cm", device = "png", dpi = 900)
