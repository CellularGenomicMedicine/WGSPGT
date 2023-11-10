###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualisation of mitochondrial genome coverage per position

# input: combined depth per mitochondrial genome position (.csv) - output from combineDepths.R, sample sheet

# output:  circos plot (.jpg), circos plot legend (.jpg)

###########################################################################################################################

# load packages

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(gridBase))

### load the data

depth = data.table::fread("dataDir/depthPerPos.csv.csv",
                          header = T, stringsAsFactors = F) %>%
  as.data.frame() %>%
  select(chr, POS)

### load the sample sheet

sampleSheet = data.table::fread("sampleSheet.csv", stringsAsFactors = F) 

### calculate the mean depth per condition
# conditions: validation (PGT-M / PGT-SR embryos), sub10x (deep-sequencing and subsampling to 10x), seq10x (direct 10x sequencing)

validation = sampleSheet %>% filter(condition == "validation") %>% select(sampleID)
sub10x = sampleSheet %>% filter(condition == "sub10x") %>% select(sampleID)
seq10x = sampleSheet %>% filter(condition == "seq10x") %>% select(sampleID)

depth$mean_validation = rowMeans(depth %>% select(all_of(validation)))
depth$mean_sub10x = rowMeans(depth %>% select(all_of(sub10x)))
depth$mean_seq10x = rowMeans(depth %>% select(all_of(seq10x)))

depth = depth %>% select(mean_validation, mean_sub10x, mean_seq10x)
                                                 
###### make a circular plot with circlize #######
  
### artificially increase the size of the segment corresponding to the indication
# the indication was MELAS (at position 3243)
repindication = depth[rep(3243, 30), ] 
depth2 = rbind(depth, repindication) %>% arrange(POS)

# set the parameters for plotting
  # split locations
split = c(3242, 3274) # split the whole replicated segment
  # number and factorise the sections
depth2$section = c(rep(1, times = 3242), rep(2, times = 31), rep(3, times = 13326))
depth2$section = factor(depth2$section, levels = c(1,2,3))
  
  # set the colour scheme
viridis = as.vector(viridis(n = 256))
fill = colorRamp2(breaks = seq(from = 0, to = 6500, by = 6500/255), colors = rev(viridis))

# initialise the output file
jpeg("outDir/coverageCircos.jpg",
       width = 3, height = 3, units = "in", res = 600)

# plot the tracks  
circos.par(start.degree = 90, gap.after = c(1.5, 1.5, 10))
set_track_gap(gap = c(0.04))
circos.heatmap(mat = depth2[ , "mean_validation"], col = fill, split = depth2$section, track.heigh = 0.12,
                 cluster = F)
circos.heatmap(mat = depth2[ , "mean_sub10x"], col = fill, split = depth2$section, track.heigh = 0.09,
                   cluster = F)
set_track_gap(gap = c(0.02))
circos.heatmap(mat = depth2[ , "mean_seq10x"], col = fill, split = depth2$section, track.heigh = 0.09,
                 cluster = F)

# end the plotting / file generation
circos.clear()
dev.off()

# separately generate the legent
jpeg("outDir/coverageCircos_legend.jpg",
       width = 3, height = 3, units = "in", res = 600)

plot.new()
circle_size = unit(1, "snpc")
lgd = Legend(at = c(0, 5000, 10000, 15000, 20000), title = "depth", col_fun = fill, title_position = "topcenter")
lgd2 = Legend(at = c(0, 5000, 10000, 15000, 20000), title = "depth", col_fun = fill, title_position = "topcenter")
lgd_list = packLegend(lgd, lgd2)
draw(lgd_list)
circos.clear()

dev.off()

