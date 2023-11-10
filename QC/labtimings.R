###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualisation of the time taken for various lab workflow steps

# input: lab step timings data (.csv)

# output: plot

###########################################################################################################################

# load libraries

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))

##########################################################################################################
############################# make the co-ordinate generating function ###################################
##########################################################################################################

# each polygon (arrow) will be made of 5 points (left_bottom, left_top, right_bottom, right_top and point)
# each point needs to have an x co-rodinate and a y co-ordinate for plotting

generate_polygon_coordinates = function(data, startcol, endcol, ymin, ymax, pointlength = 0.06){

  ## x coordinates
  # the left positions (x) represents the start time of each step 
  data$left_bottom = data[ , startcol]
  data$left_top = data[ , startcol]
  # the point should end at the end time of the step
  data$point = data[ , endcol]
  # the right start of the arrow portion should be a certain distance from the arrow point (if the step takes enough time)
  data$right_bottom = ifelse(data[ , endcol] - data[ , startcol] >= pointlength, data$point - pointlength, data$left_bottom)
  data$right_top = data$right_bottom
  
  ## re-format
  #re-arrange the data into long format for plotting (each polygon now has 5 rows)
  data_longer = data %>% pivot_longer(cols = c("left_bottom", "left_top", "point", "right_bottom", "right_top"),
                                      names_to = "position",
                                      values_to = "x") %>%
    as.data.frame()
  
  ## y coordinates
  # add ymin for bottom positions, ymax for top positions and the mid-point of these values for the arrow points
  data_longer$y = ifelse(data_longer$position %in% c("left_bottom", "right_bottom"), ymin, ymax)
  data_longer$y = ifelse(data_longer$position == "point", mean(c(ymin, ymax)), data_longer$y)
  
  ## re-order the rows
  # re-order the rows of the data frame (geom_polygon plots in the order that it finds the rows)
  data_longer$position = factor(data_longer$position, levels = c("left_bottom", "left_top", "right_top", "point", "right_bottom"))
  data_longer = data_longer %>% arrange(stepID, position)
  
  return(data_longer)
}

##########################################################################################################
############################# apply the co-ordinate generating function ##################################
##########################################################################################################

################# generate some dummy data ############################################

# each step needs to have a unique identifier (to separate the coordinates into unique polygons) and a step type identifier for the fill colour 
# method1 and method2 data will be combined into one data frame so the unique identifiers can't be the same in both methods

# load the data

timings = data.table::fread("inDir/WGS_timing.csv",
                            stringsAsFactors = F) %>%
  as.data.frame()


# filter the data and add appropriate start and end times for each step

method1 = timings %>% filter(METHOD == "GBS")

method1 = method1 %>% mutate(total = handson + incubation)
method1$start = c(0, cumsum(method1$total)[1: nrow(method1) - 1])  
  
method1_total = method1 %>% select(Step, METHOD, total, start)
method1_total$end = cumsum(method1_total$total)
method1_total$stepID = paste0("gbs", c(1:nrow(method1_total)))

method1_handson = method1 %>% select(Step, METHOD, handson, start)
method1_handson = method1_handson %>% mutate(end = start + handson,
                                             stepID = c(paste0("gbs", 1:nrow(method1_handson))))

method2 = timings %>% filter(METHOD == "WGS")

method2 = method2 %>% mutate(total = handson + incubation)
method2$start = c(0, cumsum(method2$total)[1: nrow(method2) - 1])  

method2_total = method2 %>% select(Step, METHOD, total, start)
method2_total$end = cumsum(method2_total$total)
method2_total$stepID = paste0("wgs", c(1:nrow(method2_total)))

method2_handson = method2 %>% select(Step, METHOD, handson, start)
method2_handson = method2_handson %>% mutate(end = start + handson,
                                             stepID = c(paste0("wgs", 1:nrow(method2_handson))))


##################### calculate the polygon coordinates ################################# 

## for the total time of the steps

method1_longer = generate_polygon_coordinates(method1_total, startcol = "start", endcol = "end", 
                                          ymin = 2, ymax = 3, pointlength = 8)

method2_longer = generate_polygon_coordinates(method2_total, startcol = "start", endcol = "end", 
                                          ymin = 0.5, ymax = 1.5, pointlength = 8)

## combine the data

steps = rbind(method1_longer, method2_longer)

## for the handson times

method1_handson_longer = generate_polygon_coordinates(method1_handson, startcol = "start", endcol = "end",
                                                  ymin = 2, ymax = 3, pointlength = 8)

method2_handson_longer = generate_polygon_coordinates(method2_handson, startcol = "start", endcol = "end",
                                                  ymin = 0.5, ymax = 1.5, pointlength = 8)


## combine the data

handson = rbind(method1_handson_longer, method2_handson_longer)

##########################################################################################################
############################################# plot the data ##############################################
##########################################################################################################

# set the order of the fill colours
steps$Step = factor(steps$Step, levels = c("Fragmentation", "AdapterAddition", "Ligation", 
                                           "CleanUp", "SizeSelection", "QC", "Dilution", 
                                           "SuppressionPCR", "FragmentSizeAnalysis", 
                                           "Tagmentation", "LibraryPooling", 
                                           "PrepareReagents", "DiluteMDA"))

## plot

timings_plot = 
ggplot() +
  geom_vline(xintercept = seq(from = 60, to = 540, by = 60), linetype = "dashed", colour = "grey70", size = 0.2) +
  geom_polygon(data = steps, aes(x = x, y = y, fill = Step, group = stepID)) +
  geom_polygon(data = handson, aes(x = x, y = y, group = stepID), fill = NA, colour = "black", lwd = 0.2) +
  theme_classic() +
  scale_y_continuous(limits = c(0.3, 3.2), expand = c(0,0), breaks = c(1, 2.5), labels = c("WGS", "GBS")) +
  scale_x_continuous(limits = c(0, 540), expand = c(0,0), breaks = seq(from = 0, to = 540, by = 60), labels = c(0:9)) +
  xlab("Time (hours)") +
  scale_fill_manual(values = c("Fragmentation" = "#D53E4F", 
                               "AdapterAddition" = "#E75948", 
                               "Ligation" = "#F57748", 
                               "CleanUp" = "#FA9D59", 
                               "SizeSelection" = "#FDBE6E", 
                               "QC" = "#FDDB87", 
                               "Dilution" = "#F2EA91", 
                               "SuppressionPCR" = "#E1F399", 
                               "FragmentSizeAnalysis" = "#BEE5A0", 
                               "Tagmentation" = "#99D6A4", 
                               "LibraryPooling" = "#71C6A4", 
                               "PrepareReagents" = "#50A9AF", 
                               "DiluteMDA" = "#3288BD"),
                    labels = c("AdapterAddition" = "Adapter addition", 
                               "CleanUp" = "Clean-up", 
                               "SizeSelection" = "Size selection", 
                               "SuppressionPCR" = "Suppression PCR", 
                               "FragmentSizeAnalysis" = "Fragment Size Analysis", 
                               "LibraryPooling" = "Library pooling", 
                               "PrepareReagents" = "Prepare reagents", 
                               "DiluteMDA" = "Dilute MDA")) +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 5),
        axis.text.y = element_text(size = 6, face = "bold"),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        legend.justification = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(0,0,0,0),
        plot.margin = margin(t = 1, r = 5, b = 0, l = 1),
        panel.spacing = margin(0,0,0,0)) +
  guides(fill = guide_legend(nrow = 3))

ggsave(filename = "outDir/timeline.png",
       width = 11, height = 3, units = "cm", device = "png", dpi = 900)
