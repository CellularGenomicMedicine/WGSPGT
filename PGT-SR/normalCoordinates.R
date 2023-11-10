###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: generate the chromosome ideogram plotting coordinates with a taper

# input: basic (untapered) chromosome coordinates

# output: tapered chromosome coordinates

###########################################################################################################################

# load packages

suppressPackageStartupMessages(library(dplyr))

# load the object containing the p, centromere and q arm coordinages
outlines = read.csv("chromosomeOutlines.csv")

outlines = outlines %>% arrange(chromosome, xmin)
outlines$arm = rep(c("p", "centromere", "q"), 24)

# generate the coordinates for each segment

coordinates = data.frame(chr = NA, x = NA, y = NA, arm = NA)

taper = 2000000

for(r in 1:nrow(outlines)){
  
  arm = outlines[r, "arm"]
  
  if(arm == "p"){
    
    segment = data.frame(x = c(0, taper, outlines[r, "xmax"] - taper, outlines[r, "xmax"],
                               outlines[r, "xmax"], outlines[r, "xmax"] - taper, taper, 0),
                         y = c(0.6, 1, 1, 0.6, 0.4, 0, 0, 0.4),
                         arm = arm, chr = outlines[r, "chromosome"])
  }
  
  if(arm == "centromere"){
    
    segment = data.frame(x = c(outlines[r, "xmin"], outlines[r, "xmax"], outlines[r, "xmax"], outlines[r, "xmin"]),
                         y = c(0.6, 0.6, 0.4, 0.4),
                         arm = arm, chr = outlines[r, "chromosome"])
  
  }
  
  if(arm == "q"){
    
    segment = data.frame(x = c(outlines[r, "xmin"], outlines[r, "xmin"] + taper, outlines[r, "xmax"] - taper, outlines[r, "xmax"],
                               outlines[r, "xmax"], outlines[r, "xmax"] - taper, outlines[r, "xmin"] + taper, outlines[r, "xmin"]),
                         y = c(0.6, 1, 1, 0.6, 0.4, 0, 0, 0.4),
                         arm = arm, chr = outlines[r, "chromosome"])
  }
  
  
  coordinates = rbind(coordinates, segment)
  
}

coordinates = na.omit(coordinates)

write.csv(coordinates, "/Users/user/surfdrive/PhD/WGS_PGT_trial/coveragePlotting/2mbTaperOutlines.csv", row.names = F)

####### make the "blanks"

blanks = data.frame(chr = NA, x = NA, y = NA, arm = NA)

for(r in 1:nrow(outlines)){
  
  arm = outlines[r, "arm"]
  
  if(arm == "p"){
    
    segment = data.frame(x = c(0, taper, 0, 0, taper, 0),
                         y = c(1, 1, 0.6, 0.4, 0, 0),
                         arm = arm, chr = outlines[r, "chromosome"])
  }
  
  if(arm == "centromere"){
    
    segment = data.frame(x = c(outlines[r, "xmin"] - taper, outlines[r, "xmax"] + taper, outlines[r, "xmax"], outlines[r, "xmin"],
                               outlines[r, "xmin"] - taper, outlines[r, "xmax"] + taper, outlines[r, "xmax"], outlines[r, "xmin"]),
                         y = c(1, 1, 0.6, 0.6,
                               0, 0, 0.4, 0.4),
                         arm = c(rep("cent_top", 4), rep("cent_bottom", 4)), 
                         chr = outlines[r, "chromosome"])
    
  }
  
  if(arm == "q"){
    
    segment = data.frame(x = c(outlines[r, "xmax"], outlines[r, "xmax"] - taper, outlines[r, "xmax"], outlines[r, "xmax"],
                               outlines[r, "xmax"] - taper, outlines[r, "xmax"]),
                         y = c(1, 1, 0.6, 0.4, 0, 0),
                         arm = arm, chr = outlines[r, "chromosome"])
  }
  
  
  blanks = rbind(blanks, segment)
  
}

blanks = na.omit(blanks)

write.csv(blanks, "2mbTaperBlanks.csv", row.names = F)
