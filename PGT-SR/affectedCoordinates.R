###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: generate ideogram plotting coordinates for affected PGT-SR chromosomes

# input: chromosome coordinates

# output: affected chromosome coordinates

###########################################################################################################################

# generate the chromosome outlines for Hg38
suppressPackageStartupMessages(library(dplyr))

# set the indication-specific breakpoints
affectedChrs = data.frame(Chr = c(8, 16), breakpoint = c(22092231, 89257105))

# load the object containing the p, centromere and q arm coordinages & adjust the labels
outlines = read.csv("chromosomeOutlines.csv")

outlines = outlines %>% arrange(chromosome, xmin)
outlines$arm = rep(c("p", "centromere", "q"), 24)

# set the taper size
taper = 2000000

######### for the extra chromosome segments

## chr8

outlines_8_extra = outlines %>% filter(chromosome == 8, xmin <= 22092231)
outlines_8_extra$xmax = 22092231

coordinates_8_extra = data.frame(x = c(0, taper, outlines_8_extra[1, "xmax"],
                           outlines_8_extra[1, "xmax"], taper, 0),
                     y = c(0.6, 1, 1, 0, 0, 0.4),
                     arm = "p_extra", chr = outlines_8_extra[1, "chromosome"])

blanks_8_extra = data.frame(x = c(0, taper, 0, 0, taper, 0),
                            y = c(1, 1, 0.6, 0.4, 0, 0),
                            arm = "p_extra", chr = outlines_8_extra[1, "chromosome"])

## chr16

outlines_16_extra = outlines %>% filter(chromosome == 16, arm == "q")
outlines_16_extra$xmin = 89257105

coordinates_16_extra = data.frame(x = c(outlines_16_extra[1, "xmin"], outlines_16_extra[1, "xmax"] - 1000000, outlines_16_extra[1, "xmax"],
                                       outlines_16_extra[1, "xmax"], outlines_16_extra[1, "xmax"] - 1000000, outlines_16_extra[1, "xmin"]),
                                 y = c(1, 1, 0.8, 0.2, 0, 0),
                                 arm = "q_extra", chr = outlines_16_extra[1, "chromosome"])

blanks_16_extra = data.frame(x = c(outlines_16_extra[1, "xmax"], outlines_16_extra[1, "xmax"] - 1000000, outlines_16_extra[1, "xmax"], outlines_16_extra[1, "xmax"],
                                   outlines_16_extra[1, "xmax"] - 1000000, outlines_16_extra[1, "xmax"]),
                             y = c(1, 1, 0.8, 0.2, 0, 0),
                             arm = "q_extra", chr = outlines_16_extra[1, "chromosome"])

########### for the derivative chromosomes

# calculate the new chromosome lengths

extra_8 = 22092231

extra_16 = 90338345 - 89257105

der_8 = 145138636 - extra_8 + extra_16

der_16 = 90338345 - extra_16 + extra_8

xmin_der_8 = 145138636 - der_8

## adjust the relevant chromsome coordinates

der_outlines = outlines %>% filter(chromosome %in% c(8,16))

der_outlines$xmin[4] = xmin_der_8

der_outlines$xmax[3] = der_16

der_outlines$xmax[6] = der_8

der_outlines$chromosome = paste0("der_", der_outlines$chromosome)

outlines = der_outlines

##### calculate the new coordinates for the tapered outlines & associated blanks

coordinates = data.frame(chr = NA, x = NA, y = NA, arm = NA)

for(r in 1:nrow(outlines)){
  
  arm = outlines[r, "arm"]
  
  if(arm == "p"){
    
    segment = data.frame(x = c(outlines[r, "xmin"], outlines[r, "xmin"] + taper, outlines[r, "xmax"] - taper, outlines[r, "xmax"],
                               outlines[r, "xmax"], outlines[r, "xmax"] - taper, outlines[r, "xmin"] + taper, outlines[r, "xmin"]),
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

####### make the "blanks"

blanks = data.frame(chr = NA, x = NA, y = NA, arm = NA)

for(r in 1:nrow(outlines)){
  
  arm = outlines[r, "arm"]
  
  if(arm == "p"){
    
    segment = data.frame(x = c(outlines[r, "xmin"], outlines[r, "xmin"] + taper, outlines[r, "xmin"], outlines[r, "xmin"], outlines[r, "xmin"] + taper, outlines[r, "xmin"]),
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

################################ combine all the parts and save

coords_all = rbind(coordinates, coordinates_8_extra, coordinates_16_extra)
write.csv(coords_all, file = "extraChromosomeOutlines.csv", row.names = F)

blanks_all = rbind(blanks, blanks_8_extra, blanks_16_extra)
write.csv(blanks_all, file = "extraChromosomeBlanks.csv", row.names = F)
