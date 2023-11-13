###########################################################################################################################
# Author: Anouk Janssen
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: plotting copy number variation from VeriSeq data using ggplot2
# input: _processed.xls file from VeriSeq 

#delete all objects from memory
rm(list=ls(all=TRUE))

#load packages
library(ggplot2)

# set the directories
input_dir <- "path/to/input"
output_dir <- "path/to/output"

allfiles <- list.files(input_dir, pattern = "processed", all.files = T, full.names = T, recursive = T, )

# iterate over all VeriSeq file
for (i in 1:length(allfiles)){

  file <- allfiles[i]
  
# Extract family and embryo numbers from the file path
family_number <- sub(".*/PGD_(\\d+)/.*", "\\1", file)
embryo_number <- sub(".*Emb([0-9A-Za-z]+)-.*", "\\1", file)

# Read the data
data <- fread(file, header = T, skip = 48) # skip the 48 first lines 
data$POSITION <- as.numeric(data$POSITION)

# Adjust the chromosome positions to make genome-wide scale
for (chr in unique(data$CHROMOSOME)){
  ToAdd <- as.numeric(max(data[data$CHROMOSOME==chr,]$POSITION))
  data[data$CHROMOSOME==chr+1,"POSITION"] <- data[data$CHROMOSOME==chr+1,"POSITION"] +ToAdd
  Position=max(data[data$CHROMOSOME==chr+1,"POSITION"])
  Cumlength <- cbind(chr,Position)
  if (chr == 1){
  Cumlengths <- Cumlength
  } else {
    Cumlengths <- rbind(Cumlengths, Cumlength)
  }
}

# Calculate the maximum position of each chromosome to make vertical lines in the plot
max_positions <- data %>% group_by(CHROMOSOME) %>%   summarize(Max_position = max(POSITION))
max_positions$Middle_position <- c(0, diff(max_positions$Max_position) / 2)
max_positions$Middle_position <- max_positions$Max_position - max_positions$Middle_position
max_positions$Middle_position[1] <- max_positions$Middle_position[1] / 2

Veriseq <- ggplot(data) +
  geom_vline(data = max_positions, aes(xintercept = Max_position), linetype = "dashed", color = "lightgrey") +
  geom_hline(aes(yintercept = 3), color = "darkgreen") +
  geom_hline(aes(yintercept = 1), color = "darkred") +
  geom_point(aes(x = `POSITION`, y = `BIN COPY #`), size = 0.1, col = "black") +
  geom_line(aes(x = `POSITION`, y = `BIN COPY #`), col = "darkblue", lwd = 0.1) +
  geom_line(aes(x = `POSITION`, y = `REGION COPY #`), col = "green", lwd = 1) +
  theme_classic() +
  scale_x_continuous(expand=c(0,0), breaks = max_positions$Middle_position, labels = c(1:22,"X", "Y")) +
  scale_y_continuous(expand=c(0,0), limits = c(0,4), breaks = c(0.40, 0.80, 1.20, 1.60, 2.00, 2.40, 2.80, 3.20, 3.60, 4.00)) +
  labs(x = "Chromosomal Position", y = "Copy Number") +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 7),
        axis.text.y = element_text(size =7))

ggsave(filename = file.path(output_dir, paste0("CNVPlot_PGD",family_number, "_Emb", embryo_number,".png")),
       plot = Veriseq,
       dpi = 600,
       width = 300,
       height = 100,
       units = "mm")

}

