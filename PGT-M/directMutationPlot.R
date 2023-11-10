###########################################################################################################################
# Author: Rebekka Koeck
# Lab: Cellular Genomic Medicine, Clinical Genetics, Maastricht University (Medical Centre +)

# script purpose: visualise the results from direct mutation analysis vs. expected results according to haplarithmisis

# input: variant calling information

# output: plot

###########################################################################################################################

# load libraries
library(dplyr)
library(ggplot2)
library(tidyr)

#load the data

#in this dataset the result from direct mutation analysis was reported even if it was < 5 reads at the mutation site
results <- read.table("dataDir/DirectMutationsResults.txt", 
                      header = T, sep = "\t", na.strings = c("", "NA"), stringsAsFactors = F) 

# add information about whether the point mutation detection was successful
results <- results %>%
  mutate(DirectMutationDetCor = ifelse(HaplarithmisisResult == DirectMutationResult, "TRUE", "FALSE"), # see if direct mutation was correct compared to haplotyping
         DirectMutationDetCor = ifelse(is.na(DirectMutationDetCor), "FALSE", DirectMutationDetCor), # if DirectMutationDetCor is NA then it is indicated as FALSE
         NumReads = ifelse(NumReads == 0 & !is.na(DEL), DEL, NumReads)) # when there is a deletion, the number of reads corresponds to the number of reads either side of the deletion.

# add the genotype information per sample
results <- results %>%
  mutate(Total = A + C + G + T + N + ifelse(is.na(DEL), 0, DEL),
         Percent_A = A / Total * 100,
         Percent_C = C / Total * 100,
         Percent_G = G / Total * 100,
         Percent_T = T / Total * 100,
         Percent_N = N / Total * 100,
         Percent_DEL = DEL / Total * 100) %>%
  mutate(Cell_Number = row_number())

results = results %>% mutate(ExpA_count = Exp_A * Total / 100,
                             ExpC_count = Exp_C * Total / 100,
                             ExpG_count = Exp_G * Total / 100,
                             ExpT_count = Exp_T * Total / 100)

# add the counts supporting deletion "nucleotides"
results = results %>% mutate(del_count = Total - A - C - G - T,
                                     Expdel_count = Total - ExpA_count - ExpC_count - ExpG_count - ExpT_count)

# adjust the data from wide to long format
results_longer <- pivot_longer(
  results,
  cols = c(A, C, G, T, del_count, ExpA_count, ExpC_count, ExpT_count, ExpG_count, Expdel_count),
  names_to = "Nucleotide",
  values_to = "counts") %>%
  mutate(AnalysisMethod = "DirectMutation") %>% as.data.frame()

# add columns required to distinguish samples and methods for plotting (and set orders)
results_longer$SampleGroup = paste0(results_longer$Sample, results_longer$start)

results_longer$Type = ifelse(results_longer$Nucleotide %in% c("A", "C", "T", "G", "del_count"), "obs", "exp")
results_longer$Type = factor(results_longer$Type, levels = c("exp", "obs"))

results_longer$Nucleotide2 = gsub("Exp", "", results_longer$Nucleotide)
results_longer$Nucleotide2 = gsub("_count", "", results_longer$Nucleotide2)
results_longer$Nucleotide2 = factor(results_longer$Nucleotide2, levels = c("A", "C", "G", "T", "del"))

order = results_longer %>% filter(is.na(Exp_DEL)) %>% arrange(Total) %>% .$SampleGroup %>% unique(.)
order_del = results_longer %>% filter(!is.na(Exp_DEL)) %>% arrange(Total) %>% .$SampleGroup %>% unique(.)

#stacked bar plot

plot = ggplot(results_longer, aes(x = Type, y = counts, fill = Nucleotide2, alpha = Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = c("A" = "chartreuse3", 
                               "C" = "darkorange", 
                               "G" = "cornflowerblue", 
                               "T" = "brown2",
                               "del" = "grey44")) +
  scale_alpha_discrete(range = c(0.4, 1), guide = "none") +
  scale_size_manual(values = c("exp" = 0.1, "obs" = 1)) +
  labs(x = "", y = "Read count", fill = "Nucleotide") +
  facet_grid(. ~ factor(SampleGroup, levels = c(rev(order), rev(order_del)))) +
  theme_light() +
  theme(legend.position = "bottom",
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_line(linetype = "dashed", colour = "grey70", size = 0.05),
    panel.grid.major.y = element_line(linetype = "dashed", colour = "grey60", size = 0.1),
    axis.text.x = element_text(angle = 270, size = 5, vjust = 0.5),
    axis.text.y = element_text(size = 5),
    axis.title.y = element_text(size = 6),
    legend.text = element_text(face = "plain", color = "black", size = 5, vjust = 0.5),
    legend.key.size = unit(3, "mm"),
    legend.title = element_text(size = 6),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.spacing = unit(0.3, "mm"),
    panel.border = element_rect(size = 0.3),
    plot.margin = unit(c(1, 1, 1, 1), "mm")) +
  scale_y_continuous(expand = c(0,0))

ggsave(plot, file = "outDir/directMutationsPlot.png",
       device = "png",
       width = 120, height = 70, units = "mm")

