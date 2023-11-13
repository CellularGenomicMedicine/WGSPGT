###########################################################################################################################
# author: Anouk Janssen
# research group: Cellular Genomic Medicine, clinical genetics, Maastricht University Medical Centre + (MUMC+)

# purpose script: plotting the haplotype from pilot study samples

#delete all objects from memory
rm(list=ls(all=TRUE))

#load packages
library(ggplot2)
library(dplyr)
library(ggsignif)

out.dir <- ""

pilot_concordance <- read.table("PilotConcordanceTable.txt", sep = "\t", stringsAsFactors = F, header =T)

# change commas into dots
pilot_concordance$maternal_concordance <- gsub(",",".",pilot_concordance$maternal_concordance)
pilot_concordance$paternal_concordance <- gsub(",", ".", pilot_concordance$paternal_concordance)

pilot_concordance_longer <- pilot_concordance %>%
  pivot_longer(cols = ends_with("_concordance"),
               names_to = "Parent",
               values_to = "Concordance") %>%
mutate(Parent = str_replace(Parent, "_concordance", "")) %>%
  filter(tarcov %in% c(c("5X", "10X", "20X", "30X")))

pilot_concordance_longer$tarcov <- factor(pilot_concordance_longer$tarcov, levels = c("5X", "10X", "20X", "30X"))

# Create a list of conditions to compare
conditions <- c("FULL", "30X", "20X", "10X", "5X")

plot <- ggplot(pilot_concordance_longer, aes(x = tarcov, y = as.numeric(Concordance), fill = Parent)) +
  geom_boxplot(lwd=0.4, width = 0.3, outlier.colour = NA, position = position_dodge(width = 0.7)) +
  geom_point(aes(x = as.numeric(tarcov)-0.3, y =  as.numeric(Concordance)), size = 0.1, position = position_jitter(width = 0.1)) +
  scale_fill_manual(values = c(paternal = "cornflowerblue",
                               maternal = "#ffc0cb")) +
  facet_grid(~Parent) +
  scale_x_discrete(labels = c("5X","10X", "20X", "30X")) +
  scale_y_continuous(limits = c(95,100)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 6),
    legend.position = "none",
    strip.text = element_text(colour = "white"),
    strip.background = element_rect(colour = "white", size = 0.5)) +
  labs(x = "Target coverage", y = "Haplotype concordance (%)") +
    geom_signif(data = filter(pilot_concordance_longer, Parent == "paternal"),
                comparisons = list(c("5X", "10X"), c("10X", "20X"), c("20X", "30X"), c("10X", "30X")),
                test = t.test,
                map_signif_level = T, 
                vjust = -0.1,
                y_position = c(99.1,99.25,99.4,99.75),
                tip_length = 0.01,
                textsize = 1.0,
                size = 0.2) +
    geom_signif(data = filter(pilot_concordance_longer, Parent == "maternal"),
                comparisons = list(c("5X", "10X"), c("10X", "20X"), c("20X", "30X"), c("10X", "30X")),
                map_signif_level = T, 
                vjust = -0.1,
                y_position = c(99.1,99.25,99.4,99.75),
                tip_length = 0.01,
                textsize = 1.0,
                size = 0.2)
  
ggsave(filename = file.path(out.dir, "SuppFig3.png"),
       plot = plot,
       width = 140, 
       height = 60, 
       units = "mm", 
       dpi = 600)
