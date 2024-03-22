# XYplot - Function to plot the sex chromosomes (X and Y chromosome)

XYplot <- function(script, logRs, SegLogRs, ChrsLength, ideogram, outPathXY, colID, fam_members){
  # Source script to plot ideogram for X and Y chromosome
  source(paste(script, "analyses", "functions","plotCytoXY.R", sep = "/"))
  
  # Define parameters for plot aesthetics 
  Lab	<- 1.5
  Main <- 3
  Ax <- 1.5
  Ccc	<- 1

  # Extract chromosome lengths for X and Y chromosomes
  ChrsLength_X <- subset(ChrsLength, ChrsLength[, 1] == "chrX")
  ChrsLength_Y <- subset(ChrsLength, ChrsLength[, 1] == "chrY")
  ChrsLengths	<- rbind(ChrsLength_X, ChrsLength_Y)

  # Calculate cumulative lengths and genome length
  CumLengths <- ChrsLengths
  CumLengths[, "Length"] <- cumsum(as.numeric(ChrsLengths[, 2]))
  GenomeLength <- as.numeric(CumLengths[CumLengths[, "Chromosome"] == "chrY", "Length"])
  
  # Subset logR and segmented logR data for X and Y chromosomes
  logRs_X	<- subset(logRs, Chr=="X")
  logRs_Y	<- subset(logRs, Chr=="Y")
  logRs_XY	<- rbind(logRs_X, logRs_Y)
  
  SegLogRs_X	<- subset(SegLogRs, Chr=="X")
  SegLogRs_Y	<- subset(SegLogRs, Chr=="Y")
  SegLogRs_XY	<- rbind(SegLogRs_X, SegLogRs_Y)

  print("generating XY-plot")
  
  for(chr in CumLengths[, "Chromosome"][2:nrow(CumLengths)]){
    # Adjust positions based on cumulative lengths for plotting
    ToAdd <- as.numeric(CumLengths[grep(paste(chr, "$", sep = ""), CumLengths[, "Chromosome"]) - 1, "Length"])
    logRs_XY[logRs_XY$Chr == gsub("chr","",chr),"Position"] <- logRs_XY[logRs_XY$Chr==gsub("chr","",chr),"Position"] +ToAdd
    SegLogRs_XY[SegLogRs_XY$Chr==gsub("chr","",chr),"Position"] <- SegLogRs_XY[SegLogRs_XY$Chr==gsub("chr","",chr),"Position"] +ToAdd
  }
  
  # Calculate positions for chromosome labels
  Poses <- rep(0,2)
  Poses[1] <- (as.numeric(ChrsLengths[1, "Length"]) / 2)
  
  Poses[2] <- as.numeric(CumLengths[1, "Length"]) + (as.numeric(ChrsLengths[2, "Length"]) / 2)

  # Determine plot layout based on the number of family members
  if (length(colID) == 1) {
    # Layout for a single embryo
    # Set up png file and layout matrix
    print("Layout for embryo")
    png_fn <- paste(paste(outPathXY,"",sep = "/"),"XYprofile",".png", sep="")
    png(png_fn, width = 1500, height = 350,res = 75)
    layout.matrix <- matrix(c(1, 2, 3, 4), ncol = 1)
    layout(mat = layout.matrix, heights = c(0.5, 1, 0.5, 1.5))
  }
  if (length(colID) == 3) {
    # Layout for 3 family members
    # Set up png file and layout matrix
    print("Layout for 3 fam_members")
    png_fn <- paste(paste(outPathXY, "", sep = "/"),"XYprofile",".png", sep="")
    png(png_fn,width = 1500, height = 625,res = 75)
    layout.matrix <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 1)
    layout(mat = layout.matrix, heights = c(0.5, 1, 0.5, 1.5, 1.5, 1.5))
  }
  if (length(colID) == 4) {
    # Layout for 4 family members
    # Set up png file and layout matrix
    print("Layout for 4 fam_members")
    png_fn <- paste(paste(outPathXY,"",sep = "/"),"XYprofile",".png", sep="")
    png(png_fn,width = 1500, height = 750,res = 75)
    layout.matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7), ncol = 1)
    layout(mat = layout.matrix, heights = c(0.5, 1, 0.5, 1.5, 1.5, 1.5, 1.5))
  }

  # Plot components
  #1 Title
  par(mar = c(0, 6, 2, 1))
  plot(0, xlab = "", ylab = "", type = "n", col = "white", xaxt = "n", yaxt = "n", frame = FALSE, xlim = c(0, GenomeLength))
  text((GenomeLength / 2), 0, paste("XY profile"), cex = Main, font = 2)

  #2 Ideogram
  par(mar = c(2.5, 6, 1.5, 1))
  plotCytoXY(BandName = TRUE, ChrsLength, ideogram)

  #3 X-Y chromosomes
  par(mar = c(0, 6, 0, 1))
  plot(0, xlab = "",ylab = "", type = "n", col = "white", xaxt = "n", yaxt = "n", frame = FALSE, xlim=c(0, GenomeLength))
  for(chr in CumLengths[, "Chromosome"][seq(2, 2, 2)]){
    rect((as.numeric(CumLengths[CumLengths[, "Chromosome"] == chr, "Length"]) - as.numeric(ChrsLengths[ChrsLengths[, "Chromosome"] == chr, "Length"])), -1, CumLengths[CumLengths[, "Chromosome"] == chr, "Length"], 2, col = "#84848430", border = "#84848430")
  }# end chr loop
  text(Poses, 0, c("X", "Y"), cex = Lab, font = 2)

  #4-(5-6-(7)) LogR's
  for(sample in colID){
    par(mar = c(2, 6, 0, 1))
    plot(0, ylab = "LogR", col = "white", main = "", xlab = "", frame = FALSE, ylim = c(-3, 3), xlim = c(0, GenomeLength), cex = Ccc, cex.lab = 1.2, cex.main = Main, cex.axis = Ax, xaxt = "n", yaxt = "n")
    axis(side = 2, at = c(-3, -2, -1, 0, 1, 2, 3), labels = c("", "-2", "", "0", "", "2", ""), cex.axis = Ax)
    title(xlab = paste(fam_members[fam_members[, "SampleID"] %in% sample, "Sample_MetaInfo"], sample, sep = "-"), line = 0, cex.lab = 1.5)
    for(chr in CumLengths[,"Chromosome"][seq(2, 3, 2)]){
      rect((as.numeric(CumLengths[CumLengths[, "Chromosome"] == chr, "Length"]) - as.numeric(ChrsLengths[ChrsLengths[, "Chromosome"] == chr, "Length"])), -5, CumLengths[CumLengths[, "Chromosome"] == chr, "Length"], 5, col = "#84848430", border = "#84848430") 
    }# end chr loop
    abline(h = seq(-2, 2, 1), lty = 2, ylim = c(0, 1))
    points(logRs_XY[, "Position"], logRs_XY[, colnames(logRs_XY) %in% sample], pch = 20, col = "#00000030", cex = Ccc)
    points(SegLogRs_XY[, "Position"], SegLogRs_XY[, colnames(logRs_XY) %in% sample], pch = 20, col = "orange", cex = Ccc * 0.7)
  }
 
  dev.off()

  print(paste("The plot was saved at ", outPathXY))
  
}#end function

