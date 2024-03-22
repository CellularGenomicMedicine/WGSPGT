# EmbryoTestReportPlot_genomebafplot - Function to generate a genome-wide plot

EmbryoTestReportPlot_genomebafplot <- function(script, BAFs, logRs, SegLogRs, PhBAF, ideogram, parent, Father, Mother, EmbryoID, outPathPlots){
  # Source function to plot the ideogram 
  source(paste(script, "analyses", "EmbryoTestReportPlot", "functions", "EmbryoTestReportPlot_plotCytoGenome.R", sep = "/"))

  # Define parameters for plotting aesthetics
  Lab <- 1.2
  Main <- 2.2
  Ax <- 1.2
  C = 0.4
  
  # Calculate chromosome length and genome length
  CumLengths <- ChrsLength
  CumLengths[,"Length"] <- cumsum(as.numeric(ChrsLength[, 2]))
  GenomeLength <- as.numeric(CumLengths[CumLengths[, "Chromosome"] == "chrX", "Length"])

  # Extract PhBAF data categories and for each category the segmented PhBAF data
  P1 <- PhBAF[["P1"]]
  P1Seg <- PhBAF[["P1Seg"]]
  P2 <- PhBAF[["P2"]]
  P2Seg <- PhBAF[["P2Seg"]]
  M1 <- PhBAF[["M1"]]
  M1Seg <- PhBAF[["M1Seg"]]
  M2 <- PhBAF[["M2"]]
  M2Seg <- PhBAF[["M2Seg"]]
  # Remove NA values from PhBAF data
  P1 <- na.omit(P1)
  P1Seg <- na.omit(P1Seg)
  P2 <- na.omit(P2)
  P2Seg <- na.omit(P2Seg)
  M1 <- na.omit(M1)
  M1Seg <- na.omit(M1Seg)
  M2 <- na.omit(M2)
  M2Seg <- na.omit(M2Seg)

  # Adjust positions based on chromosome length for plotting
  for(chr in CumLengths[, "Chromosome"][2:nrow(CumLengths)]){

   print(chr)
   ToAdd <- as.numeric(CumLengths[grep(paste(chr, "$", sep = ""), CumLengths[, "Chromosome"]) - 1, "Length"])

   # Adjust positions 
   BAFs[BAFs$Chr == gsub("chr", "", chr), "Position"] <- BAFs[BAFs$Chr == gsub("chr", "", chr), "Position"] +ToAdd
   logRs[logRs$Chr == gsub("chr", "", chr), "Position"] <- logRs[logRs$Chr == gsub("chr", "", chr), "Position"] +ToAdd
   SegLogRs[SegLogRs$Chr == gsub("chr", "", chr), "Position"] <- SegLogRs[SegLogRs$Chr == gsub("chr", "", chr), "Position"] +ToAdd

   P1[P1$Chr == gsub("chr", "", chr), "Position"] <- P1[P1$Chr == gsub("chr", "", chr), "Position"] +ToAdd
   P1Seg[P1Seg$Chr == gsub("chr", "", chr), "Position"] <- P1Seg[P1Seg$Chr == gsub("chr", "", chr), "Position"] +ToAdd
   P2[P2$Chr == gsub("chr", "", chr), "Position"] <- P2[P2$Chr == gsub("chr", "", chr), "Position"] +ToAdd
   P2Seg[P2Seg$Chr == gsub("chr", "", chr), "Position"] <- P2Seg[P2Seg$Chr == gsub("chr", "", chr), "Position"] +ToAdd

   M1[M1$Chr == gsub("chr", "", chr), "Position"] <- M1[M1$Chr == gsub("chr", "", chr), "Position"] +ToAdd
   M1Seg[M1Seg$Chr == gsub("chr", "", chr), "Position"] <- M1Seg[M1Seg$Chr == gsub("chr", "", chr), "Position"] +ToAdd
   M2[M2$Chr == gsub("chr", "", chr), "Position"] <- M2[M2$Chr == gsub("chr", "", chr), "Position"] +ToAdd
   M2Seg[M2Seg$Chr == gsub("chr", "", chr), "Position"] <- M2Seg[M2Seg$Chr == gsub("chr", "", chr), "Position"] +ToAdd
  }#end chr loop
  
  # Calculate positions for chromosome labels
  Poses <- rep(0, 23)
  Poses[1] <- (as.numeric(ChrsLength[1, "Length"]) / 2)

  for(i in 2:23){
   Poses[i] <- as.numeric(CumLengths[i - 1, "Length"]) + (as.numeric(ChrsLength[i, "Length"]) / 2)
  }

  # Set up png file for saving plot
  png_fn <- paste(outPathPlots, paste(parent, EmbryoID, "GenomeMultiProfile.png", sep = "_"), sep = "/")
  png(png_fn,width = 10000, height = 6500, res = 600)
  print(EmbryoID)

  # Define plot layout and plot components
  layout(rbind(matrix(1, 2, 10), matrix(2, 3.5, 10), matrix(3, 1, 10), matrix(4, 4, 10), matrix(5, 4, 10), matrix(6, 4, 10), matrix(7, 5, 10), matrix(8, 4, 10), matrix(9, 5, 10), matrix(10, 4, 10), matrix(11, 5, 10), matrix(12, 2, 10)))

  ChrsLength <- ChrsLength[ChrsLength[, "Chromosome"] != "chrY", ]
  ideogram<- ideogram[ideogram[, "Chromosome"] != "chrY", ]

  #1 - Title
  par(mar = c(0, 6, 2, 1))
  plot(0, xlab = "", ylab = "", type = "n", col = "white", xaxt = "n", yaxt = "n", frame = FALSE, xlim = c(0, GenomeLength))
  text((GenomeLength / 2), 0, paste("Whole-genome profile (", EmbryoID, ")"), cex = Main, font = 2)

  #2 - Ideogram
  par(mar = c(2.5, 6, 1.5, 1))
  EmbryoTestReportPlot_plotCytoGenome(BandName = TRUE, ChrsLength, ideogram)

  #3 - Chromosome indication background
  par(mar = c(0, 6, 0, 1))
  plot(0, xlab = "", ylab = "", type = "n", col = "white", xaxt = "n", yaxt = "n", frame = FALSE, xlim = c(0, GenomeLength))

  for(chr in CumLengths[, "Chromosome"][seq(2, 22, 2)]){
   rect((as.numeric(CumLengths[CumLengths[, "Chromosome"] == chr, "Length"]) - as.numeric(ChrsLength[ChrsLength[, "Chromosome"] == chr, "Length"])), -1, CumLengths[CumLengths[, "Chromosome"] == chr, "Length"], 2, col = "#84848430", border = "#84848430")
  }# end chromsome loop
  
  # Indicate the chromosome number
  text(Poses, 0, c(1:22, "X"), cex = Lab, font = 2)

  #4 - BAF
  par(mar = c(0, 6, 0, 1))
  plot(0, ylab = "BAF", col = "white", main = "", xlab = "", frame = FALSE, ylim = c(-0.1, 1.1), xlim = c(0, GenomeLength), cex = C, cex.lab = Lab, cex.main = Main, cex.axis = Ax, xaxt = "n", yaxt = "n")
  axis(side = 2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), cex.axis = Ax)
  for(chr in CumLengths[, "Chromosome"][seq(2, 22, 2)]){
   rect((as.numeric(CumLengths[CumLengths[, "Chromosome"] == chr, "Length"]) - as.numeric(ChrsLength[ChrsLength[, "Chromosome"] == chr, "Length"])), -1, CumLengths[CumLengths[, "Chromosome"] == chr, "Length"], 2, col = "#84848430", border = "#84848430")
  }# end chromosome loop
  abline(h = seq(0, 1, 0.25), lty = 2, ylim = c(0, 1))

  points(BAFs[, "Position"], BAFs[, EmbryoID], pch = 20, col = "#00000050", cex = C)

  #5 - PhBAF-Pat
  par(mar = c(0, 6, 0, 1))
  plot(0, ylab = "Pat-BAF", col = "white", main = "", xlab = "", frame = FALSE, ylim = c(-0.1, 1.1), xlim = c(0, GenomeLength), cex = C, cex.lab = Lab, cex.main = Main, cex.axis = Ax, xaxt = "n", yaxt = "n")
  axis(side = 2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), cex.axis = Ax)

  for(chr in CumLengths[, "Chromosome"][seq(2, 22, 2)]){
   rect((as.numeric(CumLengths[CumLengths[, "Chromosome"] == chr, "Length"]) - as.numeric(ChrsLength[ChrsLength[, "Chromosome"] == chr, "Length"])), -1, CumLengths[CumLengths[, "Chromosome"] == chr, "Length"], 2, col = "#84848430", border = "#84848430")
  }# end chromosome loop
  abline(h=seq(0, 1, 0.25), lty = 2, ylim = c(0, 1))

  points(P1[, "Position"], P1[,EmbryoID], pch = 20, col = "#0000ff10", cex = C)
  points(P2[, "Position"], P2[,EmbryoID], pch = 20, col = "#ff000010", cex = C)
  points(P1Seg[, "Position"], P1Seg[,EmbryoID], pch = 15, col = "#0000ff80", cex = C*1.5)
  points(P2Seg[, "Position"], P2Seg[,EmbryoID], pch = 15, col = "#ff000080", cex = C*1.5)

  #6 - PhBAF-Mat
  par(mar = c(0,6,0,1))
  plot(0, ylab = "Mat-BAF", col = "white", main = "", xlab = "", frame = FALSE, ylim = c(-0.1, 1.1), xlim = c(0, GenomeLength), cex = C, cex.lab = Lab, cex.main = Main, cex.axis = Ax, xaxt = "n", yaxt = "n")
  axis(side=2,at=c(0, 0.5, 1), labels = c("0", "0.5", "1"), cex.axis = Ax)

  for(chr in CumLengths[, "Chromosome"][seq(2, 22, 2)]){
   rect((as.numeric(CumLengths[CumLengths[, "Chromosome"] == chr, "Length"]) - as.numeric(ChrsLength[ChrsLength[, "Chromosome"] == chr, "Length"])), -1, CumLengths[CumLengths[, "Chromosome"] == chr, "Length"], 2, col = "#84848430", border = "#84848430")
  }# end chromosome loop
  abline(h=seq(0, 1, 0.25), lty = 2, ylim = c(0, 1))

  points(M1[, "Position"],M1[,EmbryoID], pch = 20, col = "#0000ff10", cex = C)
  points(M2[, "Position"],M2[,EmbryoID], pch = 20, col = "#ff000010", cex = C)
  points(M1Seg[, "Position"],M1Seg[,EmbryoID], pch = 15, col = "#0000ff80", cex = C*1.5)
  points(M2Seg[, "Position"],M2Seg[,EmbryoID], pch = 15, col = "#ff000080", cex = C*1.5)

  #7 - LogR Embryo
  par(mar = c(0, 6, 0, 1))
  plot(0, ylab = "logR", col = "white", main = "", xlab = "", frame = FALSE, ylim = c(-3,3), xlim = c(0,GenomeLength), cex = C, cex.lab=Lab, cex.main=Main, cex.axis=Ax,xaxt="n",yaxt="n")
  axis(side=2,at=c(-3,-2,-1,0,1,2,3),labels=c("","-2","","0","","2",""), cex.axis=Ax)

  for(chr in CumLengths[, "Chromosome"][seq(2, 22, 2)]){
   rect((as.numeric(CumLengths[CumLengths[, "Chromosome"] == chr, "Length"]) - as.numeric(ChrsLength[ChrsLength[, "Chromosome"] == chr, "Length"])), -5, CumLengths[CumLengths[, "Chromosome"] == chr, "Length"], 5, col = "#84848430",border="#84848430")
  }# end chr loop
  abline(h=seq(-2, 2, 1), lty = 2, ylim = c(0,1))

  points(logRs[, "Position"], logRs[,EmbryoID], pch =20, col = "#00000030", cex = C)
  points(SegLogRs[, "Position"], SegLogRs[,EmbryoID], pch =20, col = "orange", cex = C * 1.5)

  #8 - BAF Father
  par(mar = c(0, 6, 0, 1))
  plot(0, ylab = "Father-BAF", col = "white", main = "", xlab = "", frame = FALSE, ylim=c(-0.1, 1.1), xlim = c(0, GenomeLength), cex = C, cex.lab = Lab, cex.main = Main, cex.axis = Ax, xaxt = "n", yaxt = "n")
  axis(side = 2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), cex.axis = Ax)

  for(chr in CumLengths[, "Chromosome"][seq(2, 22, 2)]){
   rect((as.numeric(CumLengths[CumLengths[, "Chromosome"] == chr, "Length"]) - as.numeric(ChrsLength[ChrsLength[, "Chromosome"] == chr, "Length"])), -1, CumLengths[CumLengths[, "Chromosome"] == chr, "Length"], 2, col = "#84848430", border = "#84848430")
  }# end chr loop
  abline(h = seq(0, 1, 0.25), lty = 2, ylim = c(0, 1))

  points(BAFs[, "Position"], BAFs[, Father], pch =20, col = "#00000030", cex = C)
  
  #9 - LogR Father
  par(mar = c(0, 6, 0, 1))
  plot(0, ylab = "Father-logR", col = "white", main = "", xlab = "", frame = FALSE, ylim = c(-3, 3), xlim = c(0, GenomeLength), cex = C, cex.lab = Lab, cex.main = Main, cex.axis = Ax, xaxt = "n", yaxt = "n")
  axis(side = 2, at = c(-3, -2, -1, 0, 1, 2, 3), labels = c("", "-2", "", "0", "", "2", ""), cex.axis = Ax)

  for(chr in CumLengths[, "Chromosome"][seq(2,22,2)]){
   rect((as.numeric(CumLengths[CumLengths[, "Chromosome"] == chr,"Length"])-as.numeric(ChrsLength[ChrsLength[, "Chromosome"] == chr,"Length"])),-5,CumLengths[CumLengths[, "Chromosome"] == chr,"Length"],5, col = "#84848430",border="#84848430")
  }# end chr loop
  abline(h = seq(-2, 2, 1), lty = 2, ylim = c(0, 1))

  points(logRs[, "Position"],logRs[,Father], pch =20, col = "#00000030", cex = C)
  points(SegLogRs[, "Position"],SegLogRs[,Father], pch =20, col = "orange", cex = C*1.5)

  #10 - BAF Mother
  par(mar = c(0, 6, 0, 1))
  plot(0, ylab = "Mother-BAF", col = "white", main = "", xlab = "Position (Gb)",frame = FALSE, ylim = c(-0.1,1.1), xlim = c(0,GenomeLength), cex = C, cex.lab=Lab, cex.main=Main, cex.axis=Ax,yaxt="n",xaxt="n")
  axis(side = 2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), cex.axis = Ax)
  title(xlab = "Position (Gb)")
  for(chr in CumLengths[, "Chromosome"][seq(2, 22, 2)]){
   rect((as.numeric(CumLengths[CumLengths[, "Chromosome"] == chr,"Length"]) - as.numeric(ChrsLength[ChrsLength[, "Chromosome"] == chr,"Length"])), -5, CumLengths[CumLengths[, "Chromosome"] == chr, "Length"], 5, col = "#84848430", border = "#84848430")
  }# end chr loop
  abline(h = seq(0, 1, 0.25), lty = 2, ylim = c(0, 1))
  points(BAFs[, "Position"], BAFs[, Mother], pch =20, col = "#00000030", cex = C)

  #11 - LogR Mother
  par(mar = c(0, 6, 0, 1))
  plot(0, ylab = "Mother-logR", col = "white", main = "", xlab = "",frame = FALSE, ylim = c(-3,3), xlim = c(0, GenomeLength), cex = C, cex.lab = Lab, cex.main = Main, cex.axis = Ax,xaxt = "n", yaxt = "n")
  axis(side = 2, at = c(-3, -2, -1, 0, 1, 2, 3), labels = c("", "-2", "", "0", "", "2", ""), cex.axis = Ax)

  for(chr in CumLengths[, "Chromosome"][seq(2, 22, 2)]){
   rect((as.numeric(CumLengths[CumLengths[, "Chromosome"] == chr,"Length"]) - as.numeric(ChrsLength[ChrsLength[, "Chromosome"] == chr,"Length"])), -5, CumLengths[CumLengths[, "Chromosome"] == chr,"Length"], 5, col = "#84848430", border = "#84848430")
  }# end chr loop
  abline(h = seq(-2, 2, 1), lty = 2, ylim = c(0, 1))

  points(logRs[, "Position"], logRs[, Mother], pch = 20, col = "#00000030", cex = C)
  points(SegLogRs[, "Position"], SegLogRs[, Mother], pch = 20, col = "orange", cex = C * 1.5)
  axis(side = 1, at = seq(0, 3000000000, 500000000), labels = c("0.0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0"), cex.axis = Ax, xlab = "Position (Gb)")

  # Position
  par(mar = c(0.5, 6, 2, 1))
  plot(0, xlab = "",ylab = "", type = "n", col = "white", xaxt = "n", yaxt = "n", frame = FALSE, xlim = c(0, GenomeLength))
  text((GenomeLength / 2), 0, paste("Position (Gb)"), cex = Ax, font = 1)

  dev.off()

 }#end function
