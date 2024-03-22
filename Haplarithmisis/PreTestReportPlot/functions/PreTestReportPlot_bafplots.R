# PreTestReport_bafplots- Function to generate BAF plots for pre-test reports

PreTestReport_bafplots <- function(script, Chrom, upstream_boundary, downstream_boundary, ChrsLength, ideogram, outPath, colID, fam_members, bafplot, parent){
  # Load the functions for plotting the ideogram
  source(paste(script, "analyses", "functions", "plotCyto.R", sep = "/"))
  
  # Parameters used for plotting
  Main <- 3
  Ax <- 1.7

  # Determine layout based on the number of family members
  if (nrow(fam_members) == 3) {
    print("Layout for 3 fam_members")
    png_fn <- paste(outPath, "/PreTestReportPlot_BAF_Chr", Chrom, ".png", sep = "")
    png(png_fn, width = 1500, height = 625, res = 75)
    layout.matrix <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 1)
    layout(mat = layout.matrix, heights = c(0.5, 1.5, 1.5, 1.5, 1.5, 0.5))
  }
  
  if (nrow(fam_members) == 4) {
    print("Layout for 4 fam_members")
    png_fn <- paste(outPath, "/PreTestReportPlot_BAF_Chr", Chrom, ".png", sep = "")
    png(png_fn, width = 1500, height = 750, res = 75)
    layout.matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7), ncol = 1)
    layout(mat = layout.matrix, heights = c(0.5, 1.5, 1.5, 1.5, 1.5, 1.5, 0.5))
  }

  #1 - Title
  ChrLength <- as.numeric(ChrsLength[ChrsLength[, 1] == paste("chr", Chrom, sep = ""), 2])
  par(mar=c(0, 6, 2, 1))
  plot(0, xlab = "", ylab = "", type = "n", col = "white", xaxt = "n", yaxt = "n", frame = FALSE, xlim = c(0, ChrLength))
  text((ChrLength / 2), 0, paste("Chromosome", Chrom, "BAF-plot PreTestReport"), cex = Main, font = 2)

  #2 - Ideogram
  par(mar=c(4, 6, 1.5, 1))
  plotCyto(Chrom,BandName=T,ChrsLength,ideogram)
  #plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))

  #3-4-5-(6) BAF-plot
  for(sample in colID){
    par(mar = c(0, 6, 0, 1))
    plot(0, ylab = "BAF", line = 3.5, col = "white", main = "", xlab = "", frame = FALSE, ylim = c(-0.1, 1.1), xlim = c(0, ChrLength), cex = Main, cex.lab = Ax, cex.main = Main, cex.axis = Ax, xaxt = "n", yaxt = "n")
    axis(side = 2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), cex.axis = Ax)
    abline(h = seq(0, 1, 0.25), lty = 2, ylim = c(0, 1))
    points(bafplot[, "Position"], bafplot[, colnames(bafplot) %in% sample], pch = 20, col = "#00000030", cex = 1)
    upcalc <- upstream_boundary + 2000000
    downcalc <- downstream_boundary - 2000000
    if (downcalc < 0) {
      downcalc = 0
      }
    abline(v = c(downcalc, upcalc), lty = 2, lwd = 3, col = "darkorange")
    title(xlab = paste(fam_members[fam_members[, 1] %in% sample, 2], sample, sep = "-"), font = 2, line = -1, cex.lab = Ax)
  }

  # Chromosome position values
  XaxPos = round(seq(0, ChrLength, round(ChrLength / 4)) / 1000000, digits = 2)
  axis(side = 1, at = XaxPos * 1000000, labels = as.character(XaxPos), cex.axis = Ax, lwd = 0)

  #Matrix 7/8
  par(mar=c(0, 6, 2, 1))
  plot(0, xlab = "", ylab = "", type = "n", col = "white", xaxt = "n", yaxt = "n", frame = FALSE, xlim = c(0, ChrLength))
  text((ChrLength / 2), 0, paste("Position (Mb)"), cex = Ax, font = 1)

  # Save plot and print message
  dev.off()
  print(paste("The ", parent, " plot was saved at ", outPath))

} #end function
