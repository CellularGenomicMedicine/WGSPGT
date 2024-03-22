# EmbryoTestReportPlot_plotHapZoom - Function to plot the haplotype information (zoomed-in)

EmbryoTestReportPlot_plotHapZoom <- function(Chrom, data, EmbryoIDParent, opt1, opt2, colorParentopt1, colorParentopt2, Int, downcalc, upcalc){
  plot(0, xlab = "", ylab = "", type = "n", col = "white", xaxt = "n", yaxt = "n", frame = FALSE, xlim = c(downcalc, upcalc))
  abline(v = data$Position[as.character(data$Chr)==Chrom & data[,EmbryoIDParent] == opt1], col = colorParentopt1, lwd = 20)
  abline(v = data$Position[as.character(data$Chr)==Chrom & data[,EmbryoIDParent] == opt2], col =colorParentopt2, lwd = 20)#blue
  abline(v = c(Int[Int[, 1] == Chrom, 2], Int[Int[, 1] == Chrom, 3]), lty = 2, lwd = 5, col = "darkorange")
}#end Function

