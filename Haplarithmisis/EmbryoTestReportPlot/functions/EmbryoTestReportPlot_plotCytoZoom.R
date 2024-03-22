# EmbryoTestReportPlot_plotCytoZoom - Function to plot the ideogram (zoomed-in)

EmbryoTestReportPlot_plotCytoZoom <- function(Chrom, BandName = TRUE, upcalc, downcalc, ideogram){
  chromName <- paste("chr", Chrom, sep = "")
  xRange <- c(downcalc, upcalc)
  plot(0, axes = FALSE, xlab = "", ylab = "", type = "n", col = "gray", xlim = xRange)
  
  # Define upper and bottom limits for the plot
  up <- 0.5
  bottom <- -0.5
  # Draw a shaded background rectangle
  cylindrect(xRange[1],(bottom / 5), xRange[2], (up / 5), col = "#84848450", gradient = "y", nslices = 100)
  
  # Adjust band positions to fit within the plot region
  bandes <- ideogram[ideogram[, "Start"] < upcalc & ideogram[, "Stop"] > downcalc & ideogram$Chromosome == chromName, ]
  # Subset the ideogram data to the specified chromosome and region of interest
  bandes[bandes["Start"] < downcalc, "Start"] <- downcalc - 200000
  bandes[bandes["Stop"] > upcalc, "Stop"] <- upcalc + 200000
  
  # Calculate the graphical position of each band
  bandes$posGrafic <- bandes$Start + (bandes$Stop - bandes$Start) / 2
  
  # Get the number of bands for iteration
  numBandes <- dim(bandes)[1]
  gieStain <- levels(bandes$Color)
  
  # Define colors for different band types
  colorsBandes <- c("red", "gray100", "black", "gray25", "gray50", "gray75", "gray50", "gray50")

  # Loop through each band and plot them
  for(b in 1:numBandes){
    chromStart <- bandes$Start[b]
    chromEnd <- bandes$Stop[b]
    gBand <- bandes$Color[b]
    colorBanda <- colorsBandes[gBand]
    if(gBand == "gvar"){
      rect(chromStart, bottom, chromEnd, up, col = colorBanda, angle = 90, density = 0.9)
      cylindrect(chromStart, bottom, chromEnd, up, col = colorBanda, gradient = "y", nslices = 50)
     } else {
      rect(chromStart, bottom, chromEnd, up, col = colorBanda, border = "gray50")
      cylindrect(chromStart, bottom, chromEnd, up, col = colorBanda, gradient = "y", nslices = 50)
     }
  }#end b loop

  # Add band names if specified
  if (BandName){
  	text(bandes$posGrafic, rep(bottom - 0.2, numBandes), bandes[, 4], srt = -90, adj = 0, cex = 1.5, font = 2, family = "Helvetica", xpd = TRUE)
  }
}#end Function

