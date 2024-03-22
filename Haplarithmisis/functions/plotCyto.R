# plotCyto - Function to plot the cytogenetic bands on a chromosome ideogram

plotCyto <- function(Chrom, BandName = TRUE, ChrsLength, ideogram){
  # Define chromosome name
  chromName <- paste("chr", Chrom, sep = "")
  
  # Define x-axis range based on chromosome length
  xRange <- c(0, as.numeric(as.character(ChrsLength[ChrsLength[, 1] == chromName, 2])))
  
  # Create an empty plot
  plot(0, axes = FALSE, xlab = "", ylab = "", type = "n", col = "gray", xlim = xRange)
  # Define upper and bottom limits for bands
  up <- 0.5
  bottom <- -0.5
  
  # Plot background gradient
  cylindrect(xRange[1], (bottom / 5), xRange[2], (up / 5), col = "#84848450", gradient = "y", nslices = 100)
  
  # Get cytogenetic bands information for the chromosome
  bandes <- ideogram[ideogram$Chromosome == chromName, ]
  bandes$posGrafic <- bandes$Start + (bandes$Stop - bandes$Start) / 2
  numBandes <- dim(bandes)[1]
  
  # Define colors for different band types
  colorsBandes <- c("red", "gray100", "black", "gray25", "gray50", "gray75", "gray50", "gray50")
  
  # Loop through each band to plot
  for(b in 1:numBandes){
    chromStart <- bandes$Start[b]
    chromEnd <- bandes$Stop[b]
    gBand <- bandes$Color[b]
    colorBanda <- colorsBandes[gBand]
    
    # Plot bands based on their type
    if(gBand == "gvar"){
    rect(chromStart, bottom, chromEnd, up, col = colorBanda, angle = 90, density = 0.9)
    cylindrect(chromStart, bottom, chromEnd, up, col = colorBanda, gradient = "y", nslices = 50)
    } else {
      if(gBand != "acen"){
        rect(chromStart, bottom, chromEnd, up, col = colorBanda, border = "gray50")
        cylindrect(chromStart, bottom, chromEnd, up, col = colorBanda, gradient = "y", nslices = 50)
      }
    }
  } # end band  loop
  
# Add band names if BandName parameter is TRUE
if (BandName){
  	text(bandes$posGrafic, rep(bottom - 0.2, numBandes), bandes[, 4], srt = -90, adj = 0, cex = 2, font = 2, family = "Helvetica", xpd = TRUE)
}

}# end Function

