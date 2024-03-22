# plotCytoXY - Function to plot the cytogenetic bands of the sex chromosomes on a chromosome ideogram

plotCytoXY <- function(BandName, CumLengths, ideogram){
  
  # Select the sex chromosomes
  ChrsLength_X <- subset(ChrsLength, ChrsLength[, 1] == "chrX")
  ChrsLength_Y <- subset(ChrsLength, ChrsLength[, 1] == "chrY")
  ChrsLengths <- rbind(ChrsLength_X, ChrsLength_Y)
  
  # Calculate the cumulative lengths for the chromosomes
  CumLengths <- ChrsLengths
  CumLengths[, "Length"] <- cumsum(as.numeric(ChrsLengths[, 2]))
  GenomeLength <- as.numeric(CumLengths[CumLengths[, "Chromosome"] == "chrY", "Length"])
  
  # Define x-axis range based on genome length
  xRange <- c(0, GenomeLength)
  
  # Create an empty plot
  plot(0, axes = FALSE, xlab = "",ylab="", type = "n", col = "gray", xlim = xRange)
  
  # Define upper and bottom limits for bands
  up <- 0.5
  bottom <- -0.5
  
  # Plot background gradient
  cylindrect(xRange[1], (bottom / 5), xRange[2], (up / 5), col = "#84848450", gradient = "y", nslices = 100)
  
  # Subset ideogram data for sex chromosomes
  ideogramGenome <- ideogram
  ideogramGenome_X <- subset(ideogramGenome, ideogramGenome[, 1] == "chrX")
  ideogramGenome_Y <- subset(ideogramGenome, ideogramGenome[, 1] == "chrY")
  ideogramGenome <- rbind(ideogramGenome_X, ideogramGenome_Y)
  
  # Adjust band positions based on cumulative lengths
  for(chr in CumLengths[,"Chromosome"][2:nrow(CumLengths)]){
    ToAdd <- as.numeric(CumLengths[grep(paste(chr, "$", sep = ""), CumLengths[, "Chromosome"]) - 1, "Length"])
    ideogramGenome[ideogramGenome$Chromosome == chr, "Start"] <- ideogramGenome[ideogramGenome$Chromosome == chr, "Start"] + ToAdd 
    ideogramGenome[ideogramGenome$Chromosome == chr, "Stop"] <- ideogramGenome[ideogramGenome$Chromosome == chr, "Stop"] + ToAdd 
  }#end chr loop
  
  # Get cytogenetic bands information for the chromosome
  bandes <- ideogramGenome
  bandes$posGrafic <- bandes$Start + (bandes$Stop - bandes$Start) / 2
  numBandes <- dim(bandes)[1]
  
  # Define colors for different band types
  colorsBandes <- c("red", "gray100", "black", "gray25", "gray50", "gray75", "gray50", "gray50")
  
  # Plot each band with appropriate colors
  for(b in 1:numBandes){
    chromStart <- bandes$Start[b]
    chromEnd <- bandes$Stop[b]
    gBand <- bandes$Color[b]
    colorBanda <- colorsBandes[gBand]
    if(gBand == "gvar"){
      rect(chromStart,bottom,chromEnd,up,col=colorBanda,angle=90,density=0.9)
      cylindrect(chromStart,bottom,chromEnd,up,col=colorBanda,gradient="y",nslices=50)
    }else{ if(gBand != "acen"){
      rect(chromStart,bottom,chromEnd,up,col=colorBanda,border="gray50")
      cylindrect(chromStart,bottom,chromEnd,up,col=colorBanda,gradient="y",nslices=50)
     }
    }
  }#end b loop
  
  # Add band names if BandName is TRUE
  if (BandName == TRUE){
    text(bandes$posGrafic, rep(bottom - 0.2, numBandes), bandes[, 4], srt = -90, adj = 0, cex = 1, family = "Helvetica", xpd = TRUE)
  }
}#end Function



