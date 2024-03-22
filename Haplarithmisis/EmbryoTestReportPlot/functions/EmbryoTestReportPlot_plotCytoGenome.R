# EmbryoTestReportPlot_plotCytoGenome - Function to plot the ideogram

EmbryoTestReportPlot_plotCytoGenome <- function(BandName, CumLengths, ideogram){
# Copy chromosome lengths and calculate cumulative lengths
CumLengths <- ChrsLength
CumLengths[, "Length"] <- cumsum(as.numeric(ChrsLength[, 2]))
GenomeLength <- as.numeric(CumLengths[23, "Length"])

xRange <- c(0, GenomeLength)

# Plotting the basic structure of the genome plot
plot(0, axes = FALSE, xlab = "", ylab = "", type = "n", col = "gray", xlim = xRange)
	
# Draw a background rectangle for the entire genome
up <- 0.5
bottom <- -0.5

cylindrect(xRange[1], (bottom / 5), xRange[2], (up / 5), col = "#84848450", gradient = "y", nslices = 100)

ideogramGenome <- ideogram # Copy the ideogram data

# Adjust the positions of ideogram segments based on chromosome lengths
for(chr in CumLengths[,"Chromosome"][2:nrow(CumLengths)]){
  ToAdd <- as.numeric(CumLengths[grep(paste(chr, "$", sep = ""), CumLengths[, "Chromosome"]) - 1, "Length"])
  ideogramGenome[ideogramGenome$Chromosome == chr, "Start"] <- ideogramGenome[ideogramGenome$Chromosome == chr, "Start"] + ToAdd 
  ideogramGenome[ideogramGenome$Chromosome == chr, "Stop"] <- ideogramGenome[ideogramGenome$Chromosome == chr, "Stop"] + ToAdd 

	}#end chr loop

# Define the bands
bandes <- ideogramGenome
bandes$posGrafic <- bandes$Start + (bandes$Stop - bandes$Start) / 2
numBandes <- dim(bandes)[1]
gieStain <- levels(bandes$Color) 

colorsBandes <- c("red", "gray100", "black", "gray25", "gray50", "gray75", "gray50", "gray50")

# Plotting ideogram bands based on color and type
for(b in 1:numBandes){
   chromStart <- bandes$Start[b]
   chromEnd <- bandes$Stop[b]
   gBand <- bandes$Color[b]
   colorBanda <- colorsBandes[gBand]
   if(gBand == "acen"){
     #segments(chromStart,bottom,chromStart,up,col="black")
    } else { 
      if (gBand == "gvar"){
        rect(chromStart, bottom, chromEnd, up, col = colorBanda, angle = 90, density = 0.9)
        cylindrect(chromStart, bottom, chromEnd, up, col = colorBanda, gradient = "y", nslices = 50)
     } else {
      rect(chromStart, bottom, chromEnd, up, col = colorBanda, border = "gray50")
      cylindrect(chromStart, bottom, chromEnd, up, col = colorBanda, gradient = "y", nslices = 50)
     } 
    }

 }#end b loop

# Add band name if specified
if (BandName == TRUE){
  	text(bandes$posGrafic, rep(bottom - 0.2, numBandes), bandes[, 4], srt = -90, adj = 0, cex = 0.2, family = "Helvetica", xpd = TRUE)
}
}#end Function


