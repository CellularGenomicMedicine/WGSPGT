plotCytoGenome_XY <-function(BandName,CumLengths,ideogram){
  
  ChrsLength_X <- subset(ChrsLength, ChrsLength[,1]=="chrX")
  ChrsLength_Y <- subset(ChrsLength, ChrsLength[,1]=="chrY")
  ChrsLengths <- rbind(ChrsLength_X, ChrsLength_Y)
  
  CumLengths <- ChrsLengths
  CumLengths[,"Length"] <- cumsum(as.numeric(ChrsLengths[,2]))
  GenomeLength <- as.numeric(CumLengths[CumLengths[,"Chromosome"]=="chrY","Length"])
  
  xRange<-c(0,GenomeLength)
  
  
  plot(0, axes = F, xlab = "",ylab="", type = "n", col = "gray", xlim = xRange)
  
  
  up <-0.5
  bottom <- -0.5
  
  cylindrect(xRange[1],(bottom/5),xRange[2],(up/5),col="#84848450",gradient="y",nslices=100)
  
  ideogramGenome <- ideogram
  ideogramGenome_X <- subset(ideogramGenome, ideogramGenome[,1]=="chrX")
  ideogramGenome_Y <- subset(ideogramGenome, ideogramGenome[,1]=="chrY")
  ideogramGenome <- rbind(ideogramGenome_X, ideogramGenome_Y)
  
  for(chr in CumLengths[,"Chromosome"][2:nrow(CumLengths)]){
    
    ToAdd <- as.numeric(CumLengths[grep(paste(chr,"$",sep=""),CumLengths[,"Chromosome"])-1,"Length"])
    ideogramGenome[ideogramGenome$Chromosome==chr,"Start"] <- ideogramGenome[ideogramGenome$Chromosome==chr,"Start"] + ToAdd 
    ideogramGenome[ideogramGenome$Chromosome==chr,"Stop"] <- ideogramGenome[ideogramGenome$Chromosome==chr,"Stop"] + ToAdd 
    
  }#end chr loop
  bandes<-ideogramGenome
  bandes$posGrafic <- bandes$Start+(bandes$Stop-bandes$Start)/2
  numBandes<-dim(bandes)[1]
  
  colorsBandes<-c("red","gray100","black","gray25","gray50","gray75","gray50","gray50")
  
  for(b in 1:numBandes){
    chromStart<-bandes$Start[b]
    chromEnd<-bandes$Stop[b]
    gBand<-bandes$Color[b]
    colorBanda<-colorsBandes[gBand]
    if(gBand == "gvar"){
      rect(chromStart,bottom,chromEnd,up,col=colorBanda,angle=90,density=0.9)
      cylindrect(chromStart,bottom,chromEnd,up,col=colorBanda,gradient="y",nslices=50)
    }else{ if(gBand != "acen"){
      rect(chromStart,bottom,chromEnd,up,col=colorBanda,border="gray50")
      cylindrect(chromStart,bottom,chromEnd,up,col=colorBanda,gradient="y",nslices=50)
     }
    }
  }#end b loop
  
  
  if (BandName==T){
    text(bandes$posGrafic,rep(bottom-0.2,numBandes),bandes[,4],srt=-90,adj=0,cex=0.5,family="Helvetica",xpd=TRUE)
  }
  
  
}#end Function



