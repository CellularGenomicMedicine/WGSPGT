plotCyto <-function(Chrom,BandName=T,ChrsLength,ideogram){

chromName<-paste("chr",Chrom,sep="")
xRange<-c(0,as.numeric(as.character(ChrsLength[ChrsLength[,1]==chromName,2])))

plot(0, axes = F, xlab = "",ylab="", type = "n", col = "gray", xlim = xRange)

up <-0.5
bottom <- -0.5

cylindrect(xRange[1],(bottom/5),xRange[2],(up/5),col="#84848450",gradient="y",nslices=100)

bandes<-ideogram[ideogram$Chromosome==chromName,]
bandes$posGrafic <- bandes$Start+(bandes$Stop-bandes$Start)/2
numBandes<-dim(bandes)[1]

colorsBandes<-c("red","gray100","black","gray25","gray50","gray75","gray50","gray50")
#japintat<-0

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
if (BandName){
  	text(bandes$posGrafic,rep(bottom-0.2,numBandes),bandes[,4],srt=-90,adj=0,cex=2,font=2,family="Helvetica",xpd=TRUE)
}

}#end Function

