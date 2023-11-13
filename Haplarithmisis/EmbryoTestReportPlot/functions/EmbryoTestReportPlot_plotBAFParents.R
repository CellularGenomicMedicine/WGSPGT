EmbryoTestReportPlot_plotBAFParents <-function(Chrom,ChrLength,Origin1,Origin2,OriginSeg1,OriginSeg2,EmbryoID,Int,LabHead,Lab,Main,Ax,C){
  plot(0,ylab=LabHead,line=3.5,col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim=c(0,ChrLength),cex=4,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
  axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
  abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))
  points(Origin1[Origin1[,"Chr"]==Chrom,"Position"],Origin1[Origin1[,"Chr"]==Chrom,EmbryoID],pch=20,col="#0000ff10",cex=C)
  points(Origin2[Origin2[,"Chr"]==Chrom,"Position"],Origin2[Origin2[,"Chr"]==Chrom,EmbryoID],pch=20,col="#ff000010",cex=C)
  points(OriginSeg1[OriginSeg1[,"Chr"]==Chrom,"Position"],OriginSeg1[OriginSeg1[,"Chr"]==Chrom,EmbryoID],pch=15,col="#0000ff80",cex=C)
  points(OriginSeg2[OriginSeg2[,"Chr"]==Chrom,"Position"],OriginSeg2[OriginSeg2[,"Chr"]==Chrom,EmbryoID],pch=15,col="#ff000080",cex=C)
  abline(v=c(Int[Int[,1]==Chrom,2],Int[Int[,1]==Chrom,3]),lty=2,lwd=5,col="darkorange")
}#end Function

