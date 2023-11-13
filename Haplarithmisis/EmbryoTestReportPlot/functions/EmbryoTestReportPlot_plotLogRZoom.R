EmbryoTestReportPlot_plotLogRZoom <-function(Chrom,logRs,SegLogRs,EmbryoID,Int,downcalc,upcalc,LabHead,Lab,Main,Ax,Czoom,C){
  plot(0,ylab=LabHead,line=3.5,col="white",main="",xlab="",frame=F,ylim=c(-1.5,1.5),xlim = c(downcalc,upcalc),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
  axis(side=2,at=c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5),labels=c("","-2","","-1","","0","","1","","2",""),cex.axis=Ax)
  abline(h=seq(-2,2,0.5),lty=2,ylim=c(0,1))
  points(logRs[logRs[,"Chr"]==Chrom,"Position"],logRs[logRs[,"Chr"]==Chrom,EmbryoID],pch=20,col="#00000030",cex=Czoom)
  points(SegLogRs[SegLogRs[,"Chr"]==Chrom,"Position"],SegLogRs[SegLogRs[,"Chr"]==Chrom,EmbryoID],pch=20,col="#ff000080",cex=Czoom)
  abline(v=c(Int[Int[,1]==Chrom,2],Int[Int[,1]==Chrom,3]),lty=2,lwd=5,col="darkorange")
  XaxPos = c(downcalc,(downcalc+1000000),((upcalc+downcalc)/2),(upcalc-1000000),upcalc)
  XaxPosa = round((XaxPos/1000000),digit=2)
  axis(side=1,at=XaxPos,labels=as.character(XaxPosa),cex.axis=Ax,lwd=0)
}#end Function

