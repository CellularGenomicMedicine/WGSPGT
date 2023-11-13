EmbryoTestReportPlot_plotBAFZoom <-function(Chrom,BAFs,EmbryoID,Int,downcalc,upcalc,Lab,Main,Ax,Czoom,C){
  plot(0,ylab="BAF",line=3.5,col="white",main="",xlab="",frame=F,ylim=c(-0.1,1.1),xlim = c(downcalc,upcalc),cex=C,cex.lab=Lab,cex.main=Main,cex.axis=Ax,xaxt="n",yaxt="n")
  axis(side=2,at=c(0,0.5,1),labels=c("0","0.5","1"),cex.axis=Ax)
  abline(h=seq(0,1,0.25),lty=2,ylim=c(0,1))
  points(BAFs[BAFs[,"Chr"]==Chrom,"Position"],BAFs[BAFs[,"Chr"]==Chrom,EmbryoID],pch=20,col="#00000050",cex=Czoom)
  if(sum(Int[,1]==Chrom)>=1){abline(v=c(Int[Int[,1]==Chrom,2],Int[Int[,1]==Chrom,3]),lty=2,lwd=5,col="darkorange")}
}#end Function

