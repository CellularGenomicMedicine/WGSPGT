EmbryoTestReportPlot_plotPO <-function(Chrom,dataPo,EmbryoID,Int,ChrLength,Lab,Main,Ax){
  plot(dataPo$Position[dataPo$Chr == Chrom & dataPo[,EmbryoID]>0], dataPo[dataPo$Chr == Chrom & dataPo[,EmbryoID]>0,EmbryoID], "h",xaxt="n",yaxt="n",xlab="", xlim = c(0,ChrLength),ylab="PO", col = "pink",frame=F, ylim = c(-1,1),cex.lab=Lab,cex.main=Main,cex.axis=Ax)
  points(dataPo$Position[dataPo$Chr == Chrom & dataPo[,EmbryoID]<0], dataPo[dataPo$Chr == Chrom & dataPo[,EmbryoID]<0,EmbryoID], "h", col = "blue")
  points(dataPo$Position[dataPo$Chr == Chrom],rep(0,sum(dataPo[,2] == Chrom)), "p", col = "grey", pch = 8,cex=0.1)
  abline(h=seq(-1,1,by=0.5),lty=2)
  axis(side=2,at=c(-1,0,1),labels=c("-1","0","1"),cex.axis=Ax)
  abline(v=c(Int[Int[,1]==Chrom,2],Int[Int[,1]==Chrom,3]),lty=2,lwd=5,col="darkorange")
}#end Function

