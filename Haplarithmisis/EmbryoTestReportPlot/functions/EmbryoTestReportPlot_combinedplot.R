EmbryoTestReportPlot_combinedplot <- function(script,dataHap,dataHapRaw,dataPo,BAFs,logRs,SegLogRs,PhBAF,ChrsLength,ideogram,parent,EmbryoID,outPathPlots,Chrom,Int,flip){
  Lab <- 3
  Main <- 6.5
  Ax <- 2.9
  C <- 2
  Czoom <- 3
  opt1=1
  opt2=2
  if (flip == 1) {opt1=2; opt2=1}

  Time <- gsub(" ","_",gsub(":","",gsub("-","",Sys.time())))

  png_fn <- paste0(outPathPlots,paste(parent,EmbryoID,Time,paste0("Chr",Chrom,".png"),sep="_"))
  png(png_fn,width=4000,height=1500,res=75)
    #Layout for plot without SNP-table
    layout.matrix <- matrix(c(1,2,4,6,8,10,12,14,16,18,20,22,1,3,5,7,9,11,13,15,17,19,21,23),ncol=2)
    layout(mat=layout.matrix, heights = c(1,1.5,0.5,2,0.5,0.5,2,0.4,0.5,2,2,0.5))

    ChrLength <- as.numeric(ChrsLength[ChrsLength[,1]==paste("chr",Chrom,sep=""),2])

    #1 - Title
    par(mar=c(0,0,0,0))
    plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim = c(0,ChrLength))
    text((ChrLength/2),0, paste("Chromosome",Chrom," (",EmbryoID,")"),cex=Main,font=2)

    #2 - Ideogram
    par(mar=c(5,6,0,1))
    source(paste(script, "analyses", "functions","plotCyto.R",sep="/"))
    plotCyto(Chrom,BandName=T,ChrsLength,ideogram)

    #3 - Ideogram zoom
    par(mar=c(5,6,0,1))
    source(paste(script, "analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_plotCytoZoom.R",sep="/"))
    EmbryoTestReportPlot_plotCytoZoom(Chrom,BandName=T,upcalc,downcalc,ideogram)

    #4 - PO
    par(mar=c(0.7,6,0,1))
    source(paste(script, "analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_plotPO.R",sep="/"))
    EmbryoTestReportPlot_plotPO(Chrom,dataPo,EmbryoID,Int,ChrLength,Lab,Main,Ax)

    #5 - PO zoom
    par(mar=c(0.7,6,0,1))
    source(paste(script, "analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_plotPOZoom.R",sep="/"))
    EmbryoTestReportPlot_plotPOZoom(Chrom,dataPo,EmbryoID,Int,downcalc,upcalc,Lab,Main,Ax)

    #6 - BAF
    par(mar=c(0,6,0,1))
    source(paste(script, "analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_plotBAF.R",sep="/"))
    EmbryoTestReportPlot_plotBAF(Chrom,ChrLength,BAFs,EmbryoID,Int,"BAF",Lab,Main,Ax,C)

    #7 - BAF zoom
    par(mar=c(0,6,0,1))
    source(paste(script, "analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_plotBAFZoom.R",sep="/"))
    EmbryoTestReportPlot_plotBAFZoom(Chrom,BAFs,EmbryoID,Int,downcalc,upcalc,Lab,Main,Ax,Czoom,C)

    #8 - Hap-Pat-Raw
    par(mar=c(1,6,1,1))
    source(paste(script, "analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_plotHap.R",sep="/"))
    EmbryoTestReportPlot_plotHap(Chrom,ChrLength,dataHapRaw,paste(EmbryoID,"_Pat",sep=""),opt1,opt2,"blue","cornflowerblue",Int)

    #9 - Hap-Pat-Raw zoom
    par(mar=c(1,6,1,1))
    source(paste(script, "analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_plotHapZoom.R",sep="/"))
    EmbryoTestReportPlot_plotHapZoom(Chrom,dataHapRaw,paste(EmbryoID,"_Pat",sep=""),opt1,opt2,"blue","cornflowerblue",Int,downcalc,upcalc)

    #10 - Hap-Pat
    par(mar=c(0.3,6,0.3,1))
    EmbryoTestReportPlot_plotHap(Chrom,ChrLength,dataHap,paste(EmbryoID,"_Pat",sep=""),opt1,opt2,"blue","cornflowerblue",Int)

    #11 - Hap-Pat zoom
    EmbryoTestReportPlot_plotHapZoom(Chrom,dataHap,paste(EmbryoID,"_Pat",sep=""),opt1,opt2,"blue","cornflowerblue",Int,downcalc,upcalc)

    #12 - BAF-Pat
    par(mar=c(0,6,0,1))
    source(paste(script, "analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_plotBAFParents.R",sep="/"))
    EmbryoTestReportPlot_plotBAFParents(Chrom,ChrLength,PhBAF[["P1"]],PhBAF[["P2"]],PhBAF[["P1Seg"]],PhBAF[["P2Seg"]],EmbryoID,Int,"Pat-BAF",Lab,Main,Ax,C)

    #13 - BAF-Pat zoom
    par(mar=c(0,6,0,1))
    source(paste(script, "analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_plotBAFParentsZoom.R",sep="/"))
    EmbryoTestReportPlot_plotBAFParentsZoom(Chrom,PhBAF[["P1"]],PhBAF[["P2"]],PhBAF[["P1Seg"]],PhBAF[["P2Seg"]],EmbryoID,Int,"Pat-BAF",downcalc,upcalc,Lab,Main,Ax,Czoom,C)

    #14 - Hap-Mat-Raw
    par(mar=c(1,6,1,1))
    EmbryoTestReportPlot_plotHap(Chrom,ChrLength,dataHapRaw,paste(EmbryoID,"_Mat",sep=""),opt1,opt2,"red","pink",Int)

    #15 - Hap-Mat-Raw zoom
    par(mar=c(1,6,1,1))
    EmbryoTestReportPlot_plotHapZoom(Chrom,dataHapRaw,paste(EmbryoID,"_Mat",sep=""),opt1,opt2,"red","pink",Int,downcalc,upcalc)

    #16 - Hap-Mat
    par(mar=c(0.3,6,0.3,1))
    EmbryoTestReportPlot_plotHap(Chrom,ChrLength,dataHap,paste(EmbryoID,"_Mat",sep=""),opt1,opt2,"red","pink",Int)

    #17 - Hap-Mat zoom
    par(mar=c(0.3,6,0.3,1))
    EmbryoTestReportPlot_plotHapZoom(Chrom,dataHap,paste(EmbryoID,"_Mat",sep=""),opt1,opt2,"red","pink",Int,downcalc,upcalc)

    #18 BAF-Mat
    par(mar=c(0,6,0,1))
    EmbryoTestReportPlot_plotBAFParents(Chrom,ChrLength,PhBAF[["M1"]],PhBAF[["M2"]],PhBAF[["M1Seg"]],PhBAF[["M2Seg"]],EmbryoID,Int,"Mat-BAF",Lab,Main,Ax,C)

    #19 BAF-Mat zoom
    par(mar=c(0,6,0,1))
    EmbryoTestReportPlot_plotBAFParentsZoom(Chrom,PhBAF[["M1"]],PhBAF[["M2"]],PhBAF[["M1Seg"]],PhBAF[["M2Seg"]],EmbryoID,Int,"Mat-BAF",downcalc,upcalc,Lab,Main,Ax,Czoom,C)

    #20 LogR
    par(mar=c(2,6,0,1))
    source(paste(script, "analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_plotLogR.R",sep="/"))
    EmbryoTestReportPlot_plotLogR(Chrom,ChrLength,logRs,SegLogRs,EmbryoID,Int,"logR",Lab,Main,Ax,C)

    #21 LogR zoom
    par(mar=c(2,6,0,1))
    source(paste(script, "analyses","EmbryoTestReportPlot","functions","EmbryoTestReportPlot_plotLogRZoom.R",sep="/"))
    EmbryoTestReportPlot_plotLogRZoom(Chrom,logRs,SegLogRs,EmbryoID,Int,downcalc,upcalc,"logR",Lab,Main,Ax,Czoom,C)

    #22 Position
    par(mar=c(0,0,0,1))
    plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim=c(0,ChrLength))
    text((ChrLength/2),0, paste("Position (Mb)"),cex=Ax,font=1)

    #23 Position
    par(mar=c(0,0,0,1))
    plot(0, xlab = "",ylab="", type = "n" ,col = "white",xaxt="n" ,yaxt="n",frame=F,xlim = c(downcalc,upcalc))
    text(((downcalc+upcalc)/2),0, paste("Position (Mb)"),cex=Ax,font=1,line=8)
    dev.off()

    print(paste("The plots were saved for ",parent," at ",outPathPlots))
}#end function
