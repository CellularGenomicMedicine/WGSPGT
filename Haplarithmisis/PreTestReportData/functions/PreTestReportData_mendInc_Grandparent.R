PreTestReportData_mendInc_Grandparent <- function(Ref,Parent,GT,Chroms){
 for(Chr in Chroms){
   GT_chrom <- GT[GT$Chr%in%Chr,]
   GT_chrom <- GT_chrom[GT_chrom[,Parent]!="NC" & GT_chrom[,Ref]!="NC",]
   GT_mendInc <- GT_chrom[(GT_chrom[,Parent]=="AA" & GT_chrom[,Ref]=="BB") |
                          (GT_chrom[,Parent]=="BB" & GT_chrom[,Ref]=="AA"),]

   Percent_MendInc <- (nrow(GT_mendInc)/nrow(GT_chrom)*100)
   if(Chr==Chroms[1]){tot_Percent_MendInc <- Percent_MendInc} else {tot_Percent_MendInc <- rbind(tot_Percent_MendInc,Percent_MendInc)}
  }
  tot_Percent_MendInc <- cbind(Chroms,tot_Percent_MendInc)
 return(tot_Percent_MendInc)
}