PreTestReportData_mendInc_Grandparents <- function(Grandmother,Grandfather,Parent,GT,Chroms){
  for(Chr in Chroms){
    GT_chrom <- GT[GT$Chr%in%Chr,]
    GT_chrom <- GT_chrom[GT_chrom[,Parent]!="NC" & GT_chrom[,Grandmother]!="NC" & GT_chrom[,Grandfather]!="NC",]
    GT_mendInc <- GT_chrom[(GT_chrom[,Parent]=="AB" & GT_chrom[,Grandmother]=="AA" & GT_chrom[,Grandfather]=="AA") |
                         (GT_chrom[,Parent]=="AB" & GT_chrom[,Grandmother]=="BB" & GT_chrom[,Grandfather]=="BB") |

                         (GT_chrom[,Parent]=="AA" & GT_chrom[,Grandmother]=="BB" & GT_chrom[,Grandfather]=="BB") |
                         (GT_chrom[,Parent]=="AA" & GT_chrom[,Grandmother]=="AA" & GT_chrom[,Grandfather]=="BB") |
                         (GT_chrom[,Parent]=="AA" & GT_chrom[,Grandmother]=="BB" & GT_chrom[,Grandfather]=="AA") |
                         (GT_chrom[,Parent]=="AA" & GT_chrom[,Grandmother]=="AB" & GT_chrom[,Grandfather]=="BB") |
                         (GT_chrom[,Parent]=="AA" & GT_chrom[,Grandmother]=="BB" & GT_chrom[,Grandfather]=="AB") |

                         (GT_chrom[,Parent]=="BB" & GT_chrom[,Grandmother]=="AA" & GT_chrom[,Grandfather]=="AA") |
                         (GT_chrom[,Parent]=="BB" & GT_chrom[,Grandmother]=="AA" & GT_chrom[,Grandfather]=="BB") |
                         (GT_chrom[,Parent]=="BB" & GT_chrom[,Grandmother]=="BB" & GT_chrom[,Grandfather]=="AA") |
                         (GT_chrom[,Parent]=="BB" & GT_chrom[,Grandmother]=="AB" & GT_chrom[,Grandfather]=="AA") |
                         (GT_chrom[,Parent]=="BB" & GT_chrom[,Grandmother]=="AA" & GT_chrom[,Grandfather]=="AB"),]
    Percent_MendInc <- (nrow(GT_mendInc)/nrow(GT_chrom)*100)
    if(Chr==Chroms[1]){tot_Percent_MendInc <- Percent_MendInc} else {tot_Percent_MendInc <- rbind(tot_Percent_MendInc,Percent_MendInc)}
  }
  tot_Percent_MendInc <- cbind(Chroms,tot_Percent_MendInc)
 return(tot_Percent_MendInc)
}
