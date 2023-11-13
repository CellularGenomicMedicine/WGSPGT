PreTestReportData_mendInc_Sibling <- function(Father,Mother,Ref,GT,Chroms){
 for(Chr in Chroms){
  GT_chrom <- GT[GT$Chr%in%Chr,]
  GT_chrom <- GT_chrom[GT_chrom[,Ref]!="NC" & GT_chrom[,Mother]!="NC" & GT_chrom[,Father]!="NC",]
  GT_mendInc <- GT_chrom[(GT_chrom[,Ref]=="AB" & GT_chrom[,Mother]=="AA" & GT_chrom[,Father]=="AA") |
                                (GT_chrom[,Ref]=="AB" & GT_chrom[,Mother]=="BB" & GT_chrom[,Father]=="BB") |

                                (GT_chrom[,Ref]=="AA" & GT_chrom[,Mother]=="BB" & GT_chrom[,Father]=="BB") |
                                (GT_chrom[,Ref]=="AA" & GT_chrom[,Mother]=="AA" & GT_chrom[,Father]=="BB") |
                                (GT_chrom[,Ref]=="AA" & GT_chrom[,Mother]=="BB" & GT_chrom[,Father]=="AA") |
                                (GT_chrom[,Ref]=="AA" & GT_chrom[,Mother]=="AB" & GT_chrom[,Father]=="BB") |
                                (GT_chrom[,Ref]=="AA" & GT_chrom[,Mother]=="BB" & GT_chrom[,Father]=="AB") |

                                (GT_chrom[,Ref]=="BB" & GT_chrom[,Mother]=="AA" & GT_chrom[,Father]=="AA") |
                                (GT_chrom[,Ref]=="BB" & GT_chrom[,Mother]=="AA" & GT_chrom[,Father]=="BB") |
                                (GT_chrom[,Ref]=="BB" & GT_chrom[,Mother]=="BB" & GT_chrom[,Father]=="AA") |
                                (GT_chrom[,Ref]=="BB" & GT_chrom[,Mother]=="AB" & GT_chrom[,Father]=="AA") |
                                (GT_chrom[,Ref]=="BB" & GT_chrom[,Mother]=="AA" & GT_chrom[,Father]=="AB"),]

  Percent_MendInc <- (nrow(GT_mendInc)/nrow(GT_chrom)*100)
  if(Chr==Chroms[1]){tot_Percent_MendInc <- Percent_MendInc} else {tot_Percent_MendInc <- rbind(tot_Percent_MendInc,Percent_MendInc)}
  }
 tot_Percent_MendInc <- cbind(Chroms,tot_Percent_MendInc)
 return(tot_Percent_MendInc)
}
