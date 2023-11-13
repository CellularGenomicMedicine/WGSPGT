Haplarithmisis_ChrXFatherAB <-function(Gtypes,Father){
	GTFatherX <- Gtypes[Gtypes[,"Chr"]=="X",Father]
	if(sum(GTFatherX=="AB")>=1){
		warning(paste("On the chr. X  of the father",sum(GTFatherX=="AB"),"(",sum(GTFatherX=="AB")/sum(Gtypes$Chr=="X"), "%) of the SNP-calls are heterozygous!!"))
		print("These will be treated as NoCalls")
		GTFatherX[GTFatherX=="AB"] <- "NC"#There could not be htz SNPs on chromosome X of the father
	}
	return(GTFatherX)
}
