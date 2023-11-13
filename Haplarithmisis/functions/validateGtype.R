validateGtype <-function(Gtype,Chroms,outPath,Chrom){
	GtypeLess 			<- data.frame("chr" = 0, length = 0, stringsAsFactors = FALSE)
	for (chr in Chroms) {
		Gtypes <- Gtype[Gtype$Chr == chr,]
		if(nrow(Gtypes) < 25) {
			GtypeLess <- rbind(GtypeLess, c(chr, length(Gtypes)));
			write.table(chr, paste(outPath, paste("Less_then_25_Gtype_values_chr", chr, ".txt", sep = ""),sep="/"), sep = "\t")
			if (nrow(Gtypes) < 25 & chr==Chrom) {
				cat("Chromosome of interest ", chr, " has less then 25 Gtype values")
				stop(paste0("Chromosome of interest ",chr," has less then 25 Gtype"))
			}
		}
	}
}
