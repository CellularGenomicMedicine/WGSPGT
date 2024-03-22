# Haplarithmisis_sexdetermination - Function for sex determination based on haplotype scores

Haplarithmisis_sexdetermination <- function(ScSexes, ParScore) {
	for (sc in names(ParScore)){
		ScSexes[sc, 1] <- sc
		if (is.na(ParScore[[sc]]["X", "Par"])) {
		  print(paste("Sex of ", sc, "could not be determined!!"))
		} else if (ParScore[[sc]]["X", "Par"] == 12) {
		  ScSexes[sc, 2] <- "female"
		  print(paste(sc, "is female"))
			} else if (ParScore[[sc]]["X", "Par"] == 22) {
			  ScSexes[sc, 2] <- "male"
			  print(paste(sc, "is male"))
				} else if (ParScore[[sc]]["X", "Par"] == 11) {
				  ScSexes[sc, 2] <- "female"
				  print(paste(sc, "is a female but without maternal X chromosome"))
					} else if (ParScore[[sc]]["X", "Par"] == 0) {
					  ScSexes[sc, 2] <- "ND"
					  print(paste("Chr.X nullisomy of ", sc))
					} else {
						  print(paste("Sex of ", sc, "could not be determined!!"))}
	}
	return(ScSexes)
}