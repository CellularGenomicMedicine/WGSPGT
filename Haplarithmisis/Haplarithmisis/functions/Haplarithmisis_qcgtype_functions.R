# Haplarithmisis_callrate - function to compute call rates for genotype data

Haplarithmisis_callrate <- function(Gtype){
	# Calculate overall call rate
	CallRate <- (sum(Gtype != "NC") / length(Gtype)) * 100
	
	# Calculate call rates for different genotypes
	CallRateAA <- (sum(Gtype == "AA") / length(Gtype)) * 100
	CallRateBB <- (sum(Gtype == "BB") / length(Gtype)) * 100
	CallRateHet <- (sum(Gtype == "AB") / length(Gtype)) * 100
	CallRateHom <- (sum(Gtype == "AA" | Gtype == "BB") / length(Gtype)) * 100
	CallRateNoCall <- (sum(Gtype == "NC") / length(Gtype)) * 100
	
	# Combine call rates into a matrix
	CallRates <- rbind(CallRate, CallRateAA, CallRateBB, CallRateHet, CallRateHom, CallRateNoCall)
	
	# Return the call rates matrix
	CallRates

}#end Haplarithmisis_callrate function


# Haplarithmisis_mendinc - function to compute Mendelian inconsistency rate

Haplarithmisis_mendinc <- function(GTFather, GTMother, GTChild){
  # Compute Mendelian inconsistencies
	Mdc <- sum((GTChild == "BB" & (GTFather == "AB" | GTFather == "BA") & GTMother == "AA") | 
	           (GTChild == "AA" & (GTFather == "AB" | GTFather == "BA") & GTMother == "BB") | 
	           (GTChild == "BB" & (GTMother == "AB" | GTMother == "BA") & GTFather == "AA") |
	           (GTChild == "AA" & (GTMother == "AB" | GTMother == "BA") & GTFather == "BB") |
	           (GTFather == "AA" & (GTChild == "AA" | GTChild == "BB")  & GTMother == "BB") |
	           (GTFather == "BB" & (GTChild == "AA" | GTChild == "BB")  & GTMother == "AA") |
	           (GTChild == "AA" & (GTFather == "NC" | GTFather == "NoCall") & GTMother == "BB") |
	           (GTChild == "BB" & (GTFather == "NC" | GTFather == "NoCall") & GTMother == "AA") |
	           (GTChild == "AA" & (GTMother == "NC" | GTMother == "NoCall") & GTFather == "BB") |
	           (GTChild == "BB" & (GTMother == "NC" | GTMother == "NoCall") & GTFather == "AA") |
	           (GTFather == "AA" & GTMother == "AA" & GTChild != "AA" & GTChild != "NC" & GTChild != "NoCall") |
	           (GTFather == "BB" & GTMother == "BB" & GTChild != "BB" & GTChild != "NC" & GTChild != "NoCall")
	           )

	# Compute Mendelian inconsistency rate
	MdcRate <- (Mdc / sum(GTChild != "NC")) * 100
	
	# Return the Mendelian inconsistency rate
	MdcRate

} #end Haplarithmisis_mendinc function