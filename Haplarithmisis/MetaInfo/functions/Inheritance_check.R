# Inheritance_check - Function to check for inheritance patterns based on the family members (parents, referent like grandparents or siblings)

Inheritance_check <- function(fam_members, RefID, Chrom){
	Father = fam_members[grepl("Father", fam_members[, "Sample_MetaInfo"]), "Sample_Status"]
	Mother = fam_members[grepl("Mother", fam_members[, "Sample_MetaInfo"]), "Sample_Status"]
	Grandmother = fam_members[grepl("Grandmother", fam_members[, "Sample_MetaInfo"]), "Sample_Status"]
	Grandfather = fam_members[grepl("Grandfather", fam_members[, "Sample_MetaInfo"]), "Sample_Status"]
	Sibling = fam_members[grepl("Sibling", fam_members[, "Sample_MetaInfo"]), "Sample_Status"]
	
	# Checking inheritance based on reference ID
	if( RefID == "Grandparents") {
	  # Checking different combinations of grandparents and parents
		if (Father == "AF" & Mother == "UM" & Grandfather == "AGF" & Grandmother == "UGM"){
			if (Chrom == "X") {
				cat("AF/UM/AGF/UGM with X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/UM/AGF/UGM has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "AF" & Mother == "UM" & Grandfather == "UGF" & Grandmother == "AGM"){
			if (Chrom == "X") {
				cat("AF/UM/UGF/AGM with X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/UM/UGF/AGM has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "UF" & Mother == "AM" & Grandfather == "AGF" & Grandmother == "UGM"){
			if (Chrom == "X") {
				cat("UF/AM/AGF/UGM with X has been tested", "\n")
				Unknow = FALSE
			} else {
				cat("UF/AM/AGF/UGM has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "UF" & Mother == "AM" & Grandfather == "UGF" & Grandmother == "AGM"){
			if (Chrom == "X") {
				cat("UF/AM/UGF/AGM with X has been tested", "\n")
				Unknow = FALSE
			} else {
				cat("UF/AM/UGF/AGM has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "AF" & Mother == "AM" & Grandfather == "AGF" & Grandmother == "UGM"){
			if (Chrom == "X") {
				cat("AF/AM/AGF/UGM with X is not possible, no difference between parents", "\n")
				Unknow = TRUE
			} else {
				cat("AF/AM/AGF/UGM is not possible, no difference between parents", "\n")
				Unknow = TRUE
			}
		}
		if (Father == "AF" & Mother == "AM" & Grandfather == "UGF" & Grandmother == "AGM"  ){
			if (Chrom == "X") {
				cat("AF/AM/UGF/AGM with X is not possible, no difference between parents", "\n")
				Unknow = TRUE
			} else {
				cat("AF/AM/UGF/AGM is not possible, no difference between parents", "\n")
				Unknow = TRUE
			}
		}
		if (Grandmother == "AGM" & Grandfather == "AGF" ){
			if (Chrom == "X") {
				cat("Both grandparents affected on X", "\n")
				Unknow = TRUE
			} else {
				cat("Both grandparents affected", "\n")
				Unknow = TRUE
			}
		}
	}
	if( RefID == "SiblingF" ) {
		if (Father == "AF" & Mother == "UM" & Sibling == "AS"){
			if (Chrom == "X") {
				cat("AF/UM/AS X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/UM/AS has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "UF" & Mother == "AM" & Sibling == "AS"){
			if (Chrom == "X") {
				cat("UF/AM/AS X has been tested", "\n")
				Unknow = FALSE
			} else {
				cat("UF/AM/AS has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "AF" & Mother == "AM" & Sibling == "AS"){
			if (Chrom == "X") {
				cat("AF/AM/AS X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/AM/AS has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "AF" & Mother == "UM" & Sibling == "US"){
			if (Chrom == "X") {
				cat("AF/UM/US X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/UM/US has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "UF" & Mother == "AM" & Sibling == "US"){
			if (Chrom == "X") {
				cat("UF/AM/US X has been tested", "\n")
				Unknow = FALSE
			} else {
				cat("UF/AM/US has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "AF" & Mother == "AM" & Sibling == "US"){
			if (Chrom == "X") {
				cat("AF/AM/US X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/AM/US has been tested", "\n")
				Unknow = FALSE
			}
		}
	}
	if( RefID == "SiblingM") {
		if (Father == "AF" & Mother == "UM" & Sibling == "AS"){
			if (Chrom == "X") {
				cat("AF/UM/AS X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/UM/AS has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "UF" & Mother == "AM" & Sibling == "AS"){
			if (Chrom == "X") {
				cat("UF/AM/AS X has been tested", "\n")
				Unknow = FALSE
			} else {
				cat("UF/AM/AS has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "AF" & Mother == "AM" & Sibling == "AS"){
			if (Chrom == "X") {
				cat("AF/AM/AS X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/AM/AS has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "AF" & Mother == "UM" & Sibling == "US"){
			if (Chrom=="X") {
				cat("AF/UM/US X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/UM/US has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "UF" & Mother == "AM" & Sibling == "US"){
			if (Chrom == "X") {
				cat("UF/AM/US X has been tested", "\n")
				Unknow = FALSE
			} else {
				cat("UF/AM/US has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father=="AF" & Mother=="AM" & Sibling=="US"){
			if (Chrom=="X") {
				cat("AF/AM/US X is not possible", "\n")
				Unknow=TRUE
			} else {
				cat("AF/AM/US has been tested", "\n")
				Unknow=FALSE
			}
		}
	}
	if( RefID == "Grandmother") {
		if (Father == "AF" & Mother == "UM" & Grandmother == "AGM"){
			if (Chrom == "X") {
				cat("AF/UM/AGM X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/UM/AGM has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "AF" & Mother == "UM" & Grandmother == "UGM"){
			if (Chrom == "X") {
				cat("AF/UM/UGM X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/UM/UGM has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "UF" & Mother == "AM" & Grandmother == "AGM"){
			if (Chrom == "X") {
				cat("UF/AM/AGM X has been tested", "\n")
				Unknow = FALSE
			} else {
				cat("UF/AM/AGM has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "UF" & Mother == "AM" & Grandmother == "UGM"){
			if (Chrom == "X") {
				cat("UF/AM/UGM X has been tested", "\n")
				Unknow = FALSE
			} else {
				cat("UF/AM/UGM has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "AF" & Mother == "AM" & Grandmother == "AGM"){
			if (Chrom == "X") {
				cat("AF/AM/AGM X is not possible, no difference between parents", "\n")
				Unknow = TRUE
			} else {
				cat("AF/AM/AGM is not possible, no difference between parents", "\n")
				Unknow = TRUE
			}
		}
		if (Father == "AF" & Mother == "AM" & Grandmother == "UGM"){
			if (Chrom=="X") {
				cat("AF/AM/UGM X is not possible, no difference between parents", "\n")
				Unknow=TRUE
			} else {
				cat("AF/AM/UGM is not possible, no difference between parents", "\n")
				Unknow=TRUE
			}
		}
	}
	if( RefID == "Grandfather") {
		if (Father == "AF" & Mother == "UM" & Grandfather == "AGF"){
			if (Chrom == "X") {
				cat("AF/UM/AGF  X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/UM/AGF  has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "AF" & Mother == "UM" & Grandfather == "UGF"){
			if (Chrom=="X") {
				cat("AF/UM/UGF  X is not possible", "\n")
				Unknow = TRUE
			} else {
				cat("AF/UM/UGF  has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "UF" & Mother == "AM" & Grandfather == "AGF"){
			if (Chrom == "X") {
				cat("UF/AM/AGF  X has been tested", "\n")
				Unknow = FALSE
			} else {
				cat("UF/AM/AGF  has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "UF" & Mother == "AM" & Grandfather == "UGF"){
			if (Chrom == "X") {
				cat("UF/AM/UGF  X has been tested", "\n")
				Unknow = FALSE
			} else {
				cat("UF/AM/UGF  has been tested", "\n")
				Unknow = FALSE
			}
		}
		if (Father == "AF" & Mother == "AM" & Grandfather == "AGF"){
			if (Chrom=="X") {
				cat("AF/AM/AGF  X is not possible, no difference between parents", "\n")
				Unknow = TRUE
			} else {
				cat("AF/AM/AGF is not possible, no difference between parents", "\n")
				Unknow = TRUE
			}
		}
		if (Father == "AF" & Mother == "AM" & Grandfather == "UGF"){
			if (Chrom == "X") {
				cat("AF/AM/UGF  is not possible, no difference between parents", "\n")
				Unknow = TRUE
			} else {
				cat("AF/AM/UGF  is not possible, no difference between parents", "\n")
				Unknow = TRUE
			}
		}
	}
	if(Unknow == TRUE){
		cat("Inherentence has not been tested" )
		stop("Inherentence has not been tested")
	}
}
