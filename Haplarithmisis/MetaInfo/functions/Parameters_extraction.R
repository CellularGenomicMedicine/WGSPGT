Parameters_extraction <- function(outPath,fam_members){
	status_comp      <- fam_members[,"Sample_Status"]
	Fam_members_names <- fam_members[,"Sample_MetaInfo"]
	Value1_stat  <- status_comp[grepl("AF",status_comp)|grepl("AM",status_comp)]
	Value2_statGP  <- status_comp[ grepl("AGM",status_comp)| grepl("AGF",status_comp) |grepl("UGM",status_comp)| grepl("UGF",status_comp)]
	if(length(Value2_statGP)==2) {Value2 <- "Grandparents"}
	if(length(Value2_statGP)==1) {Value2 <- Fam_members_names[grepl(Value2_statGP,status_comp)]}
	Value2_statS  <- status_comp[ grepl("AS",status_comp)| grepl("US",status_comp)]
	if(length(Value2_statS)>1) {write.table(Value2_statS,paste(outPath,paste("There are multiple Sibling refs.txt",sep="_"),sep="/"),row.names=FALSE,quote=FALSE,na=''); stop("Multiple intervals")}
	if(length(Value2_statS)==1) {Value2 <- Fam_members_names[grepl(Value2_statS,status_comp)]}

	for(v1 in Value1_stat) {
		if (v1=="AF") {parent="Father"}
		if (v1=="AM") {parent="Mother"}
		parameters <- data.frame(Param=c("Parent","Seed"), Value=c(parent,Value2),stringsAsFactors=FALSE)
		write.table(parameters,paste(outPath,paste(family,parent,"parameters.txt",sep="_"),sep="/"),row.names=FALSE,quote=FALSE,na='')
	}
	return(Value2)
}