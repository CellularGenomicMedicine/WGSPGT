# PreTestReportData_mendinc - Function to calculate Mendelian Inconsistency Percentage

PreTestReportData_mendinc <- function(script, config_file_fam, Parent, REF, Father, Mother, GT, outPath, Chroms){

  # Load all the scripts needed
  source(paste(script, "analyses", "PreTestReportData", "functions", "PreTestReportData_mendInc_Grandparents.R", sep = "/"))
  source(paste(script, "analyses", "PreTestReportData", "functions", "PreTestReportData_mendInc_Sibling.R", sep = "/"))
  source(paste(script, "analyses", "PreTestReportData", "functions", "PreTestReportData_mendInc_Grandparent.R", sep = "/"))
  
  # Calculate Mendelian inconsistencies based on different references types
  if (REF == "Grandparents"){
			source(as.character(config_file_fam))
			tot_Percent_MendInc <- PreTestReportData_mendInc_Grandparents(Grandmother, Grandfather, Parent, GT, Chroms)
  }
  
  if (REF == "SiblingF" | REF == "SiblingM"){
			source(as.character(config_file_fam))
			tot_Percent_MendInc <- PreTestReportData_mendInc_Sibling(Father, Mother, RefSampleID, GT, Chroms)
  }
		
  if (REF == "Grandmother" | REF == "Grandfather"){
			source(as.character(config_file_fam))
			tot_Percent_MendInc <- PreTestReportData_mendInc_Grandparent(RefSampleID, Parent, GT, Chroms)
  }
  
  # Write Mendelian inconsistency percentages to a file
  write.table(tot_Percent_MendInc,paste(paste(outPath, "", sep = "/"), parent, "_MendInc_", REF, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

}