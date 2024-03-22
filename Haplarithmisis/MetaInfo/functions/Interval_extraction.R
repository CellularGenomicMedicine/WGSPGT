# Intervals_extraction - Function that extracts parameters from the provided family members dataframe from the samplesheet 

Intervals_extraction <- function(outPath, samplesheet){
  
  # Extracting the "Sample_Status" from the samplesheet
	status_comp <- samplesheet[, "Sample_Status"]
	
	# Extract values with 'AF' or 'AM' from the "Sample_Status" column
	Value1_stat <- status_comp[grepl("AF", status_comp) | grepl("AM", status_comp)]
	Family_Interval <- unique(samplesheet[, "Family_Interval"])
	
	# Check for multiple intervals and stop execution if found
	if(length(Family_Interval) > 1) {	
	  write.table(Family_Interval, paste(outPath, paste("Multiple intervals.txt", sep = "_"), sep = "/"), row.names = FALSE, quote = FALSE, na = ''); 
	  stop("Multiple intervals")
	}
	
	# Split the "Family_Interval" into separate columns and create a data frame
	interval_1 <- data.frame(t(unlist(strsplit(Family_Interval, "_"))), stringsAsFactors = FALSE)
	
	# Extract unique values from "Family_second_Interval column
	Family_second_Interval <- unique(samplesheet[, "Family_second_Interval"])
	
	# Check for multiple second intervals and stop execution if found
	if(length(Family_second_Interval) > 1) {	
	  write.table(Family_Interval, paste(outPath, paste("Multiple intervals.txt", sep = "_"), sep = "/"), row.names = FALSE, quote = FALSE, na = ''); 
	  stop("Multiple intervals")
	}
	
	# Split the "Family_second_Interval" into separate columns and create a data frame 
	interval_2 <- data.frame(t(unlist(strsplit(Family_second_Interval,"_"))), stringsAsFactors = FALSE)
	
	# Combine interval data frames into a single data frame "intervals"
	if(ncol(interval_2) > 1) { 
	  intervals <- rbind(interval_1, interval_2)
	  }
	if(ncol(interval_2) == 1) { 
	  intervals <- interval_1
	}
	
	# Iterate over intervals and parents to write interval information to output files
	for(interval in 1:nrow(intervals)){
		for(v1 in Value1_stat) {
			if (intervals[interval,4] == "Pat") { parent = "Father" }
			if (intervals[interval,4] == "Mat") { parent = "Mother" }
		  
		  # Create a data frame for the interval
			Interval <- data.frame(Chr = gsub("chr","",intervals[interval, 1]),
			                       Start = intervals[interval, 2],
			                       Stop = intervals[interval, 3],
			                       Origin = intervals[interval, 4],
			                       stringsAsFactors = FALSE)
			# Write interval information to output files
			
			write.table(Interval, paste(outPath, paste(family, parent, "intervals.txt", sep = "_"), sep = "/"),
			            sep="\t", row.names = FALSE, quote = FALSE, na = '')
		}
	}
	# Return the extracted intervals
	return(gsub("chr", "", intervals[interval, 1]))
}
