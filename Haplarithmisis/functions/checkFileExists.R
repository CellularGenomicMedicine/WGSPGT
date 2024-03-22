# checkFileExists - Check is a file exists at the specified filePath, and raise an error if the it does not exists

checkFileExists <- function(filePath, inputPath){
	if (!file.exists(filePath)) { # Check if the file does not exist
		write.table(filePath, paste(inputPath, paste(filePath,"has not been generated", sep = " "), sep = "/")) # Create an error message indicating the file has not been generated
		stop(paste(filePath, " has not been generated", sep = " ")) # Raise an error with the error message
	}
}
