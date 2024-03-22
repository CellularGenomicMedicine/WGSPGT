# checkFileExistsAndFread - Function to check if the file exists at the specified filePath, read the file if it exists, 
# and raises an error if it does not exist or encounters an error while reading.

checkFileExistsAndFread <- function(filePath, inputPath = NULL) {
  if (! is.null(inputPath)) {
    # Check if inputPath is provided (deprecated)
    warning(paste0("Deprecated inputPath used while reading ", filePath))
  }
  if (!file.exists(filePath)) { # Check if the file does not exist
    stop(paste(filePath, " has not been generated", sep = " ")) # Raise an error indicating that the file has not been generated
  }
  tryCatch({
	fread(filePath, header = TRUE, sep = "\t", data.table = FALSE, check.names = FALSE) # Read the file using fread with specified parameters
  }, error = function(e) {
    stop(paste0("Error reading ", filePath,": ", e)) # Raise an error if an error occurs while reading the file
  })
}
