checkFileExistsAndFread <- function(filePath, inputPath = NULL) {
  if (! is.null(inputPath)) {
    # inputPath is deprecated
    warning(paste0("Deprecated inputPath used while reading ", filePath))
  }
  if (!file.exists(filePath)) {
    stop(paste(filePath, " has not been generated", sep = " "))
  }
  tryCatch({
	fread(filePath, header = TRUE, sep = "\t", data.table = FALSE, check.names = FALSE)
  }, error = function(e) {
    stop(paste0("Error reading ", filePath,": ", e))
  })
}
