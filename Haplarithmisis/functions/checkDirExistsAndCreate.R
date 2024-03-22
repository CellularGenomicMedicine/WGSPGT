# checkDirExistsAndCreate - Function to heck if a directory exists and if not, create the directory.

checkDirExistsAndCreate <- function(dirPath){
	tryCatch({
		if (!file.exists(dirPath)){ # Cehck if directory does not exists
			dir.create(dirPath) # Create directory if it does not exist
		}
	}, error = function(e) { # Error handling for directory creation failure
		stop(paste0("Error creating ", dirPath, ": ", e)) # Raise an error message with details
	})
	dirPath # Return the directory path
}
