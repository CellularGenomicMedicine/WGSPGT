checkDirExistsAndCreate <-function(dirPath){
	tryCatch({
		if (!file.exists(dirPath)){
			dir.create(dirPath)
		}
	}, error = function(e) {
		stop(paste0("Error creating ", dirPath, ": ", e))
	})
	dirPath
}
