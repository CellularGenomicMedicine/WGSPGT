checkFileExists <-function(filePath,inputPath){
	if (!file.exists(filePath)) {
		write.table(filePath, paste(inputPath, paste(filePath,"has not been generated",sep=" "), sep = "/"))
		stop(paste(filePath," has not been generated",sep=" "))
	}

}
