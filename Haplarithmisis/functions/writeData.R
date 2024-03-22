# writeData - Function to write data to files based on family structure

writeData <- function(data, family, fam_members, parents, dataType, dataPath) {
  # Seperate family members into embryos and non-embryos
  fam_members_nonEmbryo <- fam_members[!fam_members[, "Sample_Status"] %in% "E",]
  fam_members_Embryo <- fam_members[fam_members[, "Sample_Status"] %in% "E",]
  
  # Check if there are embryos in the family
  if (nrow(fam_members_Embryo) > 0) {
    embs <- fam_members_Embryo[, 1]
    
    # Extract data for non-embryo family members
    data_non_embryo <- data.frame(data[, !names(data) %in% fam_members[, 1]], stringsAsFactors = FALSE, check.names = F)
    names(data_non_embryo) <- names(data)[!names(data) %in% fam_members[, 1]]
    
    # Merge data for non-embryo family members
    for (fam_member_nonEmbryo in fam_members_nonEmbryo[, 1]) {
      data_fam_member_nonEmbryo <- data.frame(data[, names(data) %in% fam_member_nonEmbryo], stringsAsFactors = FALSE, check.names = F)
      Names_data <- names(data_non_embryo)
      data_non_embryo <- data.frame(data_non_embryo, data_fam_member_nonEmbryo, stringsAsFactors = FALSE, check.names = F)
      names(data_non_embryo) <- c(Names_data, fam_member_nonEmbryo)
    }
    
    # Write data files for each parent-embryo combination
    for (parent in parents) {
      for (emb in embs) { # Write Gtype file per embryo
        # Exclude other embryos and include parents and Ref(s)
        data_fam_member_Embryo <- data[, names(data) %in% emb]
        data_non_embryo2 <- data.frame(data_non_embryo, data_fam_member_Embryo, stringsAsFactors = FALSE, check.names = F)
        names(data_non_embryo2) <- c(names(data_non_embryo), emb)
        parent_emb_datatype_filepath <- paste(dataPath, paste0(parent, "_", emb, "_", family, "_", dataType, ".txt"), sep = "/")
        writeDataFile(datatable = data_non_embryo2, filepath = parent_emb_datatype_filepath)
      } # End e loop
    }
  } else {
    # Write data file for the entire family if there are no embryos
    family_datatype_filepath <- paste(dataPath, paste(family, "_", dataType, ".txt", sep = ""), sep = "/")
    writeDataFile(datatable = data, filepath = family_datatype_filepath)
  }
}

# writeDataFile - Function to write data to a file

writeDataFile <- function(datatable = NULL, filepath = NULL) {
  tryCatch({
    print(paste0("Writing data to ", filepath, "..."))
    open_file <- file(description = filepath, "w")
    write.table(datatable, file = open_file, row.names = FALSE, quote = FALSE, sep = "\t")
    close(open_file)
    # Sanity check by reading back the written data
    datatable_check <- fread(filepath, header = TRUE, sep = "\t", data.table = FALSE, check.names = FALSE)
    if (nrow(datatable) != nrow(datatable_check)) {
      stop(paste0("Number of written rows does not match datatable!"))
    }
  }, error = function(e) {
    stop(paste0("Error writing ", filepath, ": ", e))
  })
}
