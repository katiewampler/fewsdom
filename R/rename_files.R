#' Rename default file names from Aqualog
#'
#' Takes the raw files downloaded from the Aqualog and renames them so file names are consistent
#' across all run types and runs for consistent analysis and processing later.
#'
#' @importFrom stringr str_detect str_sub
#'
#' @param meta the metadata table for the sample run
#' @param prjpath a string indicating the project file containing the data to process, they should not be in subfiles
#' @export

files_rename <- function(meta, prjpath){
  stopifnot(is.data.frame(meta) | is.character(prjpath) |file.exists(prjpath))

  #correct for any blanks not listed as B1
  file_correct <- list.files(prjpath)
  file_correct <- file_correct[stringr::str_detect(file_correct, ".dat")]
  correct <- file_correct[!stringr::str_detect(file_correct, "B1")]
  if(length(correct)> 0){
    correct_names <- stringr::str_sub(correct, start=3)
    correct_names <- paste("B1", correct_names, sep="")
    file.rename(from=paste(prjpath, correct, sep="/"), to=paste(prjpath, correct_names, sep="/"))
    note <- "Warning, multiple blanks were used in the sample Q, the following files were renamed for the code to run:"
    write.table(note, paste(prjpath, "/READ_ME.txt", sep=""), row.names=F, quote=F, col.names = F)
    write.table(correct, paste(prjpath, "/READ_ME.txt", sep=""), row.names=F, quote=F, col.names = F, append=T)
  }

  #determine which samples were run manually and which were run with the sample Q
  manual_rows <- which(meta$run_type == "manual" | meta$run_type == "Manual")
  sampleq_rows <- which(meta$run_type == "sampleQ" | meta$run_type == "sampleq")

  if(length(manual_rows) + length(sampleq_rows) != nrow(meta)){
    stop("Some samples couldn't be classified by run type, ensure run types in the metadata are listed as 'manual' or 'sampleQ'")
  }

  #clean manual samples
  for(f in manual_rows){
    #write names of raw data
    raw_blank <- paste(meta$data_identifier[f], " (0", meta$replicate_no[f] ,") - Waterfall Plot Blank.dat", sep="")
    raw_eem <- paste(meta$data_identifier[f], " (0", meta$replicate_no[f] ,") - Waterfall Plot Sample.dat", sep="")
    raw_abs <- paste(meta$data_identifier[f], " (0", meta$replicate_no[f] ,") - Abs Spectra Graphs.dat", sep="")

    #write names of clean data
    clean_blank <- paste(meta$unique_ID[f], "_blank.dat", sep="")
    clean_eem <- paste(meta$unique_ID[f], ".dat", sep="")
    clean_abs <- paste(meta$unique_ID[f], "_Abs.dat", sep="")

    #double check if raw data exsits
    blank_check <- file.exists(paste(prjpath, raw_blank, sep="/"))
    eem_check <- file.exists(paste(prjpath, raw_eem, sep="/"))
    abs_check <- file.exists(paste(prjpath, raw_abs, sep="/"))

    #check for missing files in metadata but not folder, if not missing, rename file
    if(blank_check == T){
      file.rename(from=paste(prjpath, raw_blank, sep="/"), to=paste(prjpath, clean_blank, sep="/"))
    }else{warning(paste("Missing sample: ", raw_blank, ". Please download from the project file.", sep=""))}

    if(eem_check == T){
      file.rename(from=paste(prjpath, raw_eem, sep="/"), to=paste(prjpath, clean_eem, sep="/"))
    }else{warning(paste("Missing sample: ", raw_eem, ". Please download from the project file.", sep=""))}

    if(abs_check == T){
      file.rename(from=paste(prjpath, raw_abs, sep="/"), to=paste(prjpath, clean_abs, sep="/"))
    }else{warning(paste("Missing sample: ", raw_abs, ". Please download from the project file.", sep=""))}
  }

  #clean sampleQ samples
  for(f in sampleq_rows){
    #write names of raw data
    raw_blank <- paste("B1", "S", meta$index[f], meta$data_identifier[f],"BEM.dat", sep="")
    raw_eem <- paste("B1", "S", meta$index[f], meta$data_identifier[f],"SEM.dat", sep="")
    raw_abs <- paste("B1", "S", meta$index[f], meta$data_identifier[f],"ABS.dat", sep="")

    #write names of clean data
    clean_blank <- paste(meta$unique_ID[f], "_blank.dat", sep="")
    clean_eem <- paste(meta$unique_ID[f], ".dat", sep="")
    clean_abs <- paste(meta$unique_ID[f], "_Abs.dat", sep="")

    #double check if raw data exsits
    blank_check <- file.exists(paste(prjpath, raw_blank, sep="/"))
    eem_check <- file.exists(paste(prjpath, raw_eem, sep="/"))
    abs_check <- file.exists(paste(prjpath, raw_abs, sep="/"))

    #check for missing files in metadata but not folder, if not missing, rename file
    if(blank_check == T){
      file.rename(from=paste(prjpath, raw_blank, sep="/"), to=paste(prjpath, clean_blank, sep="/"))
    }else{warning(paste("Missing sample: ", raw_blank, ". Please download from the project file.", sep=""))}

    if(eem_check == T){
      file.rename(from=paste(prjpath, raw_eem, sep="/"), to=paste(prjpath, clean_eem, sep="/"))
    }else{warning(paste("Missing sample: ", raw_eem, ". Please download from the project file.", sep=""))}

    if(abs_check == T){
      file.rename(from=paste(prjpath, raw_abs, sep="/"), to=paste(prjpath, clean_abs, sep="/"))
    }else{warning(paste("Missing sample: ", raw_abs, ". Please download from the project file.", sep=""))}

  }
}

#' Renames and organizes files from default Aqualog files
#'
#' Takes the default files from the Aqualog (absorbance, blank, and sample)
#' it renames the files to the sample identifier and creates folder for each
#' file type in the main file directory.
#'
#' @importFrom utils zip
#' @importFrom stringr str_detect
#'
#' @param prjpath the file path of the main file directory where data is located
#' @param meta_file file location of metadata table, used to determine how samples were run
#' @param meta_sheet optional, sheet name of the metadata only required if metadata is an .xlsx file
#' @param zip_files will create a zip file in the file directory of all the raw, non-renamed files from the Aqualog as a backup
#' @return metadata table with an added column for the unique ID
#' @export

clean_files <- function(prjpath, meta_file, meta_sheet, zip_files=T){
  stopifnot(is.character(c(prjpath, meta_sheet, meta_file))|
              is.logical(zip_files)| file.exists(prjpath))

  #Load Sample Log
  if(stringr::str_detect(meta_file, ".xlsx")){
    meta <- openxlsx::readWorkbook(meta_file, sheet=meta_sheet, detectDates = T)
  } else{
    meta <- read.csv(meta_file)
  }

  #check metadata
  stopifnot(sum(is.na(meta$RSU_area_1s)) == 0 | "data_identifier" %in% colnames(meta),
            "replicate_no" %in% colnames(meta), "integration_time_s" %in% colnames(meta))

  #check for duplicates that aren't marked as such
  meta_dup <- meta %>% group_by(data_identifier) %>% summarise(count=n())
  dups <- meta_dup$data_identifier[meta_dup$count > 1]
  for(x in dups){
    ndups <- meta_dup$count[meta_dup$data_identifier == x]
    meta$replicate_no[meta$data_identifier == x] <- 1:ndups
  }

  #get unique ID
  meta$unique_ID <- paste(meta$data_identifier, "_", meta$replicate_no, "_", meta$integration_time_s, "s", sep="")
  meta <- as.data.frame(meta)
  row.names(meta) <- meta$unique_ID

  #zip raw files to have a backup
  file_to_zip <- list.files(prjpath)[stringr::str_detect(list.files(prjpath), ".dat")]

  if(zip_files == T & length(file_to_zip) > 0){
    zip_name <- unlist(strsplit(prjpath, "/"))
    zip_name <- zip_name[length(zip_name)]
    zip(paste("rawfiles_", zip_name, ".zip", sep=""), paste(file_to_zip, sep="/"))
  }else if(length(file_to_zip) == 0 & zip_files == T){
    warning("No raw files were found to zip.")
  }

  #checks to see if there are files to rename
  if(sum(stringr::str_detect(list.files(prjpath), "SEM")) +
              sum(stringr::str_detect(list.files(prjpath), "Waterfall")) +
              sum(stringr::str_detect(list.files(prjpath), "Abs ")) +
              sum(stringr::str_detect(list.files(prjpath), "BEM")) +
              sum(stringr::str_detect(list.files(prjpath), "ABS")) > 0){
    #rename files
    files_rename(meta=meta, prjpath=prjpath)

    #identify EEMs and Absorbance files
    files <- list.files(prjpath)
    Abs <- files[str_detect(files, "_Abs.dat")]
    blank <- files[str_detect(files, "_blank.dat")]
    EEM <- files[str_detect(files, ".dat")]
    EEM <- EEM[!(EEM %in% c(Abs, blank))]

    stopifnot(length(Abs) > 0| length(blank) > 0| length(EEM) > 0)

    create_files(prjpath)

    file.copy(paste(prjpath, "/", Abs, sep=""), paste(prjpath, "/1_Absorbance/",sep=""))
    file.copy(paste(prjpath, "/", blank, sep=""), paste(prjpath, "/2_Blanks/",sep=""))
    file.copy(paste(prjpath, "/", EEM, sep=""), paste(prjpath, "/3_Samples/",sep=""))

    file.remove(paste(prjpath, "/", Abs, sep=""))
    file.remove(paste(prjpath, "/", blank, sep=""))
    file.remove(paste(prjpath, "/", EEM, sep=""))

    #check you have absorbance, blank, and sample data for each sample
    if(length(EEM) != length(blank) & length(EEM) != length(Abs)){
      stop("Warning: Your samples aren't matched, you're missing absorbance or EEM's data")
    }
  }

  return(meta)
}
