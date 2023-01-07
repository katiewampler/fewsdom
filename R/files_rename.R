#' Rename Default File Names from Aqualog
#'
#' Takes the raw files downloaded from the Aqualog and renames them so file names are consistent
#' across all run types and runs for consistent analysis and processing later.
#'
#' @importFrom stringr str_detect str_sub
#'
#' @param meta the metadata table for the sample run
#' @param prjpath a string indicating the project file containing the data to process, they should not be in subfiles
#' @export

files_rename <- function(meta=meta, prjpath){
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
