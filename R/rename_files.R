#' Rename default file names from Aqualog
#'
#' Takes the raw files downloaded from the Aqualog and renames them so file names are consistent
#' across all run types and runs for consistent analysis and processing later.
#'
#' @importFrom stringr str_detect str_sub
#' @importFrom reader get.ext
#'
#' @param meta the metadata table for the sample run
#' @param prjpath a string indicating the project file containing the data to process, they should not be in subfiles
#' @export

files_rename <- function(meta, prjpath){
  stopifnot(is.data.frame(meta) | is.character(prjpath) |file.exists(prjpath))

  #correct for any blanks not listed as B1
  file_correct <- list.files(prjpath)
  file_correct <- file_correct[reader::get.ext(file_correct)=="dat"]
  correct <- file_correct[!stringr::str_detect(file_correct, "B1") & str_sub(file_correct, 1,1) == "B"]
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
#' @importFrom tools file_ext
#' @importFrom readxl read_excel
#'
#' @param prjpath the file path of the main file directory where data is located
#' @param meta_file file location of metadata table, used to determine how samples were run
#' @param meta_sheet optional, sheet name of the metadata only required if metadata is an .xlsx file
#' @param ... arguments to pass down to functions with 'clean_files'
#' @return metadata table with an added column for the unique ID
#'
#' @export
#' @examples
#' \dontrun{
#' meta <- clean_files(prjpath=prjpath, meta_file=meta_file,
#' meta_sheet = meta_sheet)
#' }
#'

clean_files <- function(prjpath, meta_file, meta_sheet="log", ...){
  stopifnot(is.character(c(prjpath, meta_file))| file.exists(prjpath))

  #Load Sample Log
  if(tools::file_ext(meta_file) == "xlsx"){
    meta <- readxl::read_excel(path=meta_file, sheet=meta_sheet)
  } else{
    meta <- read.csv(meta_file)
  }

  #check metadata
  stopifnot(sum(is.na(meta$RSU_area_1s)) == 0 | "data_identifier" %in% colnames(meta),
            "replicate_no" %in% colnames(meta), "integration_time_s" %in% colnames(meta),
            "run_type" %in% colnames(meta))

  #check for numbers in the blanks
  num_blanks <- which(stringr::str_detect(meta$data_identifier, "#")==T) #finds blanks with # sign
  #rename files and metadata
  if(length(num_blanks) > 1){
    for(x in num_blanks){
      files <- paste("B1S", meta$index[x], meta$data_identifier[x], c(".ogw","ABS.dat","SEM.dat","BEM.dat"),sep="")
      new_files <- stringr::str_remove_all(files, "#")
      file.rename(from=paste(prjpath, files, sep="/"), to=paste(prjpath, new_files, sep="/"))
      meta$data_identifier[x] <- stringr::str_remove_all(meta$data_identifier[x], "#")
    }
  }

  #check for - in the sample names
  num_dash <- which(stringr::str_detect(meta$data_identifier, "-")==T) #finds blanks with # sign
  if(length(num_dash) >0){
    #rename files and metadata
    for(x in num_dash){
      files <- paste("B1S", meta$index[x], meta$data_identifier[x], c(".ogw","ABS.dat","SEM.dat","BEM.dat"),sep="")
      new_files <- stringr::str_remove_all(files, "-")
      file.rename(from=paste(prjpath, files, sep="/"), to=paste(prjpath, new_files, sep="/"))
      meta$data_identifier[x] <- stringr::str_remove_all(meta$data_identifier[x], "-")
    }
  }

  if((length(num_blanks) + length(num_dash)) > 0){
    #write meta again, with files renamed to remove errors
    if(stringr::str_detect(meta_file, ".xlsx")){
      openxlsx::write.xlsx(meta, meta_file, sheetName = meta_sheet,
                           colNames = TRUE, rowNames = F, append = F)} else{
                             write.csv(meta, meta_file, col.names = T, row.names = F, quote=F)}
  }else{
    if(stringr::str_detect(meta_file, ".xlsx")){
      new_meta_file <- paste(unlist(str_split(meta_file, ".xlsx"))[1], "_doc_added.xlsx", sep="")
      openxlsx::write.xlsx(meta, new_meta_file, sheetName = meta_sheet,
                           colNames = TRUE, rowNames = F, append = F)} else{
                             new_meta_file <- paste(unlist(str_split(meta_file, ".csv"))[1], "_doc_added.csv", sep="")
                             write.csv(meta, new_meta_file, col.names = T, row.names = F, quote=F)}
  }
  }


  #check for duplicates that aren't marked as such
  meta_dup <- meta %>% dplyr::filter(replicate_no == 1) %>% dplyr::group_by(data_identifier) %>% dplyr::summarise(count=n())
  dups <- meta_dup$data_identifier[meta_dup$count > 1]
  for(x in dups){
    ndups <- meta_dup$count[meta_dup$data_identifier == x]
    meta$replicate_no[meta$data_identifier == x] <- 1:ndups
  }

  #get unique ID
  meta$unique_ID <- paste(meta$data_identifier, "_", meta$replicate_no, "_", meta$integration_time_s, "s", sep="")
  meta <- as.data.frame(meta)
  row.names(meta) <- meta$unique_ID

  create_files(prjpath)

  #checks to see if there are files to rename
  if(sum(stringr::str_detect(list.files(prjpath), "SEM")) +
              sum(stringr::str_detect(list.files(prjpath), "Waterfall")) +
              sum(stringr::str_detect(list.files(prjpath), "Abs ")) +
              sum(stringr::str_detect(list.files(prjpath), "BEM")) +
              sum(stringr::str_detect(list.files(prjpath), "ABS")) > 0){
    #rename files
    files_rename(meta=meta, prjpath=prjpath)}

    #identify EEMs and Absorbance files
    files <- list.files(prjpath)
    files <- files[reader::get.ext(files)=="dat"]
    Abs <- files[stringr::str_detect(files, "_Abs.dat")]
    blank <- files[stringr::str_detect(files, "_blank.dat")]
    EEM <- files[stringr::str_detect(files, ".dat")]
    EEM <- EEM[!(EEM %in% c(Abs, blank))]

    if((length(Abs) > 0| length(blank) > 0| length(EEM) > 0)){
      file.copy(paste(prjpath, "/", Abs, sep=""), paste(prjpath, "/1_Absorbance/",sep=""))
      file.copy(paste(prjpath, "/", blank, sep=""), paste(prjpath, "/2_Blanks/",sep=""))
      file.copy(paste(prjpath, "/", EEM, sep=""), paste(prjpath, "/3_Samples/",sep=""))

      file.remove(paste(prjpath, "/", Abs, sep=""))
      file.remove(paste(prjpath, "/", blank, sep=""))
      file.remove(paste(prjpath, "/", EEM, sep=""))
    }



    #check you have absorbance, blank, and sample data for each sample
    if(length(EEM) != length(blank) & length(EEM) != length(Abs)){
      stop("Warning: Your samples aren't matched, you're missing absorbance or EEM's data")
    }

  return(meta)
}
