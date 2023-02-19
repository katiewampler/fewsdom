#' Load renamed and organized Aqualog data into R
#'
#' Takes the Aqualog files that have been renamed and put in separate folders within
#' the main directory via the 'clean_files' function and loads them into R as
#' eemlists. It will also check to see if maximum absorbance is below 0.3.
#'
#' @importFrom utils zip
#' @importFrom eemR eem_read
#'
#' @param prjpath The file path of the main file directory where data is located
#' @export
#' @examples
#' \dontrun{
#'   raw_eem <- load_eems(prjpath = prjpath)
#' }
#'  X <- raw_eem[[1]] #the raw eem samples
#'  X_blk <- raw_eem[[2]] #the matching blanks
#'  Sabs <- raw_eem[[3]] #the matching absorbance data

load_eems <- function(prjpath){
  stopifnot(is.character(prjpath) | file.exists(prjpath))
  #load in EEM's data
  X <- eemR::eem_read(paste(prjpath, "/3_Samples/",sep=""), recursive = F, import_function = "aqualog")
  X_blk <- eemR::eem_read(paste(prjpath, "/2_Blanks/",sep=""), recursive = F, import_function = "aqualog")

  #load in absorbance data
  Sabs <- absorbance_read(paste(prjpath, "/4_Clean_Absorbance", sep=""))

  #check absorbance (if over 0.3) will return warning for samples above, else will just run in background
  for(c in 2:ncol(Sabs)){
    samp <- Sabs[,c]
    max <- max(samp)
    if(max > 0.3){
      warning(paste("Sample",colnames(Sabs)[c], "has an aborbance of", round(max,3), "please dillute", sep=" "))
    }
  }

  return(list(X, X_blk, Sabs))
}

#' Reading absorbance data from txt and csv files
#'
#' Exact function from 'staRdom' package, extracted to prevent orphan package warnings.
#' Reading absorbance data from txt and csv files.
#'
#' @importFrom stringr str_extract_all regex str_remove str_extract
#' @importFrom dplyr arrange
#' @importFrom stats setNames
#'
#' @param absorbance_path directory containing absorbance data files or path to single file. See details for format of absorbance data.
#' @param order logical, data is ordered according to wavelength
#' @param recursive read files recursive, include subfolders
#' @param dec optional, either you set a decimal separator or the table is tested for . and ,
#' @param sep optional, either you set a field separator or it is tried to be determined automatically
#' @param verbose logical, provide more information
#' @param cores number of CPU cores to be used simultanuously
#' @noRd

absorbance_read <- function (absorbance_path, order = TRUE, recursive = TRUE, dec = NULL,
                         sep = NULL, verbose = FALSE, cores = parallel::detectCores(logical = FALSE)){
  if (dir.exists(absorbance_path)) {
    abs_data <- list.files(absorbance_path, full.names = TRUE,
                           recursive = recursive, no.. = TRUE, include.dirs = FALSE,
                           pattern = "*.txt|*.csv", ignore.case = TRUE)
    abs_data <- abs_data[!file.info(abs_data)$isdir]
  }
  else if (file.exists(absorbance_path)) {
    abs_data <- absorbance_path
  }
  else stop("Absorbance data was not found!")
  if (length(abs_data) < 1)
    stop("No valid files found in absorbance_path!")
  cl <- parallel::makeCluster(min(cores, length(abs_data)), type = "PSOCK")
  parallel::clusterExport(cl, c("dec", "sep", "verbose"), envir = environment())
  parallel::clusterEvalQ(cl, require(dplyr))
  parallel::clusterEvalQ(cl, require(stringr))
  abs_data <- parallel::parLapply(cl, abs_data, function(tab) {
    tryCatch({
      rawdata <- readLines(tab)
      data <- rawdata %>% sapply(stringr::str_remove, pattern = "([^0-9]*$)")
      first_number <- min(which((substr(data, 1, 1) %>%
                                   grepl("[0-9]", .))))
      last_number <- max(which((substr(data, 1, 1) %>%
                                  grepl("[0-9]", .))))
      if (is.null(sep) | is.null(dec)) {
        nsepdec <- data[first_number] %>% stringr::str_extract_all("[^-0-9eE]") %>%
          unlist()
        example_number <- data[first_number] %>% stringr::str_extract("([-]?[0-9]+[.,]?[0-9]+[eE]?[-0-9]+)$")
        if (is.null(dec) & length(nsepdec) > 1)
          dec <- example_number %>% stringr::str_replace("([-0-9eE]+)([.,]?)([-0-9eE]*)",
                                                "\\2")
        if (is.null(sep))
          sep <- gsub(pattern = dec, replacement = "",
                      x = data[first_number], fixed = TRUE) %>%
          stringr::str_extract(paste0("[^-0-9eE", dec, "]"))
        if (verbose)
          warning("processing", tab, ": using", sep,
                  "as field separator and", dec, "as decimal separator.",
                  fill = TRUE)
      }
      data <- str_split(data, sep)
      table <- data[(first_number):last_number] %>% unlist() %>%
        matrix(ncol = length(data[[first_number]]), byrow = TRUE) %>%
        data.frame(stringsAsFactors = FALSE) %>% mutate_all(gsub,
                                                            pattern = ifelse(dec != "", dec, "."), replacement = ".",
                                                            fixed = TRUE)
      table <- table %>% mutate_all(as.numeric)
      attr(table, "location") <- rep(tab, ncol(table) -
                                       1)
      if (ncol(table) == 2) {
        samples <- tab %>% basename() %>% stringr::str_replace_all(stringr::regex(".txt$|.csv$",
                                                                ignore_case = TRUE), "")
      }
      else {
        samples <- rawdata[[1]] %>% str_split(sep) %>%
          unlist() %>% matrix(ncol = length(.), byrow = TRUE) %>%
          data.frame(stringsAsFactors = FALSE) %>% .[-1]
      }
      table <- table %>% stats::setNames(c("wavelength", samples))
    }, error = function(err) {
      stop("Error while reading ", tab, ": ", err)
    })
  })
  parallel::stopCluster(cl)
  locations <- lapply(abs_data, function(tab) {
    attr(tab, "location")
  }) %>% unlist()
  if (length(abs_data) == 1)
    abs_data <- abs_data[[1]] %>% as.data.frame()
  else abs_data <- abs_data %>% list_join(by = "wavelength")
  if (order)
    abs_data <- abs_data %>% dplyr::arrange(wavelength)
  attr(abs_data, "location") <- locations
  abs_data
}


#' Full join of a list of data frames.
#'
#' Exact function from 'staRdom' package, extracted to prevent orphan package warnings.
#' Full join of a list of data frames.
#'
#' @importFrom stringr str_extract_all
#'
#' @param df_list list of data frames to by joined
#' @param by character vector containing information how to join data frames. Format to be according to by in full_join. Each data frame has to contain the column(s) used for joining.
#' @noRd
list_join <- function (df_list, by)
{
  df <- df_list[[1]]
  for (n in 2:length(df_list)) {
    df <- dplyr::full_join(df, df_list[[n]], by = by)
  }
  df
}

#' Replace matched patterns in sample names
#'
#' Exact function from 'staRdom' package, extracted to prevent orphan package warnings.
#' Replace matched patterns in sample names.
#'
#' @importFrom stringr str_extract_all
#'
#' @param eem_list data of class eemlist
#' @param pattern	 character vector containing pattern to look for.
#' @param replacement character vector of replacements. Has to have the same length as pattern
#' @noRd
eem_name_replace <- function (eem_list, pattern, replacement)
{eem_list <- lapply(eem_list, function(eem) {
    eem$sample <- eem$sample %>% stringr::str_replace_all(pattern,
                                                          replacement)
    eem
  })
  class(eem_list) <- "eemlist"
  eem_list
}


#' Checks to make sure samples and metadata match
#'
#' Takes the loaded EEMs and absorbance data and tries to match it to the 'unique_ID'
#' in the metadata, if a sample is in the metadata but not in the loaded EEMs and
#' absorbance it will remove that line in the metadata. If there are samples loaded
#' not in the metadata it will warn and remove from the loaded samples.
#'
#' @param meta metadata table with sample information, needs to have a 'unique_ID' column
#' @param eemlist object of class eemlist containing EEMs samples
#' @param blanklist object of class eemlist containing the blanks for the EEMs samples, names should be same as eemlist with "_blank" at the end
#' @param abs dataframe containing absorbance data corresponding to the EEMs samples
#' @export
#' @returns a list of metadata, eemlist, blanklist, and absorbance
check_samps <- function(meta, eemlist, blanklist, abs){
  stopifnot("unique_ID" %in% colnames(meta)| is.data.frame(meta),
            is.data.frame(abs) | .is_eem(eemlist) |
              .is_eemlist(blanklist))
    for(m in meta$index){
      #gets name of files that should be there
      x <- which(meta$index == m)
      name <- meta$unique_ID[x]

      eem_check <- name %in% eem_names(eemlist)
      blank_check <- paste(name, "_blank",sep="") %in% eem_names(blanklist)
      abs_check <- name %in% colnames(abs)

      if(eem_check == F | blank_check == F | abs_check ==F){
        warning(paste("Sample", name, "wasn't found in the loaded EEMs and absorbance. Removing from metadata"))
        meta <- meta[-x,]}
    }

    #checks if there's samples in EEMs that aren't in metadata
    if(length(eemlist) > nrow(meta)){
        missing_meta <- eem_names(eemlist)[!(eem_names(eemlist) %in% meta$unique_ID)]
        eemlist <- eem_exclude(eemlist, exclude = list(sample=missing_meta))
        blanklist <- eem_exclude(blanklist, exclude = list(sample=paste(missing_meta, "_blank", sep="")))
        abs <- abs[,-(colnames(abs) %in% missing_meta)]
        cat("Samples", paste(missing_meta, sep="\n"), "are not in metadata and will be removed", sep="\n")
    }

    #checks if absorbance matches eems
    if((ncol(abs)-1) != length(eemlist) | length(eemlist) != length(blanklist)){
      stop("mismatch occuring between absorbance and EEMs data, please check sample names")
    }

    list(meta, eemlist, blanklist, abs)
  }



