#' Load renamed and organized Aqualog data into R
#'
#' Takes the Aqualog files that have been renamed and put in separate folders within
#' the main directory via the 'clean_files' function and loads them into R as
#' eem lists. It will also modify the metadata table to link the sample, blank, and
#' absorbance data.
#'
#' @importFrom utils zip
#' @importFrom eemR eem_read
#'
#' @param prjpath The file path of the main file directory where data is located
#' @param
#' @export
#'

load_eems <- function(prjpath){
  #load in EEM's data
  X <- eemR::eem_read(paste(prjpath, "/3_Samples/",sep=""), recursive = F, import_function = "aqualog")
  X_blk <- eemR::eem_read(paste(prjpath, "/2_Blanks/",sep=""), recursive = F, import_function = "aqualog")

  #load in absorbance data
  Sabs <- absorbance_read(paste(prjpath, "/4_Clean_Absorbance", sep=""))

  #match EEM and absorbance files with information from metadata
  for(m in meta$index){
    x <- which(m == meta$index)
    absname <- paste(meta$unique_ID[x], "_Abs", sep="")
    blankname <- paste(meta$unique_ID[x], "_blank", sep="")
    abs_col <- which(meta$unique_ID[x] == colnames(Sabs))
    if(length(abs_col) == 0){
      #see if was missing 0 before site number
      unique_ID <- paste(str_sub(meta$unique_ID[x], 1,2), str_sub(meta$unique_ID[x], 4,-1), sep="")
      absname <- paste(unique_ID, "_Abs", sep="")
      blankname <- paste(unique_ID, "_blank", sep="")
      abs_col <- which(unique_ID == colnames(Sabs))
      if(length(abs_col) == 0){
        warning(paste("Sample", absname, "wasn't found in the file. Removing from metadata"))
        meta <- meta[-x,]
      } else{
        new_unique_ID <- meta$unique_ID[x]
        X <- eem_name_replace(X, unique_ID, new_unique_ID)
        X_blk <- eem_name_replace(X_blk, unique_ID, new_unique_ID)
        colnames(Sabs)[colnames(Sabs)==unique_ID] <- new_unique_ID

        df <- data.frame(absname = absname, blankname=blankname , eemname=meta$unique_ID[x], abs_col=abs_col)
        if(x==1){
          abs_index <- df
        } else{abs_index <- rbind(abs_index, df)}
      }} else{
        df <- data.frame(absname = absname, blankname=blankname , eemname=meta$unique_ID[x], abs_col=abs_col)
        if(x==1){
          abs_index <- df
        } else{abs_index <- rbind(abs_index, df)}
      }

  }

  #merge with metadata
  meta <- merge(meta, abs_index, by.x="unique_ID", by.y="eemname")
  row.names(meta) <- meta$unique_ID

  #check absorbance (if over 0.3) will return warning for samples above, else will just run in background
  for(c in 2:ncol(Sabs)){
    samp <- Sabs[,c]
    max <- max(samp)
    if(max > 0.3){
      warning(paste("Sample",colnames(Sabs)[c], "has an aborbance of", round(max,3), "please dillute", sep=" "))
    }
  }
}


#' Reading absorbance data from txt and csv files
#'
#' Exact function from 'staRdom' package, extracted to prevent orphan package warnings.
#' Reading absorbance data from txt and csv files.
#'
#' @importFrom utils zip
#'
#' @param absorbance_path directory containing absorbance data files or path to single file. See details for format of absorbance data.
#' @param order logical, data is ordered according to wavelength
#' @param recursive read files recursive, include subfolders
#' @param dec optional, either you set a decimal separator or the table is tested for . and ,
#' @param sep optional, either you set a field separator or it is tried to be determined automatically
#' @param verbose logical, provide more information
#' @param cores number of CPU cores to be used simultanuously
#'
absorbance_read <- function (absorbance_path, order = TRUE, recursive = TRUE, dec = NULL,
                         sep = NULL, verbose = FALSE, cores = parallel::detectCores(logical = FALSE),
                         ...)
{
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
  clusterExport(cl, c("dec", "sep", "verbose"), envir = environment())
  clusterEvalQ(cl, require(dplyr))
  clusterEvalQ(cl, require(stringr))
  abs_data <- parLapply(cl, abs_data, function(tab) {
    tryCatch({
      rawdata <- readLines(tab)
      data <- rawdata %>% sapply(str_remove, pattern = "([^0-9]*$)")
      first_number <- min(which((substr(data, 1, 1) %>%
                                   grepl("[0-9]", .))))
      last_number <- max(which((substr(data, 1, 1) %>%
                                  grepl("[0-9]", .))))
      if (is.null(sep) | is.null(dec)) {
        nsepdec <- data[first_number] %>% str_extract_all("[^-0-9eE]") %>%
          unlist()
        example_number <- data[first_number] %>% str_extract("([-]?[0-9]+[.,]?[0-9]+[eE]?[-0-9]+)$")
        if (is.null(dec) & length(nsepdec) > 1)
          dec <- example_number %>% str_replace("([-0-9eE]+)([.,]?)([-0-9eE]*)",
                                                "\\2")
        if (is.null(sep))
          sep <- gsub(pattern = dec, replacement = "",
                      x = data[first_number], fixed = TRUE) %>%
          str_extract(paste0("[^-0-9eE", dec, "]"))
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
        samples <- tab %>% basename() %>% str_replace_all(regex(".txt$|.csv$",
                                                                ignore_case = TRUE), "")
      }
      else {
        samples <- rawdata[[1]] %>% str_split(sep) %>%
          unlist() %>% matrix(ncol = length(.), byrow = TRUE) %>%
          data.frame(stringsAsFactors = FALSE) %>% .[-1]
      }
      table <- table %>% setNames(c("wavelength", samples))
    }, error = function(err) {
      stop("Error while reading ", tab, ": ", err)
    })
  })
  stopCluster(cl)
  locations <- lapply(abs_data, function(tab) {
    attr(tab, "location")
  }) %>% unlist()
  if (length(abs_data) == 1)
    abs_data <- abs_data[[1]] %>% as.data.frame()
  else abs_data <- abs_data %>% list_join(by = "wavelength")
  if (order)
    abs_data <- abs_data %>% arrange(wavelength)
  attr(abs_data, "location") <- locations
  abs_data
}
