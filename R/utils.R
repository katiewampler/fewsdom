utils::globalVariables(c("%dopar%","i", ".",
                         "...","wavelength", "em", "z", "ex",
                         "x","y", "data_identifier", "replicate_no", "meta","data_process",
                         "meta_clean","raw_eem"))


#' Create file structure
#'
#' Takes a file with EEMs and absorbance data and creates subfolders to store
#' raw and processed data, plots, and output tables
#'
#' @param prjpath a string indicating the project file with the data
#' @export
#'

create_files <- function(prjpath){
  stopifnot(is.character(prjpath)|file.exists(prjpath))
  #put in files
  dir.create(paste(prjpath, "/1_Absorbance", sep=""), showWarnings = F)
  dir.create(paste(prjpath, "/2_Blanks", sep=""), showWarnings = F)
  dir.create(paste(prjpath, "/3_Samples", sep=""), showWarnings = F)
  dir.create(paste(prjpath, "/4_Clean_Absorbance", sep=""), showWarnings = F) #create file for clean files if it doesn't exist yet
  dir.create(paste(prjpath, "/5_Processed", sep=""), showWarnings = F)
  dir.create(paste(prjpath,  "/5_Processed/Figures", sep=""), showWarnings = F)
}

#' Checks if object is an eemlist
#'
#' @param eem an object
#' @noRd
.is_eemlist <- function(eem) {
  ifelse(class(eem) == "eemlist", TRUE, FALSE)
}

#' Checks if object is an eem
#'
#' @param eem an object
#' @noRd
.is_eem <- function(eem) {
  ifelse(class(eem) == "eem", TRUE, FALSE)
}

#' DOC Normalizes eems data
#'
#' @param eem object of class eem or eemlist
#' @param meta dataframe with metadata, needs to have doc data
#' @noRd
#'
.eem_doc_norm <- function(eem, meta){
  #remove EEMs with no DOC data (or value of 0)
  EEM_rm <- meta$unique_ID[meta$DOC_mg_L == 0 | is.na(meta$DOC_mg_L) == T]
  eem <- eem_exclude(eem,
                   exclude=list("ex"=c(), "em"=c(),"sample"= EEM_rm ))

  #check if they've been normalized
  doc_done <- sapply(eem, function(x) attr(x, "is_doc_normalized"))
  to_norm <- which(doc_done == F)

  #normalize those that haven't
  for(x in to_norm){
    eem_name <- eem[[x]]$sample
    eem_index <- which(eem_name == meta$unique_ID)
    eem[[x]]$x <- eem[[x]]$x / as.numeric(meta$DOC_mg_L[eem_index])
    attr(eem[[x]], "is_doc_normalized") <- TRUE
  }

  #double check all are normalized
  doc_done <- sapply(eem, function(x) attr(x, "is_doc_normalized"))
  stopifnot(sum(doc_done==F)==0)
  return(eem)
}

#' Removes DOC normalization on eems data
#'
#' Removes DOC normalization if sample are DOC normalized
#' @param eem object of class eem or eemlist
#' @param meta dataframe with metadata, needs to have doc data
#' @noRd
.eem_doc_rm <- function(eem, meta){
  #check if they've been normalized
  doc_done <- sapply(eem, function(x) attr(x, "is_doc_normalized"))
  to_unnorm <- which(doc_done == T)

  #unnormalize those that haven't
  for(x in to_unnorm){
    eem_name <- eem[[x]]$sample
    eem_index <- which(eem_name == meta$unique_ID)
    eem[[x]]$x <- eem[[x]]$x  * as.numeric(meta$DOC_mg_L[eem_index])
    attr(eem[[x]], "is_doc_normalized") <- FALSE
  }

  #double check all are normalized
  doc_done <- sapply(eem, function(x) attr(x, "is_doc_normalized"))
  stopifnot(sum(doc_done==T)==0)
  return(eem)
}

#' Save files to excel
#'
#' Saves files to excel based on named sheet name
#'
#' @param file name of excel file
#' @param sheet_name sheet name of excel file
#' @param df dataframe to save to excel file
#' @param sampsascol a logical indicating how results should be oriented, TRUE puts samples as columns, FALSE puts samples as rows
#' @noRd
#'
.OSU_excel <- function(file, sheet_name, df, sampsascol){
  wb <- openxlsx::loadWorkbook(file)
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb,sheet_name, df, rowNames = sampsascol)
  openxlsx::saveWorkbook(wb,file,overwrite = TRUE)
}


#' Find difference between EEMs samples
#'
#' @param eemlist1 object of class eemlist
#' @param eemlist2 object of class eemlist
#' @param sampnum numeric, the number of the sample in the eemlist you want to compare.
#'
#' @return a dataframe with the intensity at 10 different points for both EEMs
#' @noRd
samp_dif <- function(eemlist1, eemlist2, sampnum){
  X1 <- eemlist1[[sampnum]]
  X2 <- eemlist2[[sampnum]]

  n_dif <- sum(!(X1$x == X2$x), na.rm=T)
  cat("there are", n_dif, "samples that are not the same \n")

  set.seed(9)
  ex_samp <- sample(1:length(X1$ex), 10)
  em_samp <- sample(1:length(X1$em), 10)
  ints1 <- sapply(1:10, function(x){
    ex <- ex_samp[x]
    em <- em_samp[x]
    val <- X1$x[em, ex]
    val
  })
  ints2 <- sapply(1:10, function(x){
    ex <- ex_samp[x]
    em <- em_samp[x]
    val <- X2$x[em, ex]
    val
  })

  rand_samp <- data.frame(ex = X1$ex[ex_samp], em=X1$em[em_samp],
                          int1 = ints1, int2 = ints2)

  return(rand_samp)
}


#' Check for EEMs that are empty
#'
#' @param eem object of class eem or eemlist
#' @param verbose if you want to print out which samples aren't empty
#' @returns a vector of EEMs with empty EEMs tables (eem$x)
#' @export
#'
empty_eems <- function(eem, verbose=T){
  stopifnot(.is_eem(eem) | .is_eemlist(eem))

  if (.is_eemlist(eem)) {
    res <- sapply(eem, empty_eems, verbose=verbose)
    res <- unlist(res)
    return(res)
  }
  #check if EEM has rows and columns and check's they're not all NA
  blank_check <- F
  if(is.null(dim(eem$x)) == T){blank_check <- T}
  if(nrow(eem$x) == 0 | ncol(eem$x) == 0){blank_check <- T}
  if(sum(is.na(eem$x))== (dim(eem$x)[1]*dim(eem$x)[2])){blank_check <- T}


  if(blank_check == T){
    if(verbose==T){
      cat("Sample", eem$sample, "has no EEMs data \n")
    }
    outcome <- eem$sample
  }else{
    if(verbose == T){
      cat("There were no empty EEMs samples \n")
    }
    outcome <- c()
  }
  return(outcome)
}
