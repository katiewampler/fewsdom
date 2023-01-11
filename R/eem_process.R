#' Process EEM's samples
#'
#' Takes EEMs data with blanks and absorbance data and will perform processing
#' steps including blank subtraction, removing raman and rayleigh scattering,
#' dilution correction, removing inner filtering effects, raman normalization
#'
#' @importFrom eemR eem_names
#'
#' @param prjpath The file path of the main file directory where data is located
#' @param eemlist object of class eem or eemlist containing eem's samples
#' @param blanklist object of class eem or eemlist containing the blanks for the eem's samples should be the same length as eem
#' @param abs dataframe containing absorbance data corresponding to the eem's samples
#' @param process_file logical, if TRUE it will put a text file in the processed data folder named 'processing_tracking'
#' @param meta dataframe of metadata containing unique ID's, integration time, dilutions, and raman area for each sample, see example metadata for format
#' @param replace_blank logical, if TRUE it will find the first sample labeled "blank" or "blk" and use that for blank subtraction, use when instrument blank has errors
#' @param raman logical, if TRUE will use 'raman' function to remove raman scattering
#' @param rayleigh logical, if TRUE will use 'rayleigh' function to remove rayleigh scattering
#' @param IFE logical, if TRUE will use absorbance data to remove inner filter effects
#' @param raman_norm logical, if TRUE will normalize EEM's to raman unit area
#' @param dilute logical, if TRUE will correct for dilution factors given in metadata table
#' @param ex_clip vector of length two with the excitation wavelengths to clip the EEM's to
#' @param em_clip vector of length two with the emission wavelengths to clip the EEM's to
#' @param doc_norm logical, if TRUE will normalize EEM's to 1 mg/L carbon, will remove any samples without DOC concentrations in the metadata
#' @param ... arguments passed on to scattering functions 'raman' and 'rayleigh'
#' @return an list where the first object is of class eemlist with processed EEM's samples. The
#' second object is a dataframe with the processes absorbance data.
#' If a process file is given, a file will be created in the processes folder of the file directory
#' @export

eem_proccess <- function(prjpath, eemlist, blanklist, abs,
                         meta, process_file=T, replace_blank=F,
                         raman=T, rayleigh=T, IFE=T, raman_norm=T,
                         dilute = T, em_clip = c(247,550),
                         ex_clip = c(247,450), doc_norm=T, ...){
  #create text file to track processing changes
  #create folders for finished data
  dir.create(paste(prjpath, "/5_Processed", sep=""), showWarnings = F)

  if(process_file==T){
    process_file_name <- paste(prjpath, "/5_Processed/processing_tracking.txt", sep="")
    file.create(process_file_name)
    write.table(paste("PROCESSING STEPS ON ", Sys.Date(), sep=""), process_file_name, append=F, quote=F, row.names = F, col.names = F)
  }

  #subtract blank
  X_sub <- eemlist
  for(n in 1:length(eemlist)){
    eem <- eemlist[[n]]
    if(replace_blank == T){
      blank <- eemlist[[which(stringr::str_detect(meta$unique_ID, "BLK|blk|blank|blank|BLANK") == T)[1]]]
      write.table(paste(Sys.time(), "- instrument blank was subsituted for the first run blank", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)
    } else{
      blank <- blanklist[[n]]
    }
    eem_blank <- eem$x - blank$x
    X_sub[[n]]$x <- eem_blank
  }
  write.table(paste(Sys.time(), "- blanks were subtracted from samples", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)

  #remove raman scattering
  X_mask <- X_sub
  if(raman==T){
    X_mask <- raman(X_mask, process_file=process_file_name, ...)
  }

  #remove rayleigh scattering
  if(rayleigh ==T){
    X_mask <- rayleigh(X_mask, process_file=process_file_name, ...)
  }

  #remove inner filtering effects
  X_ife <- X_mask
  if(IFE==T){
    #clip eem to be able to remove inner filtering effects
    X_ife <- eemR::eem_cut(X_mask, em=600:1000, ex=600:1000, exact=F)
    write.table(paste(Sys.time(), "- EEM's were clipped to just emission wavelengths under 600 nm", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)

    X_ife <- eem_inner_filter_effect(X_ife, absorbance=abs)
    write.table(paste(Sys.time(), "- EEM's were corrected for inner filter effects", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)
  }

  #perform raman normalization
  X_norm <- X_ife
  if(raman_norm == T){
    meta$raman_norm <- meta$integration_time_s * meta$RSU_area_1s

    for(f in 1:length(X_norm)){
      eem_name <- X_norm[[f]]$sample
      eem_index <- which(eem_name == meta$unique_ID)
      X_norm[[f]]$x <- X_ife[[f]]$x / meta$raman_norm[eem_index]}
    write.table(paste(Sys.time(), "- EEM's were normalized by raman area", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)}

  #account for dilutions
  X_dil_cor <- X_norm
  Sabs_dil_cor <- abs
  if(dilute==T & sum(meta$dilution != 1) > 0){
    dil_data <- 1/meta["dilution"] #invert to match function
    X_dil_cor <- X_dil_cor %>% eem_dilution(dil_data)
    write.table(paste(Sys.time(), "- EEM's were corrected for dilutions", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)

    #correct absorbance for dilution
    for(x in 1:nrow(meta)){
      dil_fact <- meta$dilution[x]
      Sabs_dil_cor[,meta$abs_col[x]] <- abs[,meta$abs_col[x]] / dil_fact}
    write.table(paste(Sys.time(), "- Absorbance was corrected for dilutions", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)
  }

  #normalize for DOC
  X_DOC <- X_dil_cor
  Sabs_DOC <- Sabs_dil_cor
  if(doc_norm == T){
    #remove EEMs with no DOC data (or value of 0)
    EEM_rm <- meta$unique_ID[meta$DOC_mg_L == 0 | is.na(meta$DOC_mg_L) == T]
    X_DOC <- eem_exclude(X_DOC,
                         exclude=list("ex"=c(), "em"=c(),"sample"= EEM_rm ))

    #remove absorbance with no DOC data (or value of 0)
    abs_rm <- which((colnames(Sabs_DOC) %in% EEM_rm) == T)
    if(length(abs_rm) > 0){
      Sabs_DOC <-  Sabs_DOC[,-(abs_rm)]
    }

    #EEMs
    for (x in 1:length(X_DOC)){
      eem_name <- X_DOC[[x]]$sample
      eem_index <- which(eem_name == meta$unique_ID)
      X_DOC[[x]]$x <- X_DOC[[x]]$x / as.numeric(meta$DOC_mg_L[eem_index])
    }

    #Absorbance
    for(x in colnames(Sabs_DOC)[2:ncol(Sabs_DOC)]){
      col_num <- which(x == colnames(Sabs_DOC))
      abs_index <- which(x == meta$unique_ID)
      doc <- as.numeric(meta$DOC_mg_L[abs_index])
      Sabs_DOC[,col_num] <- Sabs_DOC[,col_num] / doc}
    write.table(paste(Sys.time(), "- Absorbance and EEM's were normalized by DOC concentration", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)

  }

  #clip DOC normalized EEM's
  X_DOC_clip <- X_DOC
  X_DOC_clip<- eemR::eem_cut(X_DOC_clip, ex=c(min(X_DOC_clip[[1]]$ex):ex_clip[1], ex_clip[2]:max(X_DOC_clip[[1]]$ex),exact=F ))
  X_DOC_clip <- eemR::eem_cut(X_DOC_clip, em=c(min(X_DOC_clip[[1]]$em):em_clip[1], em_clip[2]:max(X_DOC_clip[[1]]$em),exact=F ))
  write.table(paste(Sys.time(), "- DOC Normalized EEM's were clipped to Excitation:",ex_clip[1]," to ", ex_clip[2],
                    " nm and Emission:",em_clip[1]," to ", em_clip[2], " nm", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)

 return(list(X_DOC_clip, Sabs_DOC))
}

#' Modifying fluorescence data according to dilution.
#'
#' Exact function from 'staRdom' package, extracted to prevent orphan package warnings.
#' If samples were diluted before measuring, a dilution factor has to be added to
#' the measured data. This function can do that by either multiplying each sample with
#' the same value or using a data frame with different values for each sample.
#'
#' @param data fluorescence data with class eemlist
#' @param dilution dilution factor(s), either numeric value or data frame. Row names of data frame have to be similar to sample names in eemlist.
#' @noRd
eem_dilution <- function (data, dilution = 1)
{
  if (((is.numeric(dilution) & length(dilution) == 1) | is.data.frame(dilution)) &
      inherits(data, "eemlist")) {
    res_list <- lapply(1:length(data), function(i) {
      if (is.data.frame(dilution)) {
        no_dil <- which(!(eem_names(data) %in% row.names(dilution)))
        if (length(no_dil) > 0)
          stop("Dilution factors are missing for the following samples: ",
               paste0(eem_names(data)[no_dil], collapse = ", "),
               "!")
        data[[i]]$x <- data[[i]]$x * dilution[data[[i]]$sample,
        ]
      }
      else {
        data[[i]]$x <- data[[i]]$x * dilution
      }
      data[[i]]
    })
    class(res_list) <- class(data)
  }
  else {
    stop("First argument must be of class eemlist, second argument must be a number or a data.frame to multiply with for dilution.")
  }
  res_list
}

#' Inner-filter effect correction
#'
#' Inner-filter effect correction function modified from 'eemR' to report
#' less text for putting into other functions
#'
#' @param eem An object of class eemlist.
#' @param absorbance A dataframe with aborbance data
#' @param pathlength A numeric value indicating the pathlength (in cm) of the cuvette used for absorbance measurement. Default is 1 (1cm).
#' @param verbose a logical, if TRUE will print out detailed information on filter effects, if false will only print those samples were dilution is reccomended
#' @noRd
eem_inner_filter_effect <- function (eem, absorbance, pathlength = 1, verbose=F)
{
  .is_eemlist <- function(eem) {
    ifelse(class(eem) == "eemlist", TRUE, FALSE)
  }
  .is_eem <- function(eem) {
    ifelse(class(eem) == "eem", TRUE, FALSE)
  }
  stopifnot(.is_eemlist(eem) | .is_eem(eem), is.data.frame(absorbance),
            is.numeric(pathlength))
  if (.is_eemlist(eem)) {
    res <- lapply(eem, eem_inner_filter_effect, absorbance = absorbance,
                  pathlength = pathlength, verbose=verbose)
    class(res) <- class(eem)
    return(res)
  }
  if (!any(names(absorbance) == "wavelength")) {
    stop("'wavelength' variable was not found in the data frame.",
         call. = FALSE)
  }
  wl <- absorbance[["wavelength"]]
  if (!all(is_between(range(eem$em), min(wl), max(wl)))) {
    stop("absorbance wavelengths are not in the range of\n         emission wavelengths",
         call. = FALSE)
  }
  if (!all(is_between(range(eem$ex), min(wl), max(wl)))) {
    stop("absorbance wavelengths are not in the range of\n         excitation wavelengths",
         call. = FALSE)
  }
  index <- which(names(absorbance) == eem$sample)
  if (length(index) == 0) {
    warning("Absorbance spectrum for ", eem$sample, " was not found. Returning uncorrected EEM.",
            call. = FALSE)
    return(eem)
  }
  spectra <- absorbance[[index]]
  if (attributes(eem)$is_ife_corrected) {
    return(eem)
  }
  sf <- stats::splinefun(wl, spectra)
  ex <- sf(eem$ex)
  em <- sf(eem$em)
  total_absorbance <- sapply(ex, function(x) {
    x + em
  })/pathlength
  max_abs <- max(total_absorbance)
  if(verbose == T){
    cat(eem$sample, "\n")
  }
  if (max_abs > 1.5 ) {
    cat(eem$sample, "\n")
    cat("Total absorbance is > 1.5 (Atotal = ", max_abs,
        ")\n", "A 2-fold dilution is recommended. See ?eem_inner_filter_effect.\n",
        sep = "")
  }
  ife_correction_factor <- 10^(0.5 * total_absorbance)
  if(verbose ==T | max_abs > 1.5){
    cat("Range of IFE correction factors:", round(range(ife_correction_factor),
                                                  digits = 4), "\n")
    cat("Range of total absorbance (Atotal) :", round(range(total_absorbance/pathlength),
                                                      digits = 4), "\n\n")
  }
  res <- eem
  res$x <- eem$x * ife_correction_factor

  attr(res, "is_ife_corrected") <- TRUE
  return(res)
}
