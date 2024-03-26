#' Process EEMs samples
#'
#' Takes EEMs data with blanks and absorbance data and will perform processing
#' steps including blank subtraction, removing raman and rayleigh scattering,
#' dilution correction, removing inner filtering effects, raman normalization.
#'
#' @importFrom eemR eem_names
#'
#' @param prjpath the file path of the main file directory where data is located
#' @param eemlist object of class eem or eemlist containing EEMs samples
#' @param blanklist object of class eem or eemlist containing the blanks for the EEMs samples should be the same length as eem
#' @param abs dataframe containing absorbance data corresponding to the EEMs samples
#' @param meta dataframe of metadata containing unique ID's, integration time, dilutions, and raman area for each sample, see example metadata for format
#' @param process_file logical, if TRUE it will put a text file in the processed data folder named 'processing_tracking'
#' @param replace_blank logical, if TRUE it will find the first sample labeled "blank" or "blk" and use that for blank subtraction, use when instrument blank has errors
#' @param raman logical, if TRUE will use 'raman' function to remove raman scattering
#' @param rayleigh logical, if TRUE will use 'rayleigh' function to remove rayleigh scattering
#' @param IFE logical, if TRUE will use absorbance data to remove inner filter effects
#' @param raman_norm logical, if TRUE will normalize EEMs to raman unit area
#' @param dilute logical, if TRUE will correct for dilution factors given in metadata table
#' @param ex_clip vector of length two with the excitation wavelengths to clip the EEMs to
#' @param em_clip vector of length two with the emission wavelengths to clip the EEMs to
#' @param raman_mask a vector of length 4 specifying the width of the raman line to cut, numbers 1:2 are width above and below first order line, numbers 3:4 are width above and below second order line, since you cannot use the auto method for second order raman this must be specified
#' @param raman_width either "auto" or "manual". If auto is chosen cutting widths will be found using the 'find_cut_width' function
#' @param rayleigh_mask optional if auto width method is used, a vector of length 4 specifying the width of the rayleigh line to cut, numbers 1:2 are width above and below first order line, numbers 3:4 are width above and below second order line
#' @param rayleigh_width either "auto" or "manual", if auto is chosen cutting widths will be found using the 'find_cut_width' function
#' @param ... arguments passed on to scattering functions 'raman' and 'rayleigh'
#' @return an list where the first object is of class eemlist with processed EEMs samples. The
#' second object is a dataframe with the processes absorbance data.
#' If a process file is given, a file will be created in the processes folder of the file directory
#' @export
#' @examples
#' \dontrun{
#' data_process <- eem_proccess(prjpath=prjpath, eemlist=X, blanklist=X_blk, abs=Sabs,
#' process_file=process_file, meta=meta)
#' }
#' X_clean <- data_process[[1]]
#' abs_clean <- data_process[[2]]

eem_proccess <- function(prjpath, eemlist, blanklist, abs,
                         meta, process_file=T, replace_blank=F,
                         raman=T, rayleigh=T, IFE=T, raman_norm=T,
                         dilute = T, ex_clip = c(247,450),
                         em_clip = c(247,600),
                         raman_width ="auto", raman_mask= c(8,8,1.5,1.5),
                         rayleigh_width="manual", rayleigh_mask = c(20,10,10,10),
                         ...){

  stopifnot(is.character(prjpath) | .is_eemlist(eemlist) | .is_eem(eemlist) |
              .is_eemlist(blanklist) | .is_eem(blanklist)| is.data.frame(abs)|
              is.logical(process_file)| is.logical(replace_blank)|is.logical(raman)|
              is.logical(rayleigh)|is.logical(IFE)|is.logical(raman_norm)|
              is.logical(dilute)|is.numeric(em_clip)|is.numeric(ex_clip)|file.exists(prjpath))

  #make sure file directory is good
  if(file.exists(paste(prjpath,  "/5_Processed", sep=""))==F){
    stop("Invalid file structure, please use 'create_files' function to create file for plots within file directory")
  }

  #create text file to track processing changes
  if(process_file==T){
    process_file_name <- paste(prjpath, "/5_Processed/processing_tracking.txt", sep="")
    file.create(process_file_name)
    write.table(paste("PROCESSING STEPS ON ", Sys.Date(), sep=""), process_file_name, append=F, quote=F, row.names = F, col.names = F)
  }

  #create attributes to track processing
  eemlist <- lapply(1:length(eemlist), function(i) {
    attr(eemlist[[i]], "is_doc_normalized") <- FALSE
    return(eemlist[[i]])})
  eemlist <- lapply(1:length(eemlist), function(i) {
    attr(eemlist[[i]], "is_dil_corrected") <- FALSE
    return(eemlist[[i]])})
  class(eemlist) <- class(blanklist)
  attr(abs, "is_doc_normalized") <- FALSE
  attr(abs, "is_dil_corrected") <- FALSE

  #subtract blank
  X_sub <- eemlist
  X_sub <- lapply(1:length(X_sub), function(n, replace_blank){
    eem <- eemlist[[n]]
    if(replace_blank == T){
      blank <- eemlist[[which(stringr::str_detect(meta$unique_ID, "BLK|blk|blank|blank|BLANK") == T)[1]]]
    } else{
      blank <- blanklist[[n]]
    }
    X_sub[[n]]$x <- eem$x - blank$x
    attr(X_sub[[n]], "is_blank_corrected") <- TRUE
    return(X_sub[[n]])
  }, replace_blank=replace_blank)
  class(X_sub) <- "eemlist"
  if(length(empty_eems(X_sub, verbose=F)) >0){
    stop("one or more of your EEMs has empty data after blank subtraction, use 'empty_eems' function to find out which ones")
  }
  write.table(paste(Sys.time(), "- blanks were subtracted from samples", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)

  #remove blank that was used as instrument blank so code won't break during plotting
  if(replace_blank ==T){
    write.table(paste(Sys.time(), "- instrument blank was subsituted for the first run blank", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)
    blank <- eemlist[[which(stringr::str_detect(meta$unique_ID, "BLK|blk|blank|blank|BLANK") == T)[1]]]
    X_sub <- eem_exclude(X_sub, exclude=list(ex=c(), em=c(), sample=blank$sample) )
  }

  #remove raman scattering
  X_mask <- X_sub
  if(raman==T){
    X_mask <- raman(X_mask, process_file=process_file_name, raman_width=raman_width,
                    raman_mask = raman_mask,...)
  }

  #remove rayleigh scattering
  if(rayleigh ==T){
    X_mask <- rayleigh(X_mask, process_file=process_file_name,
                       rayleigh_width = rayleigh_width, rayleigh_mask=rayleigh_mask,...)
  }
  if(length(empty_eems(X_mask, verbose=F)) >0){
    stop("one or more of your EEMs has empty data after removing scattering, try setting rayleigh_width to 'manual'")
  }
  #remove inner filtering effects
  X_ife <- X_mask
  if(IFE==T){
    #clip eem to be able to remove inner filtering effects
    X_ife <- eemR::eem_cut(X_mask, em=600:1000, ex=600:1000, exact=F)
    write.table(paste(Sys.time(), "- EEM's were clipped to just emission wavelengths under 600 nm", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)

    X_ife <- eem_inner_filter_effect(X_ife, absorbance=abs, verbose=F)
    write.table(paste(Sys.time(), "- EEM's were corrected for inner filter effects", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)
    if(length(empty_eems(X_ife, verbose=F)) >0){
      stop("one or more of your EEMs has empty data after correcting for inner filter effects, use 'empty_eems' function to find out which ones")
    }
    }

  #perform raman normalization
  X_norm <- X_ife
  if(raman_norm == T){
    meta$raman_norm <- meta$integration_time_s * meta$RSU_area_1s
    X_norm <- lapply(1:length(X_norm), function(f){
      eem_name <- X_norm[[f]]$sample
      eem_index <- which(eem_name == meta$unique_ID)
      X_norm[[f]]$x <- X_ife[[f]]$x / meta$raman_norm[eem_index]
      attr(X_norm[[f]], "is_raman_normalized") <- TRUE
      return(X_norm[[f]])
    })
    class(X_norm) <- "eemlist"
    if(length(empty_eems(X_norm, verbose=F)) >0){
      stop("one or more of your EEMs has empty data after normalizing for raman area, use 'empty_eems' function to find out which ones")
    }
    write.table(paste(Sys.time(), "- EEM's were normalized by raman area", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)}

  #account for dilutions
  X_dil_cor <- X_norm
  Sabs_dil_cor <- abs
  if(dilute == T){
    X_dil_cor <- lapply(1:length(X_dil_cor), function(x){
      attr(X_dil_cor[[x]], "is_dil_corrected") <- TRUE
      return(X_dil_cor[[x]])
    })
    class(X_dil_cor) <- "eemlist"
    attr(abs, "is_dil_corrected") <- TRUE
    write.table(paste(Sys.time(), "- EEM's were corrected for dilutions", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)
    write.table(paste(Sys.time(), "- Absorbance was corrected for dilutions", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)

    if(sum(meta$dilution != 1) > 0){
      dil_data <- meta %>% select(dilution) #invert to match function
      if(mean(unlist(dil_data), na.rm=T) < 1){dil_data <- 1 / dil_data
      warning("Inverted dilutions to be dilution factor")}

      X_dil_cor <- X_dil_cor %>% staRdom::eem_dilution(dil_data)

      #correct absorbance for dilution
      for(x in 1:nrow(meta)){
        dil_fact <- meta$dilution[x]
        if(dil_fact < 1){dil_fact <- 1 / dil_fact }
        Sabs_dil_cor[,meta$abs_col[x]] <- abs[,meta$abs_col[x]] * dil_fact}
    }
  }
  if(length(empty_eems(X_dil_cor, verbose=F)) >0){
    stop("one or more of your EEMs has empty data after correcting for dilutions, use 'empty_eems' function to find out which ones")
  }

  #clip EEM's  #clips un-normalized eem
  X_clip <- X_dil_cor
  X_clip<- eem_cut2(X_clip, ex=ex_clip, em=em_clip, exact=F)
  write.table(paste(Sys.time(), "- EEM's were clipped to Excitation:",ex_clip[1]," to ", ex_clip[2],
                    " nm and Emission:",em_clip[1]," to ", em_clip[2], " nm", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)
  if(length(empty_eems(X_clip, verbose=F)) >0){
    stop("one or more of your EEMs has empty data after clipping, use 'empty_eems' function to find out which ones")
  }
 return(list(X_clip, Sabs_dil_cor))
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
eem_inner_filter_effect <- function(eem, absorbance, pathlength = 1, verbose=F)
{
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

#' Trim EEMs to specified emission and/or excitation wavelengths
#'
#' Modified from eemR 'eem_cut' function. The specified wavelengths are now the
#' ones you keep, not the ones you remove.
#'
#' @param eem An object of class eem or eemlist.
#' @param ex numerical of length two, with range of excitation wavelengths to keep
#' @param em numerical of length two, with range of emission wavelengths to keep
#' @param exact Logical. If TRUE, only wavelengths matching em and/or ex will be removed. If FALSE, all wavelengths in the range of em and/or ex will be removed.
#' @param fill_with_na Logical. If TRUE, fluorescence values at specified wavelengths will be replaced with NA. If FALSE, these values will be removed.
#'
#' @return an object of class eem or eemlist
#' @export
#'
eem_cut2 <- function (eem, ex, em, exact = FALSE, fill_with_na = FALSE) {
  stopifnot(.is_eemlist(eem) | .is_eem(eem))
  if (.is_eemlist(eem)) {
    res <- lapply(eem, eem_cut2, ex = ex, em = em, exact = exact,
                  fill_with_na = fill_with_na)
    class(res) <- class(eem)
    return(res)
  }
  if (!missing(ex)) {
    stopifnot(is.numeric(ex), all(ex >= 0))
    ex_range <- ex[1]:ex[2]
    if (exact) {
      index <- which(!(eem$ex %in% ex_range))
      if(length(index) == length(eem$ex)){
        warning("all excitation wavelengths were removed, try 'exact=F'.")
      }
    }
    else {
      index <- which(!(is_between(eem$ex, min(ex_range, na.rm = TRUE),
                                  max(ex_range, na.rm = TRUE))))
    }
    if (length(index != 0)) {
      if (fill_with_na) {
        eem$x[, index] <- NA
      }
      else {
        eem$ex <- eem$ex[-index]
        eem$x <- eem$x[, -index]
      }
    }
  }
  if (!missing(em)) {
    stopifnot(is.numeric(em), all(em >= 0))
    em_range <- em[1]:em[2]
    if (exact) {
      index <- which(!(eem$em %in% em_range))
      if(length(index) == length(eem$ex)){
        warning("all excitation wavelengths were removed, try 'exact=F'.")
      }
    }
    else {
      index <- which(!(is_between(eem$em, min(em_range, na.rm = TRUE),
                                  max(em_range, na.rm = TRUE))))
    }
    if (length(index != 0)) {
      if (fill_with_na) {
        eem$x[index, ] <- NA
      }
      else {
        eem$em <- eem$em[-index]
        eem$x <- eem$x[-index, ]
      }
    }
  }
  return(eem)
}
