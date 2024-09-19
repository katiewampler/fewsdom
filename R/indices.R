#' Cleanly transpose data
#'
#' Takes a dataframe and transposes it (turn the rows to columns and columns to rows).
#'
#' @param df a dataframe that you want to transpose
#' @returns a transposed dataframe.
#' @noRd

#functions for getting indices
clean_transpose <- function(df){
  df_flip <- as.data.frame(t(df))
  df_names <- df_flip[1,]
  df_flip <- df_flip[2:nrow(df_flip),]
  colnames(df_flip) <- df_names
  df_flip
}

#' Get absorbance peaks and indices
#'
#' Takes a dataframe of absorbance data and calculates optical properties for
#' each sample. If data was collected at greater than 1 nm intervals, data will
#' be interpolated to 1 nm intervals.
#'
#' Absorbance indices based on Hansen et al. 2016. Measurements are defined as follows:
#'
#' SUVA254, SUVA280, SUVA350, SUVA370: SUVA at 254, 280, 350, and 370 nm.
#' Units of \eqn{\text{L mgC}^{-1} \text{m}^{-1}}.
#' Typically higher values are associated with greater aromatic content.
#'
#' SVA412, SVA440, SVA480, SVA510,
#' SVA532, SVA555: SVA at 412, 440, 480, 510, 532, 555 nm.
#' Units of \eqn{\text{L mgC}^{-1} \text{m}^{-1}}.
#' Typically higher values are associated with greater aromatic content.
#'
#' S275_295: Spectral slope between 275 to 295 nm.
#'
#' S350_400: Spectral slope between 350 to 400 nm.
#'
#' Spectral slopes are found with a nonlinear fit of an exponential function to
#' the absorption spectrum, typically higher  values are associated with
#' lower molecular weight materials and/or lower aromaticity.
#'
#' SR: Spectral slope S275_295 divided by spectral slope S350_400, negatively correlated to DOM molecular weight
#' and generally increases on irradiation.
#'
#' @importFrom stats approx
#' @param abs_data a dataframe containing the absorbance data, where each column is a sample
#' @param waves optional, a vector of wavelengths in nm to extract from the absorbance data
#' @param meta the metadata table for the sample run, should include DOC data in mg/L
#' @param keep_all a logical, TRUE will return all samples even those without DOC data, FALSE will only return samples with DOC data
#' @param cuvle the length of the cuvette in cm
#' @returns Returns a data frame of indices.
#' @references Hansen, A. M., Kraus, T. E. C., Pellerin, B. A., Fleck, J. A., Downing, B. D., & Bergamaschi, B. A. (2016). Optical properties of dissolved organic matter (DOM): Effects of biological and photolytic degradation. Limnology and Oceanography, 61(3), 1015–1032. https://doi.org/10.1002/lno.10270
#' @export

abs_parm <- function(abs_data, waves=NULL, meta, keep_all=F, cuvle = 1){
  stopifnot(is.data.frame(c(abs_data, meta)) | is.numeric(r_thresh)|is.logical(keep_all))

  #interpolate data to find values that don't fall on the excitation wavelengths used
  data <- lapply(abs_data[,-1], function(x) approx(abs_data$wavelength, x, xout=c(249:791))[[2]])
  data <- as.data.frame(data)
  data$wavelength <- 249:791

  #build data table to put in data
  abs_out <- as.data.frame(matrix(nrow=(ncol(abs_data)-1), ncol=(15 + length(waves))))
  colnam <- c("sample", paste("SUVA", c(254,280, 350, 370), sep=""), paste("SVA", c(412,440,480,510,532,555), sep=""), "S275_295", "S350_400", "SR")
  if(length(waves) > 0){
    colnam <- c(colnam, paste("a", waves, sep=""))
  }
  colnames(abs_out) <- colnam
  abs_out$sample <- colnames(abs_data)[2:ncol(abs_data)]

  #add DOC data to table for specific values
  abs_out$DOC_mgL <- NA
  for(m in meta$index){
    x <- which(m == meta$index)
    abs_col <- which(meta$unique_ID[x] == abs_out$sample)
    abs_out$DOC_mgL[abs_col] <- meta$DOC_mg_L[x]}

  #extract data for suva and sva measurements
  for(x in c(254,280, 350, 370, 412,440,480,510,532,555)){
    row_num <- which(data$wavelength == x) #row in absorbance data
    col_num <- which(stringr::str_detect(colnames(abs_out), as.character(x))==T)[1] #location in output table
    abs_out[,col_num] <- as.numeric(data[row_num, 1:(ncol(data)-1)]) / abs_out$DOC_mgL * 100 #put data in column
  }

  #extract data for added absorbance wavelengths
  for(x in waves){
    row_num <- which(data$wavelength == x) #row in absorbance data
    col_num <- which(stringr::str_detect(colnames(abs_out), paste("a",as.character(x), sep=""))==T) #location in output table
    abs_out[,col_num] <- as.numeric(data[row_num, 1:(ncol(data)-1)]) #put data in column
  }

  #get noise for ratios (will give NA if values aren't above this noise to prevent giving nonsense ratios)
  noise <- subset(data, data$wavelength >= 700)
  noise <- colSums(noise) / nrow(noise)
  noise <- mean(noise[1:(length(noise)-1)])
  noise_val <- noise * 5 #values need to be greater than 5 times the noise to be recorded

  #calculate other things (not currently included)
  #abs_out$E2_E3 <- abs_out$a250/abs_out$a365
  #abs_out$E4_E6 <- abs_out$a465/abs_out$a665

    #check if it's calculating noise
    #abs_out$E2_E3[abs_out$a250 <= noise_val] <- NA
    #abs_out$E2_E3[abs_out$a365 <= noise_val] <- NA
    #abs_out$E4_E6[abs_out$a465 <= noise_val] <- NA
    #abs_out$E4_E6[abs_out$a665 <= noise_val] <- NA

  #get ratios
  S275_295 <- sapply(lapply(data[,-ncol(data)]*log(10)*100/cuvle, staRdom:::abs_fit_slope,
                            wl=data$wavelength, lim=c(275, 295),
                            l_ref=275), function(res) res$coefficients)

  S350_400 <- sapply(lapply(data[,-ncol(data)]*log(10)*100/cuvle, staRdom:::abs_fit_slope,
                            wl=data$wavelength, lim=c(350, 400),
                            l_ref=350), function(res) res$coefficients)


  abs_out$S275_295 <- as.numeric(S275_295)
  abs_out$S350_400 <- as.numeric(S350_400)

  abs_out$SR <- abs_out$S275_295 / abs_out$S350_400

  #remove samples missing DOC
  if(keep_all == F){
    abs_out <- abs_out[is.na(abs_out$DOC_mgL) ==F,]
  }
  abs_out

}

#' Get Coble Peaks and other fluorescence indicies
#'
#' Modified from 'eem_coble_peaks' from eemR package to include some additional indices and
#' modified definitions of Coble peaks.  For indices that are ratios, the values in the ratio
#' must be 5 times (or specifed amount) greater than the area below the 1:1 line (considered noise). If the
#' desired wavelength is not explicitly measured the EEM will be interpolated using
#' the interp2 function from the pracma library.
#'
#' Coble peaks are based on Coble et al. 2014 and are defined as follows:
#'
#' Peak B (pB): ex = 270:280 nm, em = 300:320 nm, Tyrosine-like
#'
#' Peak T (pT): ex = 270:280 nm, em = 320:350 nm, Tryptophan-like
#'
#' Peak A (pA): ex = 250:260 nm, em = 380:480 nm, Humic-like.
#'
#' Peak M (pM): ex = 310:320 nm, em = 380:420 nm, Marine humic-like
#'
#' Peak C (pC): ex = 330:350 nm, em = 420:480 nm, Humic-like
#'
#' Peak D (pD): ex = 390 nm, em = 509 nm, Soil fulvic acid
#'
#' Peak E (pE): ex = 455 nm, em = 521 nm, Soil fulvic acid
#'
#' Peak N (pN): ex = 280 nm, em = 370 nm, Plankton derived
#'
#' Given that peaks A, B, C, M, and T are not defined at fixed excitation and emission wavelength, the maximum fluorescence value in the region is extracted.
#'
#'
#' Additional fluorescence indices are based on Hansen et al. 2016.
#' Measurements are defined as follows:
#'
#' rAT: The ratio of peak A to peak T, indication of the amount of humic like (recalcitrant) to fresh (liable) DOM.
#'
#' rCA: The ratio of peak C to peak A, indication of the amount of humic like to fumic like DOM.
#'
#' rCM: The ratio of peak C to peak M, indication of the amount of diagenetically altered (blueshifted) DOM.
#'
#' rCT: The ratio of peak C to peak T, indication of the amount of humic like (recalcitrant) to fresh (liable) DOM.
#'
#' Fluorescence Index (FI): Ratio of fluorescence at ex = 370 nm, em = 470 nm to em = 520 nm.
#' Identifies the relative contributions of terrestrial to microbial DOM sources.
#'
#' Max Wavelength of FI (FI_max): The location (wavelength) of the maxiumum emission peak
#' at ex= 370 nm. This is typically around 470 nm for natural materials however recent work
#' has shown that this can change considerably for pyrogenic organic matter (Egan et al. 2023).
#'
#' Humification Index (HIX): ex = 254 nm, em =\eqn{\sum}435:480 divided by em =\eqn{\sum}300:345.
#' HIX proposed by Zsolnay (1999).
#' An indication of humic substances or extent of humification. Higher values indicate an higher degree of humification.
#'
#' Humification Index (HIX_ohno): ex = 254 nm, em =\eqn{\sum}435:480 divided by em =\eqn{\sum}435:480, \eqn{\sum}300:345.
#' HIX proposed by Ohno (2002), both versions of HIX are used throughout the literature. Ohno is better when samples have
#' higher absorbance because it accounts for inner filter effects better.
#'
#' Freshness Index (\eqn{\beta:\alpha}, fresh): ex = 310 nm, ratio of em = 380 nm to max in em = 420:435 nm.
#' An indication of recently produced DOM, higher values indicate more recently produced DOM.
#'
#' Relative Fluorescence Efficiency (RFE): Ratio of fluorescence at ex = 370 nm, em = 460 nm to
#' absorbance at 370 nm. An indicator of the relative amount of algal to non-algal DOM.
#'
#' Biological Index (BIX): ex = 310 nm, ratio of em = 380 nm to em = 430 nm.
#' An indicator of autotrophic productivity, values above 1 indicate recently produced
#' autochthonous DOM.
#'
#'
#' @param eem an object of class eem or eemlist
#' @param abs_data dataframe of absorbance data matching the samples, used to calculate relative fluorescence efficiency
#' @param noise_ratio numeric indicating the noise threshold for ratio calculations
#' @param verbose logical determining if additional messages should be printed
#' @returns Returns a data frame containing fluorescence peaks and indices for each sample
#' @references Coble, P. G., Lead, J., Baker, A., Reynolds, D. M., & Spencer, R. G. M. (Eds.). (2014). Aquatic Organic Matter Fluorescence. Cambridge: Cambridge University Press. https://doi.org/10.1017/CBO9781139045452
#' @references Hansen, A. M., Kraus, T. E. C., Pellerin, B. A., Fleck, J. A., Downing, B. D., & Bergamaschi, B. A. (2016). Optical properties of dissolved organic matter (DOM): Effects of biological and photolytic degradation. Limnology and Oceanography, 61(3), 1015–1032. https://doi.org/10.1002/lno.10270
#' @references Egan, J. K., McKnight, D. M., Bowman, M. M., SanClements, M. D., Gallo, A. C., Hatten, J. A., & Matosziuk, L. M. (2023). Identifying photochemical alterations of dissolved pyrogenic organic matter using fluorescence spectroscopy. Aquatic Sciences, 85(2), 38. https://doi.org/10.1007/s00027-022-00919-7
#' @references Zsolnay, A., E. Baigar, M. Jimenez, B. Steinweg, and F. Saccomandi. 1999. Differentiating with fluorescence spectroscopy the sources of dissolved organic matter in soils subjected to drying. Chemosphere 38: 45–50. doi:10.1016/S0045-6535(98)00166-0
#' @references Ohno, T. 2002. Fluorescence inner-filtering correction for determining the humification index of dissolved organic matter. Environ Sci Technol 36: 742–746. doi:10.1021/es0155276
#' @export
#'
eem_coble_peaks2 <- function (eem, abs_data, noise_ratio = 5, verbose = FALSE){
  #subfunctions for function
  msg_warning_wavelength <- function() {
    msg <- "This metric uses either excitation or emission wavelengths that were not present in the data. Data has been interpolated to fit the requested wavelengths."
    return(msg)
  }
  #function to interpolate an eem's area over a range of excitation and emissions
  #will report the maximum value
  max_peak_val <- function(ex, em, eem){
    if(max(ex) <= max(eem$ex) & max(em) <= max(eem$em) &
       min(ex) >= min(eem$ex) & min(em) >= min(eem$em)){
      ex_p <- rep(ex, length(em))
      em_p <- rep(em, length(ex))
      em_p <- em_p[order(em_p)]
      int_res <- pracma::interp2(eem$ex, eem$em, eem$x, ex_p, em_p)
      max_res <- max(int_res, na.rm=T)
    }else{
      max_res <- NA
    }
    max_res

  }
  #gets ratio after checking for noise thresholds
  ratio_val <-function(x, y, noise){
    if(is.na(x) == F & is.na(y) == F){
      if(x >= noise & y >= noise){ratio <- x/y}else{ratio <- NA}
    }else{ratio <- NA}
    ratio
  }

  #function checks
  stopifnot(.is_eemlist(eem) | .is_eem(eem) | is.data.frame(abs_data)|
              is.numeric(noise_ratio) | is.logical(verbose))

  #run over entire eemlist is eemlist
  if (.is_eemlist(eem)) {
    res <- lapply(eem, eem_coble_peaks2, noise_ratio = noise_ratio, verbose = verbose, abs_data=abs_data)
    res <- dplyr::bind_rows(res)
    return(res)
  }

  #list peaks
  coble_em_peak <- list(pB=300:320, pT=320:350, pA=380:480, pM=380:420,
                        pC=420:480, pD=509, pE=521, pN=370, FI=c(470,520),
                        HIX=c(300:345), HIX_o=c(435:480,300:345),fresh=c(380, 420:435),
                        RFE=460, BIX=c(380, 430))

  coble_ex_peak <- list(pB=270:280, pT=270:280, pA=250:260, pM=310:320,
                        pC=330:350, pD=390, pE=455, pN=280, FI=370,
                        HIX=254, HIX_o=254, fresh=310, RFE=370, BIX=310)

  #checks if all the required peaks are in the eem dataset
  if (!all(coble_ex_peak %in% eem$ex) & verbose) {
    warning(msg_warning_wavelength(), call. = FALSE)
  }
  if (!all(coble_em_peak %in% eem$em) & verbose) {
    warning(msg_warning_wavelength(), call. = FALSE)
  }

  #get noise for checking ratios, grabs area below 1:1 line
  noise <- eem
  em <- noise$em
  ex <- noise$ex
  x <- noise$x
  ind <- mapply(function(x) em > x, ex)
  ind <- ifelse(ind == 1, NA, 1)
  x <- x * ind
  noise <- mean(x, na.rm=T)
  noise_val <- mean(noise) * noise_ratio #data must be 5 times higher than noise or it won't calculate ratios

  #get coble peak values
  pB <- max_peak_val(270:280, 300:320, eem)
  pT <- max_peak_val(270:280, 320:350, eem)
  pA <- max_peak_val(250:260, 380:480, eem)
  pM <- max_peak_val(310:320, 380:420, eem)
  pC <- max_peak_val(330:350, 420:480, eem)
  pD <- max_peak_val(390, 509, eem)
  pE <- max_peak_val(455, 521, eem)
  pN <- max_peak_val(280, 370, eem)

  #peak ratios
  rAT <- ratio_val(pA, pT, noise_val)
  rCA <- ratio_val(pC, pA, noise_val)
  rCM <- ratio_val(pC, pM, noise_val)
  rCT <- ratio_val(pC, pT, noise_val)

  #fluorescence index
  fluo_470 <- pracma::interp2(eem$ex, eem$em, eem$x, 370, 470)
  fluo_520 <- pracma::interp2(eem$ex, eem$em, eem$x, 370, 520)
  if(fluo_470 >= noise_val & fluo_520 >= noise_val){
    FI <- fluo_470/fluo_520
    #get peak FI
    fluo_370 <- pracma::interp2(eem$ex, eem$em, eem$x, rep(370, 351), 248:598)
    FI_data <- data.frame(em=248:598, intensity=fluo_370)
    FI_max <- as.numeric(na.omit(FI_data$em[FI_data$intensity == max(FI_data$intensity, na.rm=T)]))
  } else{
    FI <- NA
    FI_max <- NA
  }

  #humification index
  em_435_480 <- seq(from = 435, to = 480, by = 1)
  em_300_345 <- seq(from = 300, to = 345, by = 1)
  ex_254 <- rep(254, length(em_300_345))
  sum_em_435_480 <- sum(pracma::interp2(eem$ex, eem$em, eem$x,
                                            ex_254, em_435_480))
  sum_em_300_345 <- sum(pracma::interp2(eem$ex, eem$em, eem$x,
                                            ex_254, em_300_345))

  if(sum_em_435_480 >= noise_val & sum_em_300_345 >= noise_val){
     HIX <- sum_em_435_480/(sum_em_300_345)
  } else{HIX <- NA}

  #humification index (ohno)
  em_435_480 <- seq(from = 435, to = 480, by = 1)
  em_300_345 <- seq(from = 300, to = 345, by = 1)
  ex_254 <- rep(254, length(em_300_345))
  sum_em_435_480 <- sum(pracma::interp2(eem$ex, eem$em, eem$x,
                                        ex_254, em_435_480))
  sum_em_300_345 <- sum(pracma::interp2(eem$ex, eem$em, eem$x,
                                        ex_254, em_300_345))

  if(sum_em_435_480 >= noise_val & sum_em_300_345 >= noise_val){
    HIX_o <- sum_em_435_480/(sum_em_300_345 + sum_em_435_480)
  } else{HIX_o <- NA}

  #fresh
  fluo_380 <- pracma::interp2(eem$ex, eem$em, eem$x, 310, 380)
  fluo_420_435 <- max_peak_val(310, 420:435, eem)
  if(fluo_380 >= noise_val & fluo_420_435 >= noise_val){
    fresh <- fluo_380/fluo_420_435
  } else{fresh <- NA}

  #bix
  fluo_380 <- pracma::interp2(eem$ex, eem$em, eem$x, 310, 380)
  fluo_430 <- pracma::interp2(eem$ex, eem$em, eem$x, 310, 430)
  if(fluo_380 >= noise_val & fluo_430 >= noise_val){
    BIX <- fluo_380/fluo_430
  } else{BIX <- NA}

  #RFE
  abs_col <- which(eem$sample == colnames(abs_data)) #get sample column
  abs_370 <- stats::approx(abs_data$wavelength, abs_data[,abs_col], 370)[[2]] #interpolate at 370
  RFE <- pracma::interp2(eem$ex, eem$em, eem$x, 370, 460) / abs_370

  output <- data.frame(sample = eem$sample, pB=pB, pT=pT, pA=pA, pM=pM,
                       pC=pC, pD=pD, pE=pE, pN=pN, rAT=rAT, rCA=rCA,
                       rCM=rCM, rCT=rCT, FI=FI, FI_max = FI_max, HIX=HIX, HIX_ohno=HIX_o, fresh=fresh,
                       RFE=RFE, BIX=BIX, stringsAsFactors = FALSE)

  #clean a little
  output <- do.call(data.frame,                      # Replace Inf in data by NA
                     lapply(output,
                            function(x) replace(x, is.infinite(x), NA)))

  output <- do.call(data.frame,                      # Replace nan in data by NA
                    lapply(output,
                           function(x) replace(x, is.nan(x), NA)))
  #print data
  return(output)
}


#' Export absorbance and fluorencence indices and peaks to an excel file
#'
#' Uses 'eem_coble_peaks2' and 'abs_parm' functions to obtain absorbance and
#' fluoresence indices for an eemlist. Saves data in an excel file.
#'
#' Use 'clean_files' function before running this function to ensure file structure is correct for pre-processing.
#'
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#'
#' @param eem_list an object of class eem or eemlist
#' @param abs_data a dataframe containing the absorbance data
#' @param meta the metadata associated with the samples
#' @param prjpath location to save the exported absorbance and fluorescence metrics
#' @param doc_norm either TRUE, FALSE or "both" indicating if peaks should be normalized by DOC concentrations
#' @param sampsascol a logical indicating how results should be oriented, TRUE puts samples as columns, FALSE puts samples as rows
#' @param cuvle the length of the cuvette in cm
#' @param waves optional, a vector of wavelengths in nm to extract from the absorbance data
#' @export

get_indices <- function(eem_list, abs_data, meta, prjpath,doc_norm="both",
                        sampsascol=F, waves=NULL, cuvle=1){

  stopifnot(.is_eem(eem_list) | .is_eemlist(eem_list) | is.data.frame(c(abs_data, meta))|
              is.character(prjpath) | is.logical(sampsascol) | file.exists(prjpath))

  #information to put with peaks
  documentation <- data.frame(c("Peaks extracted using fewsdom package in R.",
  "For peak definitions see 'eem_coble_peaks2' and 'abs_parm' functions in the fewsdom package.",
  "The package can be downloaded from https://github.com/katiewampler/fewsdom"))

  #saving location
  if(file.exists(paste(prjpath,  "/5_Processed", sep=""))==F){
    stop("Invalid file structure, please use 'create_files' function to create file for plots within file directory")
  }

  file_date <- unlist(strsplit(prjpath, "/"))
  file_date <-file_date[length(file_date)]
  wb_name <- paste(prjpath, "/5_Processed/SpectralIndices_", file_date,".xlsx", sep="")

  #add documentation info
  wb <- openxlsx::createWorkbook(wb_name)
  openxlsx::addWorksheet(wb, "README")
  openxlsx::writeData(wb,"README", documentation, colNames = F)
  openxlsx::saveWorkbook(wb,wb_name,overwrite = TRUE)

  #get flourescence data
  if(doc_norm == T | doc_norm == "both"){
    eem_doc <- .eem_doc_norm(eem_list, meta)
    if(length(eem_doc) >0){
      coble <- eem_coble_peaks2(eem_doc, abs_data = abs_data, verbose = F)
      if(sampsascol == T){coble <- clean_transpose(coble)}
      .OSU_excel(wb_name, "fluor_indices_DOC", coble, sampsascol)
    } else{
      warning("No DOC data was provided, DOC normalized fluorescence metrics were not calculated")
    }

  }
  if(doc_norm == F | doc_norm == "both"){
    eem_list <- .eem_doc_rm(eem_list, meta) #remove any normalization done
    coble <- eem_coble_peaks2(eem_list, abs_data=abs_data, verbose = F)
    if(sampsascol == T){coble <- clean_transpose(coble)}
    .OSU_excel(wb_name, "fluor_indices", coble, sampsascol)
  }

  #get absorbance data
    abs_df <- abs_parm(abs_data, meta=meta, cuvle = cuvle, keep_all=F, waves=c(254, waves))

    if(sampsascol == T){abs_df <- clean_transpose(abs_df)}
    if(nrow(abs_df) > 0){
      .OSU_excel(wb_name, "abs_indices", abs_df, sampsascol)
    }else{warning("No DOC data was provided, absorbance metrics were not calculated")}
}

#' Exclude complete wavelengths or samples form data set
#'
#' Outliers in all modes should be avoided. With this functions excitation or
#' emission wavelengths as well as samples can be removed completely from your
#' sample set. Take directly from staRdom package to prevent orphaned package
#' errors
#'
#' @importFrom eemR eem_extract
#' @param eem_list object of class eemlist
#' @param exclude list of three vectors, see details
#' @param verbose stats whether additional information is given in the command line
#'
#' The argument exclude is a named list of three vectors.
#' The names must be "ex", "em" and "sample". Each element contains a vector of
#' wavelengths or sample names that are to be excluded from the data set.
#'
#' @returns object of class eemlist
#' @noRd

eem_exclude <- function (eem_list, exclude = list, verbose = FALSE)
{
  ex_exclude <- exclude[["ex"]]
  em_exclude <- exclude[["em"]]
  sample_exclude <- exclude[["sample"]]
  if (!is.null(sample_exclude)) {
    sample_exclude <- paste0("^", sample_exclude, "$")
    eem_list <- eemR::eem_extract(eem_list, sample_exclude, verbose = verbose)
  }
  eem_list <- lapply(eem_list, function(eem) {
    eem$x <- eem$x[!eem$em %in% em_exclude, !eem$ex %in%
                     ex_exclude]
    eem$ex <- eem$ex[!eem$ex %in% ex_exclude]
    eem$em <- eem$em[!eem$em %in% em_exclude]
    eem
  })
  if (!is.null(ex_exclude) & verbose)
    cat(paste0("Removed excitation wavelength(s): ", paste0(ex_exclude %>%
                                                              sort(), collapse = ", ")), fill = TRUE)
  if (!is.null(em_exclude) & verbose)
    cat(paste0("Removed emission wavelength(s): ", paste0(em_exclude %>%
                                                            sort(), collapse = ", ")), fill = TRUE)
  class(eem_list) <- "eemlist"
  eem_list
}
