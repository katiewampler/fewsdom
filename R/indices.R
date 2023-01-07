#' Cleanly Transpose Data
#'
#' Takes a dataframe and transposes it (turn the rows to columns and columns to rows).
#'
#' @param df A dataframe that you want to transpose.
#' @returns A transposed dataframe.
#' @export

#functions for getting indices
clean_transpose <- function(df){
  df_flip <- as.data.frame(t(df))
  df_names <- df_flip[1,]
  df_flip <- df_flip[2:nrow(df_flip),]
  colnames(df_flip) <- df_names
  df_flip
}

#' Get Absorbance Peaks and Indices
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
#' S290_350: Spectral slope between 290 to 350 nm.
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
#' @param abs_data A dataframe containing the absorbance data, where each column is a sample.
#' @param waves Optional, a vector of wavelengths in nm to extract from the absorbance data.
#' @param meta The metadata table for the sample run, should include DOC data in mg/L.
#' @param r_thresh The minimum \eqn{\text{R}^{2}} value for the spectral slopes for the slope to be included in the table.
#' @param keep_all A logical, TRUE will return all samples even those without DOC data, FALSE will only return samples with DOC data.
#' @returns Returns a data frame of indices.
#' @references Hansen, A. M., Kraus, T. E. C., Pellerin, B. A., Fleck, J. A., Downing, B. D., & Bergamaschi, B. A. (2016). Optical properties of dissolved organic matter (DOM): Effects of biological and photolytic degradation. Limnology and Oceanography, 61(3), 1015–1032. https://doi.org/10.1002/lno.10270
#' @export

abs_parm <- function(abs_data, waves=NULL, meta=meta,
                     r_thresh=0.8, keep_all=F){
  #interpolate data to find values that don't fall on the excitation wavelengths used
  data <- lapply(abs_data[,-1], function(x) approx(abs_data$wavelength, x, xout=c(249:791))[[2]])
  data <- as.data.frame(data)
  data$wavelength <- 249:791

  #build data table to put in data
  abs_out <- as.data.frame(matrix(nrow=(ncol(abs_data)-1), ncol=(15 + length(waves))))
  colnam <- c("sample", paste("SUVA", c(254,280, 350, 370), sep=""), paste("SVA", c(412,440,480,510,532,555), sep=""), "S275_295", "S290_350","S350_400", "SR")
  if(length(waves) > 0){
    colnam <- c(colnam, paste("a", waves, sep=""))
  }
  colnames(abs_out) <- colnam
  abs_out$sample <- colnames(abs_data)[2:ncol(abs_data)]

  #add DOC data to table for specific values
  abs_out$DOC <- NA
  for(m in meta$index){
    x <- which(m == meta$index)
    abs_col <- which(meta$unique_ID[x] == abs_out$sample)
    abs_out$DOC[abs_col] <- meta$DOC_mg_L[x]}

  #extract data for suva and sva measurements
  for(x in c(254,280, 350, 370, 412,440,480,510,532,555)){
    row_num <- which(data$wavelength == x) #row in absorbance data
    col_num <- which(stringr::str_detect(colnames(abs_out), as.character(x))==T)[1] #location in output table
    abs_out[,col_num] <- as.numeric(data[row_num, 1:(ncol(data)-1)]) / abs_out$DOC * 100 #put data in column
  }

  #extract data for added absorbance wavelengths
  for(x in waves){
    row_num <- which(data$wavelength == x) #row in absorbance data
    col_num <- which(stringr::str_detect(colnames(abs_out), as.character(x))==T) #location in output table
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
  S275_295 <- subset(data, data$wavelength >= 275 & data$wavelength <= 295)
  .get_S275_295 <- function(abs_col, waves=S275_295$wavelength){
    spec_data <- data.frame(waves = waves, abs=abs_col)
    colnames(spec_data) <- c("waves", "abs_data")
    spec_data$abs_data_m <- spec_data$abs_data * 100
    spec_data$abs_data_ln <- log(spec_data$abs_data_m)
    fit <-  lm(abs_data_m ~ waves, spec_data)
    slope <- fit$coefficients[2] * -1
    r2 <- summary(fit)$r.squared
    output <- c(slope, r2)
    output
  }

  S290_350 <- subset(data, data$wavelength >= 290 & data$wavelength <= 350)
  .get_S290_350 <- function(abs_col, waves=S290_350$wavelength){
    spec_data <- data.frame(waves = waves, abs=abs_col)
    colnames(spec_data) <- c("waves", "abs_data")
    spec_data$abs_data_m <- spec_data$abs_data * 100
    spec_data$abs_data_ln <- log(spec_data$abs_data_m)
    fit <-  lm(abs_data_m ~ waves, spec_data)
    slope <- fit$coefficients[2] * -1
    r2 <- summary(fit)$r.squared
    output <- c(slope, r2)
    output
  }

  S350_400 <- subset(data, data$wavelength >= 350 & data$wavelength <= 400)
  .get_S350_400 <- function(abs_col, waves=S350_400$wavelength){
    spec_data <- data.frame(waves = waves, abs=abs_col)
    colnames(spec_data) <- c("waves", "abs_data")
    spec_data$abs_data_m <- spec_data$abs_data * 100
    spec_data$abs_data_ln <- log(spec_data$abs_data_m)
    fit <-  lm(abs_data_m ~ waves, spec_data)
    slope <- fit$coefficients[2] * -1
    r2 <- summary(fit)$r.squared
    output <- c(slope, r2)
    output
  }

  S275_295_r <- suppressWarnings(apply(S275_295[,1:(ncol(data)-1)], 2, .get_S275_295))
  S290_350_r <- suppressWarnings(apply(S290_350[,1:(ncol(data)-1)], 2, .get_S290_350))
  S350_400_r <- suppressWarnings(apply(S350_400[,1:(ncol(data)-1)], 2, .get_S350_400))

  abs_out$S275_295 <- S275_295_r[1,]
  abs_out$S275_295[S275_295_r[2,] <= r_thresh] <- NA #remove if R2 is poor
  abs_out$S290_350 <- S290_350_r[1,]
  abs_out$S290_350[S290_350_r[2,] <= r_thresh] <- NA #remove if R2 is poor
  abs_out$S350_400 <- S350_400_r[1,]
  abs_out$S350_400[S350_400_r[2,] <= r_thresh] <- NA #remove if R2 is poor


  abs_out$SR <- abs_out$S275_295 / abs_out$S350_400
  abs_out$SR[abs_out$S275_295 <= 0.005] <- NA

  #remove samples missing DOC
  if(keep_all == F){
    abs_out <- abs_out[is.na(abs_out$DOC) ==F,]
  }
  abs_out

}

#' Get Coble Peaks and Other Fluorescence Indicies
#'
#' Modified from 'eem_coble_peaks' from eemR package to include some additional indices and
#' modified definitions of Coble peaks.  For indices that are ratios, the values in the ratio
#' must be 5 times (or specifed amount) greater than the average values above 700nm (considered noise). If the
#' desired wavelength is not explicitly measured the EEM will be interpolated using
#' the interp2 function from the pracma library.
#'
#' Coble peaks are based on Coble et al. 2014 and are defined as follows:
#'
#' Peak B (pB): ex = 270:280 nm, em = 300:320 nm, Tyrosine-like
#'
#' Peak T (pT): ex = 270:280 nm, em = 320:350 nm, Tryptophan-like
#'
#' Peak A (pA): ex = 250:260 nm, em = 380:480 nm, Humic-like
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
#' Humification Index (HIX): ex = 254 nm, em =\eqn{\sum}435:480 divided by em =\eqn{\sum}435:480, \eqn{\sum}300:345.
#' An indication of humic substances or extent of humification. Higher values indicate an higher degree of humification.
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
#' @param eem An object of class eem or eemlist.
#' @param absorbance A dataframe of absorbance data matching the samples, used to calculate relative fluorescence efficiency
#' @param verbose Logical determining if additional messages should be printed.
#' @param noise_ratio Numeric indicating the noise threshold for ratio calculations.
#' @returns Returns a data frame containing fluorescence peaks and indices for each sample.
#' @references Coble, P. G., Lead, J., Baker, A., Reynolds, D. M., & Spencer, R. G. M. (Eds.). (2014). Aquatic Organic Matter Fluorescence. Cambridge: Cambridge University Press. https://doi.org/10.1017/CBO9781139045452
#' @references Hansen, A. M., Kraus, T. E. C., Pellerin, B. A., Fleck, J. A., Downing, B. D., & Bergamaschi, B. A. (2016). Optical properties of dissolved organic matter (DOM): Effects of biological and photolytic degradation. Limnology and Oceanography, 61(3), 1015–1032. https://doi.org/10.1002/lno.10270
#'
eem_coble_peaks2 <- function (eem, absorbance, noise_ratio = 5, verbose = FALSE){
  #subfunctions for function
  .is_eemlist <- function(eem) {
    ifelse(class(eem) == "eemlist", TRUE, FALSE)
  }
  .is_eem <- function(eem) {
    ifelse(class(eem) == "eem", TRUE, FALSE)
  }
  msg_warning_wavelength <- function() {
    msg <- "This metric uses either excitation or emission wavelengths that were not present in the data. Data has been interpolated to fit the requested wavelengths."
    return(msg)
  }
  #function to interpolate an eem's area over a range of excitation and emissions
  #will report the maximum value
  max_peak_val <- function(ex, em){
    ex_p <- rep(ex, length(em))
    em_p <- rep(em, length(ex))
    em_p <- em_p[order(em_p)]
    int_res <- pracma::interp2(eem$ex, eem$em, eem$x, ex_p, em_p)
    max_res <- max(int_res, na.rm=T)
    max_res
  }
  #gets ratio after checking for noise thresholds
  ratio_val <-function(x, y, noise){
    if(x >= noise & y >= noise){ratio <- x/y}else{ratio <- NA}
    ratio
  }

  #function checks
  stopifnot(.is_eemlist(eem) | .is_eem(eem))

  #run over entire eemlist is eemlist
  if (.is_eemlist(eem)) {
    res <- lapply(eem, eem_coble_peaks2, noise_ratio = noise_ratio, verbose = verbose, absorbance=absorbance)
    res <- dplyr::bind_rows(res)
    return(res)
  }

  #list peaks
  coble_em_peak <- list(pB=300:320, pT=320:350, pA=380:480, pM=380:420,
                        pC=420:480, pD=509, pE=521, pN=370, FI=c(470,520),
                        HIX=c(435:480,300:345), fresh=c(380, 420:435),
                        RFE=460, BIX=c(380, 430))

  coble_ex_peak <- list(pB=270:280, pT=270:280, pA=250:260, pM=310:320,
                        pC=330:350, pD=390, pE=455, pN=280, FI=370,
                        HIX=254, fresh=310, RFE=370, BIX=310)

  #checks if all the required peaks are in the eem dataset
  if (!all(coble_ex_peak %in% eem$ex) & verbose) {
    warning(msg_warning_wavelength(), call. = FALSE)
  }
  if (!all(coble_em_peak %in% eem$em) & verbose) {
    warning(msg_warning_wavelength(), call. = FALSE)
  }

  #get noise for checking ratios
  noise <- eem$x[1:24,72:104]
  noise_val <- mean(noise) * noise_ratio #data must be 5 times higher than noise or it won't calculate ratios

  #get coble peak values
  pB <- max_peak_val(270:280, 300:320)
  pT <- max_peak_val(270:280, 320:350)
  pA <- max_peak_val(250:260, 380:480)
  pM <- max_peak_val(310:320, 380:420)
  pC <- max_peak_val(330:350, 420:480)
  pD <- max_peak_val(390, 509)
  pE <- max_peak_val(455, 521)
  pN <- max_peak_val(280, 370)

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
  } else{
    FI <- NA
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
    HIX <- sum_em_435_480/(sum_em_300_345 + sum_em_435_480)
  } else{HIX <- NA}


  #fresh
  fluo_380 <- pracma::interp2(eem$ex, eem$em, eem$x, 310, 380)
  fluo_420_435 <- max_peak_val(310, 420:435)
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
  abs_col <- which(eem$sample == colnames(absorbance)) #get sample column
  abs_370 <- stats::approx(abs_data$wavelength, abs_data[,abs_col], 370)[[2]] #interpolate at 370
  RFE <- pracma::interp2(eem$ex, eem$em, eem$x, 370, 460) / abs_370

  #print data
  return(data.frame(sample = eem$sample, pB=pB, pT=pT, pA=pA, pM=pM,
                    pC=pC, pD=pD, pE=pE, pN=pN, rAT=rAT, rCA=rCA,
                    rCM=rCM, rCT=rCT, FI=FI, HIX=HIX, fresh=fresh,
                    RFE=RFE, BIX=BIX, stringsAsFactors = FALSE))
}


#' Export Absorbance and Fluorencence Indices and Peaks to an Excel File
#'
#'
#'
#' @param eem_list An object of class eem or eemlist.
#' @param abs_data A dataframe containing the absorbance data.
#' @param doc_norm Either TRUE, FALSE or "both" indicating if peaks should be normalized by DOC concentrations.
#' @param meta The metadata associated with the samples.
#' @param path Location to save the exported absorbance and fluorescence metrics.

get_indices <- function(eem_list, abs_data, , meta=meta,
                        path){
  .OSU_excel <- function(file, sheet_name, df){
    wb <- openxlsx::loadWorkbook(file)
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb,sheet_name, df, rowNames = sampsascol)
    openxlsx::saveWorkbook(wb,file,overwrite = TRUE)
  }

  #copy peak file into folder
  file.copy(from="T:/LabMembers/10_Katie/11_Aqualog Processing/Peak_ReadMe.xlsx", to=paste(demopath, "/5_Processed/SpectralIndices.xlsx", sep=""), overwrite = T)
  wb_name <- paste(demopath, "/5_Processed/SpectralIndices.xlsx", sep="")
  #get flourescence data
  if(DOC_norm_index == T | DOC_norm_index == "both"){
    coble <- eem_coble_peaks2(eem_list_doc, verbose = F)
    if(sampsascol == T){coble <- clean_transpose(coble)}
    .OSU_excel(wb_name, "EEM_Indices_DOC", coble)
  }
  if(DOC_norm_index == F | DOC_norm_index == "both"){
    coble <- eem_coble_peaks2(eem_list, verbose = F)
    if(sampsascol == T){coble <- clean_transpose(coble)}
    .OSU_excel(wb_name, "EEM_Indices", coble)
  }

  #get absorbance data
  abs_wave <- c(250,254, 281,350,365,412,440, 465, 665)

  if(DOC_norm_index == T | DOC_norm_index == "both"){
    abs_df <- abs_parm(abs_data_doc, waves=abs_wave)
    if(sampsascol == T){abs_df <- clean_transpose(abs_df)}
    .OSU_excel(wb_name, "Abs_Indices_DOC", abs_df)

  }else if(DOC_norm_index == T){
    message("Error: Attempt to use DOC normalized data with missing DOC values")
  } else if(DOC_norm_index == "both"){
    message("Warning: Attempt to use DOC normalized data with missing DOC values, only unnormalized data will be reported")

  }

  #get other absorbance data
  if(DOC_norm_index == F | DOC_norm_index == "both"){
    abs_df <- abs_parm(abs_data, waves=abs_wave)
    if(sampsascol == T){abs_df <- clean_transpose(abs_df)}

    #if there's DOC, put in table and get SUVA 254
    if(sum(!is.na(meta$DOC_mg_L)) > 0){
      abs_df$DOC <- NA
      for(x in 1:nrow(abs_df)){
        val <- match(abs_df$sample[x], meta$unique_ID)
        abs_df$DOC[x] <- meta$DOC_mg_L[val]
      }
      abs_df$SUVA254 <- abs_df$a254/abs_df$DOC * 100
    }  else{
      abs_df$DOC <- NA
      abs_df$SUVA254 <- NA
    }
    .OSU_excel(wb_name, "Abs_Indices", abs_df)

  }

  cat("Spectral Indices have been saved")
}

