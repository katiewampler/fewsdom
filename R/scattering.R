#' Removes Rayleigh and Raman Scattering
#'
#' Based on 'eem_remove_scattering' function from eemR package, added functionality
#' to cut different amounts from above and below the emission line. This function will only remove
#' one kind of scattering at a time. Use on an class eem or eemlist. Should be performed
#' on data that is untransformed besides a blank subtraction.
#'
#' @importFrom stringr str_detect str_split
#' @importFrom usefun is_between
#' @importFrom pracma findpeaks pchip
#' @importFrom stats approx lm median na.omit
#' @importFrom utils read.csv read.table write.table
#' @import dplyr
#'
#' @param eem An object of class eemlist or eem.
#' @param type A string, either "raman" or "rayleigh".
#' @param order An integer number, either 1 (first order) or 2 (second order).
#' @param up_width An integer number specifying the width in nm for the cut above the emission line.
#' @param down_width An integer number specifying the width in nm for the cut below the emission line.
#' @export

eem_remove_scattering2 <- function(eem,up_width=10, down_width=10, type="raman", order=1){

  #functions to run code
  .find_raman_peaks <- function(ex) {

    # For water, the Raman peak appears at a wavenumber 3600 cm lower than the
    # incident wavenumber. For excitation at 280 nm, the Raman peak from water
    # occurs at 311 nm. Source : Principles of Fluorescence Spectroscopy (2006) -
    # Third Edition.pdf

    ## Convert wavenumber from nm to cm
    ex_wave_number <- 1 / ex

    raman_val <- 3500 / 10^7
    raman_peaks <- ex_wave_number - raman_val

    ## Bring back to nm
    raman_peaks <- 1 / raman_peaks

    # raman_peaks <- -(ex / (0.00036 * ex - 1))

    return(raman_peaks)
  }
  .is_eemlist <- function(eem) {
    ifelse(class(eem) == "eemlist", TRUE, FALSE)
  }
  .is_eem <- function(eem) {
    ifelse(class(eem) == "eem", TRUE, FALSE)
  }

  #function checks
  stopifnot(.is_eemlist(eem) | .is_eem(eem), all(type %in%
                                                   c("raman", "rayleigh")), is.numeric(order),length(order) == 1, length(type) == 1,
            usefun::is_between(order, 1, 2, include.high.value = T))

  #if its an eemlist
  if (.is_eemlist(eem)) {
    res <- lapply(eem, eem_remove_scattering2, type = type,
                  order = order, up_width = up_width, down_width=down_width)
    class(res) <- class(eem)
    return(res)
  }

  #get eems data
  x <- eem$x
  em <- eem$em
  ex <- eem$ex

  #get raman peaks if needed
  if (type == "raman") {ex <- .find_raman_peaks(eem$ex)}

  #determines which rows to cut out
  ind1 <- mapply(function(x) em <= x, order * ex - down_width)
  ind2 <- mapply(function(x) em <= x, order * ex + up_width)
  ind3 <- ifelse(ind1 + ind2 == 1, NA, 1) #if they're both true it doesn't get cut, otherwise make NA
  x <- x * ind3

  #put new eem into eem object
  res <- eem
  res$x <- x
  attributes(res) <- attributes(eem)
  attr(res, "is_scatter_corrected") <- TRUE
  class(res) <- class(eem)

  return(res)
}


#' Finds Optimal Scattering Cut Widths Based on EEM set
#'
#' Takes a group of eems and uses peak picking to determine how wide the cut should be
#' for each scattering line (does not work for 2nd order raman). Uses the exciatation at 280nm
#' and finds the start and end of the peak. If the function can't find one, it will default to NA
#' for that sample. Median of the peak widths is reported for use in the 'remove_scattering' function.
#' Should be performed on data that is untransformed besides a blank subtraction.
#'
#
#' @param eem An object of class eemlist.
#' @param type A string, either "raman" or "rayleigh".
#' @param order An integer number, either 1 (first order) or 2 (second order).
#' @export

find_cut_width <- function(eem, type="rayleigh", order=1){
  .is_eemlist <- function(eem) {
    ifelse(class(eem) == "eemlist", TRUE, FALSE)
  }
  .percent_dif <- function(start, end){
    result <- (end - start) / start * 100
    result
  }

  stopifnot(.is_eemlist(eem))
  if(type == "raman" & order == 2){
    warning("Cannot use the auto cutting method for second order raman line, providing estimate.")

  }else if(type == "raman" & order == 1){
    row_start <- 26
    row_end <- 32
    thresh <- 5
    up_cor <- 0
  }else if(type == "rayleigh" & order == 1){
    row_start <- 14
    row_end <- 18
    thresh <- 5
    up_cor <- 5
  }else if(type == "rayleigh" & order == 2){
    row_start <- 130
    row_end <- 145
    up_cor <- 0
    thresh <- 5 #smaller threshold for this one because it's right in area of interest
  }else{
    stop("Please choose either rayleigh or raman and orders 1 or 2.")
  }

  #find average width to cut
  for(n in 1:length(eem)){
    eem_ind <- eem[[n]]
    x <- eem_ind$x
    em <- eem_ind$em
    ex_280 <- x[,15]
    df <- data.frame(em = em, int=ex_280)

    #ggplot(df, aes(x=em, y=int)) + geom_line()
    peaks <- as.data.frame(pracma::findpeaks(df$int, sortstr=T, minpeakheight = 100))
    peaks <- peaks[peaks$V2 > row_start & peaks$V2 < row_end,]


    if(nrow(peaks) == 0){
      up_width <- NA
      down_width <- NA
    }else{
      peaks <- peaks[order(peaks$V1, decreasing=T),]
      peaks <- peaks[1,]

      #get percent change between points across spectra, put in df
      change <- .percent_dif(df[1:(nrow(df)-1),2],df[2:nrow(df),2])
      change <- c(change, NA)
      df$change <- change

      xp <- peaks$V2
      peak <- df[xp,]

      #get 30 nm range to look for peak
      start <- ifelse(xp - 10 < 0, 1, xp-10)
      end <- ifelse(xp + 10 < nrow(df), xp+10, nrow(df))
      range <- start:end
      clip <- df[range,]

      #identify start and end where range changes more than 20%/5 and is less than 90% of the peak
      pos_start <- clip[clip$change > thresh & clip$int < 0.9*peak$int,]
      pos_start <- pos_start[as.numeric(row.names(pos_start))< xp,]

      pos_end <- clip[clip$change > -thresh & clip$int < 0.9*peak$int,]
      pos_end <- pos_end[as.numeric(row.names(pos_end))> xp,]

      #find end points
      if((nrow(pos_end) == 0 & nrow(pos_start) == 0)){
        up_width <- NA
        down_width <- NA
      }else if(nrow(pos_end) == 0){
        start <- pos_start$em[1]
        up_width <- peak$em - start
        down_width <- NA
      }else if(nrow(pos_start) == 0){
        end <- pos_end$em[1]
        down_width <- end - peak$em
        up_width <- NA
      }else{
        start <- pos_start$em[nrow(pos_start)]
        end <- pos_end$em[1]

        #save auto found cutting points
        up_width <- peak$em - start + up_cor
        down_width <- end - peak$em
      }
    }

    df <- data.frame(eem_num = n, up_width=up_width, down_width=down_width)

    if(n == 1){
      df_test <- df
    }else{df_test <- rbind(df_test, df)}
  }
  up_width <- median(df_test$up_width, na.rm=T)
  down_width <- median(df_test$down_width, na.rm=T)

  widths <- c(round(up_width,1), round(down_width,1) )

  if(type == "raman" & order == 2){
    widths <- c(1.5, 1.5)
    }
  widths
}

#' Removes Rayleigh Scattering with Optional Interpolation
#'
#' Will remove both first and second order rayleigh scattering with the option to interpolate
#' the removed area. Also makes note that the correction was made in an optional textfile.
#'
#
#' @param eem An object of class eemlist.
#' @param rayleigh_mask Optional if auto width method is used, A vector of length 4 specifying the width of the rayleigh line to cut, numbers 1:2 are width above and below first order line, numbers 3:4 are width above and below second order line.
#' @param width_method Either "auto" or "manual". If auto is chosen cutting widths will be found using the 'find_cut_width' function.
#' @param rayleigh_interp A vector of length two, either T or F, specifying whether the first and second order lines should be interpolated. The first position refers to the first order line. The way the code is written you cannot interpolate the first order line and not the second.
#' @param process_file A file path to a .txt file, used to track processing changes to EEMs.
#' @export

rayleigh <- function(eem, rayleigh_mask=c(20,10,10,10), width_method="auto", rayleigh_interp=c(F,F), process_file=NULL){
  .is_eemlist <- function(eem) {
    ifelse(class(eem) == "eemlist", TRUE, FALSE)
  }
  #function checks
  stopifnot(.is_eemlist(eem))

  if(width_method == "auto"){
    ray1 <- find_cut_width(eem, type="rayleigh", order=1)
    ray2 <- find_cut_width(eem, type="rayleigh", order=2)
    rayleigh_mask <- c(ray1, ray2)
    cat(rayleigh_mask)
  }

  #second order rayleigh, done first in case there is interpolation
    suppressWarnings(eem_rm <- eem_remove_scattering2(eem=eem, type="rayleigh", order=2, up_width = rayleigh_mask[3], down_width = rayleigh_mask[4]))
    if(is.character(process_file)){
      write.table(paste(Sys.time(), "- Second order raleigh scattering was removed with a ", rayleigh_mask[3]," nm width above and ",
                        rayleigh_mask[4], " nm width below", sep=""), process_file, append=T, quote=F, row.names = F, col.names = F)}

  if(rayleigh_interp[2] == T){eem_rm <- eem_interp(data=eem_rm)
  if(is.character(process_file)){
    write.table(paste(Sys.time(), "- Second order raleigh scattering filled via interpolation", sep=""), process_file, append=T, quote=F, row.names = F, col.names = F)}}

  #first order rayleigh
    suppressWarnings(eem_rm <- eem_remove_scattering2(eem=eem_rm, type="rayleigh", order=1, up_width = rayleigh_mask[1], down_width = rayleigh_mask[2]))
    if(is.character(process_file)){
      write.table(paste(Sys.time(), "- First order raleigh scattering was removed with a ", rayleigh_mask[1]," nm width above and ",
                        rayleigh_mask[2], " nm width below", sep=""), process_file, append=T, quote=F, row.names = F, col.names = F)
    }

  if(rayleigh_interp[1] == T){eem_rm <- eem_interp(data=eem_rm)
  if(is.character(process_file)){
    write.table(paste(Sys.time(), "- First order raleigh scattering filled via interpolation", sep=""), process_file, append=T, quote=F, row.names = F, col.names = F)
  }}

  eem_rm
}


#' Removes Raman Scattering with Optional Interpolation
#'
#' Will remove both first and second order raman scattering with the option to interpolate
#' the removed area. Also makes note that the correct was made.
#'
#'
#' @param eem An object of class eemlist.
#' @param raman_mask A vector of length 4 specifying the width of the raman line to cut, numbers 1:2 are width above and below first order line, numbers 3:4 are width above and below second order line. Since you cannot use the auto method for second order raman this must be specified.
#' @param width_method Either "auto" or "manual". If auto is chosen cutting widths will be found using the 'find_cut_width' function.
#' @param raman_interp A vector of length two, either T or F, specifying whether the first and second order lines should be interpolated. The first position refers to the first order line. The way the code is written you cannot interpolate the first order line and not the second.
#' @param process_file A file path to a .txt file, used to track processing changes to EEMs.
#' @export

raman <- function(eem, raman_mask=c(10,10,1.5,1.5), width_method="auto", raman_interp=c(T,T), process_file=NULL){
  .is_eemlist <- function(eem) {
    ifelse(class(eem) == "eemlist", TRUE, FALSE)
  }
  #function checks
  stopifnot(.is_eemlist(eem))

  if(width_method == "auto"){
    ram1 <- find_cut_width(eem, type="raman", order=1)
    raman_mask <- c(ram1, raman_mask[3:4])
    cat(raman_mask)
  }
  #second order raman, done first for interpolation
  if(sum(raman_mask[3:4]) > 0){
    eem_rm <- eem_remove_scattering2(eem=eem, type="raman", order=2, up_width = raman_mask[3], down_width = raman_mask[4])
    if(is.character(process_file)){
      write.table(paste(Sys.time(), "- Second order raleigh scattering was removed with a ", raman_mask[3]," nm width above and ",
                        raman_mask[4], " nm width below", sep=""), process_file, append=T, quote=F, row.names = F, col.names = F)}}
  if(raman_interp[2] == T){eem_rm <- eem_interp(data=eem_rm)
  if(is.character(process_file)){
    write.table(paste(Sys.time(), "- Second order raleigh scattering filled via interpolation", sep=""), process_file, append=T, quote=F, row.names = F, col.names = F)}}

  #first order raman
    eem_rm <- eem_remove_scattering2(eem=eem_rm, type="raman", order=1, up_width = raman_mask[1], down_width = raman_mask[2])
    if(is.character(process_file)){
      write.table(paste(Sys.time(), "- First order raleigh scattering was removed with a ", raman_mask[1]," nm width above and ",
                        raman_mask[2], " nm width below", sep=""), process_file, append=T, quote=F, row.names = F, col.names = F)}

  if(raman_interp[1] == T){eem_rm <- eem_interp(data=eem_rm)
  if(is.character(process_file)){
    write.table(paste(Sys.time(), "- First order raleigh scattering filled via interpolation", sep=""), process_file, append=T, quote=F, row.names = F, col.names = F)}
  }

  eem_rm
}

#' Interpolates Missing Values in EEM Data
#'
#' 'eem_interp' function in staRdom package
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom zoo na.approx
#' @importFrom MBA mba.points
#' @importFrom tidyr gather
#' @import foreach

#'
#' @param data An object of class eemlist.
#' @param cores	 specify number of cores for parallel computation
#' @param type numeric 0 to 4 or TRUE which resembles type 1
#' @param verbose logical, whether more information on calculation should be provided
#' @param nonneg logical, whether more information on calculation should be provided
#' @param extend logical, whether data is extrapolated using type 1
#' @noRd
#'

eem_interp <- function (data, cores = parallel::detectCores(logical = FALSE),
                        type = TRUE, verbose = FALSE, nonneg = TRUE, extend = FALSE,
                        ...)
{
  cl <- makeCluster(spec = min(cores, length(data)), type = "PSOCK")
  `%dopar%` <- foreach::`%dopar%`
  doParallel::registerDoParallel(cl)
  if (verbose) {
    cat("interpolating missing data in", length(data), "EEMs",
        fill = TRUE)
  }
  eem_list <- foreach::foreach(i = 1:length(data)) %dopar% {
    eem <- data[[i]]
    if (type == 4) {
      eem$x <- cbind(zoo::na.approx(eem$x, ...), t(zoo::na.approx(t(eem$x),
                                                                  ...))) %>% array(c(nrow(eem$x), ncol(eem$x),
                                                                                     2)) %>% apply(1:2, mean, na.rm = TRUE)
    }
    if (type == 1 | type == TRUE) {
      x <- eem$x %>% data.frame() %>% `colnames<-`(eem$ex) %>%
        `rownames<-`(eem$em) %>% mutate(em = eem$em) %>%
        tidyr::gather("ex", "z", -em) %>% mutate_all(as.numeric)
      x2 <- x %>% filter(!is.na(z))
      x3 <- MBA::mba.points(xyz = x2 %>% select(em, ex,
                                                z), xy.est = expand.grid(em = eem$em, ex = eem$ex),
                            verbose = verbose, extend = extend, ...)
      eem$x[is.na(eem$x)] <- x3$xyz.est[, 3] %>% matrix(nrow = nrow(eem$x),
                                                        ncol = ncol(eem$x)) %>% .[is.na(eem$x)]
    }
    if (type == 2) {
      x1 <- try(eem$x %>% apply(1, function(row) pracma::pchip(xi = eem$ex[!is.na(row)],
                                                               yi = row %>% na.omit(), x = eem$ex), ...) %>%
                  t(), silent = TRUE)
      x2 <- try(eem$x %>% apply(2, function(col) pracma::pchip(xi = eem$em[!is.na(col)],
                                                               yi = col %>% na.omit(), x = eem$em), ...), silent = TRUE)
      if (inherits(x1, "try-error") & inherits(x2, "try-error"))
        warning(eem$sample, " could not be interpolated!")
      if (inherits(x1, "try-error"))
        x1 <- matrix(NA, nrow(eem$x), ncol(eem$x))
      if (inherits(x2, "try-error"))
        x2 <- matrix(NA, nrow(eem$x), ncol(eem$x))
      eem$x <- cbind(x1, x2) %>% array(c(nrow(x1), ncol(x2),
                                         2)) %>% apply(1:2, mean, na.rm = TRUE)
    }
    if (type == 3) {
      eem$x <- eem$x %>% apply(2, function(col) pracma::pchip(xi = eem$em[!is.na(col)],
                                                              yi = col %>% na.omit(), x = eem$em, ...))
    }
    if (type == 0) {
      eem$x[is.na(eem$x)] <- 0
    }
    if (nonneg)
      eem$x[eem$x < 0] <- 0
    eem
  }
  stopCluster(cl)
  class(eem_list) <- "eemlist"
  eem_list
}