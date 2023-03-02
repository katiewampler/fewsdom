##
#' Takes raw EEMs and absorbance data from Aqualog, returns cleaned, processed data
#' Combines all the 'fewsdom' functions to make one function to process samples.
#'
#' @importFrom stringr str_replace_all
#' @param prjpath path to folder with raw samples
#' @param meta_name the file name of the metadata with the extension
#' @param get_doc logical, if TRUE will use 'get_doc' function to match DOC data to metadata samples, only required if get_doc is TRUE
#' @param doc_file a string of the file path of the DOC file, can be .xlsx or .csv, only required if get_doc is TRUE
#' @param doc_sheet a string of the sheet name of the DOC results, only required if the DOC file is an .xlsx file, only required if get_doc is TRUE
#' @param doc_column a numeric indicating which column the DOC results are stored in in the DOC file, only required if get_doc is TRUE
#' @param name_column a numeric indicating which column the sample/site names are stored in the DOC file, only required if get_doc is TRUE
#' @param nskip a numeric indicating a number of lines to skip when reading in the DOC file, optional, only required if get_doc is TRUE
#' @param doc_delim a string indicating the separation between site name and other identifiers in the DOC data, it will only keep the first piece, only required if get_doc is TRUE
#' @param meta_sheet a string of the metadata sheet name, only required if the metadata file is an .xlsx file
#' @param site_loc a vector indicating the start and end of the site name in the metadata data identifier
#' @param process_file logical, if TRUE it will put a text file in the processed data folder named 'processing_tracking'
#' @param ... additional arguments passed to the 'get_doc', 'clean_files', 'abs_preprocess', 'load_eems', 'eem_process', 'plot_eems', 'get_indices', 'save_eems' functions
#'
#' @return saves processed EEMs as .csv and .rds files, absorbance as .csv and .rds, and metadata as .csv
#' creates plots for each sample and calculates indices.
#'
#' @details The following steps are chosen by default:
#'  -DOC is added to the metadata and the file is overwritten
#'
#'  -The file directory is created, raw files are zipped
#'
#'  -A process file is created, The instrument blank is subtracted from the
#'  samples, raman scattering is removed and interpolated, rayleigh scattering
#'  is removed but not interpolated, widths are determined automatically,
#'  samples are corrected for inner filter effects, raman normalized, and
#'  normalized for DOC, excitation is clipped from 247 to 450 nm, emission
#'  is clipped from 247 to 550 nm
#'
#'  -Sample are plotted individually and a summary plot using the parula
#'  color palette with descriptions for the plot titles, the peaks are not
#'  annotated on the plot. Plots are created for DOC normalized data (_DOC),
#'  and for non normalized data.
#'
#'  -Fluorescence indices are reported for data that has and hasn't been DOC
#'  normalized, sample are reported in rows.
#'
#'  If you want to change these defaults add the appropriate arguments into the
#'  function to change the defaults for the other functions.
#'
#' @export
#'
run_eems <- function(prjpath, meta_name, get_doc=T, doc_file, doc_sheet,
                         doc_column=7, name_column=4, nskip=3,
                         doc_delim="-", meta_sheet="log", site_loc=c(1,7),
                          process_file=T, ...){
  meta_file <- paste(prjpath,"/", meta_name, sep="")

  if(get_doc == T){
    meta <- get_doc(doc_file=doc_file,
                    doc_sheet=doc_sheet, doc_column=doc_column,
                    name_column=name_column, nskip=nskip,
                    doc_delim=doc_delim, meta_file=meta_file,
                    meta_sheet=meta_sheet, site_loc=site_loc, ...)
  }


  #rename files and create and put in folders
  cat("Renaming files and putting in files \n")
  meta <- clean_files(prjpath=prjpath, meta_file=meta_file,
                      meta_sheet = meta_sheet,...)

  #convert absorbance files from .dat to .csv
  abs_preprocess(prjpath=prjpath, "mixed", meta)

  #Load Data in R
  cat("Loading data in R \n")
  data<- load_eems(prjpath = prjpath)
  X <- data[[1]]
  X_blk <- data[[2]]
  Sabs <- data[[3]]

  #Check data with metadata, remove samples that don't have data
  test <- check_samps(meta, X, X_blk, Sabs)
  meta <- test[[1]]
  X <- test[[2]]
  X_blk <- test[[3]]
  Sabs <- test[[4]]

  ## Process the EEM's
  cat("Processing EEMs and absorbance data \n")
  data_process <- eem_proccess(prjpath=prjpath, eemlist=X, blanklist=X_blk, abs=Sabs,
                               process_file=process_file, meta=meta,...)
  X_clean <- data_process[[1]]
  abs_clean <- data_process[[2]]

  ## Report the Data
  cat("Reporting EEMs and absorbance data \n")

  #create plots
  plot_eems(prjpath = prjpath, meta=meta, eem=X_clean, doc_norm=F)

  if(sum(is.na(meta$DOC_mg_L)) < nrow(meta)){
    plot_eems(prjpath = prjpath, meta=meta, eem=X_clean, doc_norm=T,
              save_names="_DOC")
  }


  #save indices
  get_indices(X_clean, abs_clean, meta, prjpath=prjpath, doc_norm="both",
              sampsascol=F, waves=NULL)

  #save raw files (eems, abs)
  save_eems(X_clean, abs_clean, meta, prjpath=prjpath)

  cat(" \n EEM's have been successfully processed. Look in the 5_Processed folder for plots and indices \n")

}
