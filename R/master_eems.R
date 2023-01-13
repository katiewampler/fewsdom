##
#' Takes raw EEMs and absorbance data from Aqualog, returns cleaned, processed data
#' Combines all the 'fewsdom' functions to make one function to process samples.
#'
#' @importFrom stringr str_replace_all
#' @param prjpath path to folder with raw samples
#' @param run_date date the samples were run, this should correspond to the file where all the samples are stored
#' @param doc_file a string of the file path of the DOC file, can be .xlsx or .csv
#' @param doc_sheet a string of the sheet name of the DOC results, only required if the DOC file is an .xlsx file
#' @param doc_column a numeric indicating which column the DOC results are stored in in the DOC file
#' @param name_column a numeric indicating which column the sample/site names are stored in the DOC file
#' @param nskip a numeric indicating a number of lines to skip when reading in the DOC file, optional
#' @param doc_delim a string indicating the separation between site name and other identifiers in the DOC data, it will only keep the first piece
#' @param meta_sheet a string of the metadata sheet name, only required if the metadata file is an .xlsx file
#' @param site_loc a vector indicating the start and end of the site name in the metadata data identifier
#' @param zip_files will create a zip file in the file directory of all the raw, non-renamed files from the Aqualog as a backup
#' @param process_file logical, if TRUE it will put a text file in the processed data folder named 'processing_tracking'
#' @param ... additional arguments passed to the 'get_doc', 'clean_files', 'abs_preprocess', 'load_eems', 'eem_process', 'plot_eems', 'get_indices', 'save_eems' functions
#'
#' @return saves processed EEMs as .csv and .rds files, absorbance as .csv and .rds, and metadata as .csv
#' creates plots for each sample and calculates indices.
#'
#' @export
#'
process_eems <- function(prjpath, run_date, doc_file, doc_sheet,
                         doc_column=7, name_column=4, nskip=3,
                         doc_delim="-", meta_sheet="log", site_loc=c(1,7),
                         zip_files=T,process_file=T, ...){
  date_path <- stringr::str_replace_all(run_date, "-", "_")
  meta_loc <- paste(prjpath,"/metatable_manual_", date_path, ".xlsx", sep="")

  meta <- get_doc(doc_file=doc_file,
                  doc_sheet=doc_sheet, doc_column=doc_column,
                  name_column=name_column, nskip=nskip,
                  doc_delim=doc_delim, meta_file=meta_loc,
                  meta_sheet=meta_sheet, site_loc=c(1,7), ...)

  ## Section 2: Read in Sample Log and renaming files -----
  #rename files and create and put in folders
  meta <- clean_files(prjpath=prjpath, meta_file=meta_loc,
                      meta_sheet = meta_sheet, zip_files=zip_files,...)

  #convert absorbance files from .dat to .csv
  abs_preprocess(prjpath=prjpath)

  ## Section 3: Load Data in R ------
  data<- load_eems(prjpath = prjpath)
  X <- data[[1]]
  X_blk <- data[[2]]
  Sabs <- data[[3]]

  ## Section 4: Process the EEM's ------
  data_process <- eem_proccess(prjpath=prjpath, eemlist=X, blanklist=X_blk, abs=Sabs,
                               process_file=process_file, meta=meta)
  X_clean <- data_process[[1]]
  abs_clean <- data_process[[2]]

  ## Section 5: Report the Data ------
  #create plots
  plot_eems(prjpath = prjpath, meta=meta, eem=X_clean, doc_norm=F)
  plot_eems(prjpath = prjpath, meta=meta, eem=X_clean, doc_norm=T,
            save_names="_DOC")

  #save indices
  get_indices(X_clean, abs_clean, meta, prjpath=prjpath, doc_norm="both",
              sampsascol=F, waves=NULL)

  #save raw files (eems, abs)
  save_eems(X_clean, abs_clean, meta, prjpath=prjpath)
}
