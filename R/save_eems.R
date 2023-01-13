#' Save EEMs and absorbance data to file
#'
#' Saves absorbance as a .csv file and EEMs as individual .csv files and .rds file
#' for use in future analysis.
#'
#' Use 'clean_files' function before running this function to ensure file structure is correct for pre-processing.
#'
#' @param eem_list object of class eemlist or eem
#' @param abs_data a dataframe containing the absorbance data
#' @param prjpath location to save the exported EEMs and absorbance files
#' @param meta the metadata associated with the samples
#'
#' @export

save_eems <- function(eem_list, abs_data, meta, prjpath){
  #saving location
  if(file.exists(paste(prjpath,  "/5_Processed", sep=""))==F){
    stop("Invalid file structure, please use 'create_files' function to create file for plots within file directory")
  }

  save_loc <- paste(prjpath, "/5_Processed/", sep="")
  #save absorbance data
  write.csv(abs_data, paste(save_loc,"Processed_Absorbance.csv", sep=""), quote=F, row.names = F)

  #save EEMs
  lapply(1:nrow(meta), function(f){
    file_name <- meta$unique_ID[f]
    eem_file <- eem_list[[f]]
    eem_data <- eem_file$x
    eem_data <- as.data.frame(eem_data)
    row.names(eem_data) <- eem_file$em
    colnames(eem_data) <- eem_file$ex
    write.csv(eem_data, paste(save_loc, file_name, ".csv", sep=""), quote=F, row.names = T)
  })

  #save items as R data objects
  saveRDS(eem_list, paste(save_loc, "processed_eemlist.rds", sep=""))
  saveRDS(abs_data, paste(save_loc, "processed_absdata.rds", sep=""))
  saveRDS(meta, paste(save_loc, "metadata.rds", sep=""))

}
