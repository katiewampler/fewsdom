#' Convert absorbance data from a .dat file to .csv file
#'
#' Takes a folder of .dat absorbance files and converts them to .csv in the
#' correct format for reading in using the the 'absorbance_read' function from
#' the staRdom package.
#'
#' Use 'clean_files' function before running this function to ensure file structure is correct for pre-processing.
#'
#' @importFrom stringr str_split
#'
#' @param prjpath a string indicating the project file containing the absorbance data to process
#' @param runtype indicates how data was run on the Aqualog, either "manual", "sampleQ", or "mixed". If "mixed" it will get run type from the metadata
#' @param meta the metadata table for the sample run, only required if runtype is "mixed"
#' @export
#' @examples
#' \dontrun{
#'   abs_preprocess(prjpath=prjpath, "mixed", meta)
#' }

abs_preprocess <- function(prjpath, runtype="sampleQ", meta){
  stopifnot(is.character(prjpath)|runtype %in% c("manual","sampleQ", "mixed")|
              file.exists(prjpath))
  abs_files <- list.files(paste(prjpath, "/1_Absorbance", sep="")) #get files to process
  abs_files_done <- list.files(paste(prjpath, "/4_Clean_Absorbance", sep="")) #get names of finished files

  if(runtype == "mixed"){
    meta$unique_ID <- paste(meta$data_identifier, "_", meta$replicate_no, "_", meta$integration_time_s, "s", sep="")
  }

  #checks to see if all there's files that haven't been processed yet
  if(length(abs_files) > length(abs_files_done)){
    names <- sapply(stringr::str_split(abs_files, "_Abs.dat"),"[[",1) #get file names only
    if(file.exists(paste(prjpath,  "/4_Clean_Absorbance", sep=""))==F){
      stop("Invalid file structure, please use 'create_files' function to create file for plots within file directory")
    }
    #go through each sample. check, clean, and save
    for(x in names){
      #load data
      df <- read.table(paste(prjpath,"/1_Absorbance/", x, "_Abs.dat",sep=""), sep="\t", header = F, skip=3)

      if(runtype == "mixed"){
        samp_runtype <- meta$run_type[which(meta$unique_ID == x)]
      }else if(runtype %in% c("sampleQ","sampleq")){
        samp_runtype <- "sampleQ"
      }else if(runtype %in% c("manual","Manual")){
        samp_runtype <- "manual"
      }else{
        stop("Unrecognized sample type, please choose 'mixed', 'sampleQ', or 'manual'")
      }

      #if run manually it has extra rows that aren't needed
      if(samp_runtype %in% c("manual", "Manual")){
        df <- df[,c(1,10)]
      }
      if(nrow(df) > 250){
        warning(paste("Your absorbance data",  x, "has too many rows, check if transmittence data got added. Data was clipped to the correct rows.", sep=" "))
        df <- df[1:250,]
        df$V10 <- as.numeric(df$V10)
        df$V1 <- as.numeric(df$V1)
      }
      df <- subset(df, df$V1 <= 791) #it also doesn't end at the same spot, which throws an error if you have both manual and sample Q samples

      colnames(df) <- c("Wavelength", "Abs")
      write.table(df, paste(prjpath, "/4_Clean_Absorbance/", x, ".csv", sep=""), row.names = F, col.names = F, quote = F, sep=",")
    }
  }
}
