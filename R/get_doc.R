#' Input DOC data into Metadata Table
#'
#' Takes data from a separate table of DOC results and puts
#' it in the metadata table under the correct sample. The code assumes
#' DOC is reported in mg/L. See 'sample_metadata' in 'fewsdom' package
#' for formatting.
#'
#' The table with DOC data can be either a .csv or excel file.
#'
#' @importFrom stringr str_detect str_split
#' @importFrom readxl read_xlsx
#' @importFrom openxlsx write.xlsx
#'
#' @param doc_file a string of the file path of the doc file, can be .xlsx or .csv
#' @param doc_sheet a string of the sheet name of the doc results, only required if the doc file is an .xlsx file
#' @param doc_column a numeric indicating which column the doc results are stored in in the doc file
#' @param name_column a numeric indicating which column the sample/site names are stored in the doc file
#' @param nskip a numeric indicating a number of lines to skip when reading in the doc file, optional
#' @param site_delim a string indicating the separation between site name and other identifiers in the doc data, it will only keep the first piece
#'
#' @param meta_file a string indicating the file path of the metadata, can be .xlsx or .csv
#' @param meta_sheet a string of the metadata sheet name, only required if the metadata file is an .xlsx file
#' @param meta_delim a string indicating the separate between site name and other identifiers in the metadata data identifier
#' @param rewrite A logical, if TRUE original metadata will be saved over with metadata with DOC results, if FALSE it will add "_doc_added" to the end of the table
#' @export

get_doc <- function(doc_file, doc_sheet, doc_column, name_column, nskip=0, site_delim="-",
                    meta_file, meta_sheet, meta_delim, rewrite=T){
  #read in doc data
    if(stringr::str_detect(doc_file, ".xlsx")){
      doc_result <- readxl::read_xlsx(doc_file, sheet = doc_sheet, skip=nskip )
    } else{
      doc_result <- read.csv(doc_file, skip = nskip)
    }

  #condense to just DOC data
    doc_result <- doc_result[,c(name_column, doc_column)]
    colnames(doc_result) <- c("Site", "DOC")

  #clean site names, we submit samples to the lab with date and time stamps, but we match just site name (for now)
    doc_result$Site <- sapply(stringr::str_split(doc_result$Site, site_delim),"[[",1)

  #load metadata
    meta_file <- paste(prjpath,"/metatable_manual_", datepath, ".xlsx", sep="")
    if(stringr::str_detect(meta_file, ".xlsx")){
      meta <- readxl::read_xlsx(meta_file, sheet = meta_sheet)
    } else{
      meta <- read.csv(doc_file)
    }

  #clean meta names
    meta$Site <- sapply(stringr::str_split(meta$data_identifier, meta_delim),"[[",1)

  #Put DOC in metadata
    meta$DOC_mg_L <- doc_result$DOC[match(meta$Site, doc_result$Site)]

    #check it "worked"
    if(sum(is.na(meta$DOC_mg_L)) == rnow(meta)){
      warning("No matching DOC results were found, please check your file names and delimeters")
    }

  #write metafile
    if(overwrite==T){
      if(stringr::str_detect(meta_file, ".xlsx")){
        openxlsx::write.xlsx(meta, meta_file, sheetName = meta_sheet,
                   colNames = TRUE, rowNames = F, append = F)} else{
        write.csv(meta, meta_file, col.names = T, row.names = F, quote=F)}
    }else{
      if(stringr::str_detect(meta_file, ".xlsx")){
        new_meta_file <- paste(unlist(str_split(meta_file, ".xlsx"))[1], "_doc_added.xlsx", sep="")
        openxlsx::write.xlsx(meta, new_meta_file, sheetName = meta_sheet,
                   colNames = TRUE, rowNames = F, append = F)} else{
      new_meta_file <- paste(unlist(str_split(meta_file, ".csv"))[1], "_doc_added.csv", sep="")
       write.csv(meta, new_meta_file, col.names = T, row.names = F, quote=F)}
    }
}