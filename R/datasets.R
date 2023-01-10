#' An example of ideal metadata formatting for a sample set run on the Aqualog
#'
#' Many of the functions within the 'fewsdom' package require metadata for
#' smooth use with the functions it is recommended to use this sample as a
#' template for the metadata.
#'
#' @format A data frame with 22 rows and 11 columns:
#' \describe{
#'   \item{index}{The numeric order of entry (e.g. 1, 2, 3, 4, etc…)}
#'   \item{analysis_date}{The date the samples were run on instrument (not collected in the field)}
#'   \item{description}{A brief description of the sample collected 60-character max, used for plot titles}
#'   \item{data_identifier}{the data identifer for the sample that was used in the Aqualog, needs to be exact}
#'   \item{replicate_no}{A numeric denoting analytical duplicates. They must have identical Sample identifiers}
#'   \item{integration_time_s}{Integration time of the sample (e.g. 1, 2, etc…)}
#'   \item{dilution}{Include the dilution factor here as decimal format (e.g. a 2-fold dilution with 1 part sample and 1 part water will have a dilution factor of 0.5). No dilution is a 1}
#'   \item{RSU_area_1s}{The RSU Adjust Area from the Raman Test}
#'   \item{run_type}{How the sample was run, should be either "SampleQ" or "Manual", used to rename and organize samples}
#'   \item{DOC_mg_L}{DOC concentration in mg/L of the original sample (not the diluted sample)}
#'   \item{Notes}{Not required, but can be used to note abnormal sample or experimental conditions}
#' }
"sample_metadata"
