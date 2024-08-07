% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abs_preprocess.R
\name{abs_preprocess}
\alias{abs_preprocess}
\title{Convert absorbance data from a .dat file to .csv file}
\usage{
abs_preprocess(prjpath, runtype = "mixed", meta)
}
\arguments{
\item{prjpath}{a string indicating the project file containing the absorbance data to process}

\item{runtype}{indicates how data was run on the Aqualog, either "manual", "sampleQ", or "mixed". If "mixed" it will get run type from the metadata}

\item{meta}{the metadata table for the sample run, only required if runtype is "mixed"}
}
\description{
Takes a folder of .dat absorbance files and converts them to .csv in the
correct format for reading in using the the 'absorbance_read' function from
the staRdom package.
}
\details{
Use 'clean_files' function before running this function to ensure file structure is correct for pre-processing.
}
\examples{
\dontrun{
  abs_preprocess(prjpath=prjpath, "mixed", meta)
}
}
