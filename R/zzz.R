.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "The 'fewsdom' package has been deprecated and will no longer be updated.\n",
    "Please switch to the 'eemanalyzeR' package:\n",
    "https://github.com/katiewampler/eemanalyzeR"
  )
}
