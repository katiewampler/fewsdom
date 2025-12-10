# Process EEMs samples

Takes EEMs data with blanks and absorbance data and will perform
processing steps including blank subtraction, removing raman and
rayleigh scattering, dilution correction, removing inner filtering
effects, raman normalization.

## Usage

``` r
eem_proccess(
  prjpath,
  eemlist,
  blanklist,
  abs,
  meta,
  process_file = T,
  replace_blank = NULL,
  raman = T,
  rayleigh = T,
  IFE = T,
  raman_norm = T,
  dilute = T,
  ex_clip = c(247, 450),
  em_clip = c(247, 600),
  raman_width = "auto",
  raman_mask = c(8, 8, 1.5, 1.5),
  rayleigh_width = "manual",
  rayleigh_mask = c(20, 10, 10, 10),
  ...
)
```

## Arguments

- prjpath:

  the file path of the main file directory where data is located

- eemlist:

  object of class eem or eemlist containing EEMs samples

- blanklist:

  object of class eem or eemlist containing the blanks for the EEMs
  samples should be the same length as eem

- abs:

  dataframe containing absorbance data corresponding to the EEMs samples

- meta:

  dataframe of metadata containing unique ID's, integration time,
  dilutions, and raman area for each sample, see example metadata for
  format

- process_file:

  logical, if TRUE it will put a text file in the processed data folder
  named 'processing_tracking'

- replace_blank:

  a character giving the data identifier of the sample that should be
  used for blank subtraction

- raman:

  logical, if TRUE will use 'raman' function to remove raman scattering

- rayleigh:

  logical, if TRUE will use 'rayleigh' function to remove rayleigh
  scattering

- IFE:

  logical, if TRUE will use absorbance data to remove inner filter
  effects

- raman_norm:

  logical, if TRUE will normalize EEMs to raman unit area

- dilute:

  logical, if TRUE will correct for dilution factors given in metadata
  table

- ex_clip:

  vector of length two with the excitation wavelengths to clip the EEMs
  to

- em_clip:

  vector of length two with the emission wavelengths to clip the EEMs to

- raman_width:

  either "auto" or "manual". If auto is chosen cutting widths will be
  found using the 'find_cut_width' function

- raman_mask:

  a vector of length 4 specifying the width of the raman line to cut,
  numbers 1:2 are width above and below first order line, numbers 3:4
  are width above and below second order line, since you cannot use the
  auto method for second order raman this must be specified

- rayleigh_width:

  either "auto" or "manual", if auto is chosen cutting widths will be
  found using the 'find_cut_width' function

- rayleigh_mask:

  optional if auto width method is used, a vector of length 4 specifying
  the width of the rayleigh line to cut, numbers 1:2 are width above and
  below first order line, numbers 3:4 are width above and below second
  order line

- ...:

  arguments passed on to scattering functions 'raman' and 'rayleigh'

## Value

an list where the first object is of class eemlist with processed EEMs
samples. The second object is a dataframe with the processes absorbance
data. If a process file is given, a file will be created in the
processes folder of the file directory

## Examples

``` r
if (FALSE) { # \dontrun{
data_process <- eem_proccess(prjpath=prjpath, eemlist=X, blanklist=X_blk, abs=Sabs,
process_file=process_file, meta=meta)

X_clean <- data_process[[1]]
abs_clean <- data_process[[2]]} # }
```
