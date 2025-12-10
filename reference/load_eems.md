# Load renamed and organized Aqualog data into R

Takes the Aqualog files that have been renamed and put in separate
folders within the main directory via the 'clean_files' function and
loads them into R as eemlists. It will also check to see if maximum
absorbance is below 0.3.

## Usage

``` r
load_eems(prjpath)
```

## Arguments

- prjpath:

  The file path of the main file directory where data is located

## Examples

``` r
if (FALSE) { # \dontrun{
  raw_eem <- load_eems(prjpath = prjpath)
} # }
 X <- raw_eem[[1]] #the raw eem samples
 X_blk <- raw_eem[[2]] #the matching blanks
 Sabs <- raw_eem[[3]] #the matching absorbance data
```
