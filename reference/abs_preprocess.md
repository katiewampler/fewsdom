# Convert absorbance data from a .dat file to .csv file

Takes a folder of .dat absorbance files and converts them to .csv in the
correct format for reading in using the the 'absorbance_read' function
from the staRdom package.

## Usage

``` r
abs_preprocess(prjpath, runtype = "mixed", meta)
```

## Arguments

- prjpath:

  a string indicating the project file containing the absorbance data to
  process

- runtype:

  indicates how data was run on the Aqualog, either "manual", "sampleQ",
  or "mixed". If "mixed" it will get run type from the metadata

- meta:

  the metadata table for the sample run, only required if runtype is
  "mixed"

## Details

Use 'clean_files' function before running this function to ensure file
structure is correct for pre-processing.

## Examples

``` r
if (FALSE) { # \dontrun{
  abs_preprocess(prjpath=prjpath, "mixed", meta)
} # }
```
