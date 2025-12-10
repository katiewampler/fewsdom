# Rename default file names from Aqualog

Takes the raw files downloaded from the Aqualog and renames them so file
names are consistent across all run types and runs for consistent
analysis and processing later.

## Usage

``` r
files_rename(meta, prjpath)
```

## Arguments

- meta:

  the metadata table for the sample run

- prjpath:

  a string indicating the project file containing the data to process,
  they should not be in subfiles
