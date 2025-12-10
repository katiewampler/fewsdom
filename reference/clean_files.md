# Renames and organizes files from default Aqualog files

Takes the default files from the Aqualog (absorbance, blank, and sample)
it renames the files to the sample identifier and creates folder for
each file type in the main file directory.

## Usage

``` r
clean_files(prjpath, meta_file, meta_sheet = "log", ...)
```

## Arguments

- prjpath:

  the file path of the main file directory where data is located

- meta_file:

  file location of metadata table, used to determine how samples were
  run

- meta_sheet:

  optional, sheet name of the metadata only required if metadata is an
  .xlsx file

- ...:

  arguments to pass down to functions with 'clean_files'

## Value

metadata table with an added column for the unique ID

## Examples

``` r
if (FALSE) { # \dontrun{
meta <- clean_files(prjpath=prjpath, meta_file=meta_file,
meta_sheet = meta_sheet)
} # }
```
