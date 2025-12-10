# Input DOC data into metadata

Takes data from a separate table of DOC results and puts it in the
metadata table under the correct sample. The code assumes DOC is
reported in mg/L. See 'sample_metadata' in 'fewsdom' package for
formatting.

## Usage

``` r
get_doc(
  doc_file,
  doc_sheet = NULL,
  doc_column,
  name_column,
  nskip = 0,
  doc_delim = "-",
  meta_file,
  meta_sheet = NULL,
  site_loc = c(1, 7),
  rewrite = T,
  ...
)
```

## Arguments

- doc_file:

  a string of the file path of the DOC file, can be .xlsx or .csv

- doc_sheet:

  a string of the sheet name of the DOC results, only required if the
  DOC file is an .xlsx file

- doc_column:

  a numeric indicating which column the DOC results are stored in in the
  DOC file

- name_column:

  a numeric indicating which column the sample/site names are stored in
  the DOC file

- nskip:

  a numeric indicating a number of lines to skip when reading in the DOC
  file, optional

- doc_delim:

  a string indicating the separation between site name and other
  identifiers in the DOC data, it will only keep the first piece

- meta_file:

  a string indicating the file path of the metadata, can be .xlsx or
  .csv

- meta_sheet:

  a string of the metadata sheet name, only required if the metadata
  file is an .xlsx file

- site_loc:

  a vector indicating the start and end of the site name in the metadata
  data identifier

- rewrite:

  a logical, if TRUE original metadata will be saved over with metadata
  with DOC results, if FALSE it will add "\_doc_added" to the end of the
  table

- ...:

  arguments to pass down to functions with 'get_doc'

## Details

The table with DOC data can be either a .csv or excel file.
