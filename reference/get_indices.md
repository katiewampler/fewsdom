# Export absorbance and fluorencence indices and peaks to an excel file

Uses 'eem_coble_peaks2' and 'abs_parm' functions to obtain absorbance
and fluoresence indices for an eemlist. Saves data in an excel file.

## Usage

``` r
get_indices(
  eem_list,
  abs_data,
  meta,
  prjpath,
  doc_norm = "both",
  sampsascol = F,
  waves = NULL,
  cuvle = 1
)
```

## Arguments

- eem_list:

  an object of class eem or eemlist

- abs_data:

  a dataframe containing the absorbance data

- meta:

  the metadata associated with the samples

- prjpath:

  location to save the exported absorbance and fluorescence metrics

- doc_norm:

  either TRUE, FALSE or "both" indicating if peaks should be normalized
  by DOC concentrations

- sampsascol:

  a logical indicating how results should be oriented, TRUE puts samples
  as columns, FALSE puts samples as rows

- waves:

  optional, a vector of wavelengths in nm to extract from the absorbance
  data

- cuvle:

  the length of the cuvette in cm

## Details

Use 'clean_files' function before running this function to ensure file
structure is correct for pre-processing.
