# Save EEMs and absorbance data to file

Saves absorbance as a .csv file and EEMs as individual .csv files and
.rds file for use in future analysis.

## Usage

``` r
save_eems(eem_list, abs_data, meta, prjpath)
```

## Arguments

- eem_list:

  object of class eemlist or eem

- abs_data:

  a dataframe containing the absorbance data

- meta:

  the metadata associated with the samples

- prjpath:

  location to save the exported EEMs and absorbance files

## Details

Use 'clean_files' function before running this function to ensure file
structure is correct for pre-processing.
