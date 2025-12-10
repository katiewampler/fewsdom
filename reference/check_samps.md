# Checks to make sure samples and metadata match

Takes the loaded EEMs and absorbance data and tries to match it to the
'unique_ID' in the metadata, if a sample is in the metadata but not in
the loaded EEMs and absorbance it will remove that line in the metadata.
If there are samples loaded not in the metadata it will warn and remove
from the loaded samples.

## Usage

``` r
check_samps(meta, eemlist, blanklist, abs)
```

## Arguments

- meta:

  metadata table with sample information, needs to have a 'unique_ID'
  column

- eemlist:

  object of class eemlist containing EEMs samples

- blanklist:

  object of class eemlist containing the blanks for the EEMs samples,
  names should be same as eemlist with "\_blank" at the end

- abs:

  dataframe containing absorbance data corresponding to the EEMs samples

## Value

a list of metadata, eemlist, blanklist, and absorbance
