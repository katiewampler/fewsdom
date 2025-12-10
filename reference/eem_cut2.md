# Trim EEMs to specified emission and/or excitation wavelengths

Modified from eemR 'eem_cut' function. The specified wavelengths are now
the ones you keep, not the ones you remove.

## Usage

``` r
eem_cut2(eem, ex, em, exact = FALSE, fill_with_na = FALSE)
```

## Arguments

- eem:

  An object of class eem or eemlist.

- ex:

  numerical of length two, with range of excitation wavelengths to keep

- em:

  numerical of length two, with range of emission wavelengths to keep

- exact:

  Logical. If TRUE, only wavelengths matching em and/or ex will be
  removed. If FALSE, all wavelengths in the range of em and/or ex will
  be removed.

- fill_with_na:

  Logical. If TRUE, fluorescence values at specified wavelengths will be
  replaced with NA. If FALSE, these values will be removed.

## Value

an object of class eem or eemlist
