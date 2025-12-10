# Finds optimal scattering cut widths based on EEM set

Takes a group of EEMs and uses peak picking to determine how wide the
cut should be for each scattering line (does not work for 2nd order
raman). Uses the excitation at 280nm and finds the start and end of the
peak. If the function can't find one, it will default to NA for that
sample. Median of the peak widths is reported for use in the
'remove_scattering' function. Should be performed on data that is
untransformed besides a blank subtraction.

## Usage

``` r
find_cut_width(eem, type = "rayleigh", order = 1)
```

## Arguments

- eem:

  an object of class eemlist

- type:

  a string, either "raman" or "rayleigh"

- order:

  an integer number, either 1 (first order) or 2 (second order)
