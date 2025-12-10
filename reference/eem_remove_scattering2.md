# Removes rayleigh and raman scattering

Based on 'eem_remove_scattering' function from eemR package, added
functionality to cut different amounts from above and below the emission
line. This function will only remove one kind of scattering at a time.
Use on an class eem or eemlist. Should be performed on data that is
untransformed besides a blank subtraction.

## Usage

``` r
eem_remove_scattering2(
  eem,
  type = "raman",
  order = 1,
  up_width = 10,
  down_width = 10
)
```

## Arguments

- eem:

  an object of class eemlist or eem

- type:

  a string, either "raman" or "rayleigh"

- order:

  an integer number, either 1 (first order) or 2 (second order)

- up_width:

  an integer number specifying the width in nm for the cut above the
  emission line

- down_width:

  an integer number specifying the width in nm for the cut below the
  emission line
