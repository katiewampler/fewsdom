# Get absorbance peaks and indices

Takes a dataframe of absorbance data and calculates optical properties
for each sample. If data was collected at greater than 1 nm intervals,
data will be interpolated to 1 nm intervals.

## Usage

``` r
abs_parm(abs_data, waves = NULL, meta, keep_all = F, cuvle = 1)
```

## Arguments

- abs_data:

  a dataframe containing the absorbance data, where each column is a
  sample

- waves:

  optional, a vector of wavelengths in nm to extract from the absorbance
  data

- meta:

  the metadata table for the sample run, should include DOC data in mg/L

- keep_all:

  a logical, TRUE will return all samples even those without DOC data,
  FALSE will only return samples with DOC data

- cuvle:

  the length of the cuvette in cm

## Value

Returns a data frame of indices.

## Details

Absorbance indices based on Hansen et al. 2016. Measurements are defined
as follows:

SUVA254, SUVA280, SUVA350, SUVA370: SUVA at 254, 280, 350, and 370 nm.
Units of \\\text{L mgC}^{-1} \text{m}^{-1}\\. Typically higher values
are associated with greater aromatic content.

SVA412, SVA440, SVA480, SVA510, SVA532, SVA555: SVA at 412, 440, 480,
510, 532, 555 nm. Units of \\\text{L mgC}^{-1} \text{m}^{-1}\\.
Typically higher values are associated with greater aromatic content.

S275_295: Spectral slope between 275 to 295 nm.

S350_400: Spectral slope between 350 to 400 nm.

Spectral slopes are found with a nonlinear fit of an exponential
function to the absorption spectrum, typically higher values are
associated with lower molecular weight materials and/or lower
aromaticity.

SR: Spectral slope S275_295 divided by spectral slope S350_400,
negatively correlated to DOM molecular weight and generally increases on
irradiation.

## References

Hansen, A. M., Kraus, T. E. C., Pellerin, B. A., Fleck, J. A., Downing,
B. D., & Bergamaschi, B. A. (2016). Optical properties of dissolved
organic matter (DOM): Effects of biological and photolytic degradation.
Limnology and Oceanography, 61(3), 1015–1032.
https://doi.org/10.1002/lno.10270
