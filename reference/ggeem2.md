# Plot EEMs plots as contour plots

Creates a pretty EEMs contour plot with the ability to manually set the
intensity scale or use an automated one based on the minimum and maximum
of the EEM. Option to add Coble peak labels to the plot (based off peaks
from Coble et al. 2014).

## Usage

``` r
ggeem2(
  eem,
  manualscale = F,
  manualmax = 1.5,
  manualmin = 0,
  nbins = 8,
  label_peaks = F,
  prec = 2,
  palette = "parula",
  z_unit = "RU"
)
```

## Arguments

- eem:

  an object of class eem or eemlist, the data you want to plot

- manualscale:

  a logical indicating if you want to specify the minimum and maximum of
  the intensity scale manually, good for comparisons across samples

- manualmax:

  maximum value for intensity scale

- manualmin:

  minimum value used for intensity scale

- nbins:

  the number of bins (and colors) used in the contour plot, maximum of
  24

- label_peaks:

  a logical indicating if you want the main coble peaks annotated on the
  plot

- prec:

  an integer for the number of significant figures used for binning the
  intensity

- palette:

  a character with the color palette to use. Currently 'parula',
  'ocean.haline', 'cubicl', and 'kovesi.rainbow' from the 'pals' package
  are supported (and custom 'katie_pal')

- z_unit:

  either "RU" or "DOC", use "RU" for raman normalized data and "DOC" for
  raman and DOC normalized data

## Value

If class is eem, will return a single plot, if class is an eemlist will
return a list of plots
