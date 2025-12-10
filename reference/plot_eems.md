# Exports EEMs contour plots to a file

Will plot EEMs using 'ggeem2' function then will export plots to
specified file.

## Usage

``` r
plot_eems(
  prjpath,
  meta,
  eem,
  sing_plot = T,
  sum_plot = T,
  doc_norm = T,
  ncol = 4,
  title = "description",
  save_names = "",
  ...
)
```

## Arguments

- prjpath:

  a string indicating the project file

- meta:

  the metadata associated with the samples

- eem:

  object of class eemlist with EEMs to be plotted

- sing_plot:

  logical, if TRUE will create and save a separate plot for each EEMs

- sum_plot:

  logical, if TRUE will create and save a plot with all the EEMs in a
  single plot

- doc_norm:

  either TRUE, FALSE. TRUE return EEMs normalized by DOC concentration

- ncol:

  number of columns used in summary plot with all the EEMs samples

- title:

  either a dataframe with samples and titles, 'description', 'sampleID',
  or FALSE (see details)

- save_names:

  optional, a character that will be added to the unique ID for plot
  file names

- ...:

  arguments to pass to 'ggeems2' function

## Details

Use 'create_files' function before running this function to ensure file
structure is correct for pre-processing.

For title argument, you can specify the title manually, by using a
dataframe where the first column in the unique ID's for the samples and
the second column is the title. If 'description' is chosen it will use
the description column of the metadata. If 'sampleID' is chosen it will
use the unique identifier column from the metadata. Lastly if FALSE is
chosen no titles will be added.

## Examples

``` r
if (FALSE) { # \dontrun{
  plot_eems(prjpath = prjpath, meta=meta, eem=X_clean, doc_norm=F)
} # }
```
