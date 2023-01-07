#' Plot EEM's plots as contour plots
#'
#' Creates a pretty EEM's contour plot with the ability to manually set the intensity
#' scale or use an automated one based on the minimum and maximum of the EEM. Option to
#' add Coble peak labels to the plot (based off peaks from Coble et al. 2014.
#'
#' @importFrom reshape2 melt
#' @importFrom pals parula ocean.haline cubicl kovesi.rainbow
#' @importFrom plyr round_any
#' @importFrom grDevices colorRampPalette
#' @import ggplot2
#'
#' @param eem An object of class eem, the data you want to plot.
#' @param manualscale A TRUE or FALSE indicating if you want to specify the minimum and maximum of the intensity scale manually, good for comparisons across samples.
#' @param manualmax Maximum value for intensity scale.
#' @param manualmin Minimum value used for intensity scale.
#' @param nbins The number of bins (and colors) used in the contour plot. Maximum of 24.
#' @param label_peaks A TRUE or FALSE indicating if you want the main coble peaks annotated on the plot.
#' @param prec An integer for the number of significant figures used for binning the intensity.
#' @param palette A character with the color palette to use. Currently 'parula', 'ocean.haline', 'cubicl', and 'kovesi.rainbow' from the 'pals' package are supported.

#' @export

ggeem2 <- function(eem, manualscale=F, manualmax=1.5, manualmin=0,
                      nbins = 8, label_peaks=F, prec=2, palette="parula"){
  .ceiling_dec <- function(x, level=-.place_val(x)) round(x + 5.0001*10^(-level-1), level) #ceiling function with decimal
  .ceiling_int <- function(x, level=.place_val(x), precision) signif(x + 5*10^(level-precision-1), precision) #ceiling function with signif for vals above 1
  .place_val <- function(x) {
    stopifnot(class(x)=="numeric")
    place <- 0
    if(x == 1){
      place <- 1
    }
    else if(x > 1){
      while(x > 1){
        x <- x/10
        place <- place + 1
      }
    }else{
      while(x < 1){
        x <- x*10
        place <- place - 1
      }}
    place

  } #finds place value of max to determine scale
  .pretty_levels <- function(max, min){
    if(max > 0){
      max <- .ceiling_int(max, precision=prec)
      place <- .place_val(max)
      labs_df <- data.frame(start=signif(seq(min,max, length.out = nbins+1), prec))
    }else{
      max <- .ceiling_dec(max)
      place <- .place_val(max)*-1
      labs_df <- data.frame(start=sprintf(paste("%.", place+1, "f", sep=""), round(seq(min,max, length.out = nbins+1), place+1)))

    }
    labs_df
  } #make pretty even labels

  stopifnot(nbins <= 24, class(eem) == "eem")
  #get color palette
  if(palette == "parula"){
    colors <- pals::parula(n=nbins)
  }else if(palette == "ocean.haline"){
    colors <- pals::ocean.haline(n=nbins)
  }else if(palette == "cubicl"){
    colors <- pals::cubicl(n=nbins)
  }else if(palette == "kovesi.rainbow"){
    colors <- pals::kovesi.rainbow(n=nbins)
  }else if(palette == "katie_pal"){
    katie <- colorRampPalette(c("#d9ead3","#274e13"))
    colors <- katie(nbins)
  }else{
    stop("Please choose 'parula', 'ocean.haline', 'cubehelix', or 'kovesi.rainbow'")
  }
  #turn into a dataframe with x,y,z
  df_wide <- as.data.frame(eem$x)
  colnames(df_wide) <- eem$ex
  df_wide$em <- eem$em
  df <- reshape2::melt(df_wide, id.vars="em")
  df <- df[,c(2,1,3)]
  colnames(df) <- c("x", "y", "z")
  df$x <- as.numeric(as.character(df$x))

  #get things for plotting
  if(manualscale == F){
    min <- min(df$z, na.rm=T)
    max <- max(df$z, na.rm=T)
  }else{
    if(max(df$z, na.rm=T) > manualmax){
      message(paste("Warning: the maximium fill value is lower than your highest point of", round(max(df$z, na.rm=T),2)))
    }
    min <- manualmin
    max <- manualmax
  }

  #create labels
  labs_df <- .pretty_levels(max, min)
  labs_df$end <- c(labs_df$start[2:nrow(labs_df)], NA)
  if(max < 1){
    labs_df$start <- sprintf(paste("%.", (max(nchar(labs_df$start))-2), "f", sep=""), labs_df$start)
    labs_df$end <- sprintf(paste("%.", (max(nchar(labs_df$end), na.rm=T)-2), "f", sep=""), labs_df$end)
  }
  labs_df$label <- paste(as.character(labs_df$start), as.character(labs_df$end), sep="-")
  labs_df$label[nrow(labs_df)] <- paste(labs_df$start[nrow(labs_df)], "+", sep="")
  labs <- labs_df$label
  breaks <- as.numeric(labs_df$start)

  #get plot scale
  x_min <- min(df$x)
  x_max <- max(df$x)
  x_range <- plyr::round_any(((x_max - x_min) / 6), 25, f=round)
  x_min <-plyr::round_any(x_min, x_range, f = floor)
  x_max <-plyr::round_any(x_max, x_range, f = floor)


  y_min <- min(df$y)
  y_max <- max(df$y)
  y_range <- plyr::round_any(((y_max - y_min) / 6), 25, f=round)
  y_min <-plyr::round_any(y_min, y_range, f = floor)
  y_max <-plyr::round_any(y_max, y_range, f = floor)

  #create plot
  plot <- ggplot2::ggplot(df, aes(x,y, z=z)) +
    ggplot2::geom_contour_filled(aes(z = z), breaks=breaks)+
    ggplot2::coord_cartesian(expand = FALSE) + geom_contour(breaks=breaks, color=colors[1], size=0.1) +
    ggplot2::labs(x="Excitation (nm)", y="Emission (nm)")+
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = element_text(colour = 1, size = 10),
          axis.title = element_text(colour = 1, size = 12),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.position = "right")+
    ggplot2::scale_fill_manual(values = colors,  labels = labs) +
    ggplot2::guides(fill = guide_legend(title.position = "right",direction = "vertical",
                               title.theme = element_text(angle = 90, size = 12, colour = "black"),
                               barheight = .5, barwidth = .95,
                               title.hjust = 0.5, raster = FALSE,
                               title = "Intensity (R.U.)", reverse=TRUE)) +
    ggplot2::scale_x_continuous(breaks = round(seq(x_min, x_max, by = x_range),1)) +
    ggplot2::scale_y_continuous(breaks = round(seq(y_min, y_max, by = y_range),1))

  if(label_peaks == T){
    plot_pk <- data.frame(peak=c("B", "T", "A", "M", "C", "D", "E","N"),
                                x=c(275, 275, 265, 315, 340, 390, 455,280),
                                y=c(310,335, 430,400, 450,509,521,370))

    plot_pk <- subset(plot_pk, plot_pk$x > x_min & plot_pk$x < x_max & plot_pk$y > y_min & plot_pk$y < y_max)

    plot <- plot + ggplot2::annotate("text", x = plot_pk$x, y = plot_pk$y, label = plot_pk$peak, size=3.2) +
      ggplot2::annotate("rect", xmin=254, xmax = 260, ymin = 380, ymax=480, color="black", fill=NA)+
      ggplot2::annotate("rect", xmin=270, xmax = 280, ymin = 300, ymax=320, color="black", fill=NA)+
      ggplot2::annotate("rect", xmin=330, xmax = 350, ymin = 420, ymax=480, color="black", fill=NA)+
      ggplot2::annotate("rect", xmin=310, xmax = 320, ymin = 380, ymax=420, color="black", fill=NA)+
      ggplot2::annotate("rect", xmin=270, xmax = 280, ymin = 320, ymax=350, color="black", fill=NA)

  }
  plot
}
