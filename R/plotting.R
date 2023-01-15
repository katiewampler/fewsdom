#' Plot EEMs plots as contour plots
#'
#' Creates a pretty EEMs contour plot with the ability to manually set the intensity
#' scale or use an automated one based on the minimum and maximum of the EEM. Option to
#' add Coble peak labels to the plot (based off peaks from Coble et al. 2014).
#'
#' @importFrom reshape2 melt
#' @importFrom pals parula ocean.haline cubicl kovesi.rainbow
#' @importFrom plyr round_any
#' @importFrom grDevices colorRampPalette
#' @importFrom gridExtra grid.arrange
#' @importFrom grDevices dev.off png
#' @import ggplot2
#'
#' @param eem an object of class eem or eemlist, the data you want to plot
#' @param manualscale a logical indicating if you want to specify the minimum and maximum of the intensity scale manually, good for comparisons across samples
#' @param manualmax maximum value for intensity scale
#' @param manualmin minimum value used for intensity scale
#' @param nbins the number of bins (and colors) used in the contour plot, maximum of 24
#' @param label_peaks a logical indicating if you want the main coble peaks annotated on the plot
#' @param prec an integer for the number of significant figures used for binning the intensity
#' @param palette a character with the color palette to use. Currently 'parula', 'ocean.haline', 'cubicl', and 'kovesi.rainbow' from the 'pals' package are supported (and custom 'katie_pal')
#' @param z_unit either "RU" or "DOC", use "RU" for raman normalized data and "DOC" for raman and DOC normalized data
#' @returns If class is eem, will return a single plot, if class is an eemlist will return a list of plots
#' @export

ggeem2 <- function(eem, manualscale=F, manualmax=1.5, manualmin=0,
                      nbins = 8, label_peaks=F, prec=2, palette="parula",
                      z_unit="RU"){
  stopifnot(nbins <= 24 |.is_eem(eem) | .is_eemlist(eem) |
              is.logical(c(manualscale, label_peaks))|
              is.numeric(c(nbins, prec, nbins))| z_unit %in% c("RU", "DOC")|
              palette %in% c("parula", "ocean.haline","cubicl", "kovesi", "katie_pal"))

  if (.is_eemlist(eem)) {
    res <- lapply(eem, ggeem2, manualscale=manualscale,
                  manualmax=manualmax, manualmin=manualmin,
                  nbins =nbins, label_peaks=label_peaks,
                  prec=prec, palette=palette,
                  z_unit=z_unit)
    return(res)
  }
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

  #get z unit label
  if(z_unit == "RU"){
    plot_z <- expression(Intensity~(R.U.))
  }else if(z_unit == "DOC"){
    plot_z <- expression(Intensity~(R.U.~mgC~L^{-1}))
  }else{
    stop("Invalid intensity unit, please choose 'RU' or 'DOC'")
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
                               title = plot_z, reverse=TRUE)) +
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
  return(plot)
}


#' Exports EEMs contour plots to a file
#'
#' Will plot EEMs using 'ggeem2' function then will export plots to specified file.
#'
#' Use 'create_files' function before running this function to ensure file structure is correct for pre-processing.
#'
#' @param prjpath a string indicating the project file
#' @param meta the metadata associated with the samples
#' @param eem object of class eemlist with EEMs to be plotted
#' @param sing_plot logical, if TRUE will create and save a separate plot for each EEMs
#' @param sum_plot logical, if TRUE will create and save a plot with all the EEMs in a single plot
#' @param doc_norm either TRUE, FALSE. TRUE return EEMs normalized by DOC concentration
#' @param title either a dataframe with samples and titles, 'description', 'sampleID', or FALSE (see details)
#' @param ncol number of columns used in summary plot with all the EEMs samples
#' @param save_names optional, a character that will be added to the unique ID for plot file names
#' @param ... arguments to pass to 'ggeems2' function
#' @details For title argument, you can specify the title manually, by using a dataframe where the
#' first column in the unique ID's for the samples and the second column is the title. If 'description'
#' is chosen it will use the description column of the metadata. If 'sampleID' is chosen it will
#' use the unique identifier column from the metadata. Lastly if FALSE is chosen no titles will be added.
#' @export
#'

plot_eems <- function(prjpath, meta, eem, sing_plot=T, sum_plot=T, doc_norm=T,
                      ncol=4, title="description", save_names = "", ...){
  .add_title <- function(eem, title, meta){
    order <- eem_names(eem)
    if(is.data.frame(title)==T){
      index <- sapply(order, function(x) which(x == title[,1]))
      plot_title <- title[index,2]
    }else if(title =="sampleID"){
      index <- sapply(order, function(x) which(x == meta$unique_ID))
      plot_title <- meta$data_identifier[index]
    }else if(title == "description"){
      index <- sapply(order, function(x) which(x == meta$unique_ID))
      plot_title <- meta$description[index]
    }else{plot_title <- NULL}
    return(plot_title)
  }
  .save_plots <- function(n, title){
    pl <- raw_plots[[n]]
    name <- eem_names(eem)[n]
    grDevices::png(paste0(save_spot, "/", name, save_names, ".png",sep=""),width=20, height = 15, units = 'cm', res = 300)
    if(title != F){
      pl_title <- .add_title(eem[[n]], title, meta)
      pl <- pl + labs(title=pl_title) + theme(plot.title = element_text(size=10))
    }
    suppressWarnings(print(pl))
    grDevices::dev.off()
  }
  .title_plots <- function(n, title){
    pl <- raw_plots[[n]]
    if(title != F){
      pl_title <- .add_title(eem[[n]], title, meta)
      pl <- pl + labs(title=pl_title) + theme(plot.title = element_text(size=10))
    }
  }
  .is_zero <- function(eem){
    zero <- sum(eem$x, na.rm=T)
    if(zero == 0){eem$sample}
  }

  #check inputs
  stopifnot(.is_eemlist(eem) | file.exists(prjpath) | is.data.frame(meta) |
            is.logical(sing_plot) | is.logical(sum_plot) |
              doc_norm %in% c(TRUE, FALSE, "both") | file.exists(prjpath))

  #check EEM's for empty EEM's, zero EEM's or fully NA EEM's. which won't plot
   empty <- empty_eems(eem, verbose = F)
   zero <- unlist(sapply(1:length(eem), function(x){
     .is_zero(eem[[x]])}))

   eem_rm <- c(empty, zero)
   if(is.null(eem_rm) == F){
     eem <- eem_exclude(eem, exclude=list(sample=eem_rm))
     warning("The following samples were exluded from plotting because they had only missing or zero data: \n",
             paste(eem_rm, collapse ="\n"))
   }

  #have spot to put data
  if(file.exists(paste(prjpath,  "/5_Processed/Figures", sep=""))==F){
    stop("Invalid file structure, please use 'create_files' function to create file for plots within file directory")
  }
  save_spot <- paste(prjpath,  "/5_Processed/Figures", sep="")

  #if there's no DOC data in metadata, don't run DOC metrics/plots even if asked
  if(sum(is.na(meta$DOC_mg_L)) == nrow(meta)){
    doc_norm <- F
    warning("Metadata didn't contain any DOC data, unable to create DOC normalized plots.")}

  #account for doc normalization
  if(doc_norm ==T){eem <- .eem_doc_norm(eem, meta)}

  #remove doc normalization if doc_norm is false but sample are doc_normalized
  if(doc_norm ==F){eem <- .eem_doc_rm(eem, meta)}

  z_unit <- ifelse(doc_norm ==T, "DOC", "RU")
  raw_plots <- ggeem2(eem, z_unit=z_unit, ...)

  # PNGs of plots are written to output directory
  if(sing_plot == T){sapply(1:length(raw_plots), .save_plots, title=title)}
  if(sum_plot == T){
    sum_pl <- lapply(1:length(raw_plots), .title_plots, title=title)
    width <- ncol * 20
    height <- ceiling(length(sum_pl)/ncol) * 15
    grDevices::png(paste(save_spot, "/EEM_summary_plot", save_names, ".png", sep=""),width=width, height = height, units = 'cm', res = 300)
    suppressWarnings(gridExtra::grid.arrange(grobs = sum_pl, ncol=ncol))
    grDevices::dev.off()}
  }
