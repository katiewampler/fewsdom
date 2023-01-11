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
#' @param z_unit either "RU" or "DOC", use "RU" for raman normalized data and "DOC" for raman and DOC normalized data
#' @export

ggeem2 <- function(eem, manualscale=F, manualmax=1.5, manualmin=0,
                      nbins = 8, label_peaks=F, prec=2, palette="parula",
                   z_unit="RU"){
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
  plot
}


#' Plot EEM's plots as contour plots
#'
#' Creates a pretty EEM's contour plot with the ability to manually set the intensity
#' scale or use an automated one based on the minimum and maximum of the EEM. Option to
#' add Coble peak labels to the plot (based off peaks from Coble et al. 2014.
#'
#' @param prjpath a string indicating the project file
#' @param meta the metadata associated with the samples
#' @param eem object of class eemlist with EEMs to be plotted
#' @param sing_plot logical, if TRUE will create and save a separate plot for each EEMs
#' @param sum_plot logical, if TRUE will create and save a plot with all the EEMs in a single plot
#' @param doc_norm either TRUE, FALSE, or "both". TRUE or "both" will return EEMs normalized by DOC concentration
#' @param title either a vector the same length as the number of EEMs with titles, or "description" or "sampleID"
#' @param ... arguments to pass to 'ggeems2' function
#' @export
#'

plot_eems <- function(prjpath, meta, eem, sing_plot=T, sum_plot=T, doc_norm=T,
                      title="description", ...){
  #check inputs
  stopifnot(.is_eemlist(eem) | file.exists(prjpath) | is.data.frame(meta) |
            is.logical(sing_plot) | is.logical(sum_plot) |
            is.vector(title) | doc_norm %in% c(TRUE, FALSE, "both"))

  #create spot to put data
  if(file.exists(paste(prjpath,  "/5_Processed/Figures", sep=""))==F){
    stop("Invalid file structure, please use 'create_files' function to create file for plots within file directory")
  }
  save_spot <- paste(prjpath,  "/5_Processed/Figures", sep="")

  #if there's no DOC data, don't run DOC metrics/plots even if asked
  if(sum(is.na(meta$DOC_mg_L)) == nrow(meta)){
    doc_norm <- F
    warning("Metadata didn't contain any DOC data, unable to create DOC normalized plots.")}

  #account for doc normalization
  if(doc_norm ==T){
    ### stopped here ### now eems are marked with doc normalized, but need to figure out plotting
    #don't normalize again if already normalized, but do normalize if they haven't been
    if(doc_norm == T){
      #remove EEMs with no DOC data (or value of 0)
      EEM_rm <- meta$unique_ID[meta$DOC_mg_L == 0 | is.na(meta$DOC_mg_L) == T]
      X_DOC <- eem_exclude(X_DOC,
                           exclude=list("ex"=c(), "em"=c(),"sample"= EEM_rm ))

      #remove absorbance with no DOC data (or value of 0)
      abs_rm <- which((colnames(Sabs_DOC) %in% EEM_rm) == T)
      if(length(abs_rm) > 0){
        Sabs_DOC <-  Sabs_DOC[,-(abs_rm)]
      }

      #EEMs
      for (x in 1:length(X_DOC)){
        eem_name <- X_DOC[[x]]$sample
        eem_index <- which(eem_name == meta$unique_ID)
        X_DOC[[x]]$x <- X_DOC[[x]]$x / as.numeric(meta$DOC_mg_L[eem_index])
      }

      #Absorbance
      for(x in colnames(Sabs_DOC)[2:ncol(Sabs_DOC)]){
        col_num <- which(x == colnames(Sabs_DOC))
        abs_index <- which(x == meta$unique_ID)
        doc <- as.numeric(meta$DOC_mg_L[abs_index])
        Sabs_DOC[,col_num] <- Sabs_DOC[,col_num] / doc}
      write.table(paste(Sys.time(), "- Absorbance and EEM's were normalized by DOC concentration", sep=""), process_file_name, append=T, quote=F, row.names = F, col.names = F)

    }

    sing_plots_DOC <- lapply(X_DOC_clip, ggeem2, include_peaks=T)
    write.table(paste(Sys.time(), "- Plots normalized to DOC concentrations were created", sep=""), process_file, append=T, quote=F, row.names = F, col.names = F)}

  sing_plots <- lapply(X_clip, ggeem2, include_peaks=T)
  write.table(paste(Sys.time(), "- Plots NOT normalized to DOC concentrations were created", sep=""), process_file, append=T, quote=F, row.names = F, col.names = F)

  # PNGs of plots are written to output directory
  if(individual == T){
    for(n in 1:length(sing_plots)){
      pl <- sing_plots[[n]]
      name <- eem_names(X_clip)[n]
      title <- meta$description[n]
      png(paste0(name,'.png',sep=""),width=20, height = 15, units = 'cm', res = 300)
      plotp <- pl
      if(plottitle == T){
        plotp <- plotp + labs(title=title) + theme(plot.title = element_text(size=8))
      }
      suppressWarnings(print(plotp))
      dev.off()
    }

    if(DOC_norm_plot == T & DOC_norm ==T){
      for(n in 1:length(sing_plots_DOC)){
        pl <- sing_plots_DOC[[n]]
        name <- eem_names(X_DOC_clip)[n]
        title <- meta$description[which(meta$unique_ID== name)]
        png(paste0("DOC_", name,'.png',sep=""),width=20, height = 15, units = 'cm', res = 300)
        plotp <- pl
        if(plottitle == T){
          plotp <- plotp + labs(title=title) + theme(plot.title = element_text(size=8))
        }
        suppressWarnings(print(plotp))
        dev.off()
      }
    }
  }

  if(overall == T){
    overall_plots <- sing_plots
    if(plottitle == T){
      for(n in 1:length(overall_plots)){
        plotp <- overall_plots[[n]]
        name <- meta$data_identifier[n]
        title <- meta$description[n]
        if(overall_title =="description"){
          plotp <- plotp + labs(title=title) + theme(plot.title = element_text(size=9))
        }else if(overall_title == "sampleID"){
          plotp <- plotp + labs(title=name) + theme(plot.title = element_text(size=14))
        }
        overall_plots[[n]] <- plotp
      }}
    width <- overall_ncol * 20
    height <- ceiling(length(overall_plots)/overall_ncol) * 15
    png("EEM_summary_plot.png",width=width, height = height, units = 'cm', res = 300)
    suppressWarnings(grid.arrange(grobs = overall_plots, ncol=overall_ncol))
    dev.off()

    #get DOC summary
    if(DOC_norm_plot == T & DOC_norm ==T){
      overall_plots <- sing_plots_DOC
      if(plottitle == T){
        for(n in 1:length(overall_plots)){
          plotp <- overall_plots[[n]]
          name <- eem_names(X_DOC_clip)[n]
          title <- meta$description[which(meta$unique_ID==name)]
          if(overall_title =="description"){
            plotp <- plotp + labs(title=title) + theme(plot.title = element_text(size=9))
          }else if(overall_title == "sampleID"){
            plotp <- plotp + labs(title=name) + theme(plot.title = element_text(size=14))
          }
          overall_plots[[n]] <- plotp
        }}
      width <- overall_ncol * 20
      height <- ceiling(length(overall_plots)/overall_ncol) * 15
      png("EEM_summary_plot_DOC.png",width=width, height = height, units = 'cm', res = 300)
      suppressWarnings(grid.arrange(grobs = overall_plots, ncol=overall_ncol))
      dev.off()
    }

  }

}
