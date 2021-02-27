#'
#'
#'
#'
#'
#'
#' @export
rleBoxplots= function(rle_table_full, meta_qual_df, fill_column, fill_column_colors, table_name, feature){
  #' rle_table_full output of rleTableFull() above
  #' meta_qual_df is a metadata + quality assess dataframe
  #' table_name is the name of the graph
  #' part of the table name, the feature
  #' :params fill_column: column to color boxplots (make sure this column is factored)
  #' :params fill_column_colors: named vector assigning colors to factor levels of fill_column
  #' eg) c("0" = "#52854C", "1" = "#D16103")

  df = stack(rle_table_full)
  meta_qual_df$LIBRARYDATE = as.factor(meta_qual_df$LIBRARYDATE)

  stacked_rle_meta_qual_df = df %>% left_join(meta_qual_df, by=c('ind'='FASTQFILENAME'), copy=TRUE)

  ylim1 = boxplot.stats(df$values)$stats[c(1, 5)]

  # order by different variables in metadata (eg library date) <-- TODO: make this option

  # rle_boxplots = ggplot(stacked_rle_meta_qual_df, aes(reorder(ind, PROTEIN_CODING_COUNTED), values))+
  #   geom_boxplot(outlier.shape=NA)+
  #   coord_cartesian(ylim = ylim1*2.5)+
  #   theme(axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank())+
  #   xlab('Samples')+
  #   ylab('Deviation from Median')+
  #   ggtitle('90 Minute Induction RLE plots')
  #
  # plot(rle_boxplots)


  rle_boxplots = ggplot(stacked_rle_meta_qual_df, aes(ind, values))+
    geom_boxplot(outlier.shape=NA, aes(fill=!!rlang::sym(fill_column)))+
    scale_colour_manual(values = fill_column_colors, aesthetics="fill")+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    xlab('Samples')+
    ylab('Deviation from Median')+
    coord_cartesian(ylim = c(-3,3))+
    scale_y_continuous(n.breaks = 20)+
    geom_hline(yintercept = 0)+
    ggtitle(paste0(table_name, '_', feature))+theme(legend.position = "none")

  return(rle_boxplots)

} # end rleBoxPlots()

#' @export
rleScatterPlots = function(rle_summary_df, color_col, title){
  rle_scatter = ggplot(rle_summary_df, aes(ABS_SAMPLE_DEVIATION_MEDIAN, INTERQUARTILE_RANGE, color=!! rlang::sym(color_col)))+geom_point()+ggtitle(title) +xlim(0,1.5)+ylim(0,1.75)
  rle_scatter = ggMarginal(rle_scatter, type="density")
  #rle_scatter + ggMarginal(rle_scatter, type="density")
  # pdf(title)
  plot(rle_scatter)
  # dev.off()
}

#' @export
rleMedianBarPlot = function(rle_summary_df, title){
  ggplot(rle_summary_df, aes(ABS_SAMPLE_DEVIATION_MEDIAN))+
    geom_histogram()+ggtitle(title)+ylim(0,425)
  # dev.off()
}

#' @export
rleIQRBarPlot = function(rle_summary_df, title){
  ggplot(rle_summary_df, aes(INTERQUARTILE_RANGE))+
    geom_histogram()+ggtitle(title)+
    scale_x_continuous(breaks = seq(0, 2, by = .1))+ylim(0,125)
  # dev.off()
}

#' @export
createRLEPlotsByReplicateGroup = function(rle_output_dir, meta_qual_df){
  # TODO: SPECIFY NINETY MINUTE INDUCTION

  genotype_rle_full_list = Sys.glob(paste0(rle_output_dir, '/*full*'))
  genotype_rle_summary_list = Sys.glob(paste0(rle_output_dir, '/*summary*'))

  list_names = str_extract(genotype_rle_full_list, 'CNAG_[[:digit:]]+.*')
  list_names = str_remove(list_names, '_full.csv')

  names(genotype_rle_full_list) = list_names
  names(genotype_rle_summary_list) = list_names

  genotype_rle_full_df_list = lapply(genotype_rle_full_list, read_csv)

  plot_output_dir = paste(rle_output_dir, 'plots', sep='/')
  dir.create(paste(rle_output_dir, 'plots', sep='/'))
  # skip printing wildtype
  for (i in seq(1, length(genotype_rle_full_df_list))){
    # fill_column = "AUTO_AUDIT"
    # fill_column_colors = c("0" = "#52854C", "1" = "#D16103")

    fill_column = "LIBRARYDATE"
    num_colors <- length(unique(meta_qual_df$LIBRARYDATE))
    mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(num_colors)
    names(mycolors) = as.character(unique(meta_qual_df$LIBRARYDATE))
    fill_column_colors = mycolors

    rle = rleBoxplots(genotype_rle_full_df_list[[i]], meta_qual_df, fill_column, fill_column_colors, "rle_boxplots", list_names[[i]])

    pdf(paste(plot_output_dir, paste0(list_names[[i]], '.pdf'), sep='/'))
    plot(rle)
    dev.off()
  }
} # end createRLEPlotsByReplicateGroup()

#' @export
createRLEPlotsByReplicateGroup_envPert = function(rle_output_dir, meta_qual_df){


  rle_full_list = Sys.glob(paste0(rle_output_dir, '/*full*'))
  rle_summary_list = Sys.glob(paste0(rle_output_dir, '/*summary*'))

  list_names = basename(rle_full_list)
  list_names = str_remove(list_names, '_full.csv')

  names(rle_full_list) = list_names
  names(rle_summary_list) = list_names

  rle_full_list = lapply(rle_full_list, read_csv)

  plot_output_dir = paste(rle_output_dir, 'plots', sep='/')
  dir.create(paste(rle_output_dir, 'plots', sep='/'))
  # skip printing wildtype
  for (i in seq(1, length(rle_full_list))){
    # fill_column = "AUTO_AUDIT"
    # fill_column_colors = c("0" = "#52854C", "1" = "#D16103")

    fill_column = "LIBRARYDATE"
    num_colors <- length(unique(meta_qual_df$LIBRARYDATE))
    mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(num_colors)
    names(mycolors) = as.character(unique(meta_qual_df$LIBRARYDATE))
    fill_column_colors = mycolors

    rle = rleBoxplots(rle_full_list[[i]], meta_qual_df, fill_column, fill_column_colors, "rle_boxplots", list_names[[i]])

    pdf(paste(plot_output_dir, paste0(list_names[[i]], '.pdf'), sep='/'))
    plot(rle)
    dev.off()
  }
} # end createRLEPlotsByReplicateGroup()
