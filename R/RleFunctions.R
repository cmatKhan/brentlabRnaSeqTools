#' calculate RLE
#'
#' You may pass raw counts or log2 transformed counts with the logged = TRUE argument
#'
#' @param counts_df gene by samples dataframe of raw counts or logged counts (see paramter logged)
#' @param gene_id_column (Optional) if gene_id_column is not passed, a sequence from 1 to nrow(counts_df) will be used
#' @param logged (Default FALSE) set to true if log2 transformed counts are passed
#'
#' @return rle dataframe with genes x rle statistics
#'
#' @export
createFullRLETable = function(counts_df, gene_id_column = NULL, logged = FALSE){

  #if()

  # TODO: REMOVE THE GENE_ID ARGUMENT AND JUST REMOVE THE GENE_ID COLUMN IF IT EXISTS
  counts_df = as_tibble(counts_df)
  counts_df$gene_id = gene_id_column

  counts_df = counts_df %>% select(-gene_id)

  if(logged == FALSE){
    # log2 the table (add pseudo count of 1)
    log2_counts_df = log2(counts_df + 1)
  } else{
    log2_counts_df = counts_df
  }

  # calculate median expression for each gene across samples
  all_gene_medians = apply(log2_counts_df, 1, median)

  # calculate deviations from the median
  rle_table_full = sweep(log2_counts_df, 1, all_gene_medians, '-')

  return(rle_table_full)

} # end createFullRLETable()




#' rleSummary calculates summary statistics of rleFullTable
#'
#' @param rle_table_full the output of rleFullTable
#'
#' @return a dataframe sample by rle summary statistics
#'
#' @export
rleSummary = function(rle_table_full){
  # calculate median deviation by sample
  median_deviation_by_sample = apply(rle_table_full, 2, median, na.rm=TRUE)

  # calculate twenty fifth percentile of deviations from median
  q1 = apply(rle_table_full, 2, quantile, .25, na.rm=TRUE)
  # ditto seventy fifth
  q3 = apply(rle_table_full, 2, quantile, .75, na.rm=TRUE)
  # calculate interquartile range
  iqr = q3 - q1
  # min/max whisker threshold = q1/3 +/- 1.5*iqr
  min_whisker_threshold = q1 - 1.5*iqr

  min_whisker_value = c()
  for (i in seq(1,length(colnames(rle_table_full)))){
    x = min(rle_table_full[,i][rle_table_full[,i]> min_whisker_threshold[i]])
    min_whisker_value = c(min_whisker_value,x)
  }

  max_whisker_threshold = q3 + 1.5*iqr
  max_whisker_value = c()
  for (i in seq(1,length(colnames(rle_table_full)))){
    x = max(rle_table_full[,i][rle_table_full[,i]< max_whisker_threshold[i]])
    max_whisker_value = c(max_whisker_value,x)
  }
  #
  # inter_whisker_range = max_whisker_value - min_whisker_value

  # num_outliers_below =c()
  # for (i in seq(1,length(colnames(rle_table_full)))){
  #   x = length(rle_table_full[,i][rle_table_full[,i] < min_whisker_threshold[i]])
  #   num_outliers_below = c(num_outliers_below,x)
  # }

  # most_extreme_below = apply(rle_table_full,2,min)
  # for (i in seq(1,length(most_extreme_below))){
  #   if (most_extreme_below[i] >= min_whisker_threshold[i]){
  #     most_extreme_below[i] = NA
  #   }
  # }

  # num_outliers_above =c()
  # for (i in seq(1,length(colnames(rle_table_full)))){
  #   x = length(rle_table_full[,i][rle_table_full[,i] > max_whisker_threshold[i]])
  #   num_outliers_above = c(num_outliers_above,x)
  # }

  # most_extreme_above = apply(rle_table_full,2,max)
  # for (i in seq(1,length(most_extreme_above))){
  #   if (most_extreme_above[i] <= max_whisker_threshold[i]){
  #     most_extreme_above[i] = NA
  #   }
  # }

  rle_table_summary = tibble(FASTQFILENAME = colnames(rle_table_full),
                             SAMPLE_DEVIATION_MEDIAN = median_deviation_by_sample,
                             ABS_SAMPLE_DEVIATION_MEDIAN = abs(median_deviation_by_sample),
                             TWENTY_FIVE_QUANTILE = q1,
                             SEVENTY_FIVE_QUANTILE = q3,
                             INTERQUARTILE_RANGE = iqr,
                             MIN_WHISKER = min_whisker_value,
                             MAX_WHISKER = max_whisker_value)

  return (rle_table_summary)

} # end rleSummary()

#' These are messy and need to be cleaned up, but are how i calculate rle by replicate group for the 90minInduction and env_pert NOTE: only 90minInduction is currently tested with new database
#'
#' @param meta_qual_df is a dataframe pulled from the database with biosample --> quality assess (name is remnant of old system)
#' @param norm_counts from the deseq object
#' @param output_dirpath where to put the rle results. i suggest making a directory like "rle" in your working directory
#' @param protocol_selector whether or not to include library protocol in the replicate grouping
#' @param already_logged_flag are counts logged? (note: deseq norm counts are not logged unless you log them yourself)
#'
#' @return None. writes to file
#'
#' @export
extractRLEByReplicateGroup_90minInduction = function(meta_qual_df, norm_counts, output_dirpath, protocol_selector, already_logged_flag){

  dir.create(output_dirpath)

  meta_qual_df$LIBRARYDATE = as.factor(as.Date(meta_qual_df$LIBRARYDATE))
  meta_qual_df$MEDIUM = as.factor(meta_qual_df$MEDIUM)
  meta_qual_df$TEMPERATURE = as.factor(meta_qual_df$TEMPERATURE)
  meta_qual_df$ATMOSPHERE = as.factor(meta_qual_df$ATMOSPHERE)
  meta_qual_df$TREATMENT = as.factor(meta_qual_df$TREATMENT)
  meta_qual_df$OTHERCONDITIONS = as.factor(meta_qual_df$OTHERCONDITIONS)
  meta_qual_df$TIMEPOINT = as.factor(meta_qual_df$TIMEPOINT)
  meta_qual_df$LIBRARYPROTOCOL = as.factor(meta_qual_df$LIBRARYPROTOCOL)
  meta_qual_df$GENOTYPE = as.factor(meta_qual_df$GENOTYPE)

  for (genotype in unique(meta_qual_df$GENOTYPE)){
    # filter out known strain/geno problems
    genotype_filtered_df = meta_qual_df %>% filter(GENOTYPE == genotype)

    other_cond_str = 'OTHERCONDITIONS'

    if(protocol_selector){
      replicate_split_meta_qual_list = split(genotype_filtered_df, f = list(genotype_filtered_df$MEDIUM, genotype_filtered_df$TEMPERATURE, genotype_filtered_df$ATMOSPHERE, genotype_filtered_df$TREATMENT, genotype_filtered_df$OTHERCONDITIONS, genotype_filtered_df$TIMEPOINT, genotype_filtered_df$LIBRARYPROTOCOL))
    } else {
      replicate_split_meta_qual_list = split(genotype_filtered_df, f = list(genotype_filtered_df$MEDIUM, genotype_filtered_df$TEMPERATURE, genotype_filtered_df$ATMOSPHERE, genotype_filtered_df$TREATMENT, genotype_filtered_df$OTHERCONDITIONS, genotype_filtered_df$TIMEPOINT))
    }

    gene_id_col = seq(1,nrow(norm_counts))
    for (split_df in replicate_split_meta_qual_list){
      # only take genotypes with replicates > 3
      if (nrow(split_df) > 2){
        medium = as.character(unique(split_df$MEDIUM))
        temperature = as.character(unique(split_df$TEMPERATURE))
        atmosphere = as.character(unique(split_df$ATMOSPHERE))
        treatment = as.character(unique(split_df$TREATMENT))
        other_conditions = as.character(unique(split_df$OTHERCONDITIONS))
        timepoint = as.character(unique(split_df$TIMEPOINT))
        protocol = as.character(unique(split_df$LIBRARYPROTOCOL))

        ftlr_norm_counts = norm_counts[,split_df$FASTQFILENAME]
        if (protocol_selector){
          full_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint,protocol, "full.csv", sep='_'), sep='/')
          summary_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint,protocol, "summary.csv", sep='_'), sep='/')
        } else{
          full_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint, "full.csv", sep='_'), sep='/')
          summary_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint, "summary.csv", sep='_'), sep='/')
        }

        if (already_logged_flag == TRUE){
          fltr_rle_full = createFullRLETable(ftlr_norm_counts, gene_id_column=NULL, logged=TRUE)
          write_csv(as_tibble(fltr_rle_full), full_filename)

          fltr_rle_summary = rleSummary(fltr_rle_full)
          write_csv(as_tibble(fltr_rle_summary), summary_filename)

        } else{
          fltr_rle_full = createFullRLETable(ftlr_norm_counts, gene_id_column=NULL)
          write_csv(as_tibble(fltr_rle_full), full_filename)

          fltr_rle_summary = rleSummary(fltr_rle_full)
          write_csv(as_tibble(fltr_rle_summary), summary_filename)
        }
      }
    }
  }
}# end extractRLEByReplicateGroup()

#' @export
extractRLEByReplicateGroup_EnvPert = function(df, norm_counts, output_dirpath, protocol_selector, already_logged_flag){

  dir.create(output_dirpath)

  x = split(df, f = list(df$MEDIUM, df$TEMPERATURE, df$ATMOSPHERE, df$TREATMENT, df$OTHERCONDITIONS, df$TIMEPOINT, df$TREATMENT, df$TREATMENTCONC))

  gene_id_col = seq(1,nrow(norm_counts))
  for (split_df in x){
    if (nrow(split_df) != 0){
      # only take genotypes with replicates > 3
      if (nrow(split_df) > 2){
        medium = as.character(unique(split_df$MEDIUM))
        temperature = as.character(unique(split_df$TEMPERATURE))
        atmosphere = as.character(unique(split_df$ATMOSPHERE))
        treatment = as.character(unique(split_df$TREATMENT))
        other_conditions = as.character(unique(split_df$OTHERCONDITIONS))
        timepoint = as.character(unique(split_df$TIMEPOINT))
        protocol = as.character(unique(split_df$LIBRARYPROTOCOL))
        treatment = as.character(unique(split_df$TREATMENT))
        treatmentconc = as.character(unique(split_df$TREATMENTCONC))


        full_filename = paste(output_dirpath, paste(medium,temperature,atmosphere,treatment,other_conditions,timepoint, treatment, treatmentconc, "full.csv", sep='_'), sep='/')
        summary_filename = paste(output_dirpath, paste(medium,temperature,atmosphere,treatment,other_conditions,timepoint, treatment, treatmentconc, "summary.csv", sep='_'), sep='/')

        ftlr_norm_counts = norm_counts[,split_df$FASTQFILENAME]
        if (already_logged_flag){
          fltr_rle_full = createFullRLETable(ftlr_norm_counts, gene_id_column=NULL, logged=TRUE)
          write_csv(as_tibble(fltr_rle_full), full_filename)

          fltr_rle_summary = rleSummary(fltr_rle_full)
          write_csv(as_tibble(fltr_rle_summary), summary_filename)

        } else{
          fltr_rle_full = createFullRLETable(ftlr_norm_counts, seq(1,nrow(norm_counts)))
          write_csv(as_tibble(fltr_rle_full), full_filename)

          fltr_rle_summary = rleSummary(fltr_rle_full)
          write_csv(as_tibble(fltr_rle_summary), summary_filename)
        }
      }
    }
  }
}# end extractRLEByReplicateGroup()

#' @export
extractRLEByReplicateGroup_EnvPert_titrations = function(df, norm_counts, output_dirpath, protocol_selector, already_logged_flag){

  dir.create(output_dirpath)

  x = split(df, f = list(df$MEDIUM, df$TEMPERATURE, df$ATMOSPHERE, df$TREATMENT, df$TREATMENTCONC, df$OTHER_CONDITIONS, df$TIMEPOINT))

  gene_id_col = seq(1,nrow(norm_counts))
  for (split_df in x){
    if (nrow(split_df) != 0){
      # only take genotypes with replicates > 3
      if (nrow(split_df) > 2){
        medium = as.character(unique(split_df$MEDIUM))
        temperature = as.character(unique(split_df$TEMPERATURE))
        atmosphere = as.character(unique(split_df$ATMOSPHERE))
        treatment = as.character(unique(split_df$TREATMENT))
        other_conditions = as.character(unique(split_df$OTHERCONDITIONS))
        timepoint = as.character(unique(split_df$TIMEPOINT))
        protocol = as.character(unique(split_df$LIBRARYPROTOCOL))
        treatment = as.character(unique(split_df$TREATMENT))
        treatmentconc = as.character(unique(split_df$TREATMENTCONC))


        full_filename = paste(output_dirpath, paste(medium,temperature,atmosphere,treatment,other_conditions,timepoint, treatment, treatmentconc, "full.csv", sep='_'), sep='/')
        summary_filename = paste(output_dirpath, paste(medium,temperature,atmosphere,treatment,other_conditions,timepoint, treatment, treatmentconc, "summary.csv", sep='_'), sep='/')

        ftlr_norm_counts = norm_counts[,split_df$FASTQFILENAME]
        if (already_logged_flag){
          fltr_rle_full = createFullRLETable(ftlr_norm_counts, gene_id_column=NULL, logged=TRUE)
          write_csv(as_tibble(fltr_rle_full), full_filename)

          fltr_rle_summary = rleSummary(fltr_rle_full)
          write_csv(as_tibble(fltr_rle_summary), summary_filename)

        } else{
          fltr_rle_full = createFullRLETable(ftlr_norm_counts, seq(1,nrow(norm_counts)))
          write_csv(as_tibble(fltr_rle_full), full_filename)

          fltr_rle_summary = rleSummary(fltr_rle_full)
          write_csv(as_tibble(fltr_rle_summary), summary_filename)
        }
      }
    }
  }
}# end extractRLEByReplicateGroup()

#' @export
extractRLEByReplicateGroup_EnvPert = function(df, norm_counts, output_dirpath, protocol_selector, already_logged_flag){

  dir.create(output_dirpath)

  x = split(df, f = list(df$MEDIUM, df$TEMPERATURE, df$ATMOSPHERE, df$TREATMENT, df$OTHER_CONDITIONS, df$TIMEPOINT))

  gene_id_col = seq(1,nrow(norm_counts))
  for (split_df in x){
    if (nrow(split_df) != 0){
      # only take genotypes with replicates > 3
      if (nrow(split_df) > 2){
        medium = as.character(unique(split_df$MEDIUM))
        temperature = as.character(unique(split_df$TEMPERATURE))
        atmosphere = as.character(unique(split_df$ATMOSPHERE))
        treatment = as.character(unique(split_df$TREATMENT))
        other_conditions = as.character(unique(split_df$OTHERCONDITIONS))
        timepoint = as.character(unique(split_df$TIMEPOINT))
        protocol = as.character(unique(split_df$LIBRARYPROTOCOL))


        full_filename = paste(output_dirpath, paste(medium,temperature,atmosphere,treatment,other_conditions,timepoint, "full.csv", sep='_'), sep='/')
        summary_filename = paste(output_dirpath, paste(medium,temperature,atmosphere,treatment,other_conditions,timepoint, "summary.csv", sep='_'), sep='/')

        ftlr_norm_counts = norm_counts[,split_df$FASTQFILENAME]
        if (already_logged_flag){
          fltr_rle_full = createFullRLETable(ftlr_norm_counts, gene_id_column=NULL, logged=TRUE)
          write_csv(as_tibble(fltr_rle_full), full_filename)

          fltr_rle_summary = rleSummary(fltr_rle_full)
          write_csv(as_tibble(fltr_rle_summary), summary_filename)

        } else{
          fltr_rle_full = createFullRLETable(ftlr_norm_counts, seq(1,nrow(norm_counts)))
          write_csv(as_tibble(fltr_rle_full), full_filename)

          fltr_rle_summary = rleSummary(fltr_rle_full)
          write_csv(as_tibble(fltr_rle_summary), summary_filename)
        }
      }
    }
  }
}# end extractRLEByReplicateGroup()



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

#'
#' plot RLE for a given column filter (eg, metadata[metadata$MEDIUM == 'PBS']$FASTQFILENAME would give a list of fastqFileNames to filter)
#' @param deseq_object a deseq object with results from the DESeq() call
#' @param model_matrix the deseq_object model matrix
#' @param column_filter a vector of fastqFileNames (or whatever the columns -- samples -- are called)
#' @param title of the plots
#' @return list with slots norm_count_rle and effect_removed_rle
#'
#' @export
new_rlePlotFunc = function(deseq_object, model_matrix, column_filter, title){

  norm_counts = counts(deseq_object, normalize=TRUE)

  fltr_norm_counts = norm_counts[ , column_filter]

  effect_removed_counts = removeParameterEffects(deseq_object, model_matrix)

  fltr_effect_removed_counts = effect_removed_counts[ ,column_filter]

  norm_count_rle = new_rlePlot_helper(fltr_norm_counts, paste(title, 'Norm Counts', sep=" - "))
  effect_removed_rle = new_rlePlot_helper(fltr_effect_removed_counts, paste(title, 'Effect Removed', sep=' - '))

  return (list('norm_count_rle' = norm_count_rle, 'effect_removed_rle' = effect_removed_rle))


}

new_rlePlot_helper = function(count_df, title){
  rle_full_table = createFullRLETable(count_df)

  gene_id = read_tsv("~/Desktop/rnaseq_pipeline/rnaseq_pipeline/genome_files/KN99/KN99_gene_id_list.txt", col_names = 'gene_id')[1:6967,]
  rle_full_table$gene_id = gene_id

  rle_full_table %>%
    pivot_longer(!gene_id, names_to = 'FASTQFILENAME', values_to = "RLE") %>%
    ggplot()+
    geom_boxplot(aes(FASTQFILENAME, RLE), outlier.shape=NA)+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    xlab('Samples')+
    ylab('Deviation from Median')+
    coord_cartesian(ylim = c(-3,3))+
    scale_y_continuous(n.breaks = 20)+
    geom_hline(yintercept = 0)+
    ggtitle(title)
}


