#' calculate RLE of a numeric dataframe
#'
#' @param counts_df gene by samples dataframe of raw counts or logged counts (see paramter logged)
#' @param log2_transformed_flag Default FALSE set to true if log2 transformed counts are passed
#' @return rle dataframe with genes x samples. Values are the logged differences from the gene-wise medians
#'
#' @export
calculateRLE = function(counts_df, log2_transformed_flag = FALSE){

  if(!isNumeric(counts_df)){
    stop("counts_df must have all numeric columns")
  }

  counts_df = as_tibble(counts_df)

  # if log2_transformed_flag==TRUE, re-assign counts_df to log2_counts_df. else, add a pseudocount and log2
  ifelse(log2_transformed_flag, assign('log2_counts_df', counts_df),  assign('log2_counts_df', log2(counts_df + 1)))

  # calculate median expression for each gene across samples
  gene_wise_medians = apply(log2_counts_df, 1, median, na.rm=TRUE)

  # calculate deviations from the median
  rle_table_full = sweep(log2_counts_df, 1, gene_wise_medians, '-')

  return(rle_table_full)

} # end calculateRLE()

#' calculate medians across rows of dataframe
#'
#' @param count_df could be any numeric dataframe, but in this context it will typically be a count (raw or log2) df
#' @return a vector of row-wise medians (length == nrow of input df)
#'
#' @export
calculateGeneWiseMedians = function(count_df){

  # calculate median expression for each gene across samples
  all_gene_medians = apply(count_df, 1, median)

  return(all_gene_medians)

} # end calculateGeneWiseMedians

#' rleSummary calculates summary statistics of rleFullTable
#'
#' @import matrixStats
#'
#' @param rle_table_full the output of rleFullTable
#'
#' @return a dataframe sample by rle summary statistics
#'
#' @export
rleSummary = function(rle_table_full){

  # calculate median deviation by sample
  median_deviation_by_sample = apply(rle_table_full, 2, median, na.rm=TRUE)
  # calculate interquartile range
  interquartile_range = apply(rle_table_full, 2, iqr, na.rm=TRUE)
  # assemble table
  rle_table_summary = tibble(FASTQFILENAME = colnames(rle_table_full),
                             SAMPLE_DEVIATION_MEDIAN = median_deviation_by_sample,
                             ABS_SAMPLE_DEVIATION_MEDIAN = abs(median_deviation_by_sample),
                             INTERQUARTILE_RANGE = interquartile_range)

  return (rle_table_summary)

} # end rleSummary()

#'
#' calculate RLE by replicate groups
#'
#' @param replicates_vector a list of lists where each sublist represents a replicate group. Entries must be a metadata
#'                               parameter, such as fastqFileName, that corresponds to the columns of the counts.
#'                               Suggestion: use something like these dplyr functions to create the list of lists group_by() %>% group_split %>% pull(fastqFileName)
#' @param gene_quants a gene x sample dataframe with values as some sort of gene quantification (eg normed counts, or log2(norm_counts) with some effect removed), possibly already logged (@see already_logged_flag)
#' @param log2_transformed_flag a boolean where TRUE means the counts are already in log2 space
#'
#' @references rlePlotCompareEffectRemoved() to plot the norm counts and removedEffect 'counts' on the same plot
#'
#' @return a list of dataframes for each replicate group in replicateS_sample_list, each with dimensions gene x sample. values are RLE of the gene in a given sample
#'
#' @export
rleByReplicateGroup = function(replicates_vector, gene_quants, log2_transformed_flag){


  lapply(replicates_vector, function(x) calculateRLE(gene_quants[, x], log2_transformed_flag=log2_transformed_flag))

}

# TODO fix these plotting functions

#'
#' plot RLE for a given column filter (eg, metadata[metadata$MEDIUM == 'PBS']$FASTQFILENAME would give a list of fastqFileNames to filter)
#'
#' @param deseq_object a deseq object with results from the DESeq() call
#' @param model_matrix the deseq_object model matrix
#' @param column_filter a vector of fastqFileNames (or whatever the columns -- samples -- are called)
#' @param title of the plots
#' @return list with slots norm_count_rle and effect_removed_rle
#'
#' @export
rlePlot = function(deseq_object, model_matrix, column_filter, title){

  norm_counts = counts(deseq_object, normalize=TRUE)

  fltr_norm_counts = norm_counts[ , column_filter]

  effect_removed_counts = removeParameterEffects(deseq_object, model_matrix)

  fltr_effect_removed_counts = effect_removed_counts[ ,column_filter]

  norm_count_rle = rlePlot_helper(fltr_norm_counts, log2_transformed_flag=FALSE, paste(title, 'Norm Counts', sep=" - "))
  effect_removed_rle = rlePlot_helper(fltr_effect_removed_counts, log2_transformed_flag=TRUE, paste(title, 'Effect Removed', sep=' - '))

  return (list('norm_count_rle' = norm_count_rle, 'effect_removed_rle' = effect_removed_rle))


}

#' the actual plotting function for rlePlot
#'
#' @import dplyr
#' @import ggplot2
#'
#' @param count_df counts in gene x sample
#' @param log2_transformed_flag boolean where TRUE indicates the counts are in log2 space
#' @param title title of the output plot
#'
#' @return a ggplot
#'
rlePlot_helper = function(count_df, log2_transformed_flag, title){


  rle_full_table = calculateRLE(count_df, log2_transformed_flag=log2_transformed_flag)

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

#'
#' plots output of rleSummaryByReplicateGroup
#'
#' @import dplyr
#'
#' @param norm_counts_rle output of calculateRLE (maybe one of the sublists in rleByReplicateGroup())
#' @param removed_effect_rle see norm_counts_rle, but after removing some batch effects
#' @param metadata_df metadata with at least FASTQFILENAME and LIBRARYDATE
#' @param title title of the plot
#'
#' @return a ggplot with both the norm counts (more transparent) and removedEffect 'counts' on the same plot
#'
#' @export
rlePlotCompareEffectRemoved = function(norm_counts_rle, removed_effect_rle, metadata_df, title){

  norm_counts_rle %>%
    pivot_longer(colnames(.), names_to="FASTQFILENAME", values_to="RLE") %>%
    mutate(quant_type="norm_counts") %>%
    bind_rows(
      (removed_effect_rle %>%
         pivot_longer(colnames(.), names_to="FASTQFILENAME", values_to="RLE") %>%
         mutate(quant_type="removed_effect"))) %>%
    left_join(metadata_df) %>%
    mutate(LIBRARYDATE = as.factor(LIBRARYDATE)) %>%
    ggplot() +
    geom_boxplot(aes(FASTQFILENAME, RLE, fill=LIBRARYDATE, color=quant_type, alpha=quant_type), outlier.shape = NA) +
    scale_color_manual(values=c("#999999", "#000000")) +
    scale_alpha_manual(values=c(.1, 1)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ylim(-3,3)+
    ggtitle(title)

}

#'
#' composite plot of all rle_stats in one plot
#' @note no title on graphs. use the names of the returned object to title in presentation
#'
#' @param rle_df a samples x rle stats df with at minimum columns SAMPLE_DEVIATION_MEDIAN, ABS_SAMPLE_DEVIATION_MEDIAN, INTERQUARTILE_RANGE
#'
#' @return a plot with three horizontal panels for each of the rle stats
#'
# plotRLEhistograms = function(rle_df){
#   rle_df %>%
#     pivot_longer(-FASTQFILENAME, names_to="rle_stat", values_to="value") %>%
#     mutate(rle_stat = factor(rle_stat, levels=c("SAMPLE_DEVIATION_MEDIAN", "ABS_SAMPLE_DEVIATION_MEDIAN", "INTERQUARTILE_RANGE"))) %>%
#     ggplot() +
#     geom_histogram(aes(value))+
#     theme(axis.title.x=element_blank())+
#     facet_wrap('rle_stat', scales="free_x", dir='v')
# }

