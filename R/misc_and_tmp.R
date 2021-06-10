#' progressively remove max IQR sample and recalculate
#'
#' @param sample_set the metadata
#' @param logged_norm_counts counts on log2 scale and normalized
#'
#'@export
removeOneRedoIqr = function(sample_set, logged_norm_counts){

  removed_sample_list = list()
  recalculated_iqr_list = list()
  i = nrow(sample_set)

  while (i > 3){
    max_iqr = max(sample_set$INTERQUARTILE_RANGE)
    max_iqr_boolean_vector = sample_set$INTERQUARTILE_RANGE == max_iqr
    sample_number = as.character(sample_set[max_iqr_boolean_vector, "FASTQFILENUMBER"])

    removed_sample_list[[sample_number]] = max_iqr

    sample_set = sample_set[!max_iqr_boolean_vector, ]

    sample_columns = sample_set$FASTQFILENAME
    max_iqr_removed_log_norm_counts = logged_norm_counts[,sample_columns]

    rle_stats_df = rleSummary(createFullRLETable(max_iqr_removed_log_norm_counts, logged=TRUE))
    recalculated_iqr_list[[sample_number]] = rle_stats_df

    i = i-1
  }

    return (list("removed_samples" = removed_sample_list, "df_list" = recalculated_iqr_list))
  ## a recursive version -- each step works, and if you put a print statement after the df calculation then the output is correct, but the lists aren't being returned
  # if(nrow(sample_set) > 4){
  #   max_iqr = max(sample_set$INTERQUARTILE_RANGE)
  #   max_iqr_boolean_vector = sample_set$INTERQUARTILE_RANGE == max_iqr
  #   sample_number = as.character(sample_set[max_iqr_boolean_vector, "FASTQFILENUMBER"])
  #
  #   x[[sample_number]] = max_iqr
  #
  #   sample_set = sample_set[!max_iqr_boolean_vector, ]
  #
  #   sample_columns = sample_set$FASTQFILENAME
  #   max_iqr_removed_log_norm_counts = logged_norm_counts[,sample_columns]
  #
  #   rle_stats_df = rleSummary(createFullRLETable(max_iqr_removed_log_norm_counts, logged=TRUE))
  #   y[[sample_number]] = rle_stats_df
  #
  #   removeOneRedoIqr(sample_set, logged_norm_counts, x, y)
  # } # end if
}
