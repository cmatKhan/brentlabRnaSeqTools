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
