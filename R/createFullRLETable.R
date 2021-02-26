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

  if(!is.null(gene_id_column)){
    stopifnot(length(gene_id_column) == nrow(counts_df))
  } else {
    gene_id_column = seq(1:nrow(counts_df))
  }

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
