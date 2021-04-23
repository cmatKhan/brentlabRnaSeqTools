#' pull entire database (not counts) and save to output_dir for archival purposes
#'
#' saves both the individual tables and the combined_df
#'
#' @param output_dir where to deposit a subdirectory, named by todays date in this format: 20210407, with the tables and combined_df inside. eg /lts/mblab/Crypto/rnaseq_data/crypto_database_archive
#'
#' @return None, writes a directory called <today's date> with tables and combined_df as .csv to output_dir
#'
#' @export
archiveDatabase = function(output_dir){
  tablename_vector = c("BioSample",
                       "RnaSample",
                       "S1cdnaSample",
                       "S2cdnaSample",
                       "Library",
                       "FastqFiles",
                       "QualityAssess")

  table_vector = lapply(tablename_vector, function(x) joinTables(x, api_url))

  names(table_vector) = tablename_vector

  output_path = paste(output_dir, format(Sys.Date(), "%Y%m%d"), sep="/")
  dir.create(output_path)

  invisible(lapply(tablename_vector, function(x) write_csv(table_vector[[x]], paste(output_path, paste0(x, ".csv"), sep='/'))))

  # write combined_df also
  combined_df = table_vector[[1]]
  for(i in seq(2, length(table_vector))){
    combined_df = combined_df %>% left_join(table_vector[[i]])
  }
  write_csv(combined_df, paste(output_path, 'combined_df.csv', sep="/"))
}
