#' Get Combined Metadata from the database
#'
#' GET tables specified in table_vector, join on keys, return combined tables as a dataframe
#' @usage getMetadata(api_url)
#' @param api_url NOTE: api_url is a variable saved into the project environment, in addition to the parameter. You can use the usage statement directly. The argument is provided for development in the event that you want to test a local instance of the database. An example url: "http://13.59.167.2/api" No trailing /
#' @note api_url is a environmental variable saved in the package which will point to the current url (as of now, the one listed above)
#' @param table_vector list of tables to get from api_url. DEFAULT: c("BioSample", "RnaSample", "S1cdnaSample", "S2cdnaSample", "Library", "FastqFiles", "QualityAssess")
#' @return a dataframe of the joined tables in the database
#'
#' @export
getMetadata = function(api_url, table_vector = c("BioSample", "RnaSample", "S1cdnaSample", "S2cdnaSample", "Library", "FastqFiles", "QualityAssess")){

  tablename_vector = table_vector

  table_vector = lapply(tablename_vector, function(x) joinTables(x,api_url))

  names(table_vector) = tablename_vector

  combined_df = table_vector[[1]]

  for(i in seq(2, length(table_vector))){
    combined_df = combined_df %>% left_join(table_vector[[i]])
  }

  # replace empty string with NA
  combined_df[combined_df == ""] = NA

  return(combined_df)
}

#' @param table_name name of a table in the database, eg one of c("BioSample", "RnaSample", "S1cdnaSample", "S2cdnaSample", "Library", "FastqFiles", "QualityAssess")
#' @param api_url NOTE: api_url is a variable saved into the project environment, in addition to the parameter. You can use the usage statement directly. The argument is provided for development in the event that you want to test a local instance of the database. An example url: "http://13.59.167.2/api" No trailing /
#' @export
joinTables = function(tablename, api_url){
  tmp_dir = tempdir()
  tmp <- paste(tmp_dir, paste0(tablename, ".json"), sep="/")

  table_http_addr = paste(api_url, tablename, sep="/")

  curl_download(url=table_http_addr, tmp)

  response = fromJSON(tmp)

  json_file <- lapply(response, function(x) {
    x[sapply(x, is.null)] <- NA
    unlist(x)
  })

  metadata = as_tibble(do.call("cbind", json_file))

  # delete the tmp file
  unlink(tmp)

  return (metadata)
}

#' Get combined raw counts
#'
#' GET raw counts
#' @usage getRawCounts(api_url)
#' @param api_url NOTE: api_url is a variable saved into the project environment, in addition to the parameter. You can use the usage statement directly. The argument is provided for development in the event that you want to test a local instance of the database. An example url: "http://13.59.167.2/api" No trailing /
#' @return a gene by samples dataframe of all counts
#'
#' @export
getRawCounts = function(api_url){
  tablename="Counts"

  tmp_dir = tempdir()

  next_counts_url = paste(api_url, tablename, sep="/")
  page = 1
  while (length(next_counts_url > 0)){
    tmp <- paste(tmp_dir, paste0(tablename, "_", as.character(page), ".json"), sep="/")

    curl_download(url=next_counts_url, tmp)
    response = fromJSON(tmp)
    page = page+1
    next_counts_url = response[['next']]
  }

  count_paths = Sys.glob(paste(tempdir(),"Counts_*", sep="/"))
  full_df = data.frame(remove = rep(0,6967))
  for (path in count_paths){
    response = fromJSON(path)
    response_data = response$results
    count_list = lapply(response_data$rawCounts, function(x) x[-which(sapply(x, is.null))][[1]][0:6967])
    count_df = as.data.frame(count_list, check.names=FALSE)
    full_df = bind_cols(full_df, count_df)
    unlink(path)
  }
  return(full_df %>% select(-remove))
}



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
