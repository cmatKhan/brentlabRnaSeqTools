#' Get Combined Metadata from the database
#'
#' GET tables specified in table_vector, join on keys, return combined tables as a dataframe
#'
#' @param api_url is the url to the api minus any table names, eg "http://13.59.167.2/api" No trailing /
#' @param table_vector list of tables to get from api_url. DEFAULT: c("BioSample", "RnaSample", "S1cdnaSample", "S2cdnaSample", "Library", "FastqFiles", "QualityAssess")
#' @return a dataframe of the joined tables in the database
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
