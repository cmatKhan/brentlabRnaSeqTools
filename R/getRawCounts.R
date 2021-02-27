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
