#' get raw counts from database
#'
#' this downloads into $TEMPDIR a number of files (the paginations of the counts table) and then concats the counts.
#'
#' @note it only takes the first 6967 genes (hard coded currently) which are the protein coding annotations in the current gff
#' @note column names are just the fastq file name (no .fastq.gz, no _read_count.tsv or any other suffix)
#'
#' @param api_url is the url to the api minus any table names, eg "http://13.59.167.2/api" No trailing /
#'
#' @return a dataframe of gene x sample
#' @export

getRawCounts = function(api_url, tablename="Counts"){
  #' colnames are fastq filenames without any suffix (eg, no _read_counts.tsv or .fastq.gz)
  #' ONLY PROTEIN CODING 0:6967

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
