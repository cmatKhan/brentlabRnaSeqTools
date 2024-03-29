#' Get the combined metadata as a tibble from a remote database
#'
#' @description Join the biosample, rnasample, s1sample, s2sample, library, fastqFiles and qualityAssessment tables (in that order, left joins) and return the result as a tibble
#'
#' @importFrom RPostgres dbDisconnect
#' @importFrom dplyr left_join tbl
#'
#' @description Use the RPostgres package to connect to a remote postgresql database, do the table joining, and return the joined metadata as a tibble. The database connection is closed
#' @param database_host if connecting to a database hosted on AWS, it might be something like ec2-54-83-201-96.compute-1.amazonaws.com
#' @param database_name name of the database, eg for cryptococcus kn99, the database might be named kn99_database. Check with the documentation, whoever set up the database, or get into the server and check directly
#' @param database_user a user of the actual database, with some level of permissions. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param database_password password to the database user. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @note for information on using R environmental files, see \url{https://support.rstudio.com/hc/en-us/articles/360047157094-Managing-R-with-Rprofile-Renviron-Rprofile-site-Renviron-site-rsession-conf-and-repos-conf}
#' @source \url{https://rpostgres.r-dbi.org/}
#' @return A DBI connection to the remote database
#'
#' @export
getMetadata = function(database_host, database_name, database_user, database_password){

  db = connectToDatabase(database_host, database_name, database_user, database_password)

  biosample = tbl(db, 'bioSample')
  rnasample = tbl(db, 'rnaSample')
  s1sample = tbl(db, 's1cDNASample')
  s2sample = tbl(db, 's2cDNASample')
  library = tbl(db, 'library')
  fastqFiles = tbl(db, 'fastqFiles')
  quality = tbl(db, 'qualityAssessment')

  joined_meta_tables = biosample %>%
    left_join(rnasample, by = c('bioSampleNumber' = 'bioSampleNumber_id'))%>%
    left_join(s1sample, by = c('rnaSampleNumber' = 'rnaSampleNumber_id'))%>%
    left_join(s2sample, by = c('s1cDNASampleNumber' = 's1cDNASampleNumber_id'))%>%
    left_join(library, by = c('s2cDNASampleNumber' = 's2cDNASampleNumber_id'))%>%
    left_join(fastqFiles, by = c('librarySampleNumber' = 'librarySampleNumber_id'))%>%
    left_join(quality, by = c('fastqFileNumber' = 'fastqFileNumber_id'))

  metadata_df = as_tibble(joined_meta_tables)

  dbDisconnect(db)

  # cast integer64 to integers. see Utils is_integer64
  metadata_df <- metadata_df %>%
    mutate_if(is_integer64, as.integer)

  return(metadata_df)
}

#' Get combined raw counts
#'
#' @importFrom RPostgres dbGetQuery dbDisconnect
#' @importFrom dplyr bind_cols
#' @importFrom jsonlite fromJSON
#'
#' @param database_host if connecting to a database hosted on AWS,
#'                      it might be something like ec2-54-83-201-96.compute-1.amazonaws.com
#' @param database_name name of the database, eg for cryptococcus kn99, the database might be named kn99_database.
#'                      Check with the documentation, whoever set up the database, or get into the server and check
#'                      directly
#' @param database_user a user of the actual database, with some level of permissions. You'll need to check with the
#'                      database maintainer for this. It is suggested that you use a .Renviron file in your
#'                      local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param database_password password to the database user. You'll need to check with the database maintainer for this.
#'                          It is suggested that you use a .Renviron file in your local project
#'                          (make sure it is completely ignored by git, R, etc) to store this info
#'
#' @return a gene by samples dataframe of all counts
#'
#' @export
getRawCounts = function(database_host, database_name, database_user, database_password){
  db = connectToDatabase(database_host, database_name, database_user, database_password)

  counts = dbGetQuery(db, 'select "rawCounts" from counts')

  counts_df = bind_cols(lapply(seq(1, length(counts$rawCounts)), function(x) (as.data.frame(fromJSON(counts$rawCounts[x]), check.names=FALSE))))

  dbDisconnect(db)

  return (counts_df)
}


#' pull entire database (not counts) and save to output_dir for archival purposes
#'
#' @description saves both the individual tables, including counts, and the combined_df
#'
#' @importFrom RPostgres dbDisconnect
#' @importFrom readr write_csv
#' @importFrom dplyr tbl
#'
#' @param database_host if connecting to a database hosted on AWS, it might be something like ec2-54-83-201-96.compute-1.amazonaws.com
#' @param database_name name of the database, eg for cryptococcus kn99, the database might be named kn99_database. Check with the documentation, whoever set up the database, or get into the server and check directly
#' @param database_user a user of the actual database, with some level of permissions. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param database_password password to the database user. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param output_dir where to deposit a subdirectory, named by todays date in this format: 20210407, with the tables and combined_df inside. eg a mounted local directory /mnt/htcf_lts/crypto_database_archive/ --> /lts/mblab/Crypto/rnaseq_data/crypto_database_archive
#' @param archive_counts_flag boolean indicating whether or not to save the counts. default is TRUE
#' @return None, writes a directory called <today's date> with tables and combined_df as .csv to output_dir
#'
#' @export
archiveDatabase = function(database_host, database_name, database_user, database_password, output_dir, archive_counts_flag = TRUE){

  today_date = format(Sys.Date(), "%Y%m%d")
  current_output_path = file.path(output_dir, today_date)
  dir.create(current_output_path)

  db = connectToDatabase(database_host, database_name, database_user, database_password)

  tbl_list = list()

  tbl_list[['biosample']] = tbl(db, 'bioSample')
  tbl_list[['rnasample']] = tbl(db, 'rnaSample')
  tbl_list[['s1sample']] = tbl(db, 's1cDNASample')
  tbl_list[['s2sample']] = tbl(db, 's2cDNASample')
  tbl_list[['library']] = tbl(db, 'library')
  tbl_list[['fastqFiles']] = tbl(db, 'fastqFiles')
  tbl_list[['qualityAssessment']] = tbl(db, 'qualityAssessment')

  tbl_list_names = names(tbl_list)

  tbl_list = lapply(tbl_list, as_tibble)
  names(tbl_list) = tbl_list_names # maybe not necessary

  lapply(names(tbl_list), function(x) write_csv(tbl_list[[x]], file.path(current_output_path, paste0(x, ".csv") )))

  combined_df = getMetadata(database_host, database_name, database_user, database_password)
  write_csv(combined_df, file.path(current_output_path, paste0("combined_df_", today_date,'.csv')))

  if(archive_counts_flag){
    counts_df = getRawCounts(database_host, database_name, database_user, database_password)
    write_csv(counts_df, file.path(current_output_path, "counts.csv"))
  }

  dbDisconnect(db)

}

#'
#' Connect to a remote postgresql database
#'
#' @importFrom RPostgres Postgres dbConnect
#'
#' @description Use the RPostgres package to connect to a remote postgresql database
#' @param database_host if connecting to a database hosted on AWS, it might be something like ec2-54-83-201-96.compute-1.amazonaws.com
#' @param database_name name of the database, eg for cryptococcus kn99, the database might be named kn99_database. Check with the documentation, whoever set up the database, or get into the server and check directly
#' @param database_user a user of the actual database, with some level of permissions. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param database_password password to the database user. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @note for information on using R environmental files, see \url{https://support.rstudio.com/hc/en-us/articles/360047157094-Managing-R-with-Rprofile-Renviron-Rprofile-site-Renviron-site-rsession-conf-and-repos-conf}
#' @source \url{https://rpostgres.r-dbi.org/}
#' @return A DBI connection to the remote database
#'
#' @export
connectToDatabase = function(database_host, database_name, database_user, database_password){

  dbConnect(RPostgres::Postgres(),dbname = database_name,
            host = database_host, # i.e. 'ec2-54-83-201-96.compute-1.amazonaws.com'
            port = 5432, # or any other port specified by your DBA
            user = database_user,
            password = database_password)

}

#'
#' list tables in databse
#'
#' @importFrom RPostgres dbGetQuery
#'
#' @param db a connection to the database
#'
#' @seealso \url{https://www.postgresqltutorial.com/postgresql-show-tables/}
#' @return all tables in database
#'
#' @export
listTables = function(db){
  dbGetQuery(db, "SELECT *
                  FROM pg_catalog.pg_tables
                  WHERE schemaname != 'pg_catalog' AND
                  schemaname != 'information_schema';")
}

#'
#' get (via a http POST request) your user authentication token from the database
#'
#' @importFrom httr http_status content POST
#'
#' @param url check the database_info variable. for configured organisms, you can find this under database_info$organism$token_auth
#' @param username a valid username for the database. If you don't have one, then you'll need to ask for one to be created
#' @param password password associated with your username
#'
#' @note do not save your auth token in a public repository. For example, you might put it in your .Renviron and then make sure
#'       that your .Renviron is in your .gitignore. Otherwise, save it outside of a github tracked directory or otherwise ensure
#'       that it will not be pushed up to github
#'
#' @return the auth token associated with the username and password
#'
#' @export
getUserAuthToken = function(url, username, password){

  # see package httr for help
  token_response = POST(url=url,
                        body=list(username=username,
                                  password=password),
                        encode='json')

  if( http_status(token_response)$category == "Success" ){
    message("You might want to put your token in your .Renviron. If you do, please make sure the .Renviron file is in your .gitignore")
    httr::content(token_response)$token
  } else{
    message("There was a problem getting your token:")
    message(http_status(token_response)$message)
  }
}

#'
#' post new fastq sheet to database
#'
#' @importFrom httr content_type add_headers POST
#' @importFrom jsonlite toJSON
#' @importFrom dplyr select mutate
#'
#' @param database_fastq_url eg database_info$kn99_urls$FastqFiles. See see \code{\link{database_info}}
#' @param auth_token see \code{\link{getUserAuthToken}}
#' @param new_fastq_path path to new fastq sheet
#'
#' @export
postFastqSheet = function(database_fastq_url, auth_token, new_fastq_path){

  # see utils
  fastq_df = brentlabRnaSeqTools::readInData(new_fastq_path)

  fastq_colnames = colnames(fastq_df)

  augment_fastq_df = fastq_df %>%
    {if("libraryDate" %in% fastq_colnames)
      select(.,-libraryDate)
      else .} %>%
    {if("libraryPreparer" %in% fastq_colnames)
      select(.,-libraryPreparer)
      else .} %>%
    {if(!"fastqObservations" %in% fastq_colnames)
      mutate(.,fastqObservations = "")
      else .} %>%
    {if(!"laneNumber" %in% fastq_colnames)
      mutate(.,laneNumber = "")
      else .} %>%
    mutate(volumePooled = round(as.numeric(volumePooled), 15)) %>%
    mutate(tapestationConc = round(as.numeric(tapestationConc), 4))

  post_body = jsonlite::toJSON(augment_fastq_df, auto_unbox = TRUE)

  POST(url = database_fastq_url,
       add_headers(Authorization = paste("token" , auth_token, sep=" ")),
       content_type("application/json"),
       body = post_body,
       encode = 'json')
}


# TODO add dry run option

#'
#' post counts to database
#'
#' @importFrom httr content_type add_headers POST
#' @importFrom jsonlite toJSON
#' @importFrom dplyr filter pull
#' @importFrom stringr str_remove
#'
#' @description using the package httr, post the raw count .csv, which is the compiled counts for a given run, to the database
#'
#' @param database_counts_url eg database_info$kn99_urls$Counts. see \code{\link{database_info}}
#' @param run_number the run number of this counts sheet -- this is important b/c fastqFileNames aren't necessarily unique outside of their runs
#' @param auth_token see \code{\link{getUserAuthToken}}
#' @param new_counts_path path to the new counts csv
#' @param fastq_table a recent pull of the database fastq table
#' @param count_file_suffix the suffix appended to the fastqFileName in the count file column headings. default is "_read_count.tsv"
#'
#' @return a list of httr::response() objects
#'
#' @export
postCounts = function(database_counts_url, run_number, auth_token, new_counts_path, fastq_table, count_file_suffix = "_read_count.tsv"){

  # fastqFileNames may not be unique outside of their run
  fastq_table = filter(fastq_table, runNumber == run_number)
  # see utils
  count_df = brentlabRnaSeqTools::readInData(new_counts_path)

  # ensure that count_df is a dataframe of some sort
  stopifnot(sum(class(count_df) == "data.frame") > 0)

  # remove the gene_ids
  count_df = count_df[colnames(count_df) != 'gene_id']
  # remove filename suffix from colnames, leaving just the fastqFileName behind
  colnames(count_df) = str_remove(colnames(count_df), count_file_suffix)

  # remove suffixes from the fastqfiles if they exist
  fastq_table$fastqFileName = str_remove(fastq_table$fastqFileName, ".fastq.gz")

  # halt if there are sample names in the count_df that are not in the database
  stopifnot(setdiff(colnames(count_df), fastq_table$fastqFileName) == 0)

  # create named list with structure list(fastqFileName = fastqFileNumber, ...)
  fastqFileNumber_lookup_list = pull(fastq_table, fastqFileNumber)
  names(fastqFileNumber_lookup_list) = pull(fastq_table, fastqFileName)
  # filter to just those in count_df
  fastqFileNumber_lookup_list = fastqFileNumber_lookup_list[names(fastqFileNumber_lookup_list) %in% colnames(count_df)]

  # check that we still have the same number of samples
  stopifnot(length(fastqFileNumber_lookup_list) == length(colnames(count_df)))

  # send each column to the count table of the database
  res_list = list()
  for (column in colnames(count_df)){
    counts = list(as.integer(pull(count_df, column)))
    names(counts) = column

    post_body = jsonlite::toJSON(list(fastqFileNumber = fastqFileNumber_lookup_list[[column]],
                                      rawCounts = counts), auto_unbox = TRUE)

    res = POST(url=database_counts_url,
               add_headers(Authorization = paste("token" , auth_token, sep=" ")),
               content_type("application/json"),
               body=post_body,
               encode='json')

    res_list[[column]] = res

  }
  res_list
}

# TODO add dry run option

#'
#' post new qc sheet to database
#'
#' @importFrom httr content_type add_headers POST
#' @importFrom jsonlite toJSON
#' @importFrom dplyr rename filter left_join select
#' @importFrom stringr str_remove
#'
#' @description using the package httr, post the new qc sheet to the database
#'
#' @note there can be problems with dependencies and the rename function. this is working for now,
#'       but see here for more info \url{https://statisticsglobe.com/r-error-cant-rename-columns-that-dont-exist}
#'
#' @param database_qc_url eg database_info$kn99_urls$QualityAssess. see \code{\link{database_info}}.
#' @param auth_token \code{\link{getUserAuthToken}}
#' @param run_number the run number of this qc sheet -- this is important b/c fastqFileNames aren't necessarily unique
#'                   outside of their runs
#' @param new_qc_path path to the new counts csv
#' @param fastq_table a recent pull of the database fastq table
#'
#' @return a list of httr::response() objects
#'
#' @export
postQcSheet = function(database_qc_url, auth_token, run_number, new_qc_path, fastq_table) {

  # fastqFileNames may not be unique outside of their run
  fastq_table = filter(fastq_table, runNumber == run_number)

  new_qc_df = brentlabRnaSeqTools::readInData(new_qc_path)

  new_qc_df = new_qc_df %>%
    dplyr::rename(fastqFileName = FASTQFILENAME) %>%
    dplyr::rename(librarySize = LIBRARY_SIZE) %>%
    dplyr::rename(effectiveLibrarySize = EFFECTIVE_LIBRARY_SIZE) %>%
    dplyr::rename(effectiveUniqueAlignment = EFFECTIVE_UNIQUE_ALIGNMENT) %>%
    dplyr::rename(effectiveUniqueAlignmentPercent = EFFECTIVE_UNIQUE_ALIGNMENT_PERCENT) %>%
    dplyr::rename(multiMapPercent = MULTI_MAP_PERCENT) %>%
    dplyr::rename(proteinCodingTotal = PROTEIN_CODING_TOTAL) %>%
    dplyr::rename(proteinCodingTotalPercent = PROTEIN_CODING_TOTAL_PERCENT) %>%
    dplyr::rename(proteinCodingCounted = PROTEIN_CODING_COUNTED) %>%
    dplyr::rename(proteinCodingCountedPercent = PROTEIN_CODING_COUNTED_PERCENT) %>%
    dplyr::rename(ambiguousFeaturePercent = AMBIGUOUS_FEATURE_PERCENT) %>%
    dplyr::rename(noFeaturePercent = NO_FEATURE_PERCENT) %>%
    dplyr::rename(intergenicCoverage = INTERGENIC_COVERAGE) %>%
    dplyr::rename(notAlignedTotalPercent = NOT_ALIGNED_TOTAL_PERCENT) %>%
    dplyr::rename(genotype1Coverage = GENOTYPE1_COVERAGE) %>%
    dplyr::rename(genotype1Log2cpm = GENOTYPE1_LOG2CPM) %>%
    dplyr::rename(genotype2Coverage = GENOTYPE2_COVERAGE) %>%
    dplyr::rename(genotype2Log2cpm = GENOTYPE2_LOG2CPM) %>%
    dplyr::rename(overexpressionFOW = OVEREXPRESSION_FOW) %>%
    dplyr::rename(natCoverage = NAT_COVERAGE) %>%
    dplyr::rename(natLog2cpm = NAT_LOG2CPM) %>%
    dplyr::rename(g418Coverage = G418_COVERAGE) %>%
    dplyr::rename(g418Log2cpm = G418_LOG2CPM) %>%
    dplyr::rename(noMapPercent = NO_MAP_PERCENT) %>%
    dplyr::rename(homopolyFilterPercent = HOMOPOLY_FILTER_PERCENT) %>%
    dplyr::rename(readLengthFilterPercent = READ_LENGTH_FILTER_PERCENT) %>%
    dplyr::rename(tooLowAqualPercent = TOO_LOW_AQUAL_PERCENT) %>%
    dplyr::rename(rRnaPercent = rRNA_PERCENT) %>%
    dplyr::rename(nctrRnaPercent = nctrRNA_PERCENT) %>%
    dplyr::rename(autoStatus = STATUS) %>%
    dplyr::rename(autoAudit = AUTO_AUDIT) %>%
    dplyr::rename(autoStatusDecomp = STATUS_DECOMP)

  # TODO this is copied from postCounts() -- split this off into its own function to avoid repeating
  # remove suffixes from the fastqfiles if they exist
  fastq_table$fastqFileName = str_remove(fastq_table$fastqFileName, ".fastq.gz")

  new_qc_df = new_qc_df %>%
    left_join(fastq_table %>%
                filter(runNumber == run_number) %>%
                select(fastqFileName, fastqFileNumber), by="fastqFileName") %>%
    select(-fastqFileName)

  # check that we still have the same number of samples
  stopifnot(sum(is.na(new_qc_df$fastqFileNumber))==0)

  post_body = jsonlite::toJSON(new_qc_df, auto_unbox = TRUE)

  POST(url = database_qc_url,
       add_headers(Authorization = paste("token" , auth_token, sep=" ")),
       content_type("application/json"),
       body = post_body,
       encode = 'json')
}

# TODO add dry run option

#'
#' PATCH entries in database table
#'
#' @importFrom httr content_type add_headers PATCH
#' @importFrom jsonlite toJSON
#' @importFrom dplyr select
#'
#' @description using the package httr, update entries in certain fields in given rows of a table
#'
#' @param database_table_url NO TRAILING '/'. eg "http://18.224.181.136/api/v1/QualityAssess"
#' @param auth_token see brentlabRnaSeqTools::getUserAuthToken()
#' @param update_df a dataframe, preferrably a tibble, already read in, subsetted. Columns must be correct data type for db table
#' @param id_col name of the id column of the table. this number will be appended to the url to create the uri for the record
#'
#' @return a list of httr::response() objects
#'
#' @export
patchTable = function(database_table_url, auth_token, update_df, id_col){

  # send each column to the count table of the database
  res_list = list()
  for (i in seq(1,nrow(update_df))){

    id = update_df[[i,id_col]]
    row_as_list = jsonlite::toJSON(as.list(dplyr::select(update_df[i,], -id_col)), auto_unbox = TRUE, pretty=TRUE)

    base_url = str_remove(database_table_url, "/$")

    url = paste(base_url, paste0(as.character(id), "/"), sep="/")

    res = PATCH(url=url,
                add_headers(Authorization = paste("token" , auth_token, sep=" ")),
                content_type("application/json"),
                body=row_as_list,
                encode="json")

    res_list[[as.character(id)]] = res

  }
  res_list
}

#'
#' Post a table to the database
#'
#' @importFrom httr content_type add_headers POST
#' @importFrom jsonlite toJSON
#'
#' @param database_table_url see \code{\link{database_info}}. Use one of the URLS in the url slot
#' @param auth_token see \code{\link{getUserAuthToken}}
#' @param df a dataframe read in with, for example read_csv or vroom
#'
#' @return POST results object
#'
#' @export
postTable = function(database_table_url, auth_token, df){

  post_body = jsonlite::toJSON(df, auto_unbox = TRUE)

  POST(url = database_table_url,
       add_headers(Authorization = paste("token" , auth_token, sep=" ")),
       content_type("application/json"),
       body = post_body,
       encode = 'json')

}
