#' Get the combined metadata as a tibble from a remote database
#'
#' @description Join the biosample, rnasample, s1sample, s2sample, library, fastqFiles and qualityAssessment tables (in that order, left joins) and return the result as a tibble
#'
#' Use the RPostgres package to connect to a remote postgresql database, do the table joining, and return the joined metadata as a tibble. The database connection is closed
#' @usage combined_df = getMetadata(database_urls$kn99_host, database_urls$kn99_db_name, Sys.getenv("db_user"), Sys.getenv("db_password"))
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

  # TODO fix merging above so that column names aren't duplicated in qualityAssess table (fastqFileName in both)
  # deal idiosyncratically with instances in which fastqFiles and qualityAssess fastqFileName may not overlap -- this shouldn't happen, really -- need to fix in database

  # note: there are conflicts in dependencies, hence explicitely calling select and rename. see here for rename https://statisticsglobe.com/r-error-cant-rename-columns-that-dont-exist
  metadata_df = metadata_df %>%
    dplyr::select(-fastqFileName.y) %>%
    plyr::rename(c(fastqFileName.x = 'fastqFileName'))

  dbDisconnect(db)

  return(metadata_df)
}

#' Get combined raw counts
#'
#' GET raw counts
#' @usage getRawCounts(api_url)
#' @param api_url NOTE: api_url is a variable saved into the project environment, in addition to the parameter. You can use the usage statement directly. The argument is provided for development in the event that you want to test a local instance of the database. An example url: "http://13.59.167.2/api" No trailing /
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
#' saves both the individual tables, including counts, and the combined_df
#'
#' @param database_host if connecting to a database hosted on AWS, it might be something like ec2-54-83-201-96.compute-1.amazonaws.com
#' @param database_name name of the database, eg for cryptococcus kn99, the database might be named kn99_database. Check with the documentation, whoever set up the database, or get into the server and check directly
#' @param database_user a user of the actual database, with some level of permissions. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param database_password password to the database user. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param output_dir where to deposit a subdirectory, named by todays date in this format: 20210407, with the tables and combined_df inside. eg a mounted local directory /mnt/htcf_lts/crypto_database_archive/ --> /lts/mblab/Crypto/rnaseq_data/crypto_database_archive
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

#' Connect to a remote postgresql database
#'
#' Use the RPostgres package to connect to a remote postgresql database
#' @usage kn99_database = connectToDatabase(database_urls$kn99_host, database_urls$kn99_db_name, Sys.getenv("db_user"), Sys.getenv("db_password"))
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
