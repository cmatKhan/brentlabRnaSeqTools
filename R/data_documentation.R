#' URLS to active databases
#'
#' A list containing the urls to active databases. Named by organism (eg 'kn99' or 's288cr64')
#'
#' @format A list with named slots
#' \describe{
#'   \item{kn99_host}{host of the database server, eg ec2-18-224-181-136.us-east-2.compute.amazonaws.com}
#'   \item{kn99_db_name}{cryptococcus database name. probably something like kn99_database}
#'   \item{s288cr64_host}{host address of the yeast s288cr64 database, eg ec2-3-131-85-10.us-east-2.compute.amazonaws.com}
#'   \item{s288cr64_db_name}{yeast database name. probably something like yeast_database}
#'   ...
#' }
#' @source \url{https://rnaseq-databases-documentation.readthedocs.io/en/latest/}
"database_info"
