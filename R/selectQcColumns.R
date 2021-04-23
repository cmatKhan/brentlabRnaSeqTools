#' select fastqFileName, fastqFileNumber, and a pre-determined set of QC columns from a metadata df
#' @note the column names for the metadata will be cast to upper and returned in upper
#' @note must include Interquartile range (think about removing this -- user could merge with IQR df after selecting these cols)
#'
#' @param metadata a metadata df with at least the columns listed in the select statement (see source code -- notably, must include interquartile range). Column names will be cast to uppper and returned as uppers
#'
#' @export
selectQaColumns = function(metadata){

  # cast column names to upper in case they have not yet been
  colnames(metadata) = toupper(colnames(metadata))

  # todo: set the columns to vectors, cast to symbols with rlang and offer argument for custom list
  selectMetadata = metadata %>%
    select(FASTQFILENAME,
           FASTQFILENUMBER,
           AUTOSTATUSDECOMP,
           MANUALAUDIT,
           INTERQUARTILE_RANGE,
           TOTAL_DEDUPLICATED_PERCENTAGE,
           PROTEIN_CODING_TOTAL,
           LIBRARYDATE,
           NAT_COVERAGE,
           NAT_LOG2CPM,
           G418_COVERAGE,
           G418_LOG2CPM,
           NOT_ALIGNED_TOTAL_PERCENT,
           NO_FEATURE_PERCENT,
           INTERGENIC_COVERAGE,
           LIBRARY_SIZE,
           EFFECTIVE_UNIQUE_ALIGNMENT)

  return(selectMetadata)
}
