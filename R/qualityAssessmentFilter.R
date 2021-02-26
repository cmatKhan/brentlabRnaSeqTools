#' filter for manual passes (overrides auto fail) and automatic passes (unless auto failed)
#'
#' @param a metadata dataframe from the database
#' @return a metadata dataframe with column names cast to upper
#' @export
qualityAssessmentFilter = function(metadata){
  colnames(metadata) = toupper(colnames(metadata))
  passing_metadata = metadata %>% filter(MANUALAUDIT == FALSE | (is.na(MANUALAUDIT) & AUTOAUDIT == FALSE) )

  return(passing_metadata)
}
