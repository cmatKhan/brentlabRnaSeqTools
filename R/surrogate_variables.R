#'
#' run SVA
#' @param metadata metadata in the shape samples x parameters where the columns are features of the data (eg, libraryDate, RIN, etc)
#' @param raw_counts raw gene counts in the shape gene x samples where ncols matches nrow of metadata (samples == samples)
#' @param known_covariate_formula eg ~Condition+Cell, passed as an R formula
#'
#' @export
runSVA = function(raw_counts, null_model_matrix, full_model_matrix){

  # remove genes with 0 counts
  raw_counts=raw_counts[rowSums(raw_counts)>0,]

  # run sva
  sva_obj=svaseq(as.matrix(raw_counts),full_model_matrix,null_model_matrix) # run sva

  return(sva_obj)

}
