#'
#' run SVA
#'
#' @import sva
#'
#' @param raw_counts raw gene counts in the shape gene x samples where ncols matches nrow of metadata (samples == samples)
#' @param null_model_matrix a model matrix respresenting only the batch effects. Could possibly be intercept only
#' @param full_model_matrix the full model describing the experiment. Critically, this includes the parameter of interest
#'
#' @export
runSVA = function(raw_counts, null_model_matrix, full_model_matrix){

  # remove genes with 0 counts
  raw_counts=raw_counts[rowSums(raw_counts)>0,]

  # run sva
  sva_obj=svaseq(as.matrix(raw_counts),full_model_matrix,null_model_matrix) # run sva

  return(sva_obj)

}
