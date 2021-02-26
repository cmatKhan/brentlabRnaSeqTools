#' remove some effect(s) from the counts
#'
#' subtract effect from norm counts of a single factor from coef x design. coef is in normalized log space
#'
#' @param deseq_object a deseq data object
#' @param design_matrix_of_effects_to_remove: A possibly modified design matrix with any effects column that we do not want to remove set to 0
#'
#' @note the design matrix will be transposed, and the intercept removed, in the function. Don't do that before
#'
#' @return a log2 scale gene by samples matrix with desired effects removed
#'
#' @export
removeEffect = function(deseq_object, design_matrix_of_effects_to_remove_with_intercept){
  #' :params design_matrix_of_effects_to_remove: A possibly modified design matrix with any effects column that we do not want to remove set to 0
  #' NOTE: it will be transposed, and the intercept removed, in the function. Don't do that before
  #' subtract effect from norm counts of a single factor from coef x design. coef is in normalized log space

  coefficients = coef(deseq_object)
  log_norm_counts = log2(counts(deseq_object, normalized=TRUE)+1) # note psuedocount
  # coefficients is dim gene x features(including intercept, which is dropped)
  # design_matrix is dim sample x features to remove (plus one, for the intercept, which is removed below)
  return(log_norm_counts - (coefficients[,-1] %*% t(design_matrix_of_effects_to_remove_with_intercept[,-1])))

} # end removeEffect
