#' remove some effects from the counts
#'
#' subtract effect from norm counts of a single factor from coef x design. coef is in normalized log space
#' @param deseq_object a deseq data object
#' @param design_matrix_of_effects_to_remove: A possibly modified design matrix with any effects column that we do not want to remove set to 0
#' @note the design matrix will be transposed, and the intercept removed, in the function. Don't do that before
#' @return a log2 scale gene by samples matrix with desired effects removed
#'
#' @export
removeParameterEffects = function(deseq_object, design_matrix_of_effects_to_remove_with_intercept){

  coefficients = coef(deseq_object)
  log_norm_counts = log2(counts(deseq_object, normalized=TRUE)+1) # note psuedocount
  # coefficients is dim gene x features(including intercept, which is dropped)
  # design_matrix is dim sample x features to remove (plus one, for the intercept, which is removed below)
  return(log_norm_counts - (coefficients[,-1] %*% t(design_matrix_of_effects_to_remove_with_intercept[,-1])))

}



# update this function to the below. In specific sets, maybe maybe a wrapper to this that takes the dds, but this allows for more flexibility

#'
#'
#' @param effects_to_remove calculated similar to this. note -1 in the column indexing means that the intercept effect is retained in the data: (coef(deseq_object)[,-1] %*% t(design_matrix_of_effects_to_remove_with_intercept[,-1])
#'
#'
#'
# removeParameterEffects = function(counts, effects_to_remove){
#
#   coefficients = coef(deseq_object)
#   log_norm_counts = log2(counts(deseq_object, normalized=TRUE)+1) # note psuedocount
#   # coefficients is dim gene x features(including intercept, which is dropped)
#   # design_matrix is dim sample x features to remove (plus one, for the intercept, which is removed below)
#   return(log_norm_counts - effects_to_remove))
#
# }
