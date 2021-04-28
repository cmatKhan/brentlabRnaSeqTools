#'
#' run SVA
#' @param metadata metadata in the shape samples x parameters where the columns are features of the data (eg, libraryDate, RIN, etc)
#' @param raw_counts raw gene counts in the shape gene x samples where ncols matches nrow of metadata (samples == samples)
#' @param known_covariate_formula eg ~Condition+Cell, passed as an R formula
#'
#' @export
runSVA = function(metadata, raw_counts, known_covariates_formula){

  # remove genes with 0 counts
  raw_counts=raw_counts[rowSums(data)>0,]

  # create design matricies
  mod0=model.matrix(~1,data=metadata) # make intercept only model
  mod1=model.matrix(known_covariates_formula,data=metadata) # model with known covariates

  # run sva
  sva_obj=svaseq(as.matrix(data),mod1,mod0) # run sva

  return(sva_obj)

}
