#' create deseq object with protocol specific size factors
#'
#' @import DESeq2
#' @import dplyr
#'
#' @param passing_qc1_meta_qual can be any metadata df, but if you're going to run deseq you may want to filter it for passing samples first
#' @param raw_counts a dataframe of raw counts with genes in the rows and samples in the columns. sample names must be the same as the fastqFileName
#'                   column in passing_qc1_meta_qual.
#' @return a deseq data object with size factors calculated within the library protocol groups
#'
#' @export
deseqObjectWithProtocolSpecificSizeFactors = function(passing_qc1_meta_qual, raw_counts){

  colnames(passing_qc1_meta_qual) = toupper(colnames(passing_qc1_meta_qual))
  sorted_passing_meta_qual = passing_qc1_meta_qual %>% group_by(LIBRARYDATE, LIBRARYPROTOCOL) %>% arrange(.by_group = TRUE)

  sorted_passing_induction_raw_counts = raw_counts[, sorted_passing_meta_qual$FASTQFILENAME]

  old_protocol_sorted_passing_meta_qual = sorted_passing_meta_qual %>% filter(LIBRARYPROTOCOL == "SolexaPrep")
  old_protocol_sorted_passing_counts = sorted_passing_induction_raw_counts[, old_protocol_sorted_passing_meta_qual$FASTQFILENAME]

  new_protocol_sorted_passing_meta_qual = sorted_passing_meta_qual %>% filter(LIBRARYPROTOCOL == "E7420L")
  new_protocol_sorted_passing_counts = sorted_passing_induction_raw_counts[, new_protocol_sorted_passing_meta_qual$FASTQFILENAME]

  libraryprotocol_librarydate_model_matrix = createNinetyMinInductionModelMatrix(sorted_passing_meta_qual)

  old_dds = DESeqDataSetFromMatrix(colData = old_protocol_sorted_passing_meta_qual, countData = old_protocol_sorted_passing_counts, design=~1)
  new_dds = DESeqDataSetFromMatrix(colData = new_protocol_sorted_passing_meta_qual, countData = new_protocol_sorted_passing_counts, design=~1)

  old_dds = estimateSizeFactors(old_dds)

  new_dds = estimateSizeFactors(new_dds)

  size_factor_list = c(sizeFactors(old_dds), sizeFactors(new_dds))

  stopifnot(all.equal(names(size_factor_list), sorted_passing_meta_qual$FASTQFILENAME, colnames(sorted_passing_induction_raw_counts), rownames(libraryprotocol_librarydate_model_matrix)) )

  dds = DESeqDataSetFromMatrix(colData = sorted_passing_meta_qual, countData = sorted_passing_induction_raw_counts, design=libraryprotocol_librarydate_model_matrix)

  sizeFactors(dds) = size_factor_list

  stopifnot(all.equal(names(sizeFactors(dds)), as_tibble(colData(dds))$FASTQFILENAME, colnames(counts(dds)), rownames(design(dds))) )

  return(dds)
}

#'
#' temporary function to examine only PBS samples, eg
#'
#' @import DESeq2
#' @import dplyr
#'
#' @description gets all samples in qc_passing_metadata with size_factor_subset_param same as samples in row_filter
#'
#' @param qc1_passing_metadata metadata passing auto/manual flags
#' @param raw_counts must include at least the sames in the metadata
#' @param row_filter is a boolean vector, eg qc1_passing_env_pert_meta_qual$MEDIUM == "PBS"
#' @param size_factor_subset_param eg LIBRARYDATE -- the column to group the samples for size factor calculation
#' @return a list with slots metadata, raw_counts, size_factors
#' @export
examineSingleGroupWithLibDateSizeFactors = function(qc1_passing_metadata, raw_counts, row_filter, size_factor_subset_param){

  qc1_passing_metadata$LIBRARYDATE = as.factor(qc1_passing_metadata$LIBRARYDATE)

  subset_of_interest = qc1_passing_metadata[ row_filter , ]

  size_factor_grouping_vector = subset_of_interest[, size_factor_subset_param][[size_factor_subset_param]]

  all_samples_with_subset_param_metadata = qc1_passing_metadata %>%
    filter(!! rlang::sym(size_factor_subset_param) %in% size_factor_grouping_vector)

  all_samples_with_subset_param_counts = raw_counts[, all_samples_with_subset_param_metadata$FASTQFILENAME]

  if(!all.equal(colnames(all_samples_with_subset_param_counts), all_samples_with_subset_param_metadata$FASTQFILENAME)){
    stop("counts do not equal rows of metadata")
  }

  split_list = split(all_samples_with_subset_param_metadata, all_samples_with_subset_param_metadata[ ,size_factor_subset_param])

  size_factor_list = list()

  for (subset_param_level in names(split_list)){
    metadata = split_list[[subset_param_level]]
    if(nrow(metadata) == 0){
      print(paste0("does not exist in split: " ,subset_param_level))
    } else{
      counts = raw_counts[, metadata$FASTQFILENAME]
      split_dds = DESeqDataSetFromMatrix(countData = counts, colData = metadata, design=~1)
      split_dds = estimateSizeFactors(split_dds)
      split_sf = sizeFactors(split_dds)
      size_factor_list = append(size_factor_list, split_sf)
    }
  }

  if(!all.equal(colnames(all_samples_with_subset_param_counts), all_samples_with_subset_param_metadata$FASTQFILENAME, names(size_factor_list))){
    stop("columns of counts, FASTQFILENAMES, and/or names(size_factor_list) are not equal")
  }

  return (list("metadata" = all_samples_with_subset_param_metadata,
               "raw_counts" = all_samples_with_subset_param_counts,
               "size_factors" = size_factor_list))


}

#' remove some effects from the counts
#'
#' @import DESeq2
#'
#' @description subtract effect from norm counts of a single factor from coef x design. coef is in normalized log space. dds must have been created with model.matrix
#' @note works for both formula and model.matrix designs in the dds object
#'
#' @param deseq_object a deseq data object REQUIRED: the object must have been created with a model.matrix rather than a formula for the design argument
#' @param col_indicies a numeric vector corresponding to the column indicies of the batch parameters you'd like to remove
#'
#' @return a log2 scale gene by samples matrix with desired effects removed
#'
#' @export
removeParameterEffects = function(deseq_object, col_indicies){

  # if the design(dds) is a formula
  if(is_formula(design(deseq_object))){
    model_matrix = model.matrix(design(deseq_object), colData(deseq_object))
  } else if(is.matrix(design(deseq_object))){
    model_matrix = design(deseq_object)
  } else{
    stop("design(deseq_object) is not recognized as a formula or matrix")
  }

  coefficients = coef(deseq_object)[,col_indicies]
  batch_effect_matrix = model_matrix[,col_indicies]

  log_norm_counts = log2(counts(deseq_object, normalized=TRUE)+1) # note psuedocount

  # coefficients is dim gene x features
  # design_matrix is dim sample x features to remove
  return(log_norm_counts - (coefficients %*% t(batch_effect_matrix)))

}
