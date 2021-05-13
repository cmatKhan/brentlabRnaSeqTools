# group1=c("IPSC_ASDDM",
#          "IPSC_ASDDM",
#          "IPSC_ASDN",
#          "NPC_ASDDM",
#          "NPC_ASDDM",
#          "NPC_ASDN")
# group2=c("IPSC_ASDN",
#          "IPSC_TDN",
#          "IPSC_TDN",
#          "NPC_ASDN",
#          "NPC_TDN",
#          "NPC_TDN")
#
# comparisons=c("IPSC_ASDDM_ASDN",
#               "IPSC_ASDDM_TDN",
#               "IPSC_ASDN_TDN",
#               "NPC_ASDDM_ASDN",
#               "NPC_ASDDM_TDN",
#               "NPC_ASDN_TDN")
#
# ##get the results for the various contrasts
# numcomparisons=length(comparisons)
# for(i in seq(1,numcomparisons))
# {
#   res=results(dds, contrast=c("Group", group1[i],group2[i]),parallel=TRUE)
#   res$logPadj=-1*log10(res$padj)
#   res=as.data.frame(res)
#   res=na.omit(res)
#   print(comparisons[i])
#   print(res['ENSG00000196776.16-CD47',])
# }

#' Create a deseq data object with library protocol specific size factors
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
#' @description gets all samples in qc_passing_metadata with size_factor_subset_param same as samples in row_filter
#'
#' @param row_filter is a boolean vector, eg qc1_passing_env_pert_meta_qual$MEDIUM == "PBS"
#' @param size_factor_subset_param eg LIBRARYDATE -- the column to group the samples for size factor calculation
#' @return a list with slots metadata, raw_counts, size_factors
#' @export
examineSingleGroupWithLibDateSizeFactors = function(qc1_passing_metadata, raw_counts, row_filter, size_factor_subset_param){

  qc1_passing_metadata$LIBRARYDATE = as.factor(qc1_passing_metadata$LIBRARYDATE)

  subset_of_interest = qc1_passing_metadata[ row_filter , ]

  size_factor_grouping_vector = subset_of_interest[, size_factor_subset_param][[size_factor_subset_param]]

  all_samples_with_subset_param_metadata = qc1_passing_metadata %>% filter(!! rlang::sym(size_factor_subset_param) %in% size_factor_grouping_vector)

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
      sizeFactors(split_dds)
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
