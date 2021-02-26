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
  sorted_passing_meta_qual = passing_qc1_meta_qual %>% group_by(LIBRARYPROTOCOL, LIBRARYDATE) %>% arrange(LIBRARYDATE, .by_group = TRUE)

  sorted_passing_induction_raw_counts = raw_counts[, passing_qc1_meta_qual$FASTQFILENAME]

  old_protocol_sorted_passing_meta_qual = sorted_passing_meta_qual %>% filter(LIBRARYPROTOCOL == "SolexaPrep")
  old_protocol_sorted_passing_counts = raw_counts[, old_protocol_sorted_passing_meta_qual$FASTQFILENAME]

  new_protocol_sorted_passing_meta_qual = sorted_passing_meta_qual %>% filter(LIBRARYPROTOCOL == "E7420L")
  new_protocol_sorted_passing_counts = raw_counts[, new_protocol_sorted_passing_meta_qual$FASTQFILENAME]

  libraryprotocol_librarydate_model_matrix = createLibrarydateModelMatrix(sorted_passing_meta_qual)

  old_dds = DESeqDataSetFromMatrix(colData = old_protocol_sorted_passing_meta_qual, countData = old_protocol_sorted_passing_counts, design=~1)
  new_dds = DESeqDataSetFromMatrix(colData = new_protocol_sorted_passing_meta_qual, countData = new_protocol_sorted_passing_counts, design=~1)

  old_dds = estimateSizeFactors(old_dds)

  new_dds = estimateSizeFactors(new_dds)

  size_factor_list = c(sizeFactors(new_dds), sizeFactors(old_dds))

  dds = DESeqDataSetFromMatrix(colData = sorted_passing_meta_qual, countData = sorted_passing_induction_raw_counts, design=libraryprotocol_librarydate_model_matrix)

  sizeFactors(dds) = size_factor_list

  stopifnot(all.equal(names(size_factor_list), sorted_passing_meta_qual$FASTQFILENAME, colnames(sorted_passing_induction_raw_counts)))

  return(dds)
}
