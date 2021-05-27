#'
#' get total protein coding count from count dataframe
#' @descriptiongiven a count dataframe with gene_ids as rownames and quantification in a column called raw_counts, return sum of protein coding genes
#'
#' @param counts a dataframe with gene_ids in the rownames and (at minimum) a quantification column called raw_counts
#'
proteinCodingCount = function(counts, protein_coding_gene_ids){
  sum(counts[protein_coding_gene_ids,'raw_counts'])
}
