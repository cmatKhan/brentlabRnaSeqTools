#' create 90 minute induction set tally
#'
#' @import dplyr
#'
#' @param induction_meta_qual the metadata of the entire set, unfiltered
#' @param sorted_passing_induction_meta_qual metadata (with quality columns) filtered for manual/auto status
#' @param iqr_fltr_rle_summary sorted_passing_meta_qual filtered for IQR
#' @param grant_df the definition of the 90minuteInduction set. This object is available in the brentlabRnaSeqTools
#'                 package
#'
#' @export
createInductionSetTally = function(induction_meta_qual, sorted_passing_induction_meta_qual, iqr_fltr_rle_summary, grant_df){

  # add GENOTYPE column -- this is remnant of old system, but kept b/c it might help with concat double KO
  induction_meta_qual$GENOTYPE = induction_meta_qual$GENOTYPE1
  sorted_passing_induction_meta_qual$GENOTYPE = sorted_passing_induction_meta_qual$GENOTYPE1
  iqr_fltr_rle_summary$GENOTYPE = iqr_fltr_rle_summary$GENOTYPE1

  induction_samples_genotype_tally = induction_meta_qual %>% group_by(GENOTYPE) %>% tally()
  colnames(induction_samples_genotype_tally)[2] = "complete_set_no_fltr"

  qc1_passing_tally = sorted_passing_induction_meta_qual %>% group_by(GENOTYPE) %>% tally()
  colnames(qc1_passing_tally)[2] = "qc_passing"

  qc1_iqr_passing_tally = iqr_fltr_rle_summary %>% group_by(GENOTYPE) %>% tally()
  colnames(qc1_iqr_passing_tally)[2] = "qc_passing_iqr_filtered"

  # add genotypes from grant that have no replicates in the database at all
  done_grant_strains_df = grant_df %>% filter(STRAIN_STATUS !=2)
  # TODO: CHECK THAT GRANT_DF COLUMN IS GENOTYPE1
  no_geno_df = unique(done_grant_strains_df %>% filter(!GENOTYPE1 %in% unique(induction_meta_qual$GENOTYPE1)) %>% select(GENOTYPE1))
  no_replicate_genotypes = no_geno_df$GENOTYPE1

  no_replicate_genotypes_tally_df = tibble(GENOTYPE=no_replicate_genotypes,
                                           complete_set_no_fltr = rep(0, length(no_replicate_genotypes)),
                                           qc_passing = rep(0, length(no_replicate_genotypes)),
                                           qc_passing_iqr_filtered = rep(0, length(no_replicate_genotypes))
  )

  genotype_tally_summary = induction_samples_genotype_tally %>%
    # join the passing and iqr filter tables
    left_join(qc1_passing_tally) %>%
    left_join(qc1_iqr_passing_tally) %>%
    # pass values from previous column forward, if there are NA in preceeding columns
    # mutate(qc_passing = coalesce(qc_passing, complete_set_no_fltr)) %>%  # don't do this, because we don't want to pass a value from complete_set_no_fltr to qc_passing if qc_passing is 0
    mutate(qc_passing_iqr_filtered = coalesce(qc_passing_iqr_filtered, qc_passing)) %>%
    bind_rows(no_replicate_genotypes_tally_df)

  genotype_tally_summary = genotype_tally_summary %>%
    dplyr::mutate(complete_set_no_fltr = replace_na(complete_set_no_fltr, 0)) %>%
    dplyr::mutate(qc_passing = replace_na(qc_passing, 0)) %>%
    dplyr::mutate(qc_passing_iqr_filtered = replace_na(qc_passing_iqr_filtered, 0))

  return(genotype_tally_summary)

}
