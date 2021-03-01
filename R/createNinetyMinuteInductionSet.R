#' The current definition of the 90 minute induction dataset
#'
#' @param metadata is the combined tables of the metadata database
#' @param grant_df is the 2016 grant summary TODO: put this in DATA
#'
#' @return the 90minuteInduction metadata
#'
#' @export
createNinetyMinuteInductionSet = function(metadata, grant_df){
  metadata = metadata %>% dplyr::mutate(treatment = replace_na(treatment, 'noTreatment'))
  metadata = metadata %>% dplyr::mutate(otherConditions = replace_na(otherConditions, 'noOtherConditions'))
  metadata = metadata %>% dplyr::mutate(medium = replace_na(medium, 'noMedium'))
  metadata = metadata %>% dplyr::mutate(atmosphere = replace_na(atmosphere, 'noAtmosphere'))
  metadata = metadata %>% dplyr::mutate(pH = replace_na(pH, 'noPh'))

  metadata$temperature %<>% as.numeric
  metadata$timePoint %<>% as.numeric

  condition_fltr_metadata = metadata %>%
    filter(medium %in% c("DMEM"),
           temperature %in% c(37),
           atmosphere %in% c("CO2"),
           treatment %in% c("noTreatment"),
           otherConditions %in% c("noOtherConditions"),
           pH %in% c("noPh"),
           timePoint %in% c(90),
           strain != "TDY1993",
           purpose=="fullRNASeq",
           !is.na(fastqFileName),
           str_detect(genotype1, "CNAG"),
           is.na(genotype2),
           perturbation1 =="deletion" | is.na(perturbation1))

  # get wildtypes
  wt_induction_set = condition_fltr_metadata %>% filter(genotype1 == "CNAG_00000")
  # filter for genotypes in the grant summary
  perturbed_induction_set = condition_fltr_metadata %>% filter(genotype1 %in% grant_df$GENOTYPE1 & is.na(genotype2))
  # put the wt and filtered genotypes together
  induction_set = bind_rows(wt_induction_set, perturbed_induction_set)
  # set colnames to upper
  colnames(induction_set) = toupper(colnames(induction_set))
  # remove file extension
  induction_set$FASTQFILENAME = str_remove(induction_set$FASTQFILENAME, ".fastq.gz")

  return(induction_set)
}
