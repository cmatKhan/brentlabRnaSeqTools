#' filter combined_df for environmental perturbation sample set
#'
#' @param combined_df the combined tables of the database
#'
#' @return environmental pertubation set
#'
#' @export
createEnvPertSet = function(combined_df){
  env_pert_set = combined_df %>% filter((is.na(treatment) | treatment == "cAMP"),
                                        experimentDesign == 'Environmental_Perturbation',
                                        purpose == "fullRNASeq", !is.na(fastqFileName),
                                        genotype1=='CNAG_00000')
  colnames(env_pert_set) = toupper(colnames(env_pert_set))

  env_pert_set = env_pert_set %>% dplyr::mutate(TREATMENT = replace_na(TREATMENT, 'noTreatment'))
  env_pert_set = env_pert_set %>% dplyr::mutate(TREATMENTCONC = replace_na(TREATMENTCONC, 'noTreatmentConc'))
  env_pert_set = env_pert_set %>% dplyr::mutate(OTHER_CONDITIONS = replace_na(OTHERCONDITIONS, 'noOtherConditions'))
  env_pert_set = env_pert_set %>% dplyr::mutate(MEDIUM = replace_na(MEDIUM, 'noMedium'))
  env_pert_set = env_pert_set %>% dplyr::mutate(ATMOSPHERE = replace_na(ATMOSPHERE, 'noAtmosphere'))
  env_pert_set = env_pert_set %>% dplyr::mutate(TEMPERATURE = replace_na(TEMPERATURE, 'noTemperature'))
  env_pert_set = env_pert_set %>% dplyr::mutate(TIMEPOINT = replace_na(TIMEPOINT, 'noTimepoint'))

  return(env_pert_set)
}

#' The current definition of the 90 minute induction dataset, according to the 2016 grant summary (loaded into environment, see head(grant_df)) -- single KO only
#'
#' @param metadata is the combined tables of the metadata database
#' @param grant_df is the 2016 grant summary TODO: put this in DATA
#'
#' @return the set metadata -- single KO only
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

  # TODO: COMBINE THE FILTER STATEMENTS INTO SINGLE FILTER
  # get wildtypes
  wt_induction_set = condition_fltr_metadata %>% filter(genotype1 == "CNAG_00000", strain == 'TDY451')
  # filter for genotypes in the grant summary
  perturbed_induction_set = condition_fltr_metadata %>% filter(genotype1 %in% grant_df$GENOTYPE1 & is.na(genotype2))
  # put the wt and filtered genotypes together
  induction_set = bind_rows(wt_induction_set, perturbed_induction_set)
  #####

  # set colnames to upper
  colnames(induction_set) = toupper(colnames(induction_set))
  # remove file extension
  induction_set$FASTQFILENAME = str_remove(induction_set$FASTQFILENAME, ".fastq.gz")

  return(induction_set)
}

#' The current definition of the 90 minute induction dataset, according to the 2016 grant summary (loaded into environment, see head(grant_df)) -- single and double KO
#'
#' @param metadata is the combined tables of the metadata database
#' @param grant_df is the 2016 grant summary TODO: put this in DATA
#'
#' @return the set metadata
#'
#' @export
createNinetyMinuteInductionWithDoubles = function(metadata, grant_df){
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
           perturbation1 =="deletion" | is.na(perturbation1),
           perturbation2 =="deletion" | is.na(perturbation2))

  # TODO: COMBINE THE FILTER STATEMENTS INTO SINGLE FILTER
  # get wildtypes
  wt_induction_set = condition_fltr_metadata %>% filter(genotype1 == "CNAG_00000", strain == 'TDY451')
  # filter for genotypes in the grant summary
  perturbed_induction_set = condition_fltr_metadata %>% filter(genotype1 %in% grant_df$GENOTYPE1 & (genotype2 %in% grant_df$GENOTYPE1 | is.na(genotype2)))
  # put the wt and filtered genotypes together
  induction_set = bind_rows(wt_induction_set, perturbed_induction_set)
  #####

  # set colnames to upper
  colnames(induction_set) = toupper(colnames(induction_set))
  # remove file extension
  induction_set$FASTQFILENAME = str_remove(induction_set$FASTQFILENAME, ".fastq.gz")

  return(induction_set)
}


