#' filter combined_df for environmental perturbation sample set
#'
#' @param combined_df the combined tables of the database
#'
#' @return environmental pertubation set
#'
#' @export
createEnvPertSet = function(combined_df){
  env_pert_set = combined_df %>% filter((is.na(treatment) | treatment == "cAMP"), experimentDesign == 'Environmental_Perturbation', purpose == "fullRNASeq", !is.na(fastqFileName))
  colnames(env_pert_set) = toupper(colnames(env_pert_set))

  env_pert_set = env_pert_set %>% dplyr::mutate(TREATMENT = replace_na(TREATMENT, 'noTreatment'))
  env_pert_set = env_pert_set %>% dplyr::mutate(OTHER_CONDITIONS = replace_na(OTHERCONDITIONS, 'noOtherConditions'))
  env_pert_set = env_pert_set %>% dplyr::mutate(MEDIUM = replace_na(MEDIUM, 'noMedium'))
  env_pert_set = env_pert_set %>% dplyr::mutate(ATMOSPHERE = replace_na(ATMOSPHERE, 'noAtmosphere'))
  env_pert_set = env_pert_set %>% dplyr::mutate(TEMPERATURE = replace_na(TEMPERATURE, 'noTemperature'))
  env_pert_set = env_pert_set %>% dplyr::mutate(TIMEPOINT = replace_na(TIMEPOINT, 'noTimepoint'))

  return(env_pert_set)
}
