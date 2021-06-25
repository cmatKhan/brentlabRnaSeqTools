#' filter combined_df for environmental perturbation sample set
#'
#' @import magrittr
#' @importFrom dplyr mutate_if mutate filter across
#' @importFrom stringr str_remove
#' @importFrom tidyr replace_na
#'
#' @param combined_df the combined tables of the database, returned directly from getMetadata() (meaning, the df hasn't been augmented after pulling from the database)
#'
#' @return environmental pertubation set
#'
#' @export
createEnvPertSet = function(combined_df){

  # filter
  combined_df %>%
    filter((treatment == "" | treatment == "cAMP" | treatment=="noTreatment" | is.na(treatment)),
           (experimentDesign == 'Environmental_Perturbation' | experimentDesign == 'ep_cAMP_titration'),
           purpose == "fullRNASeq",
           !is.na(fastqFileName),
           genotype1=='CNAG_00000') %>%
    # cast timePoint from integer64 to integer
    mutate(timePoint = as.integer(timePoint)) %>%
    # replace empty strings in the following columns with NA
    mutate(across(c("treatment",
                    "otherConditions",
                    "medium",
                    "atmosphere",
                    "temperature",
                    "timePoint"),
                  ~ifelse(.=="", NA, .) )) %>%
    # replace NA with defined value
    mutate(treatmentConc = replace_na(treatmentConc, 'noTreatmentConc')) %>%
    mutate(pH = replace_na(pH, 'noPH')) %>%
    mutate(treatment = replace_na(treatment, 'noTreatment')) %>%
    mutate(otherConditions = replace_na(otherConditions, 'noOtherConditions')) %>%
    mutate(medium = replace_na(medium, 'noMedium')) %>%
    mutate(atmosphere = replace_na(atmosphere, 'noAtmosphere')) %>%
    mutate(temperature = replace_na(temperature, 'noTemperature')) %>%
    mutate(timePoint = replace_na(timePoint, 'noTimepoint')) %>%
    # remove extention from fastqFileName
    mutate(fastqFileName = str_remove(fastqFileName, ".fastq.gz")) %>%
    # return with uppercase column names
    rename_with(toupper)

}


#' The current definition of the 90 minute induction dataset, according to the 2016 grant summary (loaded into environment, see head(grant_df)) -- single KO only
#'
#' @import magrittr
#' @importFrom dplyr mutate_if mutate filter bind_rows
#' @importFrom stringr str_remove
#' @importFrom tidyr replace_na
#'
#' @param metadata is the combined tables of the metadata database
#' @param grant_df is the 2016 grant summary TODO: put this in DATA
#'
#' @return the set metadata -- single KO only
#'
#' @export
createNinetyMinuteInductionSet = function(metadata, grant_df){


  metadata = metadata %>%
    # replace empty strings with NA
    mutate_if(is.character, list(~na_if(.,""))) %>%
    # replace NAs with string entry
    mutate(treatment = replace_na(treatment, 'noTreatment')) %>%
    mutate(otherConditions = replace_na(otherConditions, 'noOtherConditions')) %>%
    mutate(medium = replace_na(medium, 'noMedium')) %>%
    mutate(atmosphere = replace_na(atmosphere, 'noAtmosphere')) %>%
    mutate(pH = replace_na(pH, 'noPh'))

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
  wt_induction_set = condition_fltr_metadata %>%
    filter(genotype1 == "CNAG_00000", strain == 'TDY451')
  # filter for genotypes in the grant summary
  perturbed_induction_set = condition_fltr_metadata %>%
    filter(genotype1 %in% grant_df$GENOTYPE1 & is.na(genotype2))

  # put the wt and filtered genotypes together
  induction_set = bind_rows(wt_induction_set, perturbed_induction_set)

  # set colnames to upper
  colnames(induction_set) = toupper(colnames(induction_set))
  # remove file extension
  induction_set$FASTQFILENAME = str_remove(induction_set$FASTQFILENAME, ".fastq.gz")

  return(induction_set)
}

#' The current definition of the 90 minute induction dataset, according to the 2016 grant summary (loaded into environment, see head(grant_df)) -- single and double KO
#'
#' @import magrittr
#' @importFrom dplyr mutate_if mutate filter bind_rows
#' @importFrom stringr str_remove
#' @importFrom tidyr replace_na
#'
#' @param metadata is the combined tables of the metadata database
#' @param grant_df is the 2016 grant summary TODO: put this in DATA
#'
#' @return the set metadata
#'
#' @export
createNinetyMinuteInductionWithDoubles = function(metadata, grant_df){

  metadata = metadata %>%
    mutate(treatment = replace_na(treatment, 'noTreatment')) %>%
    mutate(otherConditions = replace_na(otherConditions, 'noOtherConditions')) %>%
    mutate(medium = replace_na(medium, 'noMedium')) %>%
    mutate(atmosphere = replace_na(atmosphere, 'noAtmosphere')) %>%
    mutate(pH = replace_na(pH, 'noPh'))

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


