#'
#' create nf-co sample sheet
#'
#' @description create a samplesheet for the nf-co rnaseq pipeline
#'              \url{https://nf-co.re/rnaseq/3.2/usage#samplesheet-input}
#'
#' @note: currently hard coded to set strandedness to 'reverse' if libraryProtocol is E7420L, or 'unstranded',
#'        AND for single strand libraries (fastq2 is an empty string)
#'
#' @param metadata most likely from the database. must include runNumber, fastqFileName, libraryProtocol,
#'                 and the columns used to create the sample name (see param sample_columns)
#' @param sample_columns a character vector of columns to concat into a sample identifier.
#'                       eg c("genotype1", "treatment", "timePoint", "floodmedia", "inductionDelay")
#' @param sequence_dir_prefix path to directory where the sequence run directories are stored.
#'                            eg /mnt/htcf_lts/lts_sequence
#' @param check_files_flag default FALSE. Set to true to check that files exist
#'
#' @return a dataframe with the columns which correspond to those required for the nf-co rnaseq pipeline
#'         sample sheet. See \url{https://nf-co.re/rnaseq/3.2/usage#samplesheet-input}
#'
#' @export
createNfCorePipelineSampleSheet = function(metadata, sample_columns, sequence_dir_prefix, check_files_flag=FALSE){

  # remove trailing / if it exists
  sequence_dir_prefix = gsub("/$", "", sequence_dir_prefix)

  sample_sheet_df = metadata %>%
    mutate(strandedness = ifelse(libraryProtocol == "E7420L", "reverse", "unstranded")) %>%
    mutate(runNumber = unlist(lapply(runNumber, getRunNumberLeadingZero))) %>%
    mutate(fastq_1 = file.path(sequence_dir_prefix, paste0("run_", runNumber, "_samples"), fastqFileName)) %>%
    mutate(fastq_2 = "") %>%
    unite(sample, sample_columns, sep="_") %>%
    dplyr::select(sample, fastq_1, fastq_2, strandedness)

  if(check_files_flag){
    message("Checking whether fastq_1 files exist...")
    missing_files = unlist(lapply(sample_sheet_df$fastq_1, file.exists))

    # if there are any FALSE, return a summary
    if(sum(!missing_files)>0){
      message("Some files could not be located. See output")
      return(sample_sheet_df[['fastq_1']][!missing_files])
    }
  }

  sample_sheet_df

}
