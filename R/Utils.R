#'
#' test if argument is numeric
#' @description copied directly from the limma codebase
#' @param x any R object
#' @details copied from the limma docs: This function is used to check the validity of arguments for numeric functions. It is an attempt to emulate the behavior of internal generic math functions. IsNumeric differs from is.numeric in that data.frames with all columns numeric are accepted as numeric.
#'
#' @export
isNumeric = function (x) {
  is.numeric(x) || (is.data.frame(x) && length(x) > 0 && all(unlist(lapply(x, is.numeric))))
}

#'
#' correct run number to add leading zero where approprirate
#' @param run_number a run number, most likely from the metadata runNumber field
#'
#' @export
getRunNumberLeadingZero = function(run_number){
  run_number = as.character(run_number)
  if(!is.null(run_numbers_with_leading_zero[[run_number]])){
    run_numbers_with_leading_zero[[run_number]]
  } else{
    run_number
  }
}

#'
#' test bam path
#'
#' @export
testBamPath = function(metadata_df){
  for (i in seq(1,nrow(metadata_df))){
    run_number = getRunNumberLeadingZero(metadata_df[[1,'runNumber']])
    run_dir = paste0("run_", run_number, "_samples")

    sample_name = str_remove(metadata_df[[1,'fastqFileName']], ".fastq.gz")
    bam_name = paste0(sample_name, Sys.getenv("BAM_SUFFIX"))

    bam_path = file.path(Sys.getenv("LTS_ALIGN_EXPR_PREFIX"), run_dir, "align", bam_name)
    message(paste0("checking: ", bam_path))
    stopifnot(file.exists(bam_path))
  }
}

#'
#' read in columnar data
#'
#' @import dplyr
#'
#' @description given a csv, tsv or excel sheet, use the right function to read in the data
#'
#' @param path path to a csv, tsv or xlsx
#'
#' @export
readInData = function(path){

  switch(tools::file_ext(path),
         "csv" = read_csv(path),
         "tsv" = read_tsv(path),
         "xlsx" = readxl::read_excel(path),
         message("File extension not recognized. Must be one of {csv, tsv, xlsx}"))
}
