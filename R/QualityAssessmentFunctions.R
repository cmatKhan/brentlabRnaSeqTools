#'
#' get total protein coding count from count dataframe
#'
#' @description given a count dataframe with gene_ids as rownames and quantification in a column called raw_counts, return sum of protein coding genes
#'
#' @param counts a dataframe with gene_ids in the rownames and (at minimum) a quantification column called raw_counts
#' @param protein_coding_gene_ids a list of gene ids considered protein coding (must correspond with counts rownames)
#'
proteinCodingCount = function(counts, protein_coding_gene_ids){
  sum(counts[protein_coding_gene_ids,'raw_counts'])
}

#' filter for manual passes (overrides auto fail) and automatic passes (unless auto failed)
#'
#' @importFrom dplyr filter
#'
#' @param metadata dataframe from the database
#' @return a metadata dataframe with column names cast to upper
#'
#' @export
qualityAssessmentFilter = function(metadata){
  colnames(metadata) = toupper(colnames(metadata))
  passing_metadata = metadata %>%
    filter(MANUALAUDIT == FALSE | (is.na(MANUALAUDIT) & AUTOAUDIT == FALSE) )

  return(passing_metadata)
}

#' select fastqFileName, fastqFileNumber, and a pre-determined set of QC columns from a metadata df
#'
#' @importFrom dplyr select
#'
#' @note the column names for the metadata will be cast to upper and returned in upper
#' @note must include Interquartile range. Think about removing this -- user could merge with IQR df after selecting these cols
#'
#' @param metadata a metadata df with at least the columns listed in the select statement (see source code -- notably, must include interquartile range). Column names will be cast to uppper and returned as uppers
#'
#' @export
selectQaColumns = function(metadata){

  # cast column names to upper in case they have not yet been
  colnames(metadata) = toupper(colnames(metadata))

  # todo: set the columns to vectors, cast to symbols with rlang and offer argument for custom list
  selectMetadata = metadata %>%
    select(FASTQFILENAME,
           FASTQFILENUMBER,
           AUTOSTATUSDECOMP,
           MANUALAUDIT,
           INTERQUARTILERANGE,
           TOTALDEDUPLICATEDPERCENTAGE,
           PROTEINCODINGTOTAL,
           LIBRARYDATE,
           NATCOVERAGE,
           NATLOG2CPM,
           G418COVERAGE,
           G418LOG2CPM,
           NOTALIGNEDTOTALPERCENT,
           NOFEATUREPERCENT,
           INTERGENICCOVERAGE,
           LIBRARYSIZE,
           EFFECTIVEUNIQUE_ALIGNMENT)

  return(selectMetadata)
}

#'
#' create a sqlite database to hold the 'custom' qc data
#'
#' @importFrom stringr str_remove_all
#' @importFrom RSQLite dbConnect SQLite dbSendStatement dbClearResult dbDisconnect
#' @importFrom glue glue_sql
#'
#' @param database_dirpath path to containing directory of new qc database
#'
#' @return database path
#'
#' @export
createQCdatabase = function(database_dirpath){

  table_fields=list(fastqFileNumber='VARCHAR(150) PRIMARY KEY',
                    genotype1Coverage='INT(8)',
                    genotype1Log2cpm='FLOAT(6)',
                    genotype2Coverage='INT(8)',
                    genotype2Log2cpm='FLOAT(6)',
                    natCoverage='INT(8)',
                    natLog2cpm='FLOAT(6)',
                    g418Coverage='INT(8)',
                    g418Log2cpm='FLOAT(6)',
                    rRnaPercent='FLOAT(6)',
                    nctrRnaPercent='FLOAT(6)',
                    intergenicCoverage='FLOAT(6)')

  table_specs_vector = unlist(lapply(names(table_fields), function(x) paste(x, table_fields[[x]], sep=" ")))

  # create new database name
  db_name = paste0('qc_table_', format(Sys.time(), "%Y%m%d"), ".db")
  counter=1
  # prevent overwiting by ensuring the file does not exist
  while(file.exists(file.path(database_dirpath,db_name))){
    db_name = paste(db_name, counter, sep="_")
    counter = counter+1
  }

  db_connection = dbConnect(RSQLite::SQLite(), file.path(database_dirpath,db_name))

  sql_create_table_stmt = glue_sql("CREATE TABLE {tablename}({fields*})",
                                   tablename='qc_table',
                                   fields=table_specs_vector,
                                   .con = db_connection
                                   )
  # remove ' from entire statement
  sql_create_table_stmt = str_remove_all(sql_create_table_stmt, "\'")

  # create table
  tryCatch(
    expr = {
      res =dbSendStatement(db_connection, sql_create_table_stmt)
      message("creating qc database...")
      dbClearResult(res)
      # return the database path
      file.path(database_dirpath,db_name)
    },
    error = function(e){
      message("createQCdatabase() Error: could not create new qc_table")
      print(e)
    },
    warning = function(w){
      message("createQCdatabase() Warning: ")
      print(w)
    },
    finally = {
      dbDisconnect(db_connection)
    }
  )

}

#'
#' calcluate genotype and marker coverages for a given metadata df
#'
#' @importFrom tidyr pivot_longer drop_na fill
#' @importFrom stringr str_remove
#' @import dplyr
#'
#' @description this produces a sheet that can be used to fill the 'custom' qc table for the brentlab kn99 and yeast databases
#'
#' @param metadata_df a metadata sheet with at least the columns fastqFileName, genotype1, genotype2, runNumber, libraryProtocol
#' @param annote_db a connection to a local database -- see createQCdatabase()
#' @param bam_prefix the directory that contains the sequencing runs,, eg if locally mounted maybe "/mnt/htcf_lts/lts_align_expr"
#' @param bam_suffix the stuff appended to the end of the fastqFileName (minus the .fastq.gz). This might just be ".bam",
#'                   but could be something like "_sorted_aligned_reads_with_annote.bam"
#'
#' @note together, the bam_prefix/run_runNumber_samples/align/fastqFileName_bam_suffix form the path to the bamfile
#'
#' @return a sample sheet with the following columns: fastqFileName, fastqFileNumber, runNumber, column_param, locus, coverage
#' @export
qcGenotypeAndMarkerCoverage <- function(metadata_df, annote_db, bam_prefix, bam_suffix) {

  sample_sheet = metadata_df %>%
    dplyr::select(fastqFileName, genotype1,genotype2) %>%
    pivot_longer(!fastqFileName, names_to="column_param", values_to="locus") %>%
    mutate_all(list(~na_if(.,""))) %>%
    drop_na() %>%
    left_join(test_metadata %>%
                dplyr::select(fastqFileName, fastqFileNumber, runNumber, libraryProtocol) %>%
                mutate(libraryProtocol =ifelse(libraryProtocol == "E7420L", "reverse", "unstranded")) %>%
                dplyr::rename(strandedness = libraryProtocol)) %>%
    mutate(fastqFileName = str_remove(fastqFileName, ".fastq.gz"))

  sample_sheet = sample_sheet %>%
    mutate(fastqFileName = apply(sample_sheet, 1, function(x)
      createBamPath(
        getRunNumberLeadingZero(x[['runNumber']]),
        x[['fastqFileName']],
        bam_prefix,
        bam_suffix))) %>%
    dplyr::rename(bam_path = fastqFileName)

  # add markers rows to each group (grouping by fastqFileNumber)
  sample_sheet = sample_sheet %>%
    group_by(fastqFileNumber) %>%
    group_modify(~ add_row(.x, column_param="nat", locus="NAT"))%>%
    group_modify(~ add_row(.x, column_param="g418", locus="G418"))%>%
    arrange(fastqFileNumber)%>%
    fill(colnames(sample_sheet),.direction='down') %>%
    ungroup()

  # calculate coverage over loci, add to sample_sheet (this is returned)
  sample_sheet %>%
    mutate(coverages = apply(sample_sheet, 1, function(x) calculateCoverage(x[['bam_path']],
                                                                            annote_db,
                                                                            x[['locus']],
                                                                            x[['strandedness']])))
}

#'
#' calculate log2cpm for a given locus
#'
#' @param bam_path path to bam file. Note: the index with extension .bai also must exist
#' @param strandedness the strandedness of the library
#' @param locus_granges a granges object specifying the locus over which to count
#' @param lib_size the number of reads in the library
#'
#' @return either 0, if there are no counts, or the log2cpm of the counts over the locus. counting is done via the same
#'         method as default HTSeq
#'
#' @export
locusLog2Cpm = function(bam_path, strandedness, locus_granges, lib_size){

  locus_counts = countReadsInRanges(bam_path, locus_granges, strandedness)

  if(is.na(locus_counts)){
    stop("locusLog2Cpm Error(): counts for locus returned NA. Check input.")
  }

  ifelse(locus_counts == 0, 0, log2(locus_counts *(lib_size/1e6)) )
}

#'
#' Convenience function for edgeR's log2cpm
#' @description see \code{\link[edgeR]{cpm}}
#'
#' @param numeric_count_matrix output of \code{\link{countLibrary}} cast to data frame. See example
#' @param row_filter boolean vector length nrow(count_df) representing the rows to keep. Default is set to NA, which
#'                   returns all rows
#'
#' @return a numeric matrix in log2cpm
#'
#' @examples
#'
#' \dontrun{library(brentlabRnaSeqTools)
#' library(AnnotationDbi)
#' library(tidyverse)
#'
#' cds_count = countLibrary(bam_path, kn99_db, 'cds', 'reverse', 8)
#'
#' cds_counts_df = as_tibble(assay(cds_counts), rownames = "gene_name")
#'
#' log2cpm = edgerLog2cpm(as.matrix(dplyr::select(cds_counts_df, -gene_name)))
#' log2cpm_df = as_tibble(log2cpm) %>%
#'   mutate(gene_name = cds_counts_df$gene_name)}
#'
#' @seealso \code{\link[edgeR]{cpm}}
#'
#' @importFrom edgeR cpm DGEList
#'
#' @export
edgerLog2cpm = function(numeric_count_matrix, row_filter=NA){

  row_filter = ifelse(is.na(row_filter), rep(TRUE,nrow(numeric_count_matrix)), row_filter)

  dgelist = DGEList(numeric_count_matrix[row_filter,])
  # cpm returns the log2 of counts per million
  cpm(dgelist, log=TRUE)

}

#'
#' extract all metrics from broad rnaseqc package output
#'
#' @note this package is buggy (doesn't always get the metrics it says it will)
#'       also, all metrics are included in the nf-co/pipeline
#'
#'
#' @description extract, subset, reshape and plot metrics from the broad rnaseqc package
#'              \url{https://github.com/getzlab/rnaseqc}. Function plots Expression Profiling Efficiency (epe)
#'              vs rRNA percent and Estimated Library Complexity (elc) vs epe, both with marginal distributions.
#'
#' @note Estimated Library Complexity is a Picard metric and is similar to markDuplicates and totalDeduplicatedPercent
#'       \url{https://gatk.broadinstitute.org/hc/en-us/articles/360037591931-EstimateLibraryComplexity-Picard-}
#'
#' @importFrom dplyr select rename mutate bind_rows
#' @importFrom tidyr pivot_wider
#' @importFrom readr read_tsv
#' @import ggplot2 ggExtra
#'
#' @param rnaseqc_dir directory containing the rnaseqc output
#' @param bam_suffix the suffix to remove from the bam file sample names. Default to '.markdup.sorted.bam'
#'                   for nf-co/rnaseq_pipeline star_salmon output \url{https://nf-co.re/rnaseq}
#'
#' @return a list with items full_table, subset, rRna_vs_epe, rRna_vs_epe (see description)
#'
#' @export
parseBroadRnaseqcOutput = function(rnaseqc_dir, bam_suffix='.markdup.sorted.bam'){

  rnaseqc_summary_list = Sys.glob(file.path(rnaseqc_dir, "*metrics*"))

  message('reading in rnaseqc summaries...')
  rnaseqc_df_list = suppressMessages(lapply(rnaseqc_summary_list, read_tsv, col_names=c('metric', 'value')))

  message('re-shaping and merging rnaseqc summaries...')
  rnaseqc_df = bind_rows(lapply(rnaseqc_df_list, pivot_wider, names_from='metric', values_from='value'))

  message('subsetting rnaseqc summaries...')
  # todo: do rename in mutate, then select
  rnaseqc_df_subset = rnaseqc_df %>%
    select(Sample,
           `Expression Profiling Efficiency`,
           `Exonic Rate`,
           `Intronic Rate`,
           `Intergenic Rate`,
           `Intragenic Rate`,
           `rRNA Rate`,
           `Estimated Library Complexity`) %>%
    dplyr::rename(sample = Sample,
                  expression_profiling_efficiency = `Expression Profiling Efficiency`,
                  exonic_rate = `Exonic Rate`,
                  intronic_rate = `Intronic Rate`,
                  intergenic_rate = `Intergenic Rate`,
                  intragenic_rate = `Intragenic Rate`,
                  rRna_rate = `rRNA Rate`,
                  estimated_library_complexity = `Estimated Library Complexity`) %>%
    mutate(sample = str_remove(sample, bam_suffix),
           expression_profiling_efficiency = as.numeric(expression_profiling_efficiency),
           exonic_rate = as.numeric(exonic_rate),
           intronic_rate = as.numeric(intronic_rate),
           intergenic_rate = as.numeric(intergenic_rate),
           intragenic_rate = as.numeric(intragenic_rate),
           rRna_rate = as.numeric(rRna_rate),
           estimated_library_complexity = as.numeric(estimated_library_complexity))

  message('plotting..')
  # epe is expression profiling efficiency
  rRna_vs_epe = rnaseqc_df_subset %>%
    ggplot(aes(rRna_rate, expression_profiling_efficiency)) +
    geom_point() +
    scale_x_continuous() +
    scale_y_continuous()

  rRna_vs_epe_with_marginals = ggMarginal(rRna_vs_epe, type="histogram")

  # elc is estimated library complexity
  elc_vs_epe = rnaseqc_df_subset %>%
    ggplot(aes(estimated_library_complexity, expression_profiling_efficiency)) +
    geom_point() +
    scale_x_continuous() +
    scale_y_continuous()

  elc_vs_epe_with_marginals = ggMarginal(elc_vs_epe, type="histogram")

  results = list(rnaseqc_df = rnaseqc_df,
                 rnaseqc_df_subset = rnaseqc_df_subset,
                 rRna_vs_epe_with_marginals = rRna_vs_epe_with_marginals,
                 elc_vs_epe_with_marginals = elc_vs_epe_with_marginals)

  return(results)

}

#' decompose sums of powers of two to a list of the summed powers
#'
#' @description eg 18 = 2 + 16 decomposes to 1, 4
#'
#' @references yiming kang
#' \url{https://github.com/yiming-kang/rnaseq_pipe/blob/master/tools/utils.py}
#'
#' @param status an integer that represents the sum of powers of 2
#' @return a string of powers of 2 representing the bit, eg 1,4
#'
#' @export
decomposeStatus2Bit = function(status){

  #TODO can use negative numbers for "passing reasons"!!

  status_decomp = list()

  if(is.na(status)){
    status_decomp = "status_NA"
  } else if(status == 0){
    status_decomp = "passing_sample"
  } else if(status > 0){
    for(i in seq(floor(log2(status)),0)){
      if ((status -2**i) >= 0){
        status_decomp = append(status_decomp, i)
        status = status - 2**i
      }
    }
    status_decomp = paste(sort(unlist(status_decomp)), collapse=",")
  }
  status_decomp
}

#' calculate lower/upper fence as a default threshold on qc metrics
#'
#' @description calculate lower/upper inner and outer fences, and number of NAs
#'              in metric vector. inner fence defined as Q1/Q3 -/+ 1.5*IQR,
#'              outer fence is the same but 3*IQR
#'
#' @importFrom stats IQR quantile
#'
#' @param metric_vector numeric vector on which to calculate lower/upper fence
#' @return a list with the following slots: message, which stores a message
#'         regarding the NA count, lower_inner, lower_outer, and upper_inner,
#'         upper_outer which both store the fence values
#'
#' @export
outlierFence = function(metric_vector){
  # TODO test this
  metric_quantiles = quantile(metric_vector, c(.25, .75))

  metric_iqr = IQR(metric_vector, na.rm = TRUE)
  metric_iqr_inner_fence = 1.5*metric_iqr
  metric_iqr_outer_fence = 3*metric_iqr

  na_msg = paste0("Note: there are ",as.character(sum(is.na(metric_vector))),
                  " NAs in the argument vector.")

  list(message     = na_msg,
       lower_inner = metric_quantiles[[1]] - metric_iqr_inner_fence,
       lower_outer = metric_quantiles[[1]] - metric_iqr_outer_fence,
       upper_inner = metric_quantiles[[2]] + metric_iqr_inner_fence,
       upper_outer = metric_quantiles[[2]] + metric_iqr_outer_fence
  )
}

#'
#' Join QC tables (not qc_metrics) in qcExplorer input sqlite database
#' @description used to create the table displayed in the qcExplorer
#'
#' @importFrom purrr prepend reduce
#' @importFrom dplyr filter pull select
#' @importFrom RSQLite SQLite dbConnect dbReadTable dbDisconnect
#'
#' @param ui_selected_columns a character vector of the selected columns
#' @param id_field the common field between all QC tables
#' @param sample_df name of the table containing sample metadata
#' @param status_df name of the table containing the QC metrics (this table
#'                  does not have the id column)
#' @param qc_database_path path to the sqlite qc database
#'
#' @return the merged database based input
createMergedTable = function(ui_selected_columns, id_field, sample_df,
                             status_df, qc_database_path){
  selected_columns =
    if(is.null(ui_selected_columns)){
      id_field
    } else{
      c(id_field, ui_selected_columns)
    }

  tables_with_selected_columns = status_df %>%
    filter(metric %in% selected_columns) %>%
    pull(source) %>%
    unique()

  tables_with_selected_columns = c(tables_with_selected_columns)

  # subset the input table for only the selected fields from that table
  subsetDatabaseTables = function(db_table){

    relevant_fields =
      colnames(db_table)[colnames(db_table) %in% selected_columns]

    db_table %>% dplyr::select(relevant_fields)

  }

  # get relevant tables
  qc_database_connection = dbConnect(RSQLite::SQLite(), qc_database_path)
  selected_table_list = lapply(tables_with_selected_columns, function(x)
    dbReadTable(qc_database_connection, x,
                check.names = FALSE))
  dbDisconnect(qc_database_connection)

  # see https://stackoverflow.com/a/33177426/9708266
  selected_table_list = prepend(list(sample_df), selected_table_list, before = 1)

  selected_table_list = lapply(selected_table_list, subsetDatabaseTables)

  # cite: https://stackoverflow.com/a/34393416/9708266
  selected_table_list %>% reduce(left_join, by=id_field)
}

#'
#' create the metric status dataframe used to set thresholds and flag samples
#' in the QC explorer
#' @description columns should be source, which corresponds to the tablenames
#' in the QC explorer database, metric, which corresponds the field in the table
#' threshold, which may be set by the user and is initially NA,
#' comparative
#'
#' @importFrom dplyr bind_rows mutate
#' @importFrom stringr str_detect
#' @importFrom readr read_tsv
#'
#' @param user_metric_vector either an empty string if there are no user
#' defined metrics, or a character vector of metrics, eg c("epe", "perturbed_coverage")
#' @param multiqc_output path to the multiqc output data directory of
#' nf-co/rnaseq_pipeline
#'
#' @return a data frame with columns source, metric, alias, threshold, status
#'
#' @export
createStatusDf = function(user_metric_vector, multiqc_output){

  if(user_metric_vector != ""){
    # note: alias can be a simplified name of the metric vector. Blank here
    # just name metric something reasonable. Alias is intended for multiqc output
    custom_status_df = tibble(source = rep("user", length(user_metric_vector)),
                              metric = user_metric_vector,
                              alias = rep("", length(user_metric_vector)),
                              threshold = NA,
                              comparative = "",
                              status = seq(1, length(user_metric_vector)))
  }

  qc_data_paths = getNfCoMultiqcDataPaths(multiqc_output)

  qc_data_list = lapply(qc_data_paths, function(x)
    suppressMessages(read_tsv(x, name_repair = "minimal")))
  names(qc_data_list) = names(qc_data_paths)

  qcToStatusReformater = function(df, source_name){
    metric_vector = colnames(df)[str_detect(colnames(df), "Sample", negate = TRUE)]
    tibble(source = rep(source_name, length(metric_vector)),
           metric = metric_vector,
           alias = "",
           comparative = "",
           threshold = NA,
           status = NA)
  }

  qc_status_df = bind_rows(lapply(names(qc_data_list), function(x)
    qcToStatusReformater(qc_data_list[[x]], x) ))

  # combine user, multiqc to form the one status_df to rule them all
  last_user_status = if(user_metric_vector != ""){
    nrow(custom_status_df)
  } else{
    0
  }

  qc_status_df = qc_status_df %>%
    mutate(status = seq(1, nrow(qc_status_df))+last_user_status)

  if(user_metric_vector != ""){
    bind_rows(custom_status_df, qc_status_df)
  } else{
    qc_status_df
  }
}

#'
#' create the sqlite database as input to qcExplorer
#' @description sample_df (metadata on the samples) and any QC tables must
#' share a common ID column. qc_metrics_df contains the column names of the QC
#' metrics with columns threshold, status and comparative which store information
#' on how to flag samples which violate the threshold. Writes database to
#' specified file
#'
#' @importFrom dplyr mutate left_join select arrange
#' @importFrom readr read_tsv
#' @importFrom RSQLite dbConnect SQLite dbWriteTable dbDisconnect
#'
#' @param sample_df a tibble with an id column which corresponds to the id column
#' of the QC metrics tables and columns containing metadata on samples
#' @param status_df a tibble with columns metric, alias, threshold, status,
#' comparative for use in flagging samples
#' @param qc_df_list list of qc table dataframes -- already read in.
#' see articles in documentation
#' reports from the nf-co/rnaseq pipeline
#' @param database_output output path (full, with .sqlite extension) for the
#' sqlite database
#' @param id_column common id column between sample_df and QC tables
#'
#' @export
createDatabase = function(sample_df, status_df, qc_df_list, database_output,
                          id_column){


  # structure: Sample, id, audit, status, statusDecomp where Sample is from
  # the nextflow samplesheet, sample_id is a unique identifier for the same sample
  # which will allow join back to sample metadata in database, audit is a boolean,
  # status is the sum of bit statuses and statusDecomp is the list of powers of 2
  # to which the bit status decomposes
  #
  # there needs to be an id column in sample_df (eg make fastqFileNumber into id)

  sample_df = sample_df %>%
    mutate(audit_flag=NA, status=as.numeric(NA), statusDecomp="")

  add_id_col = function(qc_table){
    qc_colnames = c(id_column, colnames(qc_table))
    # qc_colnames = qc_colnames[!qc_colnames == "Sample"]

    qc_table %>%
      left_join(sample_df, by=id_column) %>%
      dplyr::select(qc_colnames) %>%
      arrange(!!rlang::sym(id_column))
  }

  qc_data_list_with_id = lapply(qc_df_list, add_id_col)
  sample_df = arrange(sample_df, !!rlang::sym(id_column))

  con <- dbConnect(RSQLite::SQLite(), database_output)

  dbWriteTable(con, "sample", sample_df)
  dbWriteTable(con, "qc_metrics", status_df)
  lapply(names(qc_data_list_with_id), function(x)
    dbWriteTable(con, x, qc_data_list_with_id[[x]]))

  dbDisconnect(con)
}

#'
#' get relevant QC metrics tables from nf-co/rnaseq multiqc output
#'
#' @importFrom stringr str_detect str_remove
#'
#' @param multiqc_output path to multiqc_output
#'
#' @return list of full paths
#'
#' @export
getNfCoMultiqcDataPaths = function(multiqc_output){

  qc_data_paths = list.files(str_remove(multiqc_output, "/$"),
                             full.names = TRUE)

  qc_data_paths = qc_data_paths[
    str_detect(basename(qc_data_paths), "^multiqc|.+featurecounts_biotype_.+|^tin|^qualimap")]

  qc_data_paths = qc_data_paths[
    str_detect(qc_data_paths,
               "warning|fastqc|sources|idxstats|picard|json$|log$",
               negate = TRUE)]

  names(qc_data_paths) = tools::file_path_sans_ext(
    basename(qc_data_paths))

  qc_data_paths
}

#'
#' add new_qc_database tables to old_qc_database_tables
#'
#' @importFrom RSQLite dbConnect dbDisconnect dbWriteTable dbReadTable
#'
#' @param old_qc_database_path connection to existing qc database -- the
#' one to be added to
#' @param new_qc_database_path connection to the new qc database
#' @param qc_tablenames the tables to be added
#'
#' @export
addTablesToQcDatabase = function(old_qc_database_path,
                                 new_qc_database_path,
                                 qc_tablenames){

  old_qc_database_connection = dbConnect(RSQLite::SQLite(), old_qc_database_path)

  new_qc_database_connection = dbConnect(RSQLite::SQLite(), new_qc_database_path)

  # TODO make sure names of new_qc_tables are in the old_qc_database_donnection
  #      tables

  lapply(qc_tablenames, function(x)
    dbWriteTable(old_qc_database_connection, x,
                 dbReadTable(new_qc_database_connection,x,check.names = FALSE),
                 append = TRUE))

  dbDisconnect(old_qc_database_connection)
  dbDisconnect(new_qc_database_connection)

}

#'
#' extract all metrics from the output of tin.py run on a batch of files
#'
#' @description presumably tin.py was run on a batch of samples, the output of
#'              which is in some output dir. This function compiles all the
#'              summary sheets and compiles them into a single dataframe
#' @seealso \url{http://rseqc.sourceforge.net/#tin-py}
#'
#' @note Estimated Library Complexity is a Picard metric and is similar to markDuplicates and totalDeduplicatedPercent
#'       \url{https://gatk.broadinstitute.org/hc/en-us/articles/360037591931-EstimateLibraryComplexity-Picard-}
#'
#' @importFrom dplyr select mutate select bind_rows
#' @importFrom stringr str_remove
#' @importFrom readr read_tsv
#' @import ggplot2 ggExtra
#'
#' @param tin_output_dir directory containing the rnaseqc output
#' @param bam_suffix the suffix to remove from the bam file sample names. Default to '.markdup.sorted.bam'
#'                   for nf-co/rnaseq_pipeline star_salmon output \url{https://nf-co.re/rnaseq}
#'
#' @return a list with items full_table, subset, rRna_vs_epe, rRna_vs_epe (see description)
#'
#' @export
compileTinOutput = function(tin_output_dir, bam_suffix='.markdup.sorted.bam'){

  message('reading in tin.py summaries...')
  tin_summary_list = Sys.glob(file.path(tin_output_dir, "*summary*"))
  tin_df_list = suppressMessages(lapply(tin_summary_list, read_tsv))

  message('merging tin.py summaries...')
  tin_df = bind_rows(tin_df_list)
  # todo: do rename in mutate, then select
  tin_df %>%
    mutate(Sample = str_remove(Bam_file, bam_suffix)) %>%
    select(-Bam_file)
}

#'
#' parse qualimap summary output
#' @description extract metrics from the file rnaseq_qc_results.txt in nf-co/rnaseq qualimap
#' output
#'
#' @importFrom readr read_lines
#' @importFrom tidyr separate
#' @importFrom purrr map
#' @importFrom plyr ldply
#'
#' @param sample_name name of the sample which will go into the Sample column.
#' See example
#' @param qualimap_results path to rnaseq_qc_results.txt for a given sample
#'
#' @return a dataframe with columns sample, exonic, intergenic, overlapping_exon,
#' count, percent
#'
#' @examples
#' \dontrun{
#'    pool_plate2_qualimaps_files = Sys.glob("path/to/qualimap/*/rnaseq_qc_results.txt")
#'    names(pool_plate2_qualimaps_files) =
#'        basename(dirname(pool_plate2_qualimaps_files))
#'    pool_plate2_qualimap_df = bind_rows(lapply(names(pool_plate2_qualimaps_files),
#'        function(x)
#'            compileQualimapSummaryMetrics(x, pool_plate2_qualimaps_files[[x]]))) %>%
#'            arrange(Sample)
#' }
#'
#' @export
compileQualimapSummaryMetrics = function(sample_name, qualimap_results){

  file = read_lines(qualimap_results)

  file[str_detect(file,
                  "\\s+exonic|\\s+intronic|\\s+intergenic|\\overlapping")] %>%
    str_split("=") %>%
    map(str_trim, side="both") %>%
    ldply() %>%
    separate(V2, c('count', 'percent'), sep="\\s") %>%
    mutate(count = str_remove_all(count, ","),
           percent = str_remove_all(percent, "\\(|\\)|%")) %>%
    pivot_wider(names_from=V1, values_from = c(count, percent)) %>%
    mutate(Sample = sample_name) %>%
    dplyr::select(c(Sample, colnames(.)))

}
