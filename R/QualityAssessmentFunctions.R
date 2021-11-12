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
#' @description this produces a sheet that can be used to fill the 'custom'
#' qc table for the brentlab kn99 and yeast databases
#'
#' @param metadata_df a metadata sheet with at least the columns fastqFileName,
#' genotype1, genotype2, runNumber, libraryProtocol
#' @param annote_db a connection to a local database -- see createQCdatabase()
#' @param bam_prefix the directory that contains the sequencing runs,
#' eg if locally mounted maybe "/mnt/htcf_lts/lts_align_expr"
#' @param bam_suffix the stuff appended to the end of the fastqFileName (minus the .fastq.gz). This might just be ".bam",
#'                   but could be something like "_sorted_aligned_reads_with_annote.bam"
#'
#' @note together, the bam_prefix/run_runNumber_samples/align/fastqFileName_bam_suffix form the path to the bamfile
#'
#' @return a sample sheet with the following columns: fastqFileName, fastqFileNumber, runNumber, column_param, locus, coverage
#'
# qcGenotypeAndMarkerCoverage <- function(metadata_df, annote_db, bam_prefix, bam_suffix) {
#
#   sample_sheet = metadata_df %>%
#     dplyr::select(fastqFileName, genotype1,genotype2) %>%
#     pivot_longer(!fastqFileName, names_to="column_param", values_to="locus") %>%
#     mutate_all(list(~na_if(.,""))) %>%
#     drop_na() %>%
#     left_join(test_metadata %>%
#                 dplyr::select(fastqFileName, fastqFileNumber, runNumber, libraryProtocol) %>%
#                 mutate(libraryProtocol =ifelse(libraryProtocol == "E7420L", "reverse", "unstranded")) %>%
#                 dplyr::rename(strandedness = libraryProtocol)) %>%
#     mutate(fastqFileName = str_remove(fastqFileName, ".fastq.gz"))
#
#   sample_sheet = sample_sheet %>%
#     mutate(fastqFileName = apply(sample_sheet, 1, function(x)
#       createBamPath(
#         getRunNumberLeadingZero(x[['runNumber']]),
#         x[['fastqFileName']],
#         bam_prefix,
#         bam_suffix))) %>%
#     dplyr::rename(bam_path = fastqFileName)
#
#   # add markers rows to each group (grouping by fastqFileNumber)
#   sample_sheet = sample_sheet %>%
#     group_by(fastqFileNumber) %>%
#     group_modify(~ add_row(.x, column_param="nat", locus="NAT"))%>%
#     group_modify(~ add_row(.x, column_param="g418", locus="G418"))%>%
#     arrange(fastqFileNumber)%>%
#     fill(colnames(sample_sheet),.direction='down') %>%
#     ungroup()
#
#   # calculate coverage over loci, add to sample_sheet (this is returned)
#   sample_sheet %>%
#     mutate(coverages = apply(sample_sheet, 1, function(x) calculateCoverage(x[['bam_path']],
#                                                                             annote_db,
#                                                                             x[['locus']],
#                                                                             x[['strandedness']])))
# }

