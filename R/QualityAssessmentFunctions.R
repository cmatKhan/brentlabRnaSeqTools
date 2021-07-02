#'
#' get total protein coding count from count dataframe
#'
#' @import dplyr
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
#' @import dplyr
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
           INTERQUARTILE_RANGE,
           TOTAL_DEDUPLICATED_PERCENTAGE,
           PROTEIN_CODING_TOTAL,
           LIBRARYDATE,
           NAT_COVERAGE,
           NAT_LOG2CPM,
           G418_COVERAGE,
           G418_LOG2CPM,
           NOT_ALIGNED_TOTAL_PERCENT,
           NO_FEATURE_PERCENT,
           INTERGENIC_COVERAGE,
           LIBRARY_SIZE,
           EFFECTIVE_UNIQUE_ALIGNMENT)

  return(selectMetadata)
}

#'
#' create a sqlite database to hold the 'custom' qc data
#'
#' @import dplyr
#' @import RSQLite
#' @import glue
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
#' library(brentlabRnaSeqTools)
#' library(AnnotationDbi)
#' library(tidyverse)
#'
#' cds_count = countLibrary(bam_path, kn99_db, 'cds', 'reverse', 8)
#'
#' cds_counts_df = as_tibble(assay(cds_counts), rownames = "gene_name")
#'
#' log2cpm = edgerLog2cpm(as.matrix(dplyr::select(cds_counts_df, -gene_name)))
#' log2cpm_df = as_tibble(log2cpm) %>%
#'   mutate(gene_name = cds_counts_df$gene_name)
#'
#' @seealso \code{\link[edgeR]{cpm}}
#'
#' @importFrom edgeR cpm
#'
#' @export
edgerLog2cpm = function(numeric_count_matrix, row_filter=NA){

  row_filter = ifelse(is.na(row_filter), rep(TRUE,nrow(numeric_count_matrix)), row_filter)

  dgelist = DGEList(numeric_count_matrix[row_filter,])
  # cpm returns the log2 of counts per million
  cpm(dgelist, log=TRUE)

}

