#'
#' calculate coverage over a given locus
#'
#' @import Rsamtools
#' @import GenomicFeatures
#' @import GenomicRanges
#'
#' @param bamfile_path path to a bam file @seealso brentlabRnaSeqTools::createBamPath()
#' @param annote_db a GenomicFeatures TxDb object. Maybe one made from a gtf, eg txdb = makeTxDbFromGFF("data/liftoff_h99_to_kn99.gtf", format = "gtf")
#' @param gene_id a gene_id of interest -- must be in the gene names of the annote_db object
#' @param strandedness one of c("unstranded", "reverse") indicating the strandedness of the library. note: forward not currently supported
#' @param ... additional arguments to getCoverageOverRegion()
#'
#' @return percent coverage of feature with reads above a given quality threshold and coverage depth threshold (see getCoverageOverRegion())
#'
#' @export
calculateCoverage = function(bamfile_path, annote_db, gene_id, strandedness, ...){

  coverage_df = getCoverageOverRegion(bamfile_path, annote_db, gene_id, strandedness, quality_threshold=20L,
                                      coverage_threshold=0, lts_align_expr_prefix=Sys.getenv("LTS_ALIGN_EXPR_PREFIX"),
                                      bamfile_suffix=Sys.getenv("BAM_SUFFIX"))
  # TODO this is repeated -- how to only calculate this once in getCoverageOverRegion?
  gr = featureGRanges(annote_db, gene_id, 'cds')

  length(unique(coverage_df$pos))/sum(width(gr))
}

#'
#' create a dataframe of coverage by nucleotide over a given locus
#'
#' @import Rsamtools
#' @import GenomicFeatures
#' @import GenomicRanges
#'
#' @param bamfile_path path to a bam file @seealso brentlabRnaSeqTools::createBamPath()
#' @param annote_db a GenomicFeatures TxDb object. Maybe one made from a gtf, eg txdb = makeTxDbFromGFF("data/liftoff_h99_to_kn99.gtf", format = "gtf")
#' @param gene_id a gene_id of interest -- must be in the gene names of the annote_db object
#' @param strandedness one of c("reverse", "unstranded"). NOTE: forward only strand NOT currently configured
#' @param quality_threshold quality threshold above which reads will be considered. 20l is default, which is
#'                          chosen b/c it is the default for HTSeq
#' @param coverage_threshold minimum read count above which to consider reads. Default is 0
#' @param lts_align_expr_prefix = path to the directory which stores the run_12345_samples run directories.
#'                                For example, /lts/mblab/Crypto/rnaseq_data/lts_align_expr.
#'                                By default, this looks in your .Renviron for a key LTS_ALIGN_EXPR_PREFIX
#' @param bamfile_suffix = whatever is appended after the fastqFileName (no extension).
#'                         Currently, this is "_sorted_aligned_reads_with_annote.bam". By default, this looks in your
#'                         .Renviron for a key BAM_SUFFIX
#' @references GenomicRanges, Rsamtools
#' @export
getCoverageOverRegion = function(bamfile_path, annote_db, gene_id, strandedness, quality_threshold=20L,
                                 coverage_threshold=0, lts_align_expr_prefix=Sys.getenv("LTS_ALIGN_EXPR_PREFIX"),
                                 bamfile_suffix=Sys.getenv("BAM_SUFFIX")){

  bamfile_index = getBamIndexPath(bamfile_path)

  gr = featureGRanges(annote_db, gene_id, 'cds')

  p_param <- PileupParam(min_mapq = quality_threshold,
                         min_nucleotide_depth = coverage_threshold,
                         distinguish_strands=TRUE,
                         distinguish_nucleotides=FALSE)

  # set some information for the ScanBamParam object below. gene_strand extracts the +/- strand from the GRanges object
  # and the conditional sets the minus_strand_flag used in the ScanBamParam constructor. This determines what reads are
  # returned -- either from a given strand, or from both
  # TODO add support for forward stranded libraries
  gene_strand = as.character(unique(data.frame(gr)$strand))

  minus_strand_flag = NA
  if (strandedness=='reverse' & gene_strand == '+'){
    minus_strand_flag = TRUE
  } else if(strandedness=='reverse' & gene_strand == '-'){
    minus_strand_flag = FALSE
  }

  sbp = ScanBamParam(which=gr,
                     mapqFilter=quality_threshold,
                     flag=scanBamFlag(isMinusStrand=minus_strand_flag,
                                      isSecondaryAlignment=FALSE,
                                      isNotPassingQualityControls=FALSE,
                                      isSupplementaryAlignment=FALSE,
                                      isDuplicate=FALSE,
                                      isUnmappedQuery=FALSE,
                     ))

  pileup(bamfile_path,
         index = bamfile_index,
         scanBamParam = sbp,
         pileupParam = p_param)
}

getBamIndexPath = function(bamfile_path){
  bamfile_index = paste0(bamfile_path, ".bai")

  if (file.exists(bamfile_index)==FALSE){
    stop("The bam file must be indexed (use samtools index to do this. After doing so, there should be a *.bai file in the same dir as the bam)")
  }

  return(bamfile_index)
}

#'
#' create a bam path
#'
#' @import stringr
#'
#' @description a helper function to creat a bampath from some metadata information. Also checks if index exists
#'
#' @param run_number the run_number (mind the leading zeros for old runs) of the run
#' @param fastq_filename the fastq filename, preferrably without the extension or any leading path info. However, an effort has been made to deal with full paths and extensions
#' @param lts_align_expr_prefix the location of the run directories. Eg, if you are mounted and on your local computer, it might be something like "/mnt/htcf_lts/lts_align_expr"
#' @param bam_suffix the common bam suffix for all bam files stored in /lts. Eg, it might be something like "_sorted_aligned_reads_with_annote.bam"
#' @param test boolean, default FALSE. Set to TRUE if testing this function
#'
#' @return a verified filepath to the bam file
#' @export
createBamPath = function(run_number, fastq_filename, lts_align_expr_prefix, bam_suffix="_sorted_aligned_reads_with_annote.bam", test=FALSE){

  fastqfile_basename = str_remove(basename(fastq_filename), '.fastq.gz')

  bam_path = file.path(lts_align_expr_prefix, paste0("run_", as.character(run_number),"_samples/align"),
                       paste0(fastqfile_basename, bam_suffix))

  bam_index = getBamIndexPath(bam_path)

  # if test is true, just return the bam_path
  ifelse(test,
         bam_path,
         # if the path is valid, return the bam_path. else, stop with error
         ifelse(file.exists(bam_path),
                bam_path,
                stop(paste0("The following path is invalid: ", bam_path))))
}

#' Given a GenomicFeatures annotation_db and a gene_id, extract an GRanges object of the cds
#'
#' @import GenomicFeatures
#' @import GenomicRanges
#'
#' @param annotation_db a GenomicFeatures db. You can either get this from the bioconductor resources, or create your own with a gtf
#' @param gene_id the ID of a gene in the db. Eg, for cryptococcus CKF44_05222
#' @param feature one of c("cds", "exon"), determins which feature to extract from the annotations
#' @references GenomicRanges::GRanges, GenomicFeatures
#' @return an IRanges object of the given gene's exons
#'
#' @export
featureGRanges = function(annotation_db, gene_id, feature){

  if(!feature %in% c('cds', 'exon')){
    stop("Only 'cds' and 'exon' are supported currently")
  }

  regions = tryCatch(
    expr = {
      # note, [[1]] extracts just the GRanges object
      switch(feature,
             'cds' = cdsBy(annotation_db, by="gene")[gene_id][[1]],
             'exon' = exonsBy(annotation_db, by="gene")[gene_id][[1]]
            )
    },
    error = function(e){
      message(paste0('featureGranges() Error: cannot create GRanges for gene_id: ', gene_id))
      print(e)
    },
    warning = function(w){
      message("featureGRanges() warning: ")
      print(w)
    },
    finally = {
      # none
    }
  )

  return(regions)
}

#'
#' from the metadata libraryProtocol column, determine the library strandedness.
#' @description Currently set up for cryptococcus. E7420L returns reverse, SolexaPrep returns unstranded. default return is unstranded
#'
#' @note: default is unstranded
#' @param library_protocol the library protocol of the sample (determines strandedness of the library)
#'
#' @return the strandedness of the library based on the value in the libraryProtocol column, or 'unstranded' by default
#'
#' @export
determineLibraryStrandedness = function(library_protocol){

  # stop on this condition for now. the switch statement actually defaults to 'unstranded'
  # consider allowing default to stand, but issue a warning?
  stopifnot(library_protocol %in% c('E7420L', 'SolexaPrep'))

  switch(library_protocol,
         "E7420L" = "reverse",
         "SolexaPrep" = "unstranded",
         "unstranded")

}

#'
#' given a bam file path, GRanges object, and strandedness of library, return total counts
#'
#' @import Rsamtools
#' @import GenomicRanges
#'
#' @param bamfile_path path to bam file
#' @param granges_of_interest a GRanges object
#' @param strandedness one of c("reverse", "unstranded")
#'
#' @export
countReadsInRanges = function(bamfile_path, granges_of_interest, strandedness){

  # right now, this is only configured to handle reverse and unstranded libraries
  stopifnot(strandedness %in% c('reverse', 'unstranded'))

  bamfile_index = getBamIndexPath(bamfile_path)

  # see ??scanBamParam isMinusStrand
  ignore_strand_flag = ifelse(strandedness=='unstranded', TRUE, FALSE)

  # for a stranded library, only count reads aligning to the appropriate strand
  stranded_bam_params = ScanBamParam(which=granges_of_interest,
                                     flag=scanBamFlag(
                                       isSecondaryAlignment=FALSE,
                                       isNotPassingQualityControls=FALSE,
                                       isSupplementaryAlignment=FALSE,
                                       isDuplicate=FALSE))

  ranges_hits = summarizeOverlaps(granges_of_interest,
                                  bamfile,
                                  mode="Union",
                                  param = stranded_bam_params,
                                  ignore.strand = ignore_strand_flag)
  return(ranges_hits)
}


