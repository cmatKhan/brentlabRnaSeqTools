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
#' @seealso \code{\link{getCoverageOverRegion}}, \code{\link{strandedScanBamParam}}
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

  sbp = strandedScanBamParam(gr, strandedness, quality_threshold)

  pileup(bamfile_path,
         index = bamfile_index,
         scanBamParam = sbp,
         pileupParam = p_param)
}

#'
#' create coverage scanbamparam object
#'
#' @import Rsamtools
#'
#' @description helper function to create ScanBamParam object with appropriate strandedness information
#'
#' @param locus_granges a granges object for a given gene (or some other feature on only one strand)
#' @param strandedness one of c("reverse", "same", "unstranded")
#' @param quality_threshold quality threshold above which reads will be considered. 20l is default, which is
#'                          chosen b/c it is the default for HTSeq
#'
#' @export
strandedScanBamParam = function(locus_granges, strandedness, quality_threshold=20L){

  # set some information for the ScanBamParam object below. gene_strand extracts the +/- strand from the GRanges object
  # and the conditional sets the minus_strand_flag used in the ScanBamParam constructor. This determines what reads are
  # returned -- either from a given strand, or from both
  # TODO add support for forward stranded libraries
  gene_strand = as.character(unique(data.frame(locus_granges)$strand))

  # ensure the locus is entirely on the same strand, error out if not
  stopifnot(length(gene_strand) == 1)

  minus_strand_flag = switch (paste(strandedness, gene_strand, sep="_"),
    "reverse_+" = TRUE,
    "reverse_-" = FALSE,
    "same_+" = FALSE,
    "same_-" = TRUE,
    NA
  )

  ScanBamParam(which = locus_granges,
               mapqFilter = quality_threshold,
               flag = scanBamFlag(isMinusStrand=minus_strand_flag,
                                  isSecondaryAlignment=FALSE,
                                  isNotPassingQualityControls=FALSE,
                                  isSupplementaryAlignment=FALSE,
                                  isDuplicate=FALSE,
                                  isUnmappedQuery=FALSE,
               ))
}

#'
#' helper function to add .bai to bam path
#' @param bamfile_path path to bamfile
#'
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
#' @description a helper function to create a bampath from some metadata information. Also checks if index exists
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

#'
#' Count reads in a library like HTSeq
#'
#' @description Used to generate library counts for QC purposes, eg perturbed log2cpm
#'
#' @param bam_path path to a bam file. the index must also exist
#' @param annote_db txdb object. see \code{\link[GenomicFeatures]{makeTxDbFromGFF}}
#' @param feature_type one of c("exon", "cds"). Determines what feature to count over
#' @param strandedness strandedness of library determined by protocol. One of c("reverse", "same", "unstranded")
#'                     reverse means the counted reads will be on the opposite strand of the feature, same means they
#'                     will be on the same strand. unstranded counts reads over the feature regardless of strand
#' @param num_threads number of threads available for parallelization. Default to 1
#' @param single_end_reads_flag whether or not the library is single or paired. Note: this is not written to deal
#'                              with paired end reads currently, though you can try
#' @param mapq_filter quality filter on read alignment. Default to 20, same as HTSeq
#'
#' @return A RangedSummarizedExperiment. See \code{\link[GenomicAlignment]{summarizeOverlaps}}
#'
#' @examples
#' library(brentlabRnaSeqTools)
#' library(AnnotationDbi)
#' library(tidyverse)
#'
#' kn99_db = loadDb("data/kn99_db.sqlite")
#'
#' combined_df = getMetadata(database_info$kn99_host,
#'                           database_info$kn99_db_name,
#'                           Sys.getenv("db_username"),
#'                           Sys.getenv("password"))
#'
#' run_5102 = combined_df %>%
#'  filter(runNumber == 5102,
#'         experimentDesign == '90minuteInduction',
#'         genotype1 == "CNAG_00000")
#'
#' wt_fastq = run_5102[[1, 'fastqFileName']]
#'
#' bam_path = createBamPath(5102,
#'                          wt_fastq,
#'                          lts_align_expr_prefix = "/mnt/htcf_scratch/chasem/rnaseq_pipeline/align_count_results")
#'
#' cds_count = countLibrary(bam_path, kn99_db, 'cds', 'reverse', 8)
#'
#' cds_counts_df = as_tibble(assay(cds_counts), rownames = "gene_name")
#'
#' @seealso \url{http://achri.blogspot.com/2016/02/how-not-to-use-deseq2-for-illumina.html},
#'          \url{https://htseq.readthedocs.io/en/release_0.11.1/count.html},
#'          \code{\link[GenomicAlignment]{summarizeOverlaps}}
#'
#' @import BiocParallel
#' @importFrom GenomicFeatures exonsBy
#' @importFrom GenomicFeatures cdsBy
#' @importFrom GenomicRanges invertStrand
#' @importFrom Rsamtools BamFile
#' @importFrom GenomicAlignments summarizeOverlaps
#'
#' @export
countLibrary = function(bam_path, annote_db, feature_type, strandedness, num_threads = 1, single_end_reads_flag = TRUE,
                        mapq_filter = 20L)
{

  #TODO write test for exon, cds, and strandedess. put a (subsampled, and by-hand edited) bamfile in data, write tests
  # with that

  multiparam = MulticoreParam(num_threads, progressbar = TRUE, log = TRUE)

  genome_granges = switch (feature_type,
    "exon" = exonsBy(kn99_db, by="gene"),
    "cds" = cdsBy(kn99_db, by="gene")
  )

  bam_index = getBamIndexPath(bam_path)

  bamfile = BamFile(bam_path, bam_index)

  sbp = ScanBamParam(
    mapqFilter = 20L,
    which = genome_granges@unlistData,
    flag  = scanBamFlag(isUnmappedQuery = FALSE,
                        isSecondaryAlignment = FALSE,
                        isNotPassingQualityControls = FALSE),
  )

  if (strandedness == 'reverse') {
    summarizeOverlaps(genome_granges,
                      bamfile,
                      mode="Union",
                      singleEnd=single_end_reads_flag,
                      ignore.strand=FALSE,
                      # invert GRange strand for reverse strand. See above
                      preprocess.reads=GenomicRanges::invertStrand,
                      param = sbp,
                      BPPARAM = multiparam)
  } else {
    same_strand_flag = ifelse(strandedness=='same', FALSE, TRUE)
    summarizeOverlaps(genome_granges,
                      bamfile,
                      mode="Union",
                      singleEnd=single_end_reads_flag,
                      ignore.strand=same_strand_flag,
                      param = sbp,
                      BPPARAM = multiparam)
  }

}

#'
#' plot coverage over locus
#'
#' @description ggbio plot with transcripts track and coverage track
#'
#' @import ggbio
#' @import GenomicAlignments
#'
#' @param bamfile_path path to a bam file @seealso brentlabRnaSeqTools::createBamPath()
#' @param annote_db a GenomicFeatures TxDb object. Maybe one made from a gtf, eg txdb = makeTxDbFromGFF("data/liftoff_h99_to_kn99.gtf", format = "gtf")
#' @param gene_id a gene_id of interest -- must be in the gene names of the annote_db object
#' @param strandedness one of c("reverse", "unstranded"). NOTE: forward only strand NOT currently configured
#' @param quality_threshold quality threshold above which reads will be considered. 20l is default, which is
#'                          chosen b/c it is the default for HTSeq
#'
#' @export
plotCoverageOverLocus = function(bamfile_path, annote_db, gene_id, strandedness, quality_threshold=20L){

  locus_granges = featureGRanges(annote_db, gene_id, 'exon')

  sbp = strandedScanBamParam(locus_granges, strandedness, quality_threshold)

  bamfile_index = getBamIndexPath(bamfile_path)

  parsed_alignment = readGAlignments(bamfile_path,
                                     index = bamfile_index,
                                     use.names = TRUE,
                                     param=sbp)
  # create plots
  tx = ggbio::autoplot(annote_db, which = locus_granges)
  cov = ggbio::autoplot(parsed_alignment, stat="coverage")

  # construct layers
  tracks(Reads = cov, Transcripts=tx, title = gene_id)
}
