# TODO replace coverage with coverageByTranscripts + readGAlignments over CDS
# TODO add count with Rsubreads
# TODO make the QC pipeline using the tools above

#'
#' create coverage scanbamparam object
#'
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#'
#' @description helper function to create ScanBamParam object with appropriate strandedness information
#'
#' @param locus_granges a granges object for a given gene (or some other feature on only one strand)
#' @param strandedness one of c("reverse", "same", "unstranded")
#' @param quality_threshold quality threshold above which reads will be considered. 20l is default, which is
#'                          chosen b/c it is the default for HTSeq
#'
#' @return a ScanBamParam object with certain configured options, as well as some reasonable defaults, to filter
#'         a bam file for reads of interest based on the strandedness (protocol) of the library prep.
#'
#' @seealso  \code{\link[Rsamtools]{ScanBamParam}}
#'
#' @export
strandedScanBamParam = function(locus_granges, strandedness, quality_threshold=20L){

  # set some information for the ScanBamParam object below. gene_strand extracts the +/- strand from the GRanges object
  # and the conditional sets the minus_strand_flag used in the ScanBamParam constructor. This determines what reads are
  # returned -- either from a given strand, or from both
  # TODO add support for forward stranded libraries
  gene_strand = as.character(unique(data.frame(locus_granges)$strand))

  # ensure the locus is entirely on the same strand, error out if not
  if(!length(gene_strand) == 1){
    message("The granges object contains features on both strands.
            Separate the granges object into sets so that any set is on the same strand and try again")
  }

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
#'
#' @param bamfile_path path to bamfile
#'
#' @return bamfile_path concat with '.bai'
#'
#' @export
getBamIndexPath = function(bamfile_path){

  bamfile_index = paste0(bamfile_path, ".bai")

  error_msg = paste0("The path: ", bamfile_index, " DNE.\n
  The bam file must be indexed. Use samtools index to do this.\n
  After doing so, there should be a *.bai file in the same dir as the bam.")

  if (!file.exists(bamfile_index)){
    stop(error_msg)
  }

  return(bamfile_index)
}

#'
#' create a bam path
#'
#' @importFrom stringr str_remove
#'
#' @description a helper function to create a bampath from some metadata information. Also checks if index exists
#'
#' @param run_number the run_number (mind the leading zeros for old runs) of the run
#' @param fastq_filename the fastq filename, preferrably without the extension or any leading path info. However, an effort has been made to deal with full paths and extensions
#' @param align_expr_prefix the location of the run directories. Eg, if you are mounted and on your local computer, it might be something like "/mnt/htcf_lts/lts_align_expr"
#' @param bam_suffix the common bam suffix for all bam files stored in /lts. Default is "_sorted_aligned_reads_with_annote.bam"
#' @param test boolean, default FALSE. Set to TRUE if testing this function
#'
#' @return a verified filepath to the bam file
#' @export
createBamPath = function(run_number, fastq_filename, align_expr_prefix,
                         bam_suffix="_sorted_aligned_reads_with_annote.bam", test=FALSE){

  fastqfile_basename = str_remove(basename(fastq_filename), '.fastq.gz')

  bam_path = file.path(align_expr_prefix, paste0("run_", as.character(run_number),"_samples/align"),
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
#' @importFrom GenomicFeatures cdsBy exonsBy
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
#' plot coverage over locus
#'
#' @description ggbio plot with transcripts track and coverage track
#'
#' @importFrom ggbio autoplot tracks
#' @importFrom GenomicAlignments readGAlignments
#'
#' @param bamfile_path path to a bam file
#' @param annote_db a GenomicFeatures TxDb object. Maybe one made from a gtf, eg txdb = makeTxDbFromGFF("data/liftoff_h99_to_kn99.gtf", format = "gtf")
#' @param gene_id a gene_id of interest -- must be in the gene names of the annote_db object
#' @param strandedness one of c("reverse", "unstranded"). NOTE: forward only strand NOT currently configured
#' @param quality_threshold quality threshold above which reads will be considered. 20l is default, which is
#'                          chosen b/c it is the default for HTSeq
#'
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
  cov = ggbio::autoplot(parsed_alignment, stat="coverage", geom="alignment")

  # construct layers
  tracks(Reads = cov, Transcripts=tx, title = gene_id)
}
