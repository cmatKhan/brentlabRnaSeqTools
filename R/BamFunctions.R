#'
#' @param bamfile_path path to a bam file @seealso brentlabRnaSeqTools::createBamPath()
#' @param genomic_ranges a GRanges object with the range(s) of interest
#' @param strandedness strandedness of the library (this determines how coverage is calculated)
#' @param minus_strand_flag c(TRUE, FALSE, -1). TRUE if the correctly mapped reads are on the minus strand,
#' FALSE if on the positive, -1 if unstranded (both positive and negative strand aligned reads will be returned for an unstranded library)
#'
#' @export
calculateExonicCoverage = function(bamfile_path, genomic_ranges, minus_strand_flag, strandedness, quality_threshold=20L, coverage_threshold=0, lts_align_expr_prefix=Sys.getenv("LTS_ALIGN_EXPR_PREFIX", bamfile_suffix='"_sorted_aligned_reads_with_annote.bam"')){

  bamfile_index = paste0(bamfile_path, ".bai")

  if (file.exists(bamfile_index)==FALSE){
    stop("The bam file must be indexed (use samtools index to do this. After doing so, there should be a *.bai file in the same dir as the bam)")
  }

  gr = GRanges(seqnames = as.vector(genomic_ranges@seqnames@values), ranges = genomic_ranges@ranges)

  p_param <- PileupParam(min_base_quality = quality_threshold)

  # create the scanBamParam object. If the minus_strand_flag is -1, the library is unstranded. Else, extract only the positive or negative strand reads
  sbp = ifelse(minus_strand_flag==-1,
               ScanBamParam(which=gr),
               ScanBamParam(which=gr, flag=scanBamFlag(isMinusStrand=minus_strand_flag)))

  exonic_coverage_df = pileup(bamfile_path,
                              index = bamfile_index,
                              scanBamParam = sbp,
                              pileupParam = p_param)

  # return coverage, calculated as the sum of the bp covered above a given threshold
  ifelse(minus_strand_flag=='unstranded',
         sum((res %>% distinct(pos, .keep_all=TRUE))$count > coverage_threshold)/sum(width(gr),
                                                                                     sum(exonic_coverage_df$count > coverage_threshold)/sum(width(gr))))

}

#'
#'
#' @param run_number the run_number (mind the leading zeros for old runs) of the run
#' @param fastq_filename the fastq filename, preferrably without the extension or any leading path info. However, an effort has been made to deal with full paths and extensions
#' @param lts_align_expr_prefix the location of the run directories. Eg, if you are mounted and on your local computer, it might be something like "/mnt/htcf_lts/lts_align_expr"
#' @param bam_suffix the common bam suffix for all bam files stored in /lts. Eg, it might be something like "_sorted_aligned_reads_with_annote.bam"
#'
#'
#' @return a verifiec filepath to the bam file
#' @export
createBamPath = function(run_number, fastq_filename, lts_align_expr_prefix, bam_suffix="_sorted_aligned_reads_with_annote.bam"){

  fastqfile_basename = str_remove(basename(fastq_filename), '.fastq.gz')

  bam_path = file.path(paste0(lts_align_expr_prefix, "run_", as.character(run_number),"_samples/align"),
                       paste0(fastqfile_basename, bam_suffix))

  # if the path is valid, return the bam_path. else, stop with error
  ifelse(file.exists(bam_path), bam_path, stop(paste0("The following path is invalid: ", bam_path)))
}

#' Given a GenomicFeatures annotation_db and a gene_id, extract an GRanges object of the exons
#'
#' @param annotation_db a GenomicFeatures db. You can either get this from the bioconductor resources, or create your own with a gtf
#' @param gene_id the ID of a gene in the db. Eg, for cryptococcus CKF44_05222
#'
#' @usage exon_ranges = exonRantes(kn99_db, "CKF44_05222")
#' @references GRanges, GRanges, GenomicFeatures
#' @return an IRanges object of the given gene's exons
#' @export
exonRanges = function(annotation_db, gene_id){

  # TODO add error checking -- check that gene_id is in annotation_db

  exon_regions = exonsBy(annotation_db, by="gene")[[gene_id]]

  gr = GRanges(seqnames = as.character(exon_regions@seqnames@values), ranges = exon_regions@ranges )

}

#' Given a GenomicFeatuers annotation_db and a gene_id, extract an GRanges object of the cds
#'
#' @param annotation_db a GenomicFeatures db. You can either get this from the bioconductor resources, or create your own with a gtf
#' @param gene_id the ID of a gene in the db. Eg, for cryptococcus CKF44_05222
#' @param feature one of c("cds", "exon"), determins which feature to extract from the annotations
#' @usage exon_ranges = exonRantes(kn99_db, "CKF44_05222")
#' @references GRanges, GRanges, GenomicFeatures
#' @return an IRanges object of the given gene's exons
#' @export
featureGRanges = function(annotation_db, gene_id, feature){

  if(!feature %in% c('cds', 'exon')){
    stop("Only 'cds' and 'exon' are supported currently")
  }

  regions = switch(feature,
                   'cds' = cdsBy(annotation_db, by="gene")[[gene_id]],
                   'exon' = exonsBy(annotation_db, by="gene")[[gene_id]])

  GRanges(seqnames = as.character(regions@seqnames@values), ranges = regions@ranges )
}

#'
#' from the metadata libraryProtocol column, determine the library strandedness.
#' @description Currently set up for cryptococcus. E7420L returns reverse, SolexaPrep returns unstranded. default return is unstranded
#'
#' @note: default is unstranded
#' @param metadata a metadata sheet that at leat contains the libraryProtocol column
#'
#' @return the strandedness of the library based on the value in the libraryProtocol column, or 'unstranded' by default
#' @export
determineStrandedness = function(metadata){

  switch(metadata$libraryProtocol,
         "E7420L" = "reverse",
         "SolexaPrep" = "unstranded",
         "unstranded")

}

#'
#' @param library_strandedness one of c("reverse","forward","unstranded")
#' @param gene_ranges a GRanges object, created from a GenomicFeatures db, that at least contains the gene_id of interest (maybe just the whole gene db)
#'
#' @return the strand of the gene of interest
#' @export
mappingStrand = function(library_strandedness, gene_ranges, gene_id){

  as.vector(strand(gene_ranges[names(gene_ranges) == gene_id]))

}
