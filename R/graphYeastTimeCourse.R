#' Plot time vs normalized counts of a given gene over n samples, faceted by librarydate and run
#'
#' @importFrom dplyr left_join filter
#' @import ggplot2 scales
#'
#' @description graph_output = graphTimeCourse(cst6_sample_metadata, 'CST6', 'YIL036W', cst6_norm_counts)
#'
#' @note see the vignette called yeast_timecourse_qc.Rmd. There is also sample data that comes with this
#'       package, so you can run the vignette verbatim to see how it works. The vignette includes code that
#'
#' @param metadata_df where metadata$fastqFileName is equal in format to colnames(raw_counts)
#' @param genotype_1 is the entry in the metadata genotype1 column which you would like to examine
#' @param gene_id is the gene id (the systematic name as opposed to the 'common name') most likely corresponding to genotype_1. This is how the correct row is extracted from the count data
#' @param norm_counts MUST have rownames assigned to the gene_ids includes at least all samples in metadata
#'
#' @return a ggplot graph of the gene_id in question, faceted on run and library date
#'
#' @export
graphYeastTimeCourse = function(metadata_df, genotype_1, gene_id, norm_counts){
  fltr_metadata = metadata_df %>% filter(genotype1 == genotype_1)

  norm_count_long_df = tibble(fastqFileName = fltr_metadata$fastqFileName, norm_count = norm_counts[gene_id, fltr_metadata$fastqFileName])

  timecourse_graph = norm_count_long_df %>%
    left_join(metadata_df) %>%
    ggplot(aes(timePoint, norm_count, color=treatment))+
    geom_line(aes(group=interaction(libraryDate, treatment, replicate)))+
    geom_point()+
    scale_y_continuous(trans = log2_trans(),
                       breaks = trans_breaks("log2", function(x) 2^x),
                       labels = trans_format("log2", math_format(2^.x)))+
    ggtitle(genotype_1)+
    facet_wrap(~runNumber+libraryDate)
  return(timecourse_graph)
}
