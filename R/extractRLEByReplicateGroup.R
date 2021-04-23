#' These are messy and need to be cleaned up, but are how i calculate rle by replicate group for the 90minInduction and env_pert NOTE: only 90minInduction is currently tested with new database
#'
#' @param meta_qual_df is a dataframe pulled from the database with biosample --> quality assess (name is remnant of old system)
#' @param norm_counts from the deseq object
#' @param output_dirpath where to put the rle results. i suggest making a directory like "rle" in your working directory
#' @param protocol_selector whether or not to include library protocol in the replicate grouping
#' @param already_logged_flag are counts logged? (note: deseq norm counts are not logged unless you log them yourself)
#'
#' @return None. writes to file
#'
#' @export
extractRLEByReplicateGroup_90minInduction = function(meta_qual_df, norm_counts, output_dirpath, protocol_selector, already_logged_flag){

  dir.create(output_dirpath)

  meta_qual_df$LIBRARYDATE = as.factor(as.Date(meta_qual_df$LIBRARYDATE))
  meta_qual_df$MEDIUM = as.factor(meta_qual_df$MEDIUM)
  meta_qual_df$TEMPERATURE = as.factor(meta_qual_df$TEMPERATURE)
  meta_qual_df$ATMOSPHERE = as.factor(meta_qual_df$ATMOSPHERE)
  meta_qual_df$TREATMENT = as.factor(meta_qual_df$TREATMENT)
  meta_qual_df$OTHERCONDITIONS = as.factor(meta_qual_df$OTHERCONDITIONS)
  meta_qual_df$TIMEPOINT = as.factor(meta_qual_df$TIMEPOINT)
  meta_qual_df$LIBRARYPROTOCOL = as.factor(meta_qual_df$LIBRARYPROTOCOL)
  meta_qual_df$GENOTYPE = as.factor(meta_qual_df$GENOTYPE)

  for (genotype in unique(meta_qual_df$GENOTYPE)){
    # filter out known strain/geno problems
    genotype_filtered_df = meta_qual_df %>% filter(GENOTYPE == genotype)

    other_cond_str = 'OTHERCONDITIONS'

      if(protocol_selector){
        replicate_split_meta_qual_list = split(genotype_filtered_df, f = list(genotype_filtered_df$MEDIUM, genotype_filtered_df$TEMPERATURE, genotype_filtered_df$ATMOSPHERE, genotype_filtered_df$TREATMENT, genotype_filtered_df$OTHERCONDITIONS, genotype_filtered_df$TIMEPOINT, genotype_filtered_df$LIBRARYPROTOCOL))
      } else {
        replicate_split_meta_qual_list = split(genotype_filtered_df, f = list(genotype_filtered_df$MEDIUM, genotype_filtered_df$TEMPERATURE, genotype_filtered_df$ATMOSPHERE, genotype_filtered_df$TREATMENT, genotype_filtered_df$OTHERCONDITIONS, genotype_filtered_df$TIMEPOINT))
      }

      gene_id_col = seq(1,nrow(norm_counts))
      for (split_df in replicate_split_meta_qual_list){
        # only take genotypes with replicates > 3
        if (nrow(split_df) > 2){
          medium = as.character(unique(split_df$MEDIUM))
          temperature = as.character(unique(split_df$TEMPERATURE))
          atmosphere = as.character(unique(split_df$ATMOSPHERE))
          treatment = as.character(unique(split_df$TREATMENT))
          other_conditions = as.character(unique(split_df$OTHERCONDITIONS))
          timepoint = as.character(unique(split_df$TIMEPOINT))
          protocol = as.character(unique(split_df$LIBRARYPROTOCOL))

          ftlr_norm_counts = norm_counts[,split_df$FASTQFILENAME]
          if (protocol_selector){
            full_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint,protocol, "full.csv", sep='_'), sep='/')
            summary_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint,protocol, "summary.csv", sep='_'), sep='/')
          } else{
            full_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint, "full.csv", sep='_'), sep='/')
            summary_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint, "summary.csv", sep='_'), sep='/')
          }

          if (already_logged_flag == TRUE){
            fltr_rle_full = createFullRLETable(ftlr_norm_counts, gene_id_column=NULL, logged=TRUE)
            write_csv(as_tibble(fltr_rle_full), full_filename)

            fltr_rle_summary = rleSummary(fltr_rle_full)
            write_csv(as_tibble(fltr_rle_summary), summary_filename)

          } else{
            fltr_rle_full = createFullRLETable(ftlr_norm_counts, gene_id_column=NULL)
            write_csv(as_tibble(fltr_rle_full), full_filename)

            fltr_rle_summary = rleSummary(fltr_rle_full)
            write_csv(as_tibble(fltr_rle_summary), summary_filename)
          }
      }
    }
  }
}# end extractRLEByReplicateGroup()

#' @export
extractRLEByReplicateGroup_EnvPert = function(df, norm_counts, output_dirpath, protocol_selector, already_logged_flag){

  dir.create(output_dirpath)

  x = split(df, f = list(df$MEDIUM, df$TEMPERATURE, df$ATMOSPHERE, df$TREATMENT, df$OTHER_CONDITIONS, df$TIMEPOINT))

  gene_id_col = seq(1,nrow(norm_counts))
  for (split_df in x){
    if (nrow(split_df) != 0){
      # only take genotypes with replicates > 3
      if (nrow(split_df) > 2){
        medium = as.character(unique(split_df$MEDIUM))
        temperature = as.character(unique(split_df$TEMPERATURE))
        atmosphere = as.character(unique(split_df$ATMOSPHERE))
        treatment = as.character(unique(split_df$TREATMENT))
        other_conditions = as.character(unique(split_df$OTHERCONDITIONS))
        timepoint = as.character(unique(split_df$TIMEPOINT))
        protocol = as.character(unique(split_df$LIBRARYPROTOCOL))


        full_filename = paste(output_dirpath, paste(medium,temperature,atmosphere,treatment,other_conditions,timepoint, "full.csv", sep='_'), sep='/')
        summary_filename = paste(output_dirpath, paste(medium,temperature,atmosphere,treatment,other_conditions,timepoint, "summary.csv", sep='_'), sep='/')

        ftlr_norm_counts = norm_counts[,split_df$FASTQFILENAME]
        if (already_logged_flag){
          fltr_rle_full = createFullRLETable(ftlr_norm_counts, gene_id_column=NULL, logged=TRUE)
          write_csv(as_tibble(fltr_rle_full), full_filename)

          fltr_rle_summary = rleSummary(fltr_rle_full)
          write_csv(as_tibble(fltr_rle_summary), summary_filename)

        } else{
          fltr_rle_full = createFullRLETable(ftlr_norm_counts, seq(1,nrow(norm_counts)))
          write_csv(as_tibble(fltr_rle_full), full_filename)

          fltr_rle_summary = rleSummary(fltr_rle_full)
          write_csv(as_tibble(fltr_rle_summary), summary_filename)
        }
      }
    }
  }
}# end extractRLEByReplicateGroup()
