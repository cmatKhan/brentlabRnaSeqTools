#'
#'
#'
#'
#'
#' @export
extractRLEByReplicateGroup_90minInduction = function(meta_qual_df, norm_counts, output_dirpath, protocol_selector, already_logged_flag){

  dir.create(output_dirpath)

  for (genotype in unique(meta_qual_df$GENOTYPE)){
    print(genotype)
    # filter out known strain/geno problems
    if (!genotype %in% c('CNAG_01523.CNAG_05431', 'CNAG_02153.CNAG_05431', 'CNAG_05222.CNAG_02566', 'CNAG_02566.CNAG_01438')){
      df = meta_qual_df %>% filter(GENOTYPE == genotype)

      if(protocol_selector){
        x = split(df, f = list(df$MEDIUM, df$TEMPERATURE, df$ATMOSPHERE, df$TREATMENT, df$OTHER_CONDITIONS, df$TIMEPOINT, df$LIBRARYPROTOCOL))
      } else {
        x = split(df, f = list(df$MEDIUM, df$TEMPERATURE, df$ATMOSPHERE, df$TREATMENT, df$OTHER_CONDITIONS, df$TIMEPOINT))
      }

      gene_id_col = seq(1,nrow(norm_counts))
      for (split_df in x){
        # only take genotypes with replicates > 3
        if (nrow(split_df) > 2){
          medium = as.character(unique(split_df$MEDIUM))
          temperature = as.character(unique(split_df$TEMPERATURE))
          atmosphere = as.character(unique(split_df$ATMOSPHERE))
          treatment = as.character(unique(split_df$TREATMENT))
          other_conditions = as.character(unique(split_df$OTHER_CONDITIONS))
          timepoint = as.character(unique(split_df$TIMEPOINT))
          protocol = as.character(unique(split_df$LIBRARYPROTOCOL))

          ftlr_norm_counts = norm_counts[,split_df$FASTQFILENAME]
          if(protocol_selector){
            full_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint,protocol, "full.csv", sep='_'), sep='/')
            summary_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint,protocol, "summary.csv", sep='_'), sep='/')
          }else{
            full_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint, "full.csv", sep='_'), sep='/')
            summary_filename = paste(output_dirpath, paste(genotype,medium,temperature,atmosphere,treatment,other_conditions,timepoint, "summary.csv", sep='_'), sep='/')
          }

          if (already_logged_flag){
            fltr_rle_full = createFullRLETableAlreadyLogged(ftlr_norm_counts)
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
        other_conditions = as.character(unique(split_df$OTHER_CONDITIONS))
        timepoint = as.character(unique(split_df$TIMEPOINT))
        protocol = as.character(unique(split_df$LIBRARYPROTOCOL))


        full_filename = paste(output_dirpath, paste(medium,temperature,atmosphere,treatment,other_conditions,timepoint, "full.csv", sep='_'), sep='/')
        summary_filename = paste(output_dirpath, paste(medium,temperature,atmosphere,treatment,other_conditions,timepoint, "summary.csv", sep='_'), sep='/')

        ftlr_norm_counts = norm_counts[,split_df$FASTQFILENAME]
        if (already_logged_flag){
          fltr_rle_full = createFullRLETableAlreadyLogged(ftlr_norm_counts)
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
