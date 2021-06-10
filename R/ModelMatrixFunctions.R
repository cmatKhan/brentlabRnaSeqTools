# TODO CLEAN THIS UP

#' Create the libraryProtocol + libraryDate model matrix with the earliest date of each library protocol dropped
#'
#' @param metadata_df the joined tables of the database (biosample to quality assess)
#' @return a model matrix constructed as specified in the description
#'
#' @export
createNinetyMinInductionModelMatrix = function(metadata_df){
  wt_genotype = "CNAG_00000"
  colnames(metadata_df) = toupper(colnames(metadata_df))
  # cast librarydate to datetime object
  metadata_df$LIBRARYDATE = as.Date(metadata_df$LIBRARYDATE)
  # create model.frame
  x = model.frame(~LIBRARYPROTOCOL+LIBRARYDATE+GENOTYPE1, metadata_df)
  # x(the model.frame) and metadata are in the same order -- label rows by metadata_df$FASTQFILENAME column
  rownames(x) = metadata_df$FASTQFILENAME
  # extract min old protocol and min new protocol dates
  min_old_protocol_date = as.character(min(x[x$LIBRARYPROTOCOL=="SolexaPrep", ]$LIBRARYDATE))
  print(min_old_protocol_date)
  min_new_protocol_date = as.character(min(x[x$LIBRARYPROTOCOL=="E7420L", ]$LIBRARYDATE))
  print(min_new_protocol_date)
  # remove the min_old_protocol_date from the list of unique library dates
  librarydate_columns = unique(as.character(x$LIBRARYDATE))[unique(as.character(x$LIBRARYDATE)) != min_old_protocol_date] # combine these lines with %in%
  librarydate_columns = librarydate_columns[librarydate_columns != min_new_protocol_date]
  # create genotype columns
  genotype_columns = unique(as.character(x$GENOTYPE1))[unique(as.character(x$GENOTYPE1)) != wt_genotype]
  # create a model matrix with a column for the intercept, the protocol, and the number of unique library dates - 1
  model_matrix = matrix(0L, nrow=nrow(x), ncol=length(librarydate_columns)+2+length(genotype_columns))
  # label rows/columns
  rownames(model_matrix) = rownames(x)
  colnames(model_matrix) = c("(Intercept)", "LIBRARYPROTOCOL", librarydate_columns, genotype_columns)
  # set intercept to 1
  model_matrix[,1] = 1

  for (i in 1:nrow(x)){
    # extract info of interest for a given sample
    fastqfilename = rownames(x)[i]
    librarydate = as.character(x$LIBRARYDATE[i])
    protocol = x$LIBRARYPROTOCOL[i]
    genotype1 = x$GENOTYPE1[i]

    if(protocol == "E7420L"){
      model_matrix[fastqfilename, "LIBRARYPROTOCOL"] = 1
    }
    if(librarydate != min_old_protocol_date & librarydate != min_new_protocol_date){
      model_matrix[fastqfilename, librarydate] = 1
    }
    if(genotype1 != wt_genotype){
      model_matrix[fastqfilename, genotype1] = 1
    }
  }

  return(model_matrix)

} # end createLibrarydateModelMatrix()


#' filter low replicate parameters from metadata
#'
#' @import dplyr
#'
#' @description given a model formula, remove samples with less than a specified number of replicates from the metadata
#'
#' @param metadata_df a data frame that contains at least the model paramters of interest
#' @param design_formula an R formula, eg ~libraryDate+treatment, of parameters contained in the metadata_df
#' @param replicate_threshold the number of replicates below which samples will be removed. Default is 2
#'
#' @return the input metadata with samples in replicate groups with less than the specified thershold filtered out
#'
#' @export
fltrLowReplicateParams = function(metadata_df, design_formula, replicate_threshold=2){

  # TODO write test and error handling for this function

  # create a character vector of the formula parameters
  # eg 'librarydate', 'treatment' if the original formula is ~librarydate+treatment
  formula_vector = str_trim(unlist(str_split(as.character(design_formula)[[2]], "\\+")))

  low_rep_params = -1
  while(length(low_rep_params != 0)){

    # create a model matrix and summary of number of replicates in each of the parameters
    model_matrix = model.matrix(design_formula, metadata_df)
    mm_summary_df = tibble(model_params = colnames(model_matrix), replicate_tally = colSums(model_matrix))

    # get the model parameters' factor levels with samples that fall in replicate groups below a specified number
    low_rep_params = mm_summary_df %>%
      filter(replicate_tally < replicate_threshold) %>%
      pull(model_params)

    low_rep_params = str_remove(low_rep_params, formula_vector)

    # remove samples with the factor levels in the parameters specified in low_rep_params
    samples_to_remove = unlist(lapply(formula_vector, function(x)
      filter(metadata_df, !!rlang::sym(x) %in% low_rep_params) %>% pull(FASTQFILENAME)))

    metadata_df = droplevels(filter(metadata_df, !FASTQFILENAME %in% samples_to_remove))
  }

  return(metadata_df)

} # end fltrLowReplicateParams()
