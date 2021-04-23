#' Create the libraryProtocol + libraryDate model matrix with the earliest date of each library protocol dropped
#'
#' @param metadata_df the joined tables of the database (biosample to quality assess)
#' @return a model matrix constructed as specified in the description
#'
#' @export
createLibraryProtocolLibrarydateModelMatrix = function(metadata_df){
  colnames(metadata_df) = toupper(colnames(metadata_df))
  # cast librarydate to datetime object
  metadata_df$LIBRARYDATE = as.Date(metadata_df$LIBRARYDATE)
  # create model.frame
  x = model.frame(~LIBRARYPROTOCOL+LIBRARYDATE, metadata_df)
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
  # create a model matrix with a column for the intercept, the protocol, and the number of unique library dates - 1
  model_matrix = matrix(0L, nrow=nrow(x), ncol=length(librarydate_columns)+2)
  # label rows/columns
  rownames(model_matrix) = rownames(x)
  colnames(model_matrix) = c("(Intercept)", "LIBRARYPROTOCOL", librarydate_columns)
  # set intercept to 1
  model_matrix[,1] = 1

  for (i in 1:nrow(x)){
    # extract info of interest for a given sample
    fastqfilename = rownames(x)[i]
    librarydate = as.character(x$LIBRARYDATE[i])
    protocol = x$LIBRARYPROTOCOL[i]

    if(protocol == "E7420L"){
      model_matrix[fastqfilename, "LIBRARYPROTOCOL"] = 1
    }
    if(librarydate != min_old_protocol_date & librarydate != min_new_protocol_date){
      model_matrix[fastqfilename, librarydate] = 1
    }
  }
  return(model_matrix)
} # end createLibrarydateModelMatrix()
