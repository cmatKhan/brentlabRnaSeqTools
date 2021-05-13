#'
#' test if argument is numeric
#' @description copied directly from the limma codebase
#' @param any_object
#' @details copied from the limma docs: This function is used to check the validity of arguments for numeric functions. It is an attempt to emulate the behavior of internal generic math functions. IsNumeric differs from is.numeric in that data.frames with all columns numeric are accepted as numeric.
#'
#' @export
isNumeric = function (x) {
  is.numeric(x) || (is.data.frame(x) && length(x) > 0 && all(unlist(lapply(x, is.numeric))))
}
