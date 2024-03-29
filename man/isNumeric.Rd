% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{isNumeric}
\alias{isNumeric}
\title{test if argument is numeric}
\usage{
isNumeric(x)
}
\arguments{
\item{x}{any R object}
}
\description{
copied directly from the limma codebase
}
\details{
copied from the limma docs: This function is used to check the validity of arguments for numeric functions. It is an attempt to emulate the behavior of internal generic math functions. IsNumeric differs from is.numeric in that data.frames with all columns numeric are accepted as numeric.
}
