\name{convertToFactors}
\alias{convertToFactors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
converts variables of the data frame to factors
}
\description{
the user specifies which variables he/she would like to consider as factors, the functions converts them to factors
}
\usage{
convertToFactors(data, facs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  data frame }

 \item{facs}{
  vector with names of variables that the user would like to convert to factors }
}

\value{
returns the same data frame as in the input but  with the specified variables converted to factors
}

\author{
Alexandra Kuznetsova
}
\examples{
library(SensMixed)

str(ham)

ham <- convertToFactors(ham, c("Gender"))

str(ham)
}