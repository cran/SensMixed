\name{TVbo}
\alias{TVbo}
\title{
  TV dataset
}
\description{
  The \code{TVbo} dataset comes from Bang and Olufsen company. The main purpose was to test products, specified by two attributes Picture and TVset. 15 different response variables (characteristics of the product) were assessed by trained panel list.
  
}
\usage{
TVbo
}
\format{
\describe{
  \item{\code{Assessor}}{
    factor: numbering identifying assessors
  }
  \item{\code{TVset}}{
    factor: attribute of the product
  }
  \item{\code{Picture}}{
    factor: attribute of the product
  }
  \item{\code{15 Characteristics of the product}}{
    numeric variables: Coloursaturation, Colourbalance, Noise, Depth, Sharpness, Lightlevel, Contrast, Sharpnessofmovement, Flickeringstationary, Flickeringmovement,
 Distortion, Dimglasseffect, Cutting, Flossyedges, Elasticeffect   
  }
  
}}
\source{
  Bang and Olufsen company
}



\keyword{datasets}

\examples{
## import SensMixed package
library(SensMixed)

## convert some variables to factors in TVbo
TVbo <- convertToFactors(TVbo, c("Assessor", "Repeat", "Picture"))

## run automated selection process
res <- sensmixed(c("Coloursaturation", "Colourbalance"),
prod_effects = c("TVset", "Picture"), 
assessor="Assessor", data=TVbo, MAM=FALSE)
res

}