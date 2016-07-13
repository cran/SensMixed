\name{sensmixed}
\alias{sensmixed}

\title{Analysis of sensory data within a mixed effects model framework}

\description{
Constructs mixed effects models for each of the selected by a user attributes. By default the largest possible models (that contain all possible interactions in fixed and random parts) are fitted. The complexity of the fitted models can be changed.  Non-significant random effects are eliminated (by default). The likelihood ratio test (LRT) is used for testing the random terms, F-type hypothesis test is used for testing the terms. The type of the model and the type of the analysis can be changed with the \code{control} argument (see \code{sensmixedControl()})
}
\usage{
sensmixed(attributes, prod_effects, assessor, 
          replication = NULL, data, product_structure = 3,
          error_structure ="ASS-REP", MAM = TRUE,
          control = sensmixedControl())
}
\arguments{
  \item{attributes}{
    
    a vector with names of sensory attributes
}

\item{prod_effects}{
    
     a vector with the names of the variables related to the product
}

\item{replication}{
    
     a character with the name of the replication column in the data, if present
}

\item{assessor}{
    
    a characthe with the name of the column in the data that represents assessors
}

\item{data}{
    a data frame (data from sensory studies)
}

\item{product_structure}{
    numeric, takes values in  in c(1, 2, 3). Specifies the complexity of the fixed part (product effects) of the mixed effects models for all attributes. 
    \describe{
      \item{product_structure = 1}{only main effects}
      \item{product_structure = 2}{main effects and 2-way interactions}
      \item{product_structure = 3}{all main effects and all possible interaction}
    }
}
\item{error_structure}{
    character, takes values in c("ONLY-ASS", "ASS-REP"). Specifies the complexity of the random part of the mixed effects models for all attributes. 
    \describe{
      \item{error_structure = "ONLY-ASS"}{assessor effect and all possible interactions between assessor and product effects}
      \item{error_structure = "ASS-REP"}{assessor and replicate effect and interaction 
    between them and interaction between them and product effects}
    }
}

\item{MAM}{
    logical. if TRUE then mixed assessor models (MAM) are fitted for the selected attributes  (see Brockhoff, P. B., Schlich, P., & Skovgaard, I. (2015))
}
\item{control}{
    a list (of class \code{sensmixedControl()}) containing control parameters. See the \code{sensmixedControl()} documentation for details.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{FCHi}{matrix with Chi square values from LRT test and F values 
 form F-type test for the selected attributes}
  \item{pvalue}{matrix withp-values for all effects for the selected attributes}
  }
  
\seealso{\link{sensmixedControl}, \link{conjoint}, \link{SensMixedUI}}
  
\author{
Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}

\references{
Brockhoff, P. B., Schlich, P., & Skovgaard, I. (2015). Taking individual
scaling differences into account by analyzing profile data with the mixed
assessor model. Food Quality and Preference, 39 , 156-166.

Kuznetsova, A., Christensen, R. H., Bavay, C., & Brockhoff,
P. B. (2015). Automated mixed ANOVA modeling of sensory and consumer data. 
Food Quality and Preference, 40,
Part A, 31 { 38. URL: http://www.sciencedirect.com/science/
article/pii/S0950329314001724. doi:http://dx.doi.org/10.1016/j.
foodqual.2014.08.004.
}
}

\examples{


## import SensMixed package
library(SensMixed)

## convert some variables to factors in TVbo
TVbo <- convertToFactors(TVbo, c("Assessor", "Repeat", "Picture"))

## run automated selection process
res <- sensmixed(c("Coloursaturation", "Colourbalance"),
prod_effects = c("TVset", "Picture"), 
assessor="Assessor", data=TVbo, MAM=TRUE)

res


## run MAManalysis function
res_MAM <- sensmixed(c("Coloursaturation", "Colourbalance"),
                 prod_effects=c("TVset"), replication="Repeat", 
                 assessor="Assessor", data=TVbo, control = list(MAM_balanced=TRUE))
## print is not yet implemented
## get anova part
res_MAM[[3]][,,1]

## compare with the general implementation
res <- sensmixed(c("Coloursaturation", "Colourbalance"),
                  prod_effects=c("TVset"), 
                  assessor="Assessor", data=TVbo, MAM=TRUE, 
                  control = list(reduce.random=FALSE))
res$fixed      

\dontrun{
res <- sensmixed(names(TVbo)[5:(ncol(TVbo) - 1)],
                  prod_effects=c("TVset", "Picture"), 
                  assessor="Assessor", 
                  data=TVbo)


plot F and Chi square values
plot(res)
}
}