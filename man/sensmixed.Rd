\name{sensmixed}
\alias{sensmixed}

\title{Analysis of sensory data within a mixed effects model framework}

\description{
Constructs a mixed effects model for each of the selected by user attributes 
according to the specified by the user structure arguments. If required, then
the random structures are reduced by eliminating NS random effects.
The likelihood ratio test (LRT) is used for testing random terms, F-type 
hypothesis test is used for testing fixed terms
}
\usage{
sensmixed(attributes=NULL, Prod_effects, replication = NULL, 
                              individual, data, product_structure = 3, 
                              error_structure ="No_Rep", MAM = FALSE,
                              mult.scaling = FALSE, oneway_rand = TRUE,
                              MAM_PER = FALSE, adjustedMAM = FALSE, 
                              alpha_conditionalMAM = 1,
                              calc_post_hoc = FALSE, parallel = FALSE, 
                              reduce.random=TRUE, alpha.random = 0.1, 
                              alpha.fixed = 0.05, interact.symbol = ":", 
                              keep.effs = NULL,  ...)
}
\arguments{
  \item{attributes}{
    
     vector with names of sensory attributes
}

\item{Prod_effects}{
    
     names of the variables related to the product
}

\item{replication}{
    
     names of the replication column in the data, if present
}

\item{individual}{
    
    name of the column in the data that represent assessors
}

\item{data}{
    
     data frame (data from sensory studies)
}

\item{product_structure}{
    
    one of the values in c(1, 2, 3). 1: only main effects will enter the initial
    biggest model. 2: main effects and 2-way interaction. 
    3: all main effects and all possible interaction 
}

\item{error_structure}{
    
    one of the values in c("No_Rep", "2-WAY", "3-WAY"). "No_Rep" and "2-WAY" - 
    assessor effect and all possible interactions between assessor and 
    Product_effects. "3-WAY" - assessor and replicate effect and interaction 
    between them and interaction between them and Product_effects
}

\item{MAM}{
    
    logical. if MAM model should be constructed (scaling correction)
}
\item{mult.scaling}{
    
    logical. Whether multiple scaling should be used
}
\item{oneway_rand}{
    
    logical. Whether there should be just prod effect as part of the random part in MAM
}
\item{MAM_PER}{
    
    logical. if MAManalysis function should be called (scaling correction)
}
\item{adjustedMAM}{
    
    logical. should MAM be adjusted for the scaling
}
\item{alpha_conditionalMAM}{
    
    logical. scaling should be part of the model in case its p-value 
    is less than alpha_conditionalMAM
}
\item{calc_post_hoc}{
    
    logical. Should the post hoc analysis be performed on the final resuced 
    models for all the attributes
}
\item{parallel}{
    logical. Should the computation be done in parallel. the default is FALSE   
}

\item{reduce.random}{
    
     logical. Eliminate non-significant random effects according to 
     alpha.random or not. The default is TRUE
}

\item{alpha.random}{
    
     significance level for elimination of the random part (for LRT test)
}

\item{alpha.fixed}{
    
     significance level for elimination of the fixed part (for F test)
}
\item{interact.symbol}{
    
     symbol for the indication of the interaction between effects. the default 
     one is ":".
}
\item{keep.effs}{
    
     which effects should be kept in a model.
}

  \item{\dots}{
     other potential arguments.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{FCHi}{matrix with Chi square values from LRT test and F values 
 form F-type test for the selected attributes}
  \item{pvalue}{matrix withp-values for all effects for the selected attributes}
  }
  
\author{
Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}

\examples{


#import SensMixed package
library(SensMixed)

#convert some variables to factors in TVbo
TVbo <- convertToFactors(TVbo, c("Assessor", "Repeat", "Picture"))

#run automated selection process
res <- sensmixed(c("Coloursaturation", "Colourbalance"),
Prod_effects = c("TVset", "Picture"), replication="Repeat", 
individual="Assessor", data=TVbo, MAM=TRUE)


## run MAManalysis function
res_MAM <- sensmixed(c("Coloursaturation", "Colourbalance"),
                 Prod_effects=c("TVset"), replication="Repeat", 
                  individual="Assessor", data=TVbo, MAM_PER=TRUE)
## print is not yet implemented
## get anova part
res_MAM[[3]][,,1]

## compare with the general implementation
res <- sensmixed(c("Coloursaturation", "Colourbalance"),
                  Prod_effects=c("TVset"), 
                  individual="Assessor", data=TVbo, MAM=TRUE, 
                  reduce.random=FALSE)
res$fixed      

\dontrun{
plot F and Chi square values
plot(result)
}



}