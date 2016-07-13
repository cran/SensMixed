\name{sensmixedControl}
\alias{sensmixedControl}
\title{Control of sensmixed function}
\description{
  Construct control structures for \code{sensmixed} function.. All arguments
  have defaults, and can be grouped into
  
  \itemize{
    \item mixed assessor model (MAM) parameters  (are only used if MAM = TRUE in the \code{sensmixed} function)
    \item model analysis specifications
  }
}
\usage{
sensmixedControl(MAM_mult_scaling = FALSE, MAM_oneway_rand = FALSE,
                              MAM_balanced = FALSE, MAM_adjusted = FALSE, 
                              MAM_alpha_conditional = 1,
                              calc_post_hoc = FALSE, parallel = FALSE, 
                              reduce.random=TRUE, alpha.random = 0.1, 
                              alpha.fixed = 0.05, interact.symbol = ":", 
                              keep.effs = NULL)
}
\arguments{
\item{MAM_mult_scaling}{
    logical. Whether multiple scaling should be used. This option is not fully investigated (recommened is FALSE)
}
\item{MAM_oneway_rand}{
    logical. Whether there should be just one \code{prod_effect} (in a multi-factorial product structure setting this means the highest order interaction) as part of the random part in MAM
}
\item{MAM_balanced}{
    logical. if \code{MAManalysis} function should be called (scaling correction). Can only be used for one product effect, no replications and balanced data. 
}
\item{MAM_adjusted}{
    logical. should MAM be adjusted for the scaling (is only used if \code{MAM_balanced = TRUE})
}
\item{MAM_alpha_conditional}{
    logical. scaling should be part of the model in case its p-value 
    is less than \code{MAM_alpha_conditional} (is only used if \code{MAM_balanced = TRUE})
}
\item{calc_post_hoc}{
    logical. Should the post hoc analysis be performed  for the attributes
}

\item{reduce.random}{
     logical. Eliminate non-significant random effects according to 
     alpha.random or not. The default is TRUE
}

\item{keep.effs}{
     a vector with the names which effects should be kept in a model (even if non-significant). By default the assessor effect and interaction between assessor and highest order product effects are always kept.
}

\item{alpha.random}{
     significance level for elimination of the random part (for LRT test)
}

\item{alpha.fixed}{
     significance level for elimination of the fixed part (for F test)
}
\item{interact.symbol}{
     symbol for the indication of the interaction between effects (for presenting in tables and plots). the default is ":".
}


\item{parallel}{
    logical. Should the computation be done in parallel for the attributes (saves computational time). the default is FALSE   
}

}
\value{
  The \code{sensmixedControl} function returns a list containing
  \enumerate{
      \item  mixed assessor model fitting related control parameters (are used only if 
    \code{MAM=TRUE} in the \code{sensmixed} function).
    \item  control parameters for the analysis.
  }
}
\seealso{\link{sensmixed}}

\examples{
str(sensmixedControl())
}
