\name{plotLSMEANS}
\alias{plotLSMEANS}


\title{plots bars for the LSMEANS or differences of LSMEANS}

\description{
plots  bars for the LSMEANS or differences of LSMEANS for product factors
and confidence intervals
}
\usage{
plotLSMEANS(table, response, 
                        which.plot=c("LSMEANS", "DIFF of LSMEANS"), 
                        main = NULL, cex = 1.4, effs = NULL, mult = TRUE)
}
\arguments{
  \item{table}{    
     data table containing LSMEANS/ DIFFLSMEANS table from the step function of the lmerTest package}

\item{response}{
    
    vector with the name of the attribute, for which the LSMEANS / DIFFLSMEANS are calculated
}

\item{which.plot}{
    
    name, indicating the type of plot to generate.
}

\item{main}{
    
    name of the title for the plot
}

\item{cex}{
    
     cex for representing the plot (UNUSED?..)
}

\item{effs}{
    
   vector with the names for the effects, for which to plot the LSMEANS / DIFFLSMEANS
}

\item{mult}{
    
     logical. TRUE means plot LSMEANS for all effects in one layout
}
\item{\dots}{
     other potential arguments.
}
}

\value{
 barplots created via ggplot2 package

  }
  
\author{
Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}

