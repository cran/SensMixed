\name{plot}
\alias{plot.sensmixed}

\title{
function creates plots for the sensmixed object 
}
\description{
function creates barplots for the  square roots of F statistics and 
square roots of chi square values for all attributes
}
\usage{
\method{plot}{sensmixed} (x, mult = FALSE, dprime = FALSE, sep = FALSE, 
                             cex = 2,  interact.symbol = ":",
                             isRand = TRUE, isScaling = TRUE, stacked = TRUE, ...)
}

\arguments{
  \item{x}{
object of class sensmixed
}
 \item{mult}{
logical. Should multiple plots be plotted, that is barplots 
for each effect in a separate plot
}
 \item{dprime}{
logical. Should multiattribute plot for product effects use average squared dprimes 
instead of square root of F statistics
}
\item{sep}{
logical. If TRUE then separate plot is plotted for each effect 
(mult argument should be then also TRUE)
}
\item{cex}{
The magnification to be used
}
\item{interact.symbol}{
The symbol to be used for the interaction effects
}
\item{isRand}{
logical. Whether to plot tests of the random effects
}

\item{isScaling}{
logical. Whether to plot the scaling factor if present
}

\item{stacked}{
logical. Whether  bars should be stacked
}

 \item{\dots}{
     other potential arguments.
}
}

\value{
 NULL is returned
}

\author{
Alexandra Kuznetsova
}




\examples{
res <- sensmixed(c("Coloursaturation", "Colourbalance", "Noise"),
                  prod_effects=c("TVset", "Picture"), 
                  assessor="Assessor", data=TVbo, MAM=TRUE, 
                  control=list(reduce.random=FALSE))
plot(res, isRand = TRUE)
plot(res, isRand = FALSE)
plot(res, stacked = FALSE, mult = TRUE)
plot(res, isRand = FALSE, stacked = FALSE, interact.symbol = " x ")
}
