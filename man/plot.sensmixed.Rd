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
                             cex = 2,  interact.symbol = ":", isFixed = TRUE, 
                             isRand = TRUE, isScaling = TRUE, ...)
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
\item{isFixed}{
logical. Whether to plot tests of the fixed effects
}
\item{isRand}{
logical. Whether to plot tests of the random effects
}

\item{isScaling}{
logical. Whether to plot the scaling factor if present
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
res <- sensmixed(c("Coloursaturation", "Colourbalance"),
                  Prod_effects=c("TVset"), 
                  individual="Assessor", data=TVbo, MAM=TRUE, 
                  reduce.random=FALSE)
plot(res)
plot(res, mult = TRUE)
plot(res, interact.symbol = " x ")
}
