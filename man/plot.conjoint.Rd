\name{plot.conjoint}
\alias{plot.conjoint}

\title{
plots the post-hoc for the conjoint object
}
\description{
plots the least squares means and differences of least squares means 
together with the confidence intervals for the fixed effects
}
\usage{
\method{plot}{conjoint} (x, main = NULL, cex = 1.4, 
                           which.plot = c("LSMEANS", "DIFF of LSMEANS"),
                           test.effs = NULL, ...)
}

\arguments{
  \item{x}{
 object of class conjoint
}
\item{main}{
 string. Title for the plots
}
\item{cex}{
 A numerical value giving the amount by which plotting text and symbols 
 should be magnified relative to the default
}
\item{which.plot}{
 type of plot to be drawn
}
\item{test.effs}{
 name of the effext for which to draw the plots
}
 \item{\dots}{
     other potential arguments.
}
}

\value{
returns NULL
}

\author{
Alexandra Kuznetsova 
}


\examples{

#convert some variables to factors in Tham
response <- c("Liking")
fixed <- list(Product=c("Product", "Information"), Consumer="Gender")
random <- c("Consumer")
facs <- c("Consumer", "Product", "Information", "Gender")

res.ham <- conjoint(structure=3, ham, response, fixed, random, facs)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
