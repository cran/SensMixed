\name{plot.consmixed}
\alias{plot.consmixed}

\title{
plots the post-hoc for the consmixed object
}
\description{
plots the least squares means and differences of least squares means 
together with the confidence intervals for the fixed effects
}
\usage{
\method{plot}{consmixed} (x, main = NULL, cex = 1.4, 
                           which.plot = c("LSMEANS", "DIFF of LSMEANS"),
                           effs = NULL, ...)
}

\arguments{
  \item{x}{
 object of class consmixed
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
\item{effs}{
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
ham <- convertToFactors(ham, c("Consumer", "Product", "Information", "Gender"))

res <- consmixed(response="Liking", 
Prod_effects= c("Product","Information"), 
Cons_effects=c("Gender","Age"), Cons = "Consumer", data =ham, structure=1)

plot(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
