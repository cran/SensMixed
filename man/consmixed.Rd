\name{consmixed}
\alias{consmixed}


\title{Automated model selection process for the Consumer data}

\description{
Constructs the biggest possible model and reduces it to the best by principle of parcimony. First elimination of random effects is performed following by elimination of fixed effects.
The LRT test is used for testing random terms, F-type hypothesis test is used for testing fixed terms. The post-hoc and plots are provided
}
\usage{
consmixed(response, Prod_effects, Cons_effects=NULL,
Cons, data, structure = 3, alpha.random = 0.1, alpha.fixed = 0.05, ...)
}
\arguments{
  \item{response}{    
     name of the liking variable in the Consumer data}

\item{Prod_effects}{
    
    vector with names of the variables associated with products
}

\item{Cons_effects}{
    
    vector with names of the effects associated with consumers
}

\item{Cons}{
    
    name of the column in the data that represents consumers
}

\item{data}{
    
     data frame (data from consumer studies)
}

\item{structure}{
    
    one of the values in c(1,2,3). 1:Analysis of main effects, Random consumer effect AND interaction  between consumer and the main effects(Automized reduction in random part, NO reduction in fixed part). 2: Main effects AND all 2-factor interactions. Random consumer effect AND interaction between consumer and  all fixed effects (both main and interaction ones). (Automized reduction in random part, NO reduction in fixed part). 3: Full factorial model with ALL possible fixed and random effects.  (Automized reduction in random part, AND automized reduction in fixed part).
}

\item{alpha.random}{
    
     significance level for elimination of the random part (for LRT test)
}

\item{alpha.fixed}{
    
     significance level for elimination of the fixed part (for F test)
}

  \item{\dots}{
     other potential arguments.
}
}

\value{
\item{rand.table}{table with value of Chi square test, p-values e t.c. for the random effects}
  \item{anova.table }{table which tests whether the model fixed terms are significant (Analysis of Variance)}
  \item{model}{Final model - object of class lmer or gls (after all the required reduction has been performed)}

  }
  
\author{
Alexandra Kuznetsova, Per Bruun Brockhoff, Rune Haubo Bojesen Christensen
}

\examples{
library(SensMixed)

#convert some variables to factors in Tham
ham <- convertToFactors(ham, c("Consumer", "Product", "Information", "Gender"))


consmixed(response="Liking", 
Prod_effects= c("Product","Information"), 
Cons_effects=c("Gender","Age"), Cons = "Consumer", data =ham, structure=1)

}