\name{saveToDoc}
\alias{saveToDoc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
save the result in tables into a doc file for sensmixed or conjoint objects
}
\description{
save the tests for the random and fixed effects into a doc file 
for sensmixed or conjoint objects
}
\usage{
saveToDoc(x, file = NA, bold = FALSE, append = TRUE, type = "html", typeEffs = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
object of class sensmixed or conjoint.
}
  \item{file}{
a character string naming the file to write to, or NULL to stop sink-ing.
}
  \item{bold}{
logical. Should the significance be in bold text instead of the stars. 
The default is FALSE
}
  \item{append}{
  logical. If TRUE, output will be appended to file; 
  otherwise, it will overwrite the contents of file.
  }
  \item{type}{
  type of output as in xtable. Either "html" or "latex"
  }
  \item{typeEffs}{
  one of the numbers in c(1,2,3,4). 
  1: save results for the random effects in a table
  2: save results for the fixed effects in a table
  3: save results for the scaling effects in a table
  4: save results for all effects in a table  
  }
}



\author{
Alexandra Kuznetsova
}



\examples{
\dontrun{
res <- sensmixed(c("Coloursaturation", "Colourbalance"),
                  prod_effects=c("TVset"), 
                  assessor="Assessor", data=TVbo)

saveToDoc(res, file = "C:/Desktop/output.doc")
}
}

