\name{conjoint}
\alias{conjoint}
\title{
Conjoint analysis within a mixed effects model framework.
}
\description{
Performs conjoint analysis within a mixed effects model framework. This function is then used in the conjoint analysis of the ConsumerCheck software tool.
}
\usage{
conjoint(structure = 1, data, response, fixed, random, facs, 
         corr = FALSE, alpha.random = 0.1, alpha.fixed = 0.05)

}

\arguments{
  \item{structure}{
     numerical, takes values in c(1,2,3). Specifies the complexity of the conjoint model. 
    \describe{
      \item{structure = 1}{mixed effects model includes fixed main effects. Random effects consist of random consumer effect and interaction between consumer and the main effects.}
      \item{structure = 2}{mixed effects model includes main effects and all 2-factor interactions. Random effects consist of consumer effect and interaction between consumer and all fixed effects (both main and interaction ones).}
	    \item{structure = 3}{
	    This is a full factorial model with all possible fixed and random effects (i.e. including all main effects and all higher-way interactions). The automated reduction in
random part is followed by an automated reduction in fixed part. }
    }
     }
  \item{data}{
   a data frame.
  }
  \item{response}{
   a vector. Names of the variables that correspond to the liking scores.
  }
  \item{fixed}{
  a list. first element is a vector with a name Product includes names of the design variables. second element with a name Consumer includes names of the consumer characteristics variables
  }
  \item{random}{
  a string with the name for a consumer variable
  }
  \item{facs}{
  a vector with the names of the variables that need to be considered as factors
  }
  \item{corr}{
  a logical vaue. TRUE if the correlations between random effects are included in the model. FALSE if correlations between random effects are set to 0. 
  }
  \item{alpha.random}{
     significance level for elimination of the random part (for LRT test)
  }

  \item{alpha.fixed}{
      significance level for elimination of the fixed part (for F test)
  }

}
\details{
Conjoint analysis (Green and Rao 1971; Green and Srinivasan 1978) is a method for analysing
the effects of design factors and consumer characteristics on consumer likings. A
common approach is to analyse it in a mixed effects model framework, where random effects
consist of consumer effect and interactions between consumer effects and design factors, and
fixed effects consist of design factors and consumer characteristics and possibly interactions
between them.}
\value{
 \item{rand.table}{data frame with value of Chi square statistics, p-values  for the likelihood ratio test for random effects. The order of elimination of non-significant random effects for strucutre = 3}
  \item{anova.table }{data frame with tests for whether the model fixed terms are significant (Analysis of Variance)}
  \item{lsmeans.table}{Least Squares Means data frame with p-values and Confidence intervals}
  \item{diffs.lsmeans.table}{Differences of Least Squares Means data frame with p-values and Confidence intervals. The multiplicity correction for multiple comparisons tests uses Bonferroni method}
  \item{residuals}{These are the residuals of the final model (for structure = 3, the reduced one)
  }
  \item{residualsFixed}{
    These are the residuals that are extracted from a fixed effects model with a saturated fixed structure (main effects and all higher-way interactions form fixed part) and one fixed
Consumer effect.
  }
}

\note{
The tests for the random
effects use likelihood ratio tests while the tests for the fixed effects use the F-test
with Satterthwaites approximation to degrees of freedom. The automated reduction
in the fixed part uses the principle of marginality, i.e. the highest order interactions
are tested first: if they are significant, the lower order effects are not eliminated even
if being non-significant.
}


\author{
Alexandra Kuznetsova, Per Bruun Brockhoff
}


\references{
Green P, Rao V (1971). Conjoint measurement for quantifying judgemental data. Journal
of Marketing Research, 8, 355 - 363

Green P, Srinivasan V (1978). Conjoint analysis in consumer research: Issues and outlook.
Journal of Consumer Research, 5, 103 - 123

}

\seealso{\link{sensmixed}, \link{SensMixedUI}}

\examples{
## check with the  ham
response <- c("Liking")
fixed <- list(Product=c("Product", "Information"), Consumer="Gender")
random <- c("Consumer")
facs <- c("Consumer", "Product", "Information", "Gender")

res.ham <- conjoint(structure=3, ham, response, fixed, random, facs)
res.ham
res.ham$Liking
}