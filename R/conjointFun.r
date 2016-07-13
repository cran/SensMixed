##############################################################################
### performs  Conjoint Analysis for lmer object
##############################################################################


conjointFun <- function(structure = 1, data, response, fixed, random, facs, 
                        corr = FALSE, alpha.random = 0.1, alpha.fixed = 0.05)
{
  ## structure=1  (default structure) : mixed effects model includes fixed main effects. 
  ##             Random effects consist of random consumer effect and 
  ##             interaction between consumer and the main effects.
  ## structure=2 : mixed effects model includes main effects and all 2-factor interactions. 
  ##             Random effects consist of consumer effect and interaction between 
  ##             consumer and all fixed effects (both main and interaction ones).
  ## structure=3 : This is a full factorial model with all possible fixed and 
  ##             random effects (i.e. including all main effects and 
  ##             all higher-way interactions). The automated reduction in random part 
  ##             is followed by an automated reduction in fixed part. 
  

  ## the result that will be returned
  resultFULL <- vector("list",length(response))
  
  #convert some of the variables specified by user to factors
  data  <-  convertToFactors(data, facs)
  
  
  for(ind.resp in 1:length(response))
  {
    
    print(paste("Calculating ", response[ind.resp],"...", sep=" "))
    
    fmodelfull <- createLMERmodelConjoint(3, data, response[ind.resp], fixed, 
                                          random, corr, isFormula = TRUE)
    model <- createLMERmodelConjoint(structure, data, response[ind.resp], fixed, 
                                     random, corr, isFormula = FALSE)
    
    
    isRandReduce <- TRUE
    isFixReduce <- TRUE
    isLsmeans <- TRUE
    ## check if reduction of the fixed/random part is required
    if(structure == 1 || structure == 2){  
      isFixReduce <- FALSE
      isRandReduce <- FALSE
    }
    
    ## check if there are correlations between intercepts and slopes
    checkCorr <- function(model)
    {
      corr.intsl <- FALSE
      lnST <- length(getME(model, "ST"))
      for(i in 1:lnST)
      {    
        if(nrow(getME(model, "ST")[[i]])>1)
          corr.intsl <- TRUE
      } 
      return(corr.intsl) 
    }
    
    if(checkCorr(model))
      isRandReduce <- FALSE
    
    suppressMessages(t <- step(model, reduce.fixed=isFixReduce, 
                               reduce.random=isRandReduce, alpha.random=alpha.random, 
                               alpha.fixed=alpha.fixed))
    
    
    fillresult <- function(t, fmodelfull)
    {
      result  <-  NULL
      result$rand.table <- t$rand.table
      result$anova.table <- t$anova.table
      result$lsmeans.table <- t$lsmeans.table
      result$diffs.lsmeans.table <- t$diffs.lsmeans.table
      ### calculate p adjusted  
      if("elim.num" %in% colnames(t$anova.table))
        final.facs <- rownames(t$anova.table)[t$anova.table[,"elim.num"]=="kept"]
      else
        final.facs <- rownames(t$anova.table)
      rnames <- rownames(t$diffs.lsmeans.table)
      diffs.facs <- sapply(rnames, function(x) 
        substring(x, 1, substring.location(x, " ")$first[1]-1), USE.NAMES = FALSE)
      pv.adjust <- rep(0, length(diffs.facs))
      for(i in 1:length(final.facs))
      {
        find.fac <- diffs.facs %in% final.facs[i]
        pv.adjust[find.fac] <- p.adjust(t$diffs.lsmeans.table$"p-value"[find.fac], 
                                        method="bonferroni")
      }
      result$diffs.lsmeans.table$"p-value.adjust" <- pv.adjust
      
      res <- residuals(t$model)
      result$residuals <- res
      
      result$residualsFixed <- calcResidSaturFixedCons(data, fmodelfull, random)

      #format p-values
      if(!is.null(result$rand.table))
             result$rand.table[,which(colnames(result$rand.table)=="p.value")] <- 
        format.pval(result$rand.table[,which(colnames(result$rand.table)=="p.value")], 
                    digits=3, eps=1e-3)
      if(class(t$model)!="lm")
        result$anova.table[,which(colnames(result$anova.table)=="Pr(>F)")] <- 
        format.pval(result$anova.table[,which(colnames(result$anova.table)=="Pr(>F)")], 
                    digits=3, eps=1e-3)
      
      return(result)
    }
    
    resultFULL[[ind.resp]] <- fillresult(t, fmodelfull)
  }
  
  names(resultFULL) <- response
  return(resultFULL)
}
