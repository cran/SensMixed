################################################################################
## LSMEANS and DIFFLSMEANS related functions
################################################################################

lsmeans.calc <- function(model, alpha, test.effs = NULL, 
                         lsmeansORdiff = TRUE, ddf = "Satterthwaite") {
  rho <- list() ## environment containing info about model
  rho <- rhoInit(rho, model, TRUE) ## save lmer outcome in rho envir variable
  rho$A <- calcApvar(rho) ## asymptotic variance-covariance matrix for theta and sigma
  result <- list(summ.data = NA, response = NA)
  result$summ.data <- calcLSMEANS_LT(rho, alpha, 
                             test.effs,
                             lsmeansORdiff, ddf)$summ.data
  result$response <- rownames(attr(terms(rho$model),"factors"))[1]
  result
}

################################################################################
## calculate LSMEANS DIFFS and CI for all effects
################################################################################
calcLSMEANS_LT <- function(rho, alpha, test.effs = NULL, 
                        lsmeansORdiff = TRUE, ddf = "Satterthwaite") {  
  m <- refitLM_LT(rho$model)
  effs <- attr(terms(m),"term.labels")
  if(!is.null(test.effs))
    effs <- effs[effs %in% test.effs]
  dclass <- attr(terms(m),"dataClasses")
  facs <- names(dclass[which(dclass=="factor")])
  ## get standard deviation of random parameters from model
  std.rand <- c(unlist(lapply(VarCorr(rho$model), function(x) attr(x,"stddev"))), 
                attr(VarCorr(rho$model), "sc"))^2 #as.numeric(rho$s@REmat[,3])
  
  dd <- model.frame(rho$model) 
  
  ## init lsmeans summary
  if(lsmeansORdiff) {
    lsmeans.summ <-  matrix(ncol=length(facs)+7,nrow=0)
    colnames(lsmeans.summ) <- c(facs,"Estimate","Standard Error", "DF", 
                                "t-value", "Lower CI", "Upper CI", "p-value")
    summ.data <- as.data.frame(lsmeans.summ)
  } else {
    ## init diff summary
    diff.summ <-  matrix(ncol=7,nrow=0)
    colnames(diff.summ) <- c("Estimate","Standard Error", "DF", "t-value", 
                             "Lower CI", "Upper CI", "p-value")
    summ.data <- as.data.frame(diff.summ)
  }
  
  for(eff in effs) {
    split.eff  <-  unlist(strsplit(eff,":"))
    if(checkAllCov(split.eff, dd))
      next
    mat  <-  popMatrix(m, split.eff)
    fac.comb <- getFacCombForLSMEANS(split.eff, dd)  
    if(!lsmeansORdiff)
      summ.data <- rbind(summ.data, calcDiffsForEff_LT(
        facs, fac.comb, split.eff, eff, effs, rho, alpha, mat, ddf))
    else
      summ.data <- rbind(summ.data, calcLsmeansForEff_LT(
        lsmeans.summ, fac.comb, eff, split.eff, alpha, mat, rho, facs, ddf))
  }
  return(list(summ.data = summ.data))
}


## calculate DIFFERENCES OF LSMEANS and STDERR for effect
calcDiffsForEff_LT <- function(facs, fac.comb, split.eff, eff, effs, rho, 
                            alpha, mat, ddf) {
  ## calculating diffs for 2 way interaction
  if(length(split.eff) >= 1 && length(split.eff) <= 2) {   
    if(length(split.eff) == 2) {            
      fac.comb.names <- concatLevs(fac.comb)
      main.eff <- effs[effs %in% split.eff]
      
      mat.names.diffs <- combn(fac.comb.names,2)
      mat.nums.diffs <- apply(mat.names.diffs, c(1, 2), 
                              function(x) which(fac.comb.names == x))
    } else {
      mat.names.diffs <- combn(fac.comb, 2)
      mat.nums.diffs <- apply(mat.names.diffs, c(1, 2), 
                              function(x) which(fac.comb == x))
    }
    mat.diffs <- matrix(0, nrow=ncol(mat.nums.diffs), ncol=ncol(mat))
    colnames(mat.diffs) <- colnames(mat)
    for(ind.diffs in 1:ncol(mat.nums.diffs)) {
      mat.diffs[ind.diffs,] <- mat[mat.nums.diffs[1, ind.diffs], ] - 
        mat[mat.nums.diffs[2, ind.diffs], ]
    }
    names.combn <- apply(mat.names.diffs, 2, 
                         function(x) paste(x[1], x[2], sep=" - "))
    rownames(mat.diffs) <- paste(eff, names.combn)
    
    diffs.summ <-  matrix(NA, ncol=7, nrow=nrow(mat.diffs))
    colnames(diffs.summ) <- c("Estimate","Standard Error", "DF", "t-value", 
                              "Lower CI", "Upper CI", "p-value")
    rownames(diffs.summ) <- rownames(mat.diffs)
    
    diffs.summ <- as.data.frame(fillLSMEANStab_LT(mat.diffs, rho, diffs.summ, 0, 
                                               alpha, ddf))
    return(roundLSMEANStab(diffs.summ, 0))
  }
}


## calculate LSMEANS and STDERR for effect
calcLsmeansForEff_LT <- function(lsmeans.summ, fac.comb, eff, split.eff, alpha, mat, 
                              rho, facs, ddf) {
  summ.eff <- matrix(NA, ncol=ncol(lsmeans.summ), nrow=nrow(fac.comb))
  colnames(summ.eff) <- colnames(lsmeans.summ)
  #rownames(summ.eff) <- rep(eff, nrow(fac.comb))
  summ.eff[,split.eff] <- fac.comb
  names.arg <- concatLevs(summ.eff[,split.eff])
  summ.eff <- as.data.frame(fillLSMEANStab_LT(mat, rho, summ.eff, length(facs), 
                                           alpha, ddf))
  summ.eff <- convertFacsToNum(summ.eff, length(facs)+1, ncol(summ.eff))
  summ.eff <- roundLSMEANStab(summ.eff, length(facs))
  rownames(summ.eff) <- paste(rep(eff, nrow(fac.comb)), names.arg)
  return(summ.eff) 
}

###################################################################
#fills the LSMEANS and DIFF summary matrices
###################################################################
fillLSMEANStab_LT <- function(mat, rho, summ.eff, nfacs, alpha, 
                           ddf = "Satterthwaite") {
  newcln <- colnames(mat)[colnames(mat) %in% names(rho$fixEffs)]
  ## check estimability
  if(sum(!colnames(mat) %in% names(rho$fixEffs)) > 0){
    ids.est <- checkForEstim(mat, rho)
    mat <- mat[ids.est, , drop = FALSE]
  } else {
    ids.est <- 1:nrow(mat)
  }
  if(length(ids.est) == 0)
    return(summ.eff)
  mat <- matrix(mat[, colnames(mat) %in% names(rho$fixEffs)], nrow=nrow(mat), 
                ncol=length(newcln), dimnames=list(rownames(mat), newcln))
  estim.lsmeans <- mat %*% rho$fixEffs
  summ.eff[ids.est, nfacs+1] <- estim.lsmeans  
  
  ttest.res <- if(ddf == "Satterthwaite") {
    do.call(rbind, lapply(1:nrow(mat), function(i)
      calcSatterth1DF(rho=rho, L=mat[i, ], isF=FALSE)))
  } else {
    ttest.res <- do.call(rbind, lapply(1:nrow(mat), function(i)
      calcKR1DF(rho=rho, L=mat[i, ])))
  }
  if(is.vector(ttest.res))
    ttest.res <- t(as.matrix(ttest.res))
  
  summ.eff[ids.est, nfacs+2] <- ttest.res[, 4]#stdErrLSMEANS(rho, std.rand, mat)
  #df
  summ.eff[ids.est, (nfacs+3)] <- ttest.res[, 1]
  #t values
  summ.eff[ids.est, (nfacs+4)] <- ttest.res[, 2]
  #p values
  summ.eff[ids.est, (nfacs+7)] <- ttest.res[, 3]
  # CIs
  summ.eff[ids.est, nfacs+5] <- estim.lsmeans - 
    abs(qt(alpha/2, ttest.res[,1])) * ttest.res[, 4]
  summ.eff[ids.est ,nfacs+6] <- estim.lsmeans + 
    abs(qt(alpha/2, ttest.res[,1])) * ttest.res[, 4]
  return(summ.eff)
}

## refit model to lm in order to use popMatrix function
refitLM_LT <- function(obj) {
  mm <- model.matrix(obj)
  l.lmerTest.private.contrast <- attr(mm,"contrasts")
  contr <- l.lmerTest.private.contrast
  mm2 <- model.frame(obj)
  colnames(mm2)[1] <- "y"
  fo <- getFormula(obj, withRand = FALSE) ## formula(obj,fixed.only=TRUE)
  #l.lmerTest.private.contrast <- attr(mm, "contrasts")
  #contr <- l.lmerTest.private.contrast
  
  ## change contrasts for F tests calculations
  ## list of contrasts for factors
  if(length(which(unlist(contr)!="contr.SAS")) > 0)
  {
    names.facs <- names(contr)
    l.lmerTest.private.contrast <- as.list(rep("contr.SAS", length(names.facs)))
    names(l.lmerTest.private.contrast) <- names(contr)
  }
  # if(fo != as.formula(.~.))
  # {
  #   inds <-  names(l.lmerTest.private.contrast) %in% attr(terms(fo), "term.labels")
  #   ## update contrast l
  #   l.lmerTest.private.contrast <- l.lmerTest.private.contrast[inds]
  # }
  fo <- update(fo, y ~ .)
  lm(fo, data=mm2, contrasts = l.lmerTest.private.contrast)
}
