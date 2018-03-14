testOST <- FALSE

library(SensMixed)
load(system.file("testdata","bugResiddat.RData",package="SensMixed"))
load(system.file("testdata","bb.RData",package="SensMixed"))
load(system.file("testdata","ost.RData",package="SensMixed"))

response <- c("liking")
fixed <- list(Product = c("Sugar", "Acid"))
random <- c("consumer")
facs <- c("consumer", "Sugar", "Acid")

res <- conjoint(structure = 3, dat, response, fixed, random, facs)
names(res)
names(res$liking)
length(res$liking$residualsFixed)
dim(dat)
dat.resid <- cbind(dat, res$liking$residualsFixed)
colnames(dat.resid)[6] <- "doubresid"
dat.resid <- dat.resid[,-c(3,4,5)]
str(dat.resid)

if(requireNamespace("reshape2", quietly = TRUE)) {
  dr <- reshape2::dcast(dat.resid, prod ~ consumer)

  ## check 0 mean across columns and rows (residuals double-centered)
  stopifnot(all.equal(sum(rowSums(dr[,-1])),0)) ## across products OK
  stopifnot(all.equal(sum(colSums(dr[,-1])),0))
}


#########################################################################
## check BB data
response <- c("Liking")
fixed <- list(Product = c("Barley", "Salt"), Consumer = c("Sex", "Age"))
random <- c("Consumer")
facs<-c("Consumer", "Barley", "Salt", "Sex", "Age")

bb.N$Barley <- as.factor(bb.N$Barley)
bb.N$Age<- as.factor(bb.N$Age)
bb.N$Salt<- as.factor(bb.N$Salt)
bb.N$Sex<- as.factor(bb.N$Sex)
bb.N$Consumer<- as.factor(bb.N$Consumer)

res <- conjoint(structure = 1, bb.N, response, fixed, random, facs)

# m <- lmer(Liking ~ Barley + Salt + Sex + Age + (1|Barley:Consumer) + (1|Salt:Consumer)
#           + (1|Consumer), data = bb.N)
# 
# an <- anova(m)


TOL <- 1e-6
TOL2 <- 1e-2
## check the structure parameter
stopifnot(all.equal(rownames(res[[1]][[1]]),
                    c("Barley:Consumer", "Salt:Consumer", "Consumer")))

## check analysis of random effects
stopifnot(all.equal(res[[1]][[1]][,3], c("<0.001","0.45","1.00")))


# ## check anova table
# stopifnot(all.equal(res[[1]][[2]][,5], an[,5], tol=TOL))

## check degrees of freedom
stopifnot(all.equal(round(res[[1]][[2]][,4]), c(170,111,182,182)))


## check bonferroni correction
stopifnot(all.equal(res[[1]][[4]][,8], c(0.0000,0.0004,0.0061,
                                         1,1,1,0.5720,1,1,1,1,0.8270,1)))



## check that covariate is converted to a factor (CC 1.3.3 cannot handle covariates)
bb.N$Age <- as.numeric(as.character(bb.N$Age))

res2 <- conjoint(structure = 1, bb.N, response, fixed, random, facs)

stopifnot(all.equal(res[[1]][[2]], res2[[1]][[2]]))

## check structure 2
bb.N$Age <- as.factor(bb.N$Age)
res2 <- conjoint(structure = 2, bb.N, response, fixed, random, facs)

# m <- lmer(Liking ~ (Barley + Salt + Sex + Age)^2 - Sex:Age + (1|Barley:Consumer) + (1|Salt:Consumer)
#           + (1|Consumer),
#           data = bb.N)
# 
# m2 <- lmer(Liking ~ Barley + Salt + Sex + Age + Barley:Salt + Barley:Sex + Barley:Age + Salt:Sex + Salt:Age  + (1|Barley:Consumer) + (1|Salt:Consumer)
#            + (1|Consumer),
#            data = bb.N)
# 
# m3 <- lmer(Liking ~ Barley + (1 | Barley:Consumer) + Salt + (1 | Salt:Consumer) +
#              Sex + Age + Barley:Salt + Barley:Sex + Barley:Age + Salt:Sex +
#              Salt:Age + (1 | Consumer), data = bb.N)
# 
# stopifnot(all.equal(anova(m), anova(m3), tol=TOL))
# stopifnot(all.equal(anova(m), anova(m2), tol=TOL))
# 
# ## checked with SAS (testConjoint.r)
# stopifnot(all.equal(res2[[1]][[2]][-6], anova(m)[,-6], check.attributes = FALSE))

## check structure 3
res <- conjoint(structure = 3, bb.N, response, fixed, random, facs)


stopifnot(all.equal(res[[1]][[1]][,3], c("1","2","kept")))
stopifnot(all.equal(res[[1]][[1]][,4], c("0.864","0.388","<0.001")))

stopifnot(all.equal(rownames(res[[1]][[2]]), c("Barley:Salt:Sex","Salt:Sex","Barley:Salt:Age",
                                               "Salt:Age","Barley:Salt","Barley:Sex","Barley:Age",
                                               "Age","Barley","Salt","Sex")))


######################################################################################
## check ost data (checks also the lack of consumer data)

if(testOST){
  response <- c("liking")
  fixed <- list(Product = c("Pasteur", "Package", "Organic", "Omega", "Price"))
  random <- c("consumer")
  facs <- c("consumer", "Pasteur", "Package", "Organic", "Omega", "Price")

  res <- conjoint(structure = 2, ost, response, fixed, random, facs)


  (randTab <- res[[1]][1])
  (anovaTab <- res[[1]][2])
  (lsmeansTab <- res[[1]][3])
  (diffLsTab <- res[[1]][4])
  resid <- res[[1]][5]
  resid.ind <- res[[1]][6]

  # ## check ANOVA with SAS
  # m.ost <- lmer(liking ~ (Pasteur + Package + Organic + Omega + Price)^2 + (1|Pasteur:consumer)+
  #                 (1|Package:consumer) + (1|Organic:consumer) + (1|Omega:consumer) + (1|Price:consumer) +
  #                 (1|Organic:Omega:consumer), data = ost)
  # 
  # an <- anova(m.ost)
  # stopifnot(all.equal(an[,5], c(15.22, 0.05, 0.17, 1.45, 87.49, NA, NA,
  #                               17.45, NA,  3.19, 34.855, 3.31), tol=TOL2))
}


######################################################################################
## check with the  ham
response <- c("Liking")
fixed <- list(Product=c("Product", "Information"), Consumer="Gender")
random <- c("Consumer")
facs <- c("Consumer", "Product", "Information", "Gender")

#######################
## structure 1
res.ham <- conjoint(structure=1, ham, response, fixed, random, facs)

## check nemas of the random effects
stopifnot(all.equal(rownames(res.ham[[1]][[1]]), c("Product:Consumer",
                                                   "Information:Consumer", "Consumer")))
## check names of the fixed effects
stopifnot(all.equal(rownames(res.ham[[1]][[2]]), c("Product", "Information",
                                                   "Gender")))
## check F test with SAS
stopifnot(all.equal(res.ham[[1]][[2]][,5], c(3.83, 3.34,0.88), tol=TOL2))

## check ddf with SAS
stopifnot(all.equal(round(res.ham[[1]][[2]][,4]), c(240,80,79)))

## check difflsmeans for the product effect
stopifnot(all.equal(res.ham[[1]][[4]][1:6,1], c(0.7037, -0.2840, -0.1173, -0.9877,
                                                -0.8210, 0.1667), tol=TOL))
## check bonferroni correction for the product effect
stopifnot(all.equal(res.ham[[1]][[4]][1:6,8], c(0.156,1,1,0.0112,0.0571,1), tol=TOL2))



## check which LRT test (REML or ML)
ham$Consumer <-as.factor(ham$Consumer)
ham$Product <-as.factor(ham$Product)
ham$Information <-as.factor(ham$Information)
ham$Gender <-as.factor(ham$Gender)


# m.ham <- lmer(Liking ~ Product + Information + Gender + (1|Consumer) + (1|Information:Consumer)
#               + (1|Product:Consumer), data=ham)
# m.ham.red <- lmer(Liking ~ Product + Information + Gender  + (1|Information:Consumer)
#                   + (1|Product:Consumer), data=ham)
# 
# ## test Consumer effect
# anova(m.ham, m.ham.red, refit = FALSE) ## conjoint now usis LRT based on REML models!!!





########################
## check structure=2
res.ham2 <- conjoint(structure=2, ham, response, fixed, random, facs)

## check nemas of the random effects
stopifnot(all.equal(rownames(res.ham2[[1]][[1]]), c("Product:Consumer",
                                                    "Information:Consumer", "Consumer")))
## check names of the fixed effects
stopifnot(all.equal(rownames(res.ham2[[1]][[2]]), c("Product", "Information",
                                                    "Gender", "Product:Information",
                                                    "Product:Gender","Information:Gender")))
## check F test with SAS
stopifnot(all.equal(res.ham2[[1]][[2]][,5], c(3.82, 3.29,0.88,2.2,0.35,0.72), tol=TOL2))


########################
## check structure=3
res.ham3 <- conjoint(structure=3, ham, response, fixed, random, facs)


stopifnot(all.equal(res.ham3[[1]][[2]][,6], c("1","2","3","4","5","kept","kept")))


##########################################################################################
#check with the  TVbo data (multiple responses)
response <- c("Coloursaturation", "Colourbalance")
fixed <- list(Product=c("TVset", "Picture"))
random <- c("Assessor")
facs <- c("Assessor", "TVset", "Picture")

res.tv <- conjoint(structure=3, TVbo, response, fixed, random, facs)


stopifnot(all.equal(names(res.tv), c("Coloursaturation","Colourbalance")))

