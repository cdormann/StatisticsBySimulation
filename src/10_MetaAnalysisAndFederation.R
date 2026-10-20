## Meta-analysis, federated analysis and its context

#### study set up ####
# 20 predictors, 10 relevant, with some important interactions thrown in, N=20000, correlation in [-0.5, 0.5] (or rnorm of that); 20 studies; sample size U(100, 1000), without replacement; in local model, only provide p predictors, sampled with replacement, with 50-60\% overlap (or full overlap, for comparison); stepwise regression model ($X_1$ in all models, not selected out; vs. pure dredging, and reports without $X_1$ are not in the meta-analysis data; vs. all negative effects are excluded from publications/funding); none/overlap x $X_1$ not/kept/non-negative; plot with truth as reference line; give values for SE of meta-analysis result; use random effect MA

set.seed(123)
S <- matrix(rnorm(20*20, 0, 0.21), 20, 20)
Sigma <- S %*% t(S) # to ensure positive definite covariance matrix
range(cov2cor(Sigma)[lower.tri(Sigma)])
library(mvtnorm)
X <- rmvnorm(20000, mean=rep(0, 20), sigma=Sigma)
X <- scale(X) # scale to get standardised regression coefficients as effect size
range(cor(X)[lower.tri(cor(X))])
mean(cor(X)[lower.tri(cor(X))])
sd(cor(X)[lower.tri(cor(X))])
slopes <- runif(20, -2, 2)
Y <- X %*% slopes + rnorm(nrow(X), 0, 5)

fullData <- cbind.data.frame(Y, X) 
colnames(fullData) <- c("Y", paste0("X", 1:20))

# summary(lm(Y ~ X))
coef(lm(Y ~ X))[-1] - slopes


Nstudies <- 21
Ns <- sample(10:100, Nstudies) * 10 # between 100 and 1000
index <- sample(20000, sum(Ns)) # pick data from full simulation

studyList <- list()
for (i in 1:Nstudies){
  studyList[[i]] <- fullData[index[1:Ns[i]], sort(c(1, 2, 3, sample(4:21, 8)))] # pick 10 predictors, 8 randomly
  index <- index[-c(1:Ns[i])] # remove data already assigned
}


#### meta-analysis ####
# local analyses, always keeping X1 and X2
fmList <- lapply(studyList, \(x) step(lm(Y ~ ., data=x), 
                                      scope=list("lower"=Y~X1+X2, "upper"=Y ~ .), 
                                      k=log(nrow(x)), trace=F))

# compute effect size = standardised regression coefficient
effectSizes <- matrix(0, Nstudies, 20)
colnames(effectSizes) <- colnames(fullData)[-1]
for (i in 1:Nstudies){
  cs <- summary(fmList[[i]])$coefficients[-1, 1] 
  effectSizes[i, names(cs)] <- cs
}
# effectSizes


# using only those studies that actually report on an effect:
effectSizesNAs <- effectSizes
effectSizesNAs[effectSizesNAs==0] <- NA        
effectSizesNAs.wmeans <- apply(effectSizesNAs, 2, \(x) weighted.mean(x, w=Ns, na.rm=T)) # sample-size weighted mean
effectSizes.means <- apply(effectSizes, 2, \(x) weighted.mean(x, w=Ns, na.rm=T)) # weighted mean
#plot(colMeans(effectSizesNAs, na.rm=T), slopes)
#abline(0,1)
#summary(lm(slopes ~ colMeans(effectSizesNAs, na.rm=T) ) )  # -0.05, 0.54 : underestimationg of effects when using missing as NA


# using all studies, and setting missing predictors to 0:
plot(colMeans(effectSizes), slopes)
abline(0,1)
summary(lm(slopes ~ colMeans(effectSizes) ) ) # -0.01, 1.35 : overestimated with 0s for missing

biasNAs <- effectSizesNAs.wmeans - slopes
biasAll <- effectSizes.means - slopes

relbiasNAs <- (effectSizesNAs.wmeans - slopes) / slopes
relbiasAll <- (effectSizes.means - slopes)/ slopes

mean(abs(biasNAs))
mean(abs(biasAll))

dats <- data.frame(estimate=c(biasNAs, biasAll), 
                   inclusion=rep(c("reported", "all"), each=20),
                   coefficient=factor(rep(colnames(effectSizes), times=2), levels=paste0("X", 1:20)))
library(ggplot2)
pdf(file="../figures/metaAnalysis1.pdf", width=7, height=5)
ggplot(data=dats, aes(coefficient, estimate, fill=inclusion)) + 
  geom_col(position=position_dodge()) + 
  ylab("bias of estimate") +
  scale_fill_manual(values=c("grey20", "grey70")) + 
  theme_minimal() + theme(legend.position = c(0.8, 0.9))
dev.off()



## repeat many times ##
library(mvtnorm)
Nstudies <- 21

set.seed(123)
effectSizes.list <- list()
for (k in 1:100){ # repetitions
  S <- matrix(rnorm(20*20, 0, 0.21), 20, 20)
  Sigma <- S %*% t(S) # to ensure positive definite covariance matrix
  X <- rmvnorm(20000, mean=rep(0, 20), sigma=Sigma)
  X <- scale(X) # scale to get standardised regression coefficients as effect size
  slopes <- runif(20, -2, 2)
  Y <- X %*% slopes + rnorm(nrow(X), 0, 5)
  
  fullData <- cbind.data.frame(Y, X) 
  colnames(fullData) <- c("Y", paste0("X", 1:20))
  
  Ns <- sample(10:100, Nstudies) * 10 # between 100 and 1000
  index <- sample(20000, sum(Ns)) # pick data from full simulation
  
  studyList <- list()
  for (i in 1:Nstudies){
    studyList[[i]] <- fullData[index[1:Ns[i]], sort(c(1, 2, 3, sample(3:21, 8)))] # pick 10 predictors, 8 randomly
    index <- index[-c(1:Ns[i])] # remove data already assigned
  }
  
  # local analyses, always keeping X1 and X2
  fmList <- lapply(studyList, \(x) step(lm(Y ~ ., data=x), 
                                        scope=list("lower"=Y~X1+X2, "upper"=Y ~ .), 
                                        k=log(nrow(x)), trace=F))
  
  # compute effect size = standardised regression coefficient
  effectSizes <- matrix(0, Nstudies, 20)
  colnames(effectSizes) <- colnames(fullData)[-1]
  for (i in 1:Nstudies){
    cs <- summary(fmList[[i]])$coefficients[-1, 1] 
    effectSizes[i, names(cs)] <- cs
  }
  
  effectSizesNAs <- effectSizes
  effectSizesNAs[effectSizesNAs==0] <- NA 
  effectSizesNAs.wmeans <- apply(effectSizesNAs, 2, \(x) weighted.mean(x, w=Ns, na.rm=T)) # sample-size weighted mean
  effectSizes.means <- apply(effectSizes, 2, \(x) weighted.mean(x, w=Ns, na.rm=T)) # weighted mean
  
  
  relbiasNAs <- (effectSizesNAs.wmeans  - slopes) / slopes
  relbiasAll <- (effectSizes.means - slopes) / slopes
  
  effectSizes.list[[k]] <- list(effectSizes, relbiasNAs, relbiasAll)
  
  print(k)
  
}
effectSizes.list
# save(effectSizes.list, file="effectSizeslist.Rdata")

relbiasNAs100 <- sapply(effectSizes.list, \(x) x[[2]])
relbiasAll100 <- sapply(effectSizes.list, \(x) x[[3]])

dats <- cbind(rbind(stack(as.data.frame(t(relbiasNAs100))), 
                    stack(as.data.frame(t(relbiasAll100)))), 
              "NAhandling"=rep(c("missing=NA", "missing=0"), each=20*100))
#summary(dats)

library(ggplot2)
pdf(file="../figures/metaAnalysis_relbias.pdf", width=7, height=5)
ggplot(dats, aes(y=values, x=ind, fill=NAhandling)) + 
  geom_hline(yintercept=0, col="grey") +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values=c("grey40", "grey70")) + 
  coord_cartesian(ylim = c(-3, 5)  ) + xlab(label="predictor") + ylab(label="relative bias") +
  theme_minimal() + theme(legend.title=element_blank())
dev.off()       


##### MRP excursion #####
str(fmList, 1)
fmList[[1]]
# focus: coefficient for X1:
ds <- sapply(fmList, \(x) summary(x)$coefficients[2,1])
ses <- sapply(fmList, \(x) summary(x)$coefficients[2,2])
# plot(ds/ses ~ I(1/ses))
lm(ds/ses ~ I(1/ses) - 1)
UWLS <- unname(coef(lm(ds/ses ~ I(1/ses) - 1)))
(1 - pnorm(1.96 - abs(UWLS)/median(ses))) *100
# 52.5

# focus: X2:
ds <- sapply(fmList, \(x) summary(x)$coefficients[3,1])
ses <- sapply(fmList, \(x) summary(x)$coefficients[3,2])
UWLS <- unname(coef(lm(ds/ses ~ I(1/ses) - 1)))
(1 - pnorm(1.96 - abs(UWLS)/median(ses))) *100
#99.99




#### IPL-meta-analysis ####
# re-run lines 1-35 (data simulation)

#allStudies <- matrix(NA, nrow=sum(Ns), ncol=22)
#colnames(allStudies) <- c(colnames(fullData), "ID")
#Nscumsum <- c(0, cumsum(Ns))
#for (i in 1:21){
#  allStudies[(Nscumsum[i]+1):Nscumsum[i+1], pmatch(colnames(studyList[[i]]), colnames(allStudies))] <- as.matrix(studyList[[i]])
#}
#allStudies <- as.data.frame(allStudies, stringsAsFactor=T)
#allStudies$ID <- rep(LETTERS[1:21], times=Ns)
#summary(allStudies)

allStudies <- fullData[1:sum(Ns), ]
allStudies$ID <- rep(LETTERS[1:21], times=Ns)

library(lmerTest)
f <- as.formula(paste("Y~", paste0("X", 1:20, collapse=" + "), "+(1|ID)"))
fmer <- lmer(f, data=allStudies)
summary(fmer)
summary(fmer)$coefficients[-1, 1]
plot(summary(fmer)$coefficients[-1, 1], slopes)
abline(0,1)

# in the spirit of previous step-selection, here the step version for lmer:
fmer.step <- step(fmer, reduce.random=FALSE)
fmer.step
# well-hidden formula:
fs <- attr(attr(fmer.step, "drop1"), "heading")[3]
summary(lmer(fs, data=allStudies))


###### many replications ######
set.seed(123)
Nreps <- 100
IPLMAsimu <- replicate(Nreps, {
  S <- matrix(rnorm(20*20, 0, 0.21), 20, 20)
  Sigma <- S %*% t(S) # to ensure positive definite covariance matrix
  library(mvtnorm)
  X <- rmvnorm(20000, mean=rep(0, 20), sigma=Sigma)
  X <- scale(X) # scale to get standardised regression coefficients as effect size
  slopes <- runif(20, -2, 2)
  Y <- X %*% slopes + rnorm(nrow(X), 0, 5)
  
  fullData <- cbind.data.frame(Y, X) 
  colnames(fullData) <- c("Y", paste0("X", 1:20))
  
  Nstudies <- 21
  Ns <- sample(10:100, Nstudies) * 10 # between 100 and 1000
  allStudies <- fullData[1:sum(Ns), ]
  allStudies$ID <- rep(LETTERS[1:21], times=Ns)
  f <- as.formula(paste("Y~", paste0("X", 1:20, collapse=" + "), "+(1|ID)"))
  fmer <- lmer(f, data=allStudies)
  (summary(fmer)$coefficients[-1, 1] - slopes) / slopes
})
#save(IPLMAsimu, file="IPLMAsimu.Rdata")
boxplot(t(IPLMAsimu), ylim=c(-2, 2))
abline(h=0)


# plot: MA vs IPL-MA
load(file="IPLMAsimu.Rdata")
load(file="effectSizeslist.Rdata") # (effectSizes, relbiasNAs, relbiasAll)
relbiasNAs100 <- sapply(effectSizes.list, \(x) x[[2]])

dats <- cbind(rbind(stack(as.data.frame(t(relbiasNAs100))), 
                    stack(as.data.frame(t(IPLMAsimu)))), 
              "MAtype"=rep(c("meta analysis", "IPL-MA"), each=20*100))
#summary(dats)

library(ggplot2)
pdf(file="../figures/metaAnalysisVsIPLMA_relbias.pdf", width=7, height=5)
ggplot(dats, aes(y=values, x=ind, fill=MAtype)) + 
  geom_hline(yintercept=0, col="grey") +
  geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(-3, 4)  ) + 
  scale_fill_manual(values=c("grey40", "grey70")) + 
  xlab(label="predictor") + ylab(label="relative bias") +
  theme_minimal() + theme(legend.title=element_blank())
dev.off()       



#### one-step federated analysis ####
# like meta-analysis, but imposing a single model proposed by one of the groups based on their results
# same split as in meta-analysis; but same model for all, more complicated model than in meta-analysis; pick largest (or random?) site as model generator, use that model for everybody else (dropping missing predictors from the model); average across model estimates, using 0 for missing predictors in a model (as a shrinkage; possibly try NA otherwise); use random effect

##### identify common model #####
# re-run lines 1-35;

set.seed(123)
coordinatingSite <- sample(length(studyList), 1)

fmCommon <- step(lm(Y ~ (.)^2, data=studyList[[coordinatingSite]]), 
                scope=list("lower"=Y~X1+X2, "upper"=Y ~ (.)^2), k=log(nrow(studyList[[coordinatingSite]])), trace=F)

formulaCommon <- fmCommon$call$formula
predictorsCommon <- attr(terms(fmCommon) , "term.labels")
#predictorsCommon <- predictorsCommon[-grep(":", predictorsCommon)] # remove interactions

fullPredictorSetAvailable <- rep(FALSE, length(studyList)) # list of studies that contain all predictors of the common model
onestep.coefficients.mat <- matrix(NA, nrow=length(studyList), ncol=length(predictorsCommon))
colnames(onestep.coefficients.mat) <- predictorsCommon
for (i in seq_along(studyList)){
  onOffer <- colnames(model.matrix(~(.)^2, data=studyList[[i]][1:2,-1]))
  studyPredictors <- intersect(onOffer, predictorsCommon)
  if (length(studyPredictors) == length(predictorsCommon))  fullPredictorSetAvailable[i] <- TRUE
  fmone <- lm(formula(paste("Y~", paste(studyPredictors, collapse="+"))), data=studyList[[i]])
  onestep.coefficients.mat[i, match(studyPredictors, colnames(onestep.coefficients.mat))] <- coef(fmone)[-1]
  rm(fmone, onOffer, studyPredictors)
}

onestep.coefficients.mat
# final result: same figure as before, but comparing the few studies with full set of predictors vs. all
# use all studies:
apply(onestep.coefficients.mat, 2, \(x) weighted.mean(x, w=Ns, na.rm=T)) 
# use only studies with full set of predictors:
apply(onestep.coefficients.mat[fullPredictorSetAvailable,, drop=F], 2, \(x) weighted.mean(x, w=Ns[fullPredictorSetAvailable], na.rm=T)) 


##### now do this many times: #####
library(mvtnorm)
Nstudies <- 21
MM <- model.matrix(~(.)^2, data=fullData[1:2,-1]) # for names

reps <- 100
onestep.simu <- array(NA, dim=c(2, 20, reps)) # all/full x coefs x reps
set.seed(12)
for (j in 1:reps){
  ###### set up the data and split over studies: ######
  S <- matrix(rnorm(20*20, 0, 0.21), 20, 20)
  Sigma <- S %*% t(S) # to ensure positive definite covariance matrix
  X <- rmvnorm(20000, mean=rep(0, 20), sigma=Sigma)
  X <- scale(X) # scale to get standardised regression coefficients as effect size
  slopes <- runif(20, -2, 2)
  Y <- X %*% slopes + rnorm(nrow(X), 0, 5)
  fullData <- cbind.data.frame(Y, X) 
  colnames(fullData) <- c("Y", paste0("X", 1:20))
  Ns <- sample(10:100, Nstudies) * 10 # between 100 and 1000
  index <- sample(20000, sum(Ns)) # pick data from full simulation
  
  studyList <- list()
  for (i in 1:Nstudies){
    studyList[[i]] <- fullData[index[1:Ns[i]], sort(c(1, 2, 3, sample(4:21, 8)))] # pick 10 predictors, 8 randomly
    index <- index[-c(1:Ns[i])] # remove data already assigned
  }
  
  ###### identify common model based on one site ######
  coordinatingSite <- sample(length(studyList), 1)
  fmCommon <- step(lm(Y ~ (.)^2, data=studyList[[coordinatingSite]]), 
                   scope=list("lower"=Y~X1+X2, "upper"=Y ~ . + (.)^2), 
                   k=log(nrow(studyList[[coordinatingSite]])), trace=F)
  formulaCommon <- fmCommon$call$formula
  predictorsCommon <- attr(terms(fmCommon) , "term.labels")
  
  ###### run common model (or subset) in all sites ######
  fullPredictorSetAvailable <- rep(FALSE, length(studyList)) # list of studies that contain all predictors of the common model
  
  onestep.coefficients.mat <- matrix(NA, nrow=length(studyList), ncol=ncol(MM)-1) #without intercept
  colnames(onestep.coefficients.mat) <- colnames(MM)[-1]
  
  for (i in seq_along(studyList)){
    onOffer <- colnames(model.matrix(~(.)^2, data=studyList[[i]][1:2,-1]))
    studyPredictors <- intersect(onOffer, predictorsCommon)
    if (length(studyPredictors) == length(predictorsCommon))  fullPredictorSetAvailable[i] <- TRUE
    fmone <- lm(formula(paste("Y~", paste(studyPredictors, collapse="+"))), data=studyList[[i]])
    # print(coef(fmone))
    onestep.coefficients.mat[i, match(studyPredictors, colnames(onestep.coefficients.mat))] <- coef(fmone)[-1] #without intercept
    rm(fmone, onOffer, studyPredictors)
  }
  
  ###### return relbias of estimates ######
  onestep.simu[,,j] <- rbind( 
    (apply(onestep.coefficients.mat[, 1:20], 2, \(x) weighted.mean(x, w=Ns, na.rm=T)) - slopes)/slopes, 
    (apply(onestep.coefficients.mat[fullPredictorSetAvailable, 1:20, drop=F], 2, \(x) weighted.mean(x, w=Ns[fullPredictorSetAvailable, drop=F], na.rm=T)) - slopes)/slopes)
  
  cat("Replicate number", j, "done!\n")
} # end for loop 

# or using replicate:
#onestep.simu <- replicate(reps, {
  ###### set up the data and split over studies: ######
  S <- matrix(rnorm(20*20, 0, 0.21), 20, 20)
  Sigma <- S %*% t(S) # to ensure positive definite covariance matrix
  X <- rmvnorm(20000, mean=rep(0, 20), sigma=Sigma)
  X <- scale(X) # scale to get standardised regression coefficients as effect size
  slopes <- runif(20, -2, 2)
  Y <- X %*% slopes + rnorm(nrow(X), 0, 5)
  fullData <- cbind.data.frame(Y, X) 
  colnames(fullData) <- c("Y", paste0("X", 1:20))
  Ns <- sample(10:100, Nstudies) * 10 # between 100 and 1000
  index <- sample(20000, sum(Ns)) # pick data from full simulation
  
  studyList <- list()
  for (i in 1:Nstudies){
    studyList[[i]] <- fullData[index[1:Ns[i]], sort(c(1, 2, 3, sample(4:21, 8)))] # pick 10 predictors, 8 randomly
    index <- index[-c(1:Ns[i])] # remove data already assigned
  }
  
  ###### identify common model based on one site ######
  coordinatingSite <- sample(length(studyList), 1)
  fmCommon <- step(lm(Y ~ (.)^2, data=studyList[[coordinatingSite]]), 
                   scope=list("lower"=Y~X1+X2, "upper"=Y ~ . + (.)^2), 
                   k=log(nrow(studyList[[coordinatingSite]])), trace=F)
  formulaCommon <- fmCommon$call$formula
  predictorsCommon <- attr(terms(fmCommon) , "term.labels")
  
  ###### run common model (or subset) in all sites ######
  fullPredictorSetAvailable <- rep(FALSE, length(studyList)) # list of studies that contain all predictors of the common model

  onestep.coefficients.mat <- matrix(NA, nrow=length(studyList), ncol=ncol(MM)-1) #without intercept
  colnames(onestep.coefficients.mat) <- colnames(MM)[-1]
  
  for (i in seq_along(studyList)){
    onOffer <- colnames(model.matrix(~(.)^2, data=studyList[[i]][1:2,-1]))
    studyPredictors <- intersect(onOffer, predictorsCommon)
    if (length(studyPredictors) == length(predictorsCommon))  fullPredictorSetAvailable[i] <- TRUE
    fmone <- lm(formula(paste("Y~", paste(studyPredictors, collapse="+"))), data=studyList[[i]])
    # print(coef(fmone))
    onestep.coefficients.mat[i, match(studyPredictors, colnames(onestep.coefficients.mat))] <- coef(fmone)[-1] #without intercept
    rm(fmone, onOffer, studyPredictors)
  }
  
  ###### return relbias of estimates ######
  rbind( 
    (apply(onestep.coefficients.mat[, 1:20], 2, \(x) weighted.mean(x, w=Ns, na.rm=T)) - slopes)/slopes, 
   (apply(onestep.coefficients.mat[fullPredictorSetAvailable, 1:20, drop=F], 2, \(x) weighted.mean(x, w=Ns[fullPredictorSetAvailable, drop=F], na.rm=T)) - slopes)/slopes)
  
})
onestep.simu  # all/full x coefs x reps
#save(onestep.simu, file="onestep.simu.RData")

rowMeans(onestep.simu[1, 1:10,], na.rm=T)
apply(onestep.simu[1, 1:10,], 1, median, na.rm=T)
apply(onestep.simu[2, 1:10,], 1, median, na.rm=T)

median(onestep.simu[1, 1:20,], na.rm=T)
median(onestep.simu[2, 1:20,], na.rm=T)

dats <- as.data.frame(t(cbind(onestep.simu[1,,], onestep.simu[2,,])))
dats$type <- rep(c("all", "full"), each=dim(onestep.simu)[1])
dats <- reshape(dats, direction="long", varying=1:20, times=colnames(dats)[1:20], v.names="relBias", timevar="predictor")
dats$predictor <- gsub("V","X", dats$predictor)
dats$predictor <- factor(dats$predictor, levels=paste0("X", 1:20))

#tail(dats)
library(ggplot2)
pdf(file="../figures/onestep_relbias.pdf", width=7, height=5)
ggplot(dats, aes(predictor, relBias, fill=type)) + 
  geom_hline(yintercept=0, col="grey") +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values=c("grey40", "grey70")) + 
  coord_cartesian(ylim = c(-3, 5)  ) + xlab(label="predictor") + ylab(label="relative bias") +
  theme_minimal() + theme(legend.title=element_blank())
dev.off()   





#### two-step federated analysis ####
# re-start R, then:
# re-run lines 1-35 (data simulation) !!
# yielding: studyList
# local analyses, always keeping X1 and X2; allow for interactions (which is not the case for the meta-analysis!)

fmList <- lapply(studyList, \(x) step(lm(Y ~ (.)^2, data=x), 
                                      scope=list("lower"=Y~X1+X2, "upper"=Y ~ .), 
                                      k=log(nrow(x)), trace=F))


##### variables present in ANY study #####
(frequencyOfVariables <- sort(table(unlist(sapply(fmList, \(x) names(x$coefficients)[-1])))))
commonThird <- which(frequencyOfVariables > 7)
predictorsCommon <- sort(names(frequencyOfVariables)[commonThird])
commonModel <- as.formula(paste("Y ~", paste(predictorsCommon, collapse="+")))
# # collect means across all studies that provide this variable:
# XColMeans <- lapply(studyList, colMeans)
# XMeans <- matrix(NA, ncol=21, nrow=21)
# colnames(XMeans) <- colnames(fullData)
# for (i in 1:21){
#   XMeans[i, pmatch(names(XColMeans[[i]]), colnames(XMeans))] <- XColMeans[[i]]
# }
# colMeans(XMeans, na.rm=T)

# run this common model on all studies (akin to one-step, just different model):
MM <- model.matrix(~(.)^2, data=fullData[1:2,-1]) # for names
twostep.coefficients.mat <- matrix(NA, nrow=length(studyList), ncol=ncol(MM)-1)
colnames(twostep.coefficients.mat) <- colnames(MM)[-1]
for (i in seq_along(studyList)){
  onOffer <- colnames(model.matrix(~(.)^2, data=studyList[[i]][1:2,-1]))
  studyPredictors <- intersect(onOffer, predictorsCommon)
  fmone <- lm(formula(paste("Y~", paste(studyPredictors, collapse="+"))), data=studyList[[i]])
  twostep.coefficients.mat[i, match(studyPredictors, colnames(twostep.coefficients.mat))] <- coef(fmone)[-1]
  rm(fmone, onOffer, studyPredictors)
}
twostep.coefficients.mat
apply(twostep.coefficients.mat, 2, \(x) weighted.mean(x, w=Ns, na.rm=T)) # not shown in book


##### twostep: repeat many times #####
reps <- 100
set.seed(123)
twostep.simu <- replicate(reps, {  # takes hours!!
  ###### set up the data and split over studies: ######
  S <- matrix(rnorm(20*20, 0, 0.21), 20, 20)
  Sigma <- S %*% t(S) # to ensure positive definite covariance matrix
  X <- rmvnorm(20000, mean=rep(0, 20), sigma=Sigma)
  X <- scale(X) # scale to get standardised regression coefficients as effect size
  slopes <- runif(20, -2, 2)
  Y <- X %*% slopes + rnorm(nrow(X), 0, 5)
  fullData <- cbind.data.frame(Y, X) 
  colnames(fullData) <- c("Y", paste0("X", 1:20))
  Ns <- sample(10:100, Nstudies) * 10 # between 100 and 1000
  index <- sample(20000, sum(Ns)) # pick data from full simulation
  
  studyList <- list()
  for (i in 1:Nstudies){
    studyList[[i]] <- fullData[index[1:Ns[i]], sort(c(1, 2, 3, sample(3:21, 8)))] # pick 10 predictors, 8 randomly
    index <- index[-c(1:Ns[i])] # remove data already assigned
  }
  
  ###### run local analyses ######
  fmList <- lapply(studyList, \(x) step(lm(Y ~ (.)^2, data=x), 
                                      scope=list("lower"=Y~X1+X2, "upper"=Y ~ .), 
                                      k=log(nrow(x)), trace=F))

  ###### variables present in ANY study ######
  (frequencyOfVariables <- sort(table(unlist(sapply(fmList, \(x) names(x$coefficients)[-1])))))
  commonThird <- which(frequencyOfVariables > 7)
  predictorsCommon <- sort(names(frequencyOfVariables)[commonThird])
  commonModel <- as.formula(paste("Y ~", paste(predictorsCommon, collapse="+")))

  # run this common model on all studies (akin to one-step, just different model):
  twostep.coefficients.mat <- matrix(NA, nrow=length(studyList), ncol=ncol(MM)-1)
  colnames(twostep.coefficients.mat) <- colnames(MM)[-1]
  for (i in seq_along(studyList)){
    onOffer <- colnames(model.matrix(~(.)^2, data=studyList[[i]][1:2,-1]))
    studyPredictors <- intersect(onOffer, predictorsCommon)
    fmone <- lm(formula(paste("Y~", paste(studyPredictors, collapse="+"))), data=studyList[[i]])
    twostep.coefficients.mat[i, match(studyPredictors, colnames(twostep.coefficients.mat))] <- coef(fmone)[-1]
    rm(fmone, onOffer, studyPredictors)
  }
  #twostep.coefficients.mat
  
  ###### return relbias of estimates ######
  (apply(twostep.coefficients.mat[, 1:20], 2, \(x) weighted.mean(x, w=Ns, na.rm=T)) - slopes)/slopes

})
twostep.simu  # all/full x coefs x reps
rowMeans(twostep.simu, na.rm=T)
apply(twostep.simu, 1, median, na.rm=T)

#save(twostep.simu, file="twostep.simu.Rdata")





#### full model analysis ####
# simulate publication bias!

set.seed(123)
# make 40 studies, remove those with non-significant effects for X1 or X2
Nstudies <- 40
Ns <- sample(10:100, Nstudies) * 10 # between 100 and 1000

S <- matrix(rnorm(20*20, 0, 0.21), 20, 20)
Sigma <- S %*% t(S) # to ensure positive definite covariance matrix
library(mvtnorm)
X <- rmvnorm(sum(Ns), mean=rep(0, 20), sigma=Sigma)
X <- scale(X) # scale to get standardised regression coefficients as effect size
slopes <- runif(20, -2, 2)
Y <- X %*% slopes + rnorm(nrow(X), 0, 5)

fullData <- cbind.data.frame(Y, X) 
colnames(fullData) <- c("Y", paste0("X", 1:20))

fullData$ID <- rep(paste0("study", 1:40), times=Ns)

index <- sample(sum(Ns), sum(Ns)) # pick data from full simulation

studyList <- list()
for (i in 1:Nstudies){
  studyList[[i]] <- fullData[index[1:Ns[i]], sort(c(1, 2, 3, sample(4:21, 8)))] # pick 10 predictors, 8 randomly
  index <- index[-c(1:Ns[i])] # remove data already assigned
}
sapply(studyList, dim)

# now find out which substudy would NOT have been published because of insignificant X1 or X2, hence representing our publication bias:
fmList <- lapply(studyList, \(x) step(lm(Y ~ ., data=x), 
                                      scope=list("lower"=Y~X1+X2, "upper"=Y ~ .), 
                                      k=log(nrow(x)), trace=F))
notSignif <- which(sapply(fmList, \(x) (any(summary(x)$coefficients[2:3, 4] > 0.05))))

# produce a second data set with only published studies:
publishedData <- fullData[!(fullData$ID %in% paste0("study", notSignif)), ]

# now two simple LMs, one for full, one for published:
flmfull <- step(lm(Y ~ ., data=fullData), 
            scope=list("lower"=Y~X1+X2, "upper"=Y ~ .), 
            k=log(nrow(fullData)), trace=F)                 

flmpublished <- step(lm(Y ~ ., data=publishedData), 
                scope=list("lower"=Y~X1+X2, "upper"=Y ~ .), 
                k=log(nrow(publishedData)), trace=F)                 

coef(flmfull)
coef(flmpublished)

mean( (coef(flmfull)[-1] - slopes[as.numeric(gsub("X", "", names(coef(flmfull)[-1])))])/slopes[as.numeric(gsub("X", "", names(coef(flmfull)[-1])))] )
mean( (coef(flmpublished)[-1] - slopes[as.numeric(gsub("X", "", names(coef(flmpublished)[-1])))])/slopes[as.numeric(gsub("X", "", names(coef(flmpublished)[-1])))] )

##### repeat many times #####
reps <- 100

set.seed(123)
full.simu <- replicate(reps, {
  # make 40 studies, remove those with non-significant effects for X1 or X2
  Nstudies <- 40
  Ns <- sample(10:100, Nstudies) * 10 # between 100 and 1000
  
  S <- matrix(rnorm(20*20, 0, 0.21), 20, 20)
  Sigma <- S %*% t(S) # to ensure positive definite covariance matrix
  X <- rmvnorm(sum(Ns), mean=rep(0, 20), sigma=Sigma)
  X <- scale(X) # scale to get standardised regression coefficients as effect size
  slopes <- runif(20, -2, 2)
  Y <- X %*% slopes + rnorm(nrow(X), 0, 5)
  
  fullData <- cbind.data.frame(Y, X) 
  colnames(fullData) <- c("Y", paste0("X", 1:20))
  fullData$ID <- rep(paste0("study", 1:40), times=Ns)
  
  index <- sample(sum(Ns), sum(Ns)) # pick data from full simulation
  
  studyList <- list()
  for (i in 1:Nstudies){
    studyList[[i]] <- fullData[index[1:Ns[i]], sort(c(1, 2, 3, sample(4:21, 8)))] # pick 10 predictors, 8 randomly
    index <- index[-c(1:Ns[i])] # remove data already assigned
  }

  # now find out which substudy would NOT have been published because of insignificant X1 or X2, hence representing our publication bias:
  fmList <- lapply(studyList, \(x) step(lm(Y ~ ., data=x), 
                                        scope=list("lower"=Y~X1+X2, "upper"=Y ~ .), 
                                        k=log(nrow(x)), trace=F))
  notSignif <- which(sapply(fmList, \(x) (any(summary(x)$coefficients[2:3, 4] > 0.05))))
  
  # produce a second data set with only published studies:
  publishedData <- fullData[!(fullData$ID %in% paste0("study", notSignif)), ]
  
  # now two simple LMs, one for full, one for published:
  flmfull <- step(lm(Y ~ ., data=fullData), 
                  scope=list("lower"=Y~X1+X2, "upper"=Y ~ .), 
                  k=log(nrow(fullData)), trace=F)                 
  
  flmpublished <- step(lm(Y ~ ., data=publishedData), 
                       scope=list("lower"=Y~X1+X2, "upper"=Y ~ .), 
                       k=log(nrow(publishedData)), trace=F)                 
  
  list("full"=(coef(flmfull)[-1] - slopes[as.numeric(gsub("X", "", names(coef(flmfull)[-1])))])/slopes[as.numeric(gsub("X", "", names(coef(flmfull)[-1])))],
  "published"=(coef(flmpublished)[-1] - slopes[as.numeric(gsub("X", "", names(coef(flmpublished)[-1])))])/slopes[as.numeric(gsub("X", "", names(coef(flmpublished)[-1])))] )
})  
# save(full.simu, file="full.simu.RData")

dim(full.simu)
full.simu[2,2]

dats <- data.frame()
for (i in 1:reps){
  onedat <- data.frame( # all
    "coef"=names(full.simu[1, i][[1]]),
    "estimate"=full.simu[1, i][[1]],
    "type"=rep(names(full.simu[1, i]), length(full.simu[1, i][[1]]))
  )
  twodat <- data.frame( # published
    "coef"=names(full.simu[2, i][[1]]),
    "estimate"=full.simu[2, i][[1]],
    "type"=rep(names(full.simu[2, i]), length(full.simu[2, i][[1]]))
  )
  dats <- rbind(dats, onedat, twodat)
}
summary(dats)
dats$coef <- factor(dats$coef, levels=paste0("X", 1:20))

library(ggplot2)
pdf(file="../figures/fullAnalysis.pdf", width=7, height=5)
ggplot(dats, aes(coef, estimate, fill=type)) +
  geom_hline(yintercept=0, col="grey") +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values=c("grey40", "grey70")) + 
  coord_cartesian(ylim = c(-1.5, 1.5)  ) + xlab(label="predictor") + ylab(label="relative bias") +
  theme_minimal() + theme(legend.title=element_blank())
dev.off()



#### final comparison of bias across the 20 main effects ####
# MA
load(file="effectSizeslist.Rdata") # (effectSizes, relbiasNAs, relbiasAll)
relbiasNAs100 <- sapply(effectSizes.list, \(x) x[[2]])
MAdats <- stack(as.data.frame(relbiasNAs100))
MAdats$type <- rep("meta-analysis", nrow(MAdats))
head(MAdats)

# IPL MA
load(file="IPLMAsimu.Rdata")
head(IPLMAsimu)
IPLMAdats <- stack(as.data.frame(IPLMAsimu))
IPLMAdats$type <- rep("IPL-MA", nrow(MAdats))
head(IPLMAdats)

# onestep
load(file="onestep.simu.RData")
onestep.simu[1,,]
onestepdats <- stack(as.data.frame(onestep.simu[1,,]))
onestepdats$type <- rep("one-step", nrow(onestepdats))
onestepdats <- onestepdats[!is.na(onestepdats$values), ]
head(onestepdats)

# twostep
load(file="twostep.simu.Rdata")
head(twostep.simu)
twostepdats <- stack(as.data.frame(twostep.simu))
twostepdats$type <- rep("two-step", nrow(twostepdats))
twostepdats <- twostepdats[!is.na(twostepdats$values), ]
head(twostepdats)

# full
load(file="full.simu.RData")
fullsimudats <- data.frame()
for (i in 1:100){
  onedat <- data.frame( # all
    "values"=full.simu[1, i][[1]],
    "ind"=names(full.simu[1, i][[1]]),
    "type"=rep(names(full.simu[1, i]), length(full.simu[1, i][[1]]))
  )
  fullsimudats <- rbind(fullsimudats, onedat)
}
head(fullsimudats)

allMA <- rbind.data.frame(MAdats, IPLMAdats, onestepdats, twostepdats, fullsimudats)
allMA$ind <- gsub("V", "X", allMA$ind)
allMA$ind <- factor(allMA$ind, levels=paste0("X", 1:10))
allMA$type <- factor(allMA$type, levels=c("meta-analysis", "IPL-MA", "one-step", "two-step", "full"))
summary(allMA)

allMA <- allMA[allMA$ind %in% paste0("X", 1:20), ]
head(allMA)
library(ggplot2)
pdf(file="../figures/MetaAnalysis_final.pdf", width=7, height=5)
ggplot(allMA, aes(ind, values, fill=type)) +
  geom_hline(yintercept=0, col="grey") +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_grey(start=0.8, end=0.2) +
  coord_cartesian(ylim = c(-3.5, 3.5)  ) + xlab(label="predictor") + ylab(label="relative bias") +
  theme_minimal() + theme(legend.title=element_blank())
dev.off()
