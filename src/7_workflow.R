
1-(1-0.02)^20

#### imputation with pima data: demonstrate that type of imputation has an effect ####
library(faraway)
summary(pima)
# insulin is the level of insulin 2 hours after a glucose administration; in non-diabetic people this value is below 50 or so (https://www.quora.com/Is-2-hour-after-glucose-administration-test-the-same-as-2-hour-serum-insulin-mu-U-ml-test-What-is-the-normal-range-for-the-2-hour-serum-insulin-mu-U-ml-test); this is a clinical test

for (i in 3:7) pima[which(pima[, i] == 0), i] <- NA
summary(pima)

# library(mi)
# pima.mdf <- missing_data.frame(pima)
# show(pima.mdf)
# image(pima.mdf)
# # imputation uses bayesglm (see ?01missing_variable) of all other predictors and its posterior predictive density (ppd)
# set.seed(123)
# options(mc.cores=4)
# imputations <- mi(pima.mdf, n.iter=30, n.chains=4, verbose=F)
# round(mipply(imputations, mean, to.matrix = TRUE), 3) # all chains should be similar
# Rhats(imputations) 
# imputations <- mi(imputations, n.iter=100, verbose=F) # run more (for insulin)
# round(mipply(imputations, mean, to.matrix = TRUE), 3)
# Rhats(imputations) # better!
# 
# hist(imputations)
# #for GLM analysis:
# analysis <- pool(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=imputations, family=binomial, m=10) 
# summary(analysis)
# str(complete(imputations, m=2)) # export a list of two imputed data sets

# so clearly there is potential!


# for imputation approaches use MICE
#see help file for mice::pool for nice example, which allows any kind of model, it seems.

library(mice)
# https://data.library.virginia.edu/getting-started-with-multiple-imputation-in-r/
md.pattern(pima)
# 393 complete, 144 with 1 missing, 192 with 2 missing asf.
# !! check in Ecol Stats book for other plot of missingness pattern, e.g. values of bmi for which insulin is missing
with(pima, plot(as.factor(is.na(insulin)) ~ bmi))
# missingness decreases with BMI (except at extremely high BMI)

# analysis with complete cases only ("listwise deletion"):
summary(fm.raw <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=pima, family=binomial))

# impute using different approaches:
# mean imputation, with cirucularity!! (using y to infer x to infer y):
#pima.impMean <- complete(mice(pima, method="mean", m=5)) 
#pima.impPMM <- complete(mice(pima, method="pmm", m=5)) 
#pima.impRF <- complete(mice(pima, method="rf", m=5)) 

# show what imputation does
pima.impMice.pmm <- mice(pima, method="pmm", print=F)
densityplot(pima.impMice.pmm) # not shown

# imputation using only predictors:
pima.impMean <- cbind.data.frame(complete(mice(pima[,-9], method="mean", m=5, print=F)), "test"=pima$test ) # exclude "y"
pima.impPMM <- cbind.data.frame(complete(mice(pima[,-9], method="pmm", m=5, print=F)), "test"=pima$test ) # exclude "y"
pima.impRF <- cbind.data.frame(complete(mice(pima[,-9], method="rf", m=5, print=F)), "test"=pima$test ) # exclude "y"

summary(fm.mean <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=pima.impMean, family=binomial))
summary(fm.pmm <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=pima.impPMM, family=binomial))
summary(fm.rf <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=pima.impRF, family=binomial))

barplot(t(cbind(coef(fm.raw), coef(fm.mean), coef(fm.pmm), coef(fm.rf))[-1,]), beside=T)
# having some value is better than none, it seems; insulin not affected, despite having most holes: --> mcar!!

# plot coefs with error bars:
library(ggplot2)
# estimates without intercept:
dats <- data.frame("estimate"=c(coef(fm.raw)[-1], coef(fm.mean)[-1], coef(fm.pmm)[-1], coef(fm.rf)[-1]), "type"=factor(rep(c("raw", "mean", "PMM", "RF"), each=6), levels=c("raw", "mean", "PMM", "RF")), "effect"=rep(c("pregnant", "diastolic", "triceps", "bmi", "age", "insulin"), times=4), "SE"=c(summary(fm.raw)$coefficients[-1,2], summary(fm.mean)$coefficients[-1,2], summary(fm.pmm)$coefficients[-1,2], summary(fm.rf)$coefficients[-1,2]))
#pdf(file="../figures/imputation1.pdf", width=8, height=5)
ggplot(dats, aes(x=effect, y=estimate, fill=type)) + geom_bar(stat="identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=estimate-2*SE, ymax=estimate+2*SE), width=.2, position=position_dodge(.9)) +
  scale_fill_grey() + theme_minimal()
#dev.off()
# Caption: Imputation can affect model estimates and indeed significance. Here all three imputation methods render the predictor ``pregnang'' significant, which it is not in the raw data. Estimates are different for imputation approaches, but significances are not.


##### But: is this correct? SIMULATION! MCAR #####
# idea: use imputed pima data set, with set test-relationship, and thin in simulations!
coef(fm.raw)
sim <- pima.impPMM
dim(sim)

##### once #####
sim <- pima.impPMM
set.seed(123)
sim$test <- with(pima.impPMM,  rbinom(nrow(pima), size=1, prob=plogis(-6 + 0.06*pregnant + 0.004*diastolic + 0.015*triceps + 0.06*bmi + 0.05*age + 0.004*insulin)))
# remove randomly 30% of predictor values:
# it is simpler to shoot holes into a matrix than into a data.frame:
holes <- matrix(1, ncol=8, nrow=738)
holes[sample(round(738*8), size=0.3*8*738)] <- NA
# now we identify position of holes, using row/column indexing:
NAindex <- which(is.na(holes), arr.ind=T)
# transfer of holes:
sim[NAindex] <- NA
summary(sim)
md.pattern(sim) # very different to original data!

# now impute:
sim.impMean <- cbind.data.frame(complete(mice(sim[,-9], method="mean", m=5, print=F)), "test"=sim$test )
sim.impPMM <- cbind.data.frame(complete(mice(sim[,-9], method="pmm", m=5, print=F)), "test"=sim$test )

# now run models and extract coefficients into results matrix
summary(fm.raw <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=sim, family=binomial))
summary(fm.mean <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=sim.impMean, family=binomial))
summary(fm.pmm <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=sim.impPMM, family=binomial))

cbind(coef(fm.raw), coef(fm.mean), coef(fm.pmm))

##### repeat many times #####
set.seed(123)
res <- replicate(5, { # 500
  sim <- pima.impPMM
  sim$test <- with(pima.impPMM, rbinom(nrow(pima), size=1, prob=plogis(-6 + 0.06*pregnant + 0.004*diastolic + 0.015*triceps + 0.06*bmi + 0.05*age + 0.004*insulin)))
  # remove randomly 30% of predictor values:
  # it is simpler to shoot holes into a matrix than into a data.frame:
  holes <- matrix(1, ncol=8, nrow=738)
  holes[sample(round(738*8), size=0.3*8*738)] <- NA
  # now we identify position of holes, using row/column indexing:
  NAindex <- which(is.na(holes), arr.ind=T)
  # transfer of holes:
  sim[NAindex] <- NA

  # now impute:
  sim.impMean <- cbind.data.frame(complete(mice(sim[,-9], method="mean", m=5, print=F)), "test"=sim$test )
  sim.impPMM <- cbind.data.frame(complete(mice(sim[,-9], method="pmm", m=5, print=F)), "test"=sim$test )
  
  # now run models and extract coefficients into results matrix
  fm.raw <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=sim, family=binomial)
  fm.mean <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=sim.impMean, family=binomial)
  fm.pmm <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=sim.impPMM, family=binomial)
  
  cbind("none"=coef(fm.raw), "mean"=coef(fm.mean), "PMM"=coef(fm.pmm))
})
#save(res, file="pimaImputationMCAR.Rdata")
load(file="pimaImputationMCAR.Rdata")
#dimnames(res)[[2]] <- c("none", "mean", "PMM")
signif(apply(res, c(1,2), mean), 2)
signif(apply(res, c(1,2), sd), 3)

dats <- as.data.frame.table(res[-1,,], responseName="estimate") # omit intercept
colnames(dats) <- c("effect", "imputation", "replicate", "estimate")
summary(dats)
#pdf(file="../figures/imputation2MCAR.pdf", width=8, height=5)
ggplot(dats, aes(x=effect, y=estimate, fill=imputation)) + geom_boxplot() + scale_fill_grey() + 
  geom_segment(aes(x = .5, y = 0.06, xend = 1.5, yend = 0.06)) + 
  geom_segment(aes(x = 1.5, y = 0.004, xend = 2.5, yend = 0.004)) + 
  geom_segment(aes(x = 2.5, y = 0.015, xend = 3.5, yend = 0.015)) + 
  geom_segment(aes(x = 3.5, y = 0.06, xend = 4.5, yend = 0.06)) + 
  geom_segment(aes(x = 4.5, y = 0.05, xend = 5.5, yend = 0.05)) + 
  geom_segment(aes(x = 5.5, y = 0.004, xend = 6.5, yend = 0.004)) + 
  guides(fill=guide_legend(title="")) +
  theme_minimal() +
  theme(legend.position = c(.8, .25), axis.title.x=element_blank())
#dev.off()


##### But: is this correct? SIMULATION! MAR #####
# make probability of missing a function of response
##### once #####
set.seed(123)
simMAR <- pima.impPMM
simMAR$prob <- with(pima.impPMM,  plogis(-6 + 0.06*pregnant + 0.004*diastolic + 0.015*triceps + 0.06*bmi + 0.05*age + 0.004*insulin))
simMAR$test <- rbinom(nrow(pima), size=1, prob=simMAR$prob)

# remove 30% of two variables with probability proportional to simMAR$prob:
holes1 <- rbinom(nrow(simMAR), prob=simMAR$prob, size=1)
holes2 <- rbinom(nrow(simMAR), prob=simMAR$prob, size=1)
# transfer of holes:
simMAR$triceps[which(holes1==1)] <- NA
simMAR$insulin[which(holes2==1)] <- NA
summary(simMAR)
md.pattern(simMAR)

# now rerun imputation as for MCAR
##### repeat many times #####
set.seed(123)
resMAR <- replicate(5, { # 500
  simMAR <- pima.impPMM
  simMAR$prob <- with(pima.impPMM,  plogis(-6 + 0.06*pregnant + 0.004*diastolic + 0.015*triceps + 0.06*bmi + 0.05*age + 0.004*insulin))
  simMAR$test <- rbinom(nrow(pima), size=1, prob=simMAR$prob)
  # remove 30% of two variables with probability proportional to simMAR$prob:
  holes1 <- rbinom(nrow(simMAR), prob=simMAR$prob, size=1)
  holes2 <- rbinom(nrow(simMAR), prob=simMAR$prob, size=1)
  # transfer of holes:
  simMAR$triceps[which(holes1==1)] <- NA
  simMAR$insulin[which(holes2==1)] <- NA

  # now impute:
  simMAR.impMean <- cbind.data.frame(complete(mice(simMAR[,-9], method="mean", m=5, print=F)), "test"=simMAR$test )
  simMAR.impPMM <- cbind.data.frame(complete(mice(simMAR[,-9], method="pmm", m=5, print=F)), "test"=simMAR$test )
  
  # now run models and extract coefficients into results matrix
  fm.raw <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=simMAR, family=binomial)
  fm.mean <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=simMAR.impMean, family=binomial)
  fm.pmm <- glm(test ~ pregnant + diastolic + triceps + bmi + age + insulin, data=simMAR.impPMM, family=binomial)
  
  cbind("none"=coef(fm.raw), "mean"=coef(fm.mean), "PMM"=coef(fm.pmm))
})
#save(resMAR, file="pimaImputationMAR.Rdata")
load(file="pimaImputationMAR.Rdata")
signif(apply(resMAR, c(1,2), mean), 2)
signif(apply(resMAR, c(1,2), sd), 3)


dats <- as.data.frame.table(resMAR[-1,,], responseName="estimate") # omit intercept
colnames(dats) <- c("effect", "imputation", "replicate", "estimate")
summary(dats)
#pdf(file="../figures/imputation3MAR.pdf", width=8, height=5)
ggplot(dats, aes(x=effect, y=estimate, fill=imputation)) + geom_boxplot() + scale_fill_grey() + 
  geom_segment(aes(x = .5, y = 0.06, xend = 1.5, yend = 0.06)) + 
  geom_segment(aes(x = 1.5, y = 0.004, xend = 2.5, yend = 0.004)) + 
  geom_segment(aes(x = 2.5, y = 0.015, xend = 3.5, yend = 0.015)) + 
  geom_segment(aes(x = 3.5, y = 0.06, xend = 4.5, yend = 0.06)) + 
  geom_segment(aes(x = 4.5, y = 0.05, xend = 5.5, yend = 0.05)) + 
  geom_segment(aes(x = 5.5, y = 0.004, xend = 6.5, yend = 0.004)) + 
  guides(fill=guide_legend(title="")) +
  theme_minimal() +
  theme(legend.position = c(.8, .25), axis.title.x=element_blank())
#dev.off()









#### Missing Data  #### 
# ill advice on the internet, particular from computer scientists: https://www.askpython.com/python/examples/impute-missing-data-values
# goal of simulation: demonstrate problem with mean/median vs kNN/randomForest-based replacement for MCAR and MAR (vs MNAR)
# simu with mice: https://statisticsglobe.com/regression-imputation-stochastic-vs-deterministic/
# titanic imputation: https://jyu-theartofml.github.io/posts/imputation

library(effects)
data("TitanicSurvival")
attach(TitanicSurvival)
plot(survived ~ is.na(age)) # more dead without age
plot(sex ~ is.na(age)) # more male without age
plot(passengerClass ~ is.na(age)) # more 3rd without age; largest effect on missingness!
detach(TitanicSurvival)

# ? Does it matter whether I impute?
# ? Does it matter how I impute?
# ?? Why use simulation, not directly compare effect on data?? --> differentiate between MCAR, MAR and MNAR effects
# --> ? Why does it matter how NAs emerged?
#   
# 1. Simulate complete data set. 
# 2. Then introduce NAs as MCAR, MAR and MNAR. 
# 3. Then impute by kNN and estimate effects: do they differ depending on how NAs were produced?
# 4. Repeat 100 times.

# Ad 1.:
Tita <- na.omit(TitanicSurvival)
fm <- glm(survived ~ (age + sex + passengerClass)^2, data=Tita, family=binomial)
summary(fm)
# take similar values for simulating data:
beta <- c(3, 0, -2.5, 0, -3, -0.03, -0.05, 0, -1, 2)
X <- model.matrix(fm)
Yhat <- plogis(as.matrix(X) %*% beta)
set.seed(1234)
Tita$survived <- rbinom(length(Yhat), prob=Yhat, size=1)
fmS <- glm(survived ~ (age + sex + passengerClass)^2, data=Tita, family=binomial)
#summary(fmS)

# Ad 2:
Tita.mnar <- Tita.mar <- Tita.mcar <- Tita
NAprop <- 0.2 # proportion of NAs
N <- length(Yhat)
set.seed(123)
# MCAR:
Tita.mcar$age[sample(N, size=floor(NAprop*N), replace=F)] <- NA # MCAR
# MAR: 80% in males, 20% in females
table(Tita$sex)
set.seed(123)
NAmales <- which(Tita$sex=="male")[sample(658, size=floor(NAprop*658*0.8))]
NAfemales <- which(Tita$sex=="female")[sample(388, size=floor(NAprop*388*0.2))]
Tita.mar$age[c(NAmales, NAfemales)] <- NA
summary(Tita.mar)
# MNAR: 80% of NAs onto dead, 20% onto survived
table(Tita$survived)
set.seed(123)
NAdead <- which(Tita$survived==0)[sample(619, size=floor(NAprop*619*0.8))]
NAalive <- which(Tita$survived==1)[sample(427, size=floor(NAprop*427*0.2))]
Tita.mnar$age[c(NAdead, NAalive)] <- NA
summary(Tita.mnar)

coef(fmS)
coef(update(fmS, .~., data=Tita.mcar))
coef(update(fmS, .~., data=Tita.mar))
coef(update(fmS, .~., data=Tita.mnar))
# so yes, coef #6 is rather different (age:sexmale), as is age:passengerClass3rd

# Ad 3:
library(mice)
?mice
# pmm similar to kNN
# mean imputation:
Tita.mcar.impMean <- complete(mice(Tita.mcar, method="mean", m=1)) 
#summary(Tita.mcar.impMean)
Tita.mar.impMean <- complete(mice(Tita.mar, method="mean", m=1))
Tita.mnar.impMean <- complete(mice(Tita.mnar, method="mean", m=1))
# randomForest imputation, using only predictors:
Tita.mcar.impRF <- complete(mice(Tita.mcar[,-1], method="rf", m=1))   
Tita.mar.impRF <- complete(mice(Tita.mar[,-1], method="rf", m=1))   
Tita.mnar.impRF <- complete(mice(Tita.mnar[,-1], method="rf", m=1))   
Tita.mcar.impRF$survived <- Tita.mar.impRF$survived <- Tita.mnar.impRF$survived <- Tita$survived


coef.mat <- rbind(
  # true
  c(3, 0, -2.5, 0, -3, -0.03, -0.05, 0, -1, 2),
  # with NAs
  coef(update(fmS, .~., data=Tita.mcar)),
  coef(update(fmS, .~., data=Tita.mar)),
  coef(update(fmS, .~., data=Tita.mnar)),
  #mean  
  coef(update(fmS, .~., data=Tita.mcar.impMean)),
  coef(update(fmS, .~., data=Tita.mar.impMean)),
  coef(update(fmS, .~., data=Tita.mnar.impMean)),
  #rf
  coef(update(fmS, .~., data=Tita.mcar.impRF)),
  coef(update(fmS, .~., data=Tita.mar.impRF)),
  coef(update(fmS, .~., data=Tita.mnar.impRF))
)
# subtract truth from all estimates:
diff.mat <- t(t(coef.mat) - coef.mat[1,])
barplot(diff.mat[,c(2, 6, 7, 8)], beside=T)
# reorder, with missing data types of a kind next to each other:
barplot(diff.mat[c(2,5,8, 3,6,9, 4,7,10),c(2, 6, 7, 8)], beside=T, col=c("grey30", "grey50", "grey70"), las=0, ylab="bias of estimate", space=c(0, 1), border=NA)
legend("bottomright", legend=c("with NAs", "mean imputed", "rf imputed"), bty="n", col=c("grey30", "grey50", "grey70"), pch=15, cex=1.5)
text(c(2.5, 5.5, 8.5), c(0.001, 0.001, 0.001), c("MCAR", "MAR", "MNAR"), srt=0)
text(c(12.5, 15.5, 18.5), -c(0.001, 0.001, 0.001), c("MCAR", "MAR", "MNAR"), srt=0)
text(c(22.5, 25.5, 28.5), -c(0.001, 0.001, 0.001), c("MCAR", "MAR", "MNAR"), srt=0)
text(c(32.5, 35.5, 38.5), -c(0.001, 0.001, 0.001), c("MCAR", "MAR", "MNAR"), srt=0)


# Ad 4: Repeat 100 times ...
beta <- c(3, 0, -2.5, 0, -3, -0.03, -0.05, 0, -1, 2)
X <- model.matrix(fm)
Yhat <- plogis(as.matrix(X) %*% beta)
Tita$survived <- rbinom(length(Yhat), prob=Yhat, size=1)
fmS <- glm(survived ~ (age + sex + passengerClass)^2, data=Tita, family=binomial)
NAprop <- 0.4 # proportion of NAs
N <- length(Yhat)

reps <- 50
diff.arr <- array(NA, dim=c(10, 10, reps))
set.seed(123)
for (i in 1:reps){
  # Ad 2:
  Tita.mnar <- Tita.mar <- Tita.mcar <- Tita
  # MCAR:
  Tita.mcar$age[sample(N, size=floor(NAprop*N), replace=F)] <- NA # MCAR
  # MAR: 80% in males, 20% in females
  NAmales <- which(Tita$sex=="male")[sample(658, size=floor(NAprop*658*0.8))]
  NAfemales <- which(Tita$sex=="female")[sample(388, size=floor(NAprop*388*0.2))]
  Tita.mar$age[c(NAmales, NAfemales)] <- NA
  # MNAR: 80% of NAs onto dead, 20% onto survived
  NAdead <- which(Tita$survived==0)[sample(619, size=floor(NAprop*619*0.8))]
  NAalive <- which(Tita$survived==1)[sample(427, size=floor(NAprop*427*0.2))]
  Tita.mnar$age[c(NAdead, NAalive)] <- NA
  # Ad 3:
  # mean imputation:
  Tita.mcar.impMean <- complete(mice(Tita.mcar, method="mean", m=1, print=FALSE )) 
  Tita.mar.impMean  <- complete(mice(Tita.mar,  method="mean", m=1, print=FALSE ))
  Tita.mnar.impMean <- complete(mice(Tita.mnar, method="mean", m=1, print=FALSE ))
  # randomForest imputation, using only predictors:
  Tita.mcar.impRF <- complete(mice(Tita.mcar[,-1], method="rf", m=1, print=FALSE ))   
  Tita.mar.impRF  <- complete(mice(Tita.mar[,-1],  method="rf", m=1, print=FALSE ))   
  Tita.mnar.impRF <- complete(mice(Tita.mnar[,-1], method="rf", m=1, print=FALSE ))   
  Tita.mcar.impRF$survived <- Tita.mar.impRF$survived <- Tita.mnar.impRF$survived <- Tita$survived
  
  coef.mat <- rbind(
    # true
    c(3, 0, -2.5, 0, -3, -0.03, -0.05, 0, -1, 2),
    # with NAs
    coef(update(fmS, .~., data=Tita.mcar)),
    coef(update(fmS, .~., data=Tita.mar)),
    coef(update(fmS, .~., data=Tita.mnar)),
    #mean  
    coef(update(fmS, .~., data=Tita.mcar.impMean)),
    coef(update(fmS, .~., data=Tita.mar.impMean)),
    coef(update(fmS, .~., data=Tita.mnar.impMean)),
    #rf
    coef(update(fmS, .~., data=Tita.mcar.impRF)),
    coef(update(fmS, .~., data=Tita.mar.impRF)),
    coef(update(fmS, .~., data=Tita.mnar.impRF))
  )
  # subtract truth from all estimates:
  diff.mat <- t(t(coef.mat) - coef.mat[1,])
  
  diff.arr[,,i] <- diff.mat
  cat(i, " ")
}  
diff.arr

diff.means <- apply(diff.arr, c(1,2), mean)
diff.sds <- apply(diff.arr, c(1,2), sd)
barplot(diff.means[c(2,5,8, 3,6,9, 4,7,10),c(2, 6, 7, 8)], beside=T, col=c("grey30", "grey50", "grey70"), las=0, ylab="bias of estimate", space=c(0, 1), border=NA)
points(0.5+c(1:9, 11:19, 21:29, 31:39), as.vector(diff.means[c(2,5,8, 3,6,9, 4,7,10),c(2, 6, 7, 8)]) + as.vector(diff.sds[c(2,5,8, 3,6,9, 4,7,10),c(2, 6, 7, 8)]), pch="_")
points(0.5+c(1:9, 11:19, 21:29, 31:39), as.vector(diff.means[c(2,5,8, 3,6,9, 4,7,10),c(2, 6, 7, 8)]) - as.vector(diff.sds[c(2,5,8, 3,6,9, 4,7,10),c(2, 6, 7, 8)]), pch="_")

# plot reduction in bias for rf:
#rf.bias.correction <- coef.mat[2:4,] - coef.mat[8:10,]
#barplot(rf.bias.correction[,c(2, 6, 7, 8)], beside=T)




#### Censoring ####
library(NADA)
data(Recon) # regression and/or categorical!
# see here for a more complete analysis of this data set: 
# https://cran.r-project.org/web/packages/NADA2/vignettes/Censored_Data_Analysis.html
summary(Recon)
Recon$sqrtArea <- sqrt(Recon$Area)
hist(log10(Recon$AtraConc))
#to deal with censoring miscoding
Recon$AtraPHLevel <- Recon$AtraCen #Put existing AtraCen into AtraPHLevel
Recon$AtraCen <- FALSE
Recon$AtraCen[Recon$AtraConc<=.05]<-TRUE #Change AtraCen to TRUE if \leq 0.05

attach(Recon)

plot(AtraConc ~ sqrt(Area), pch=ifelse(Recon$AtraCen==T, 1, 16), log="y", las=1)
plot(AtraConc ~ Temp, pch=ifelse(Recon$AtraCen==T, 1, 16), log="y", las=1)
plot(AtraConc ~ Dyplant, pch=ifelse(Recon$AtraCen==T, 1, 16), log="y", las=1)
plot(AtraConc ~ Pctl, pch=ifelse(Recon$AtraCen==T, 1, 16), log="y", las=1)
plot(AtraConc ~ PctCorn, pch=ifelse(Recon$AtraCen==T, 1, 16), log="y", las=1)
plot(AtraConc ~ Applic, pch=ifelse(Recon$AtraCen==T, 1, 16), log="y", las=1)
# AtracConc: Atrazine concentrations in midwestern US streams
# Temp: mean annual temperature
# PctCorn: percentage corn grown in catchment
# Dyplant: days since planting
# Pctl: flow percentile of river

# 1. analysis without censored data:
fmNoCens <- lm(log10(AtraConc) ~ Temp+PctCorn+Dyplant+Pctl+Precip+sqrt(Area)+Applic, data=Recon[Recon$AtraCen==FALSE,])
summary(fmNoCens)
# only area:
fmNoCens <- lm(log10(AtraConc) ~ sqrt(Area), data=Recon[Recon$AtraCen==FALSE,])
summary(fmNoCens)
#range(sqrt(Recon$Area)) # [1]  13.45362 221.69348
10^(coef(fmNoCens)[1] + range(sqrt(Recon$Area)) * coef(fmNoCens)[2]) #y

#now, correctly, versions 2 and 3 return the same results!
# 2. analysis with censored data set to lowest measured value:
# (we cannot use 0 and then log-transform the data ...)
Recon$lowlim <- ifelse(Recon$AtraCen==FALSE, Recon$AtraConc, min(Recon$AtraConc))
fmLowlim <- lm(log10(lowlim) ~ sqrt(Area), data=Recon)
fmlowlim <- lm(log10(AtraConc)~ sqrt(Area), data=Recon)

summary(fmLowlim)
10^(coef(fmLowlim)[1] + range(sqrt(Recon$Area)) * coef(fmLowlim)[2]) #y


# 3. use all data, replacing censored by detection limit:
Recon$highlim <- ifelse(Recon$AtraCen==FALSE, Recon$AtraConc, .05)
fmHighlim <- lm(log10(highlim) ~ sqrt(Area), data=Recon)
summary(fmHighlim)
10^(coef(fmHighlim)[1] + range(sqrt(Recon$Area)) * coef(fmHighlim)[2]) #y

4. analysis with censor model:
library(NADA2)
fmcen <- cencorreg(Recon$AtraConc, Recon$AtraCen, sqrt(Recon$Area), LOG=T)
summary(fmcen)
10^(coef(fmcen)[1] + range(Recon$sqrtArea) * coef(fmcen)[2]) #y


#Fails to converge; survival library not designed for environmental data: USE NADA2
#library(survival)
#summary(fmsurv <- survreg(Surv(AtraConc, AtraCen, type="left") ~ sqrtArea,
#	dist="loggaussian", data=Recon ))
#exp(coef(fmsurv)[1] + range(Recon$sqrtArea) * coef(fmsurv)[2]) #y


# notice change in p-value for sqrt(Area) from ignore to cen


# make ggplot: data, plus line for only black, all, and handling censoring
library(ggplot2)
Recon$sqrtArea <- sqrt(Recon$Area)

#pdf(file="../figures/AtrazinePlot.pdf", width=5, height=5)
ggplot(Recon, aes(y=AtraConc, x=sqrtArea, color=AtraCen)) + 
  geom_point() + 
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "grey")) +
  scale_y_continuous(trans="log10") +
  labs(y=expression(paste('Atrazine concentration [',mu,'g/L]')),
  	x=expression(paste(sqrt('catchment area'),' [',km^2,']')) ) +
  labs(color="Below detection limit:") +
  geom_hline(yintercept=.051, linetype=2, color="grey") + 
  annotate("text", 225, .06, label=expression(paste('Detection limit (0.05', mu,'g/ml)')), hjust=1, color="grey") +
  geom_segment(aes(x = 13.45, y = 1.17, xend = 221.69, yend = 0.70), colour = "black") + #only values > .05
  annotate("text", 225, 1, label="Using only values > 0.05", color="black", hjust=1) +
  geom_segment(aes(x = 13.45, y = 0.77, xend = 221.69, yend = 0.69), colour = "grey70") + #
  annotate("text", 225, .6, label="Using all values as reported", color="grey70", hjust=1) +
  geom_segment(aes(x = 13.45, y = 0.39, xend = 221.69, yend = 0.40), colour = "grey40") +
  annotate("text", 225, .35, label="Using censored regression", color="grey40", hjust=1) +
  theme_minimal() +
  theme(legend.position = c(.8, .85))
#dev.off()

# But: which is correct? And is there a best way to fill in (if next steps in the pipeline cannot handle censoring regression (e.g. ML, spatial/temporal models) and we don't want to go Bayesian)

##### simu censoring #####
# set one detection limit
# use one of the Recon models as rough truth
# try: replace by lowlim, median of cens class, uplim, omit, ignore
# evaluate distance to true parameter, no model selection!
attach(Recon)
set.seed(123)
#simulate the response
y <- 10^(-2 + .15*Temp + 0.015*PctCorn -0.004*Dyplant + 0.01*Pctl + 0.0003*Precip + 0.0015*sqrt(Area) -0.007*Applic + rnorm(423, 0, 0.5))
#hist(log10(y))
# detection limit: 0.05 (as in Recon) 
detectlim <- 0.05 #so we can play with multiple detection limits
# --> set all values below to 0!
y <- ifelse(y<detectlim, 0, y)
Recon$AtCo_NA <- ifelse(y < detectlim, NA, y)
Recon$AtCo_half <- ifelse(y < detectlim, detectlim*0.5, y) # half of detection limit
Recon$AtCo_uplim <- ifelse(y < detectlim, detectlim, y)
library(mice)
Recon$AtCo_impute <- complete(mice(data=Recon[, c(1:8, 11)], m=5, method="pmm", printFlag=F))$AtCo_NA

# analysis:
coefs <- rbind(
#  "ignore"=coef(lm(log10(AtCo_raw) ~ Temp+PctCorn+Dyplant+Pctl+Precip+sqrt(Area)+Applic, data=Recon)),
  "omit"=coef(lm(log10(AtCo_NA) ~ Temp+PctCorn+Dyplant+Pctl+Precip+sqrt(Area)+Applic, data=Recon)),
  "half"=coef(lm(log10(AtCo_half) ~ Temp+PctCorn+Dyplant+Pctl+Precip+sqrt(Area)+Applic, data=Recon)),
  "upper"=coef(lm(log10(AtCo_uplim) ~ Temp+PctCorn+Dyplant+Pctl+Precip+sqrt(Area)+Applic, data=Recon)),
  "impute"=coef(lm(log10(AtCo_impute) ~ Temp+PctCorn+Dyplant+Pctl+Precip+sqrt(Area)+Applic, data=Recon)),
  "true"=c(-2, .15, .015, -.004, .01, .0003, .0015, -.007)
)
round(coefs, 4)

## many runs ##
set.seed(1234)
detectlim<-.05
coefsMany <- replicate(100, {
  y <- 10^(-2 + .15*Temp + 0.015*PctCorn -0.004*Dyplant + 0.01*Pctl + 0.0003*Precip + 0.0015*sqrt(Area) -0.007*Applic + rnorm(423, 0, 0.5))
  # detection limit: detectlim #as in single runs
  Recon$AtCo_raw <- y # use uncensored data for comparison
  y <- ifelse(y<detectlim, 0, y)  
  Recon$AtCo_NA <- ifelse(y < detectlim, NA, y)
  Recon$AtCo_half <-  ifelse(y < detectlim, detectlim*0.5, y)
  Recon$AtCo_uplim <- ifelse(y < detectlim, detectlim, y)
  Recon$AtCo_impute <- complete(mice(data=Recon[, c(1:8, 11)], m=5, method="pmm", printFlag=F))$AtCo_NA

  coefs <- rbind(
    "omit"=coef(lm(log10(AtCo_NA) ~ Temp+PctCorn+Dyplant+Pctl+Precip+sqrt(Area)+Applic, data=Recon)),
    "half"=coef(lm(log10(AtCo_half) ~ Temp+PctCorn+Dyplant+Pctl+Precip+sqrt(Area)+Applic, data=Recon)),
    "upper"=coef(lm(log10(AtCo_uplim) ~ Temp+PctCorn+Dyplant+Pctl+Precip+sqrt(Area)+Applic, data=Recon)),
    "impute"=coef(lm(log10(AtCo_impute) ~ Temp+PctCorn+Dyplant+Pctl+Precip+sqrt(Area)+Applic, data=Recon)),
    "raw"=coef(lm(log10(AtCo_raw) ~ Temp+PctCorn+Dyplant+Pctl+Precip+sqrt(Area)+Applic, data=Recon))
    )
})
round(apply(coefsMany, c(1,2), mean), 4)
# plot for Area:
boxplot(t(coefsMany[,7,]), las=1)
abline(h=0.0015)

dats <- stack(as.data.frame(t(coefsMany[,7,])))
colnames(dats) <- c("estimate", "method")

library(ggplot2)
#pdf(file="../figures/censoringFillins.pdf", width=5, height=5)
ggplot(dats, aes(y=estimate, x=method))  + 
  ylim(-0.00175,0.004)+
  geom_hline(yintercept=0.0015, linetype=2, color="grey") + 
  geom_boxplot(fill=c(rep("gray", 4), "white"))  + geom_jitter(width=0.2) +
  theme_minimal()
#dev.off()



# do the same without censoring to show that otherwise this would work!
# outcomment: 1. the ifelse-line after creating y in the replicate expression; 2. the line re-defining the NA-column

detach(Recon)







data(Atra)
summary(Atra)
AtraCen <- Atra
AtraCen$June <- ifelse(Atra$JuneCen, 0, Atra$June)
AtraCen$Sept <- ifelse(Atra$SeptCen, 0, Atra$Sept)
AtraCen
boxplot(Atra[,1:2], log="y")
abline(h=0.01)

# raw:
boxplot(AtraCen[,1:2]+.001, log="y")
abline(h=0.01)
t.test(log(AtraCen$June+.001), log(AtraCen$Sept+.001)) # *
# 0 omit:
boxplot(AtraCen[,1:2], log="y", ylim=c(0.005, .1))
abline(h=0.01)
t.test(log(AtraCen$June[!AtraCen$JuneCen]), log(AtraCen$Sept[!AtraCen$SeptCen])) # **
# 0 replace by LDT:
t.test(log(Atra$June), log(Atra$Sept)) # **




#### Regression of residual #### 
#ACSWR::usc: effect of unemployment on crime, after correcting for white/non-white?
#  or: wealth -> police -> reporting

#carData:Ginzberg: depression ~ simplicity and fatalism
#::Highway1: accident rate ~ trks after correcting for adt
#::Salaries: salary ~ sex, after accounting for rank
#::UN: lifeExpF ~ ppgdp, correcting for pctUrban; see also UN98

# Use simulation of Freckleton 2002 as motivation: reproduce them!

library(mvtnorm)
set.seed(123)
##### For Fig. 1a, use N=50 and beta=0.2 #####
N <- 50
X <- rmvnorm(N, mean=c(0,0), sigma=matrix(0.05*c(1, 0.7, .7, 1), 2, 2))
cor(X)
x1 <- X[,1]
x2 <- X[,2]
y = 1 + .2*x1 + .2*x2 + runif(N)
fm1 <- lm(y ~ x1)
fm2 <- lm(residuals(fm1) ~ x2)
summary(fm2)
fm3 <- lm(y ~ x1 + x2)
summary(fm3)

# simulation for 10 values of covariance from 0 to 1 
cova <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9, .99)
coefs <- matrix(NA, ncol=2, nrow=10)
colnames(coefs) <- c("resReg", "multReg")
set.seed(123)
for (i in 1:10){
  X <- rmvnorm(N, mean=c(0,0), sigma=matrix(c(1, cova[i], cova[i], 1), 2, 2))
  x1 <- X[,1]
  x2 <- X[,2]
  y = 1 + .2*x1 + .2*x2 + runif(N)
  fm1 <- lm(y ~ x1)
  coefs[i, 1] <- coef(lm(residuals(fm1) ~ x2))[2]
  coefs[i, 2] <- coef(lm(y ~ x1 + x2))[3]
}
coefs
plot(cova^2, coefs[,1])


set.seed(123)
coef.arr <- replicate(1000, {
  coefs <- matrix(NA, ncol=2, nrow=10)
  colnames(coefs) <- c("resReg", "multReg")
  for (i in 1:10){
    X <- rmvnorm(N, mean=c(0,0), sigma=matrix(0.05*c(1, cova[i], cova[i], 1), 2, 2)) # multiply c(.) by 0.05 or so to get Freckleton's spread
    x1 <- X[,1]
    x2 <- X[,2]
    y = 1 + .2*x1 + .2*x2 + runif(N)
    fm1 <- lm(y ~ x1)
    coefs[i, 1] <- coef(lm(residuals(fm1) ~ x2))[2]
    coefs[i, 2] <- coef(lm(y ~ x1 + x2))[3]
  }
  coefs
})
#str(coef.arr)
apply(coef.arr, c(1,2), mean)
apply(coef.arr, c(1,2), quantile, c(0.025, 0.975))
rowMeans(coef.arr[, 1, ])

#plot(cova^2, rowMeans(coef.arr[, 1, ]), ylim=c(0, 0.5), type="o", pch=16, las=1, xlab="squared relationship between independent variables", ylab="estimate of slope")
#abline(h=0.2, col="grey")
#lines(cova^2, apply(coef.arr[, 1, ], 1, quantile, 0.025), lty=2)
#lines(cova^2, apply(coef.arr[, 1, ], 1, quantile, 0.975), lty=2)
#
#lines(cova^2, apply(coef.arr[, 2, ], 1, mean), lty=1, col="green", type="o", pch=16)
#lines(cova^2, apply(coef.arr[, 2, ], 1, quantile, 0.025), lty=2, col="green")
#lines(cova^2, apply(coef.arr[, 2, ], 1, quantile, 0.975), lty=2, col="green")

library(ggplot2)
dats <- data.frame("cor2"=cova^2, "coefResid"=rowMeans(coef.arr[, 1, ]), "CIlowResid"=apply(coef.arr[, 1, ], 1, quantile, 0.025), "CIhighResid"=apply(coef.arr[, 1, ], 1, quantile, 0.975),
                   "coefMult"=rowMeans(coef.arr[, 2, ]), "CIlowMult"=apply(coef.arr[, 2, ], 1, quantile, 0.025), "CIhighMult"=apply(coef.arr[, 2, ], 1, quantile, 0.975))
#pdf(file="../figures/Freckleton1a.pdf", width=5, height=5)
ggplot(dats, aes(x=cor2, y=coefResid)) + 
  xlab("squared correlation between independent variables") +
  ylab("estimate of slope") + 
  geom_ribbon(aes(ymin=CIlowResid, ymax=CIhighResid), fill="grey50", alpha=0.4) + 
  geom_line() +
  geom_line(aes(y=coefMult), col="grey40") +
  geom_ribbon(aes(x=cor2, ymin=CIlowMult, ymax=CIhighMult), fill="grey50", alpha=0.4) + 
  geom_hline(yintercept=0.2, linetype=2, color="grey") + 
  theme_minimal() +
  annotate("text", .85, .3, label="multiple regression", color="grey40") +
  annotate("text", .6, 0, label="regression on residuals") +
  coord_cartesian(ylim = c(-1, 1), xlim = c(0,1)) 
#dev.off()

##### now the same for bottom-right figure, with beta=0.7 and N=200: #####
N <- 200
set.seed(123)
coef.arr <- replicate(1000, {
  coefs <- matrix(NA, ncol=2, nrow=10)
  colnames(coefs) <- c("resReg", "multReg")
  for (i in 1:10){
    X <- rmvnorm(N, mean=c(0,0), sigma=matrix(0.05*c(1, cova[i], cova[i], 1), 2, 2)) # multiply c(.) by 0.05 or so to get Freckleton's spread
    x1 <- X[,1]
    x2 <- X[,2]
    y = 1 + .7*x1 + .7*x2 + runif(N)
    fm1 <- lm(y ~ x1)
    coefs[i, 1] <- coef(lm(residuals(fm1) ~ x2))[2]
    coefs[i, 2] <- coef(lm(y ~ x1 + x2))[3]
  }
  coefs
})

library(ggplot2)
dats <- data.frame("cor2"=cova^2, "coefResid"=rowMeans(coef.arr[, 1, ]), "CIlowResid"=apply(coef.arr[, 1, ], 1, quantile, 0.025), "CIhighResid"=apply(coef.arr[, 1, ], 1, quantile, 0.975),
                   "coefMult"=rowMeans(coef.arr[, 2, ]), "CIlowMult"=apply(coef.arr[, 2, ], 1, quantile, 0.025), "CIhighMult"=apply(coef.arr[, 2, ], 1, quantile, 0.975))
#pdf(file="../figures/Freckleton1d.pdf", width=5, height=5)
ggplot(dats, aes(x=cor2, y=coefResid)) + 
  xlab("squared correlation between independent variables") +
  ylab("estimate of slope") + 
  geom_ribbon(aes(ymin=CIlowResid, ymax=CIhighResid), fill="grey50", alpha=0.4) + 
  geom_line() +
  geom_line(aes(y=coefMult), col="grey40") +
  geom_ribbon(aes(x=cor2, ymin=CIlowMult, ymax=CIhighMult), fill="grey50", alpha=0.4) + 
  geom_hline(yintercept=0.7, linetype=2, color="grey") + 
  theme_minimal() +
  annotate("text", .85, .75, label="multiple regression", color="grey40") +
  annotate("text", .6, 0.1, label="regression on residuals") +
  coord_cartesian(ylim = c(-.5, 1), xlim = c(0,1)) 
#dev.off()


#### preselection of predictors ####

# https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic: see there for reference and DOI
# colnames copy-pasted from webpage, not part of the data!!!
#Idea: preselect predictors with randomForest, then use GLM for modelling
# (optional) first show that pre-selection makes a difference
# then evaluate the actual difference to truth using simu
# simu from multiple logistic regression (at least third-order effects!)
# evaluate: 1: estimation of coefficients as relative bias; 2: classification accuracy
# because we need the relative bias, we have to use standardised regression coefficients (so that they are comparable)
# thus, we scale the data before analysis and simulations!


#BCWD <- read.csv("../Data/breast+cancer+wisconsin+diagnostic/wdbc.data", stringsAsFactors = T, head=F)
#colnames(BCWD) <- c("ID", "Diagnosis", "radius1", "texture1", "perimeter1", "area1", "smoothness1", "compactness1", "concavity1", "concave_points1", "symmetry1", "fractal_dimension1", "radius2", "texture2", "perimeter2", "area2", "smoothness2", "compactness2", "concavity2", "concave_points2", "symmetry2", "fractal_dimension2", "radius3", "texture3", "perimeter3", "area3", "smoothness3", "compactness3", "concavity3", "concave_points3", "symmetry3", "fractal_dimension3")
#summary(BCWD)
#write.csv(BCWD, file="../Data/breast+cancer+wisconsin+diagnostic/BCWD.csv", row.names = F)


BCWD <- read.csv(file="../Data/breast+cancer+wisconsin+diagnostic/BCWD.csv", stringsAsFactors = T, row.names = 1)
summary(BCWD)
BCWD$Diagnosis <- ifelse(BCWD$Diagnosis == "B", 0, 1) # recode as 0 (benign) and 1 (malign)

# round(cor(BCWD[,-1]), 2) # lots of high correlation, making importance values dubious (good!)
# 
# library(ranger)
# fr <- ranger(Diagnosis ~., data=BCWD, importance="impurity")
# sort(importance(fr), decr=T)
# 
# library(rpart)
# frpart <- rpart(Diagnosis ~., data=BCWD, control=rpart.control(minsplit=5))
# plot(frpart)
# text(frpart)
BCWDs <- as.data.frame(scale(BCWD))
BCWDs[,1] <- BCWD[,1]

fglm <- glm(Diagnosis ~ (radius3+concave_points3+texture3+texture1+concave_points1)^3 + I(area3^2), family=binomial, data=BCWDs)
summary(fglm)
# plot(Diagnosis ~ area3, data=BCWDs)

fstep <- step(fglm, k=2, trace=F) # log(nrow(BCWD))
summary(fstep)
# add another 3-way interaction:
fmmanu1 <- update(fstep, . ~ . + radius3:texture1:concave_points1 )
summary(fmmanu1)

# use this a data-generating model?

# Alternative: just randomly pick some combination of predictors, fit them to get good values, then run with that!
set.seed(123)
picks <- sample(30, size=6)
colnames(BCWD)[picks]

fDGM <- glm(Diagnosis ~ (area2+concave_points2+perimeter2+texture1+symmetry1)^3 + I(concavity2^2), family=binomial, data=BCWDs)
summary(fDGM)

beta <- coef(fDGM)
X <- model.matrix(fDGM)
# hist(plogis(X %*% beta))

set.seed(123)
table(rbinom(569, prob=plogis(X %*% beta), size=1))

sim <- cbind.data.frame(Diagnosis=rbinom(569, prob=plogis(X %*% beta), size=1), BCWDs[,-1])

##### simu analysis #####
# two-step GLM: always pick k best predictors, then step:
k <- 5
## screen for correlation with response (see Murray & Conner 2009 for justification)
sort(abs(cor(sim)[,1]))
varsC <- names(sort(abs(cor(sim)[,1]), decr=T))[2:(k+1)]
RHS1C <- paste0(varsC, collapse=" + ") # first-order interactions
RHS2C <- paste0(varsC, collapse="^2) + I(") # quadratic terms
fC <- formula(paste0("Diagnosis ~ (", RHS1C, ")^2 + I(", RHS2C, "^2)")) # as formula for GLM
# fC <- formula(paste0("Diagnosis ~ ", RHS1C, "+ I(", RHS2C, "^2)")) # as formula for GLM
fglmC <- glm(fC, family=binomial, data=sim)
fstepC <- step(fglmC, k=log(569), trace=F)
summary(fstepC)

# set up a matrix with true and estimated coefficients:
paraEsts <- t(beta)
paraEsts <- rbind(paraEsts, rep(0, ncol(paraEsts)))
for (i in 1:length(coef(fstepC))){
  putthere <- which(colnames(paraEsts) == names(coef(fstepC))[i])
  if (length(putthere) == 0) { # if parameter name does not exist yet
    paraEsts <- cbind(paraEsts, c(0, coef(fstepC)[i]))
    colnames(paraEsts)[ncol(paraEsts)] <- names(coef(fstepC))[i]
  } else { # if parameter name does exist
    paraEsts[2, putthere] <- coef(fstepC)[i]
  }
}
paraEsts

## two-step ML:
library(ranger)
fr <- ranger(Diagnosis ~., data=sim, importance="impurity")
sort(importance(fr), decr=T)[1:k]
varsR <- names(sort(importance(fr), decr=T)[1:k])

RHS1R <- paste0(varsR, collapse=" + ") # first-order interactions
RHS2R <- paste0(varsR, collapse="^2) + I(") # quadratic terms
fR <- formula(paste0("Diagnosis ~ ", RHS1R, "+ I(", RHS2R, "^2)")) # as formula for GLM
fglmR <- glm(fR, family=binomial, data=sim)
system.time(fstepR <- step(fglmR, k=log(569), trace=F) )
summary(fstepR)
# add coefs to matrix:
paraEsts <- rbind(paraEsts, rep(0, ncol(paraEsts)))
for (i in 1:length(coef(fstepR))){
  putthere <- which(colnames(paraEsts) == names(coef(fstepR))[i])
  if (length(putthere) == 0) { # if parameter name does not exist yet
    paraEsts <- cbind(paraEsts, c(0, 0, coef(fstepR)[i]))
    colnames(paraEsts)[ncol(paraEsts)] <- names(coef(fstepR))[i]
  } else { # if parameter name does exist
    paraEsts[3, putthere] <- coef(fstepR)[i]
  }
}
rownames(paraEsts) <- c("true", "cor", "rf")
paraEsts
#relative bias:
(paraEsts[2,] - paraEsts[1,]) / pmax(abs(paraEsts[2,]), abs(paraEsts[1,]))
(paraEsts[3,] - paraEsts[1,]) / pmax(abs(paraEsts[3,]), abs(paraEsts[1,]))

mean(abs((paraEsts[2,] - paraEsts[1,]) / pmax(abs(paraEsts[2,]), abs(paraEsts[1,]))), na.rm=T)
mean(abs((paraEsts[3,] - paraEsts[1,]) / pmax(abs(paraEsts[3,]), abs(paraEsts[1,]))), na.rm=T)

### repeat ###
set.seed(123)
relBiasSim <- replicate(100, {
  sim <- cbind.data.frame(Diagnosis=rbinom(569, prob=plogis(X %*% beta), size=1), BCWDs[,-1])
  # two-step GLM: always pick k best predictors, then step:
  k <- 5
  ## screen for correlation with response (see Murray & Conner 2009 for justification)
  varsC <- names(sort(abs(cor(sim)[,1]), decr=T))[2:(k+1)]
  RHS1C <- paste0(varsC, collapse=" + ") # first-order interactions
  RHS2C <- paste0(varsC, collapse="^2) + I(") # quadratic terms
  fC <- formula(paste0("Diagnosis ~ (", RHS1C, ")^2 + I(", RHS2C, "^2)")) # as formula for GLM
  fglmC <- glm(fC, family=binomial, data=sim)
  fstepC <- step(fglmC, k=log(569), trace=F) 

  # set up a matrix with true and estimated coefficients:
  paraEsts <- t(beta)
  paraEsts <- rbind(paraEsts, rep(0, ncol(paraEsts)))
  for (i in 1:length(coef(fstepC))){
    putthere <- which(colnames(paraEsts) == names(coef(fstepC))[i])
    if (length(putthere) == 0) { # if parameter name does not exist yet
      paraEsts <- cbind(paraEsts, c(0, coef(fstepC)[i]))
      colnames(paraEsts)[ncol(paraEsts)] <- names(coef(fstepC))[i]
    } else { # if parameter name does exist
      paraEsts[2, putthere] <- coef(fstepC)[i]
    }
  }

  ## two-step ML:
  fr <- ranger(Diagnosis ~., data=sim, importance="impurity")
  varsR <- names(sort(importance(fr), decr=T)[1:k])
  RHS1R <- paste0(varsR, collapse=" + ") # first-order interactions
  RHS2R <- paste0(varsR, collapse="^2) + I(") # quadratic terms
  fR <- formula(paste0("Diagnosis ~ ", RHS1R, "+ I(", RHS2R, "^2)")) # as formula for GLM
  fglmR <- glm(fR, family=binomial, data=sim)
  fstepR <- step(fglmR, k=log(569), trace=F)
  # add coefs to matrix:
  paraEsts <- rbind(paraEsts, rep(0, ncol(paraEsts)))
  for (i in 1:length(coef(fstepR))){
    putthere <- which(colnames(paraEsts) == names(coef(fstepR))[i])
    if (length(putthere) == 0) { # if parameter name does not exist yet
      paraEsts <- cbind(paraEsts, c(0, 0, coef(fstepR)[i]))
      colnames(paraEsts)[ncol(paraEsts)] <- names(coef(fstepR))[i]
    } else { # if parameter name does exist
      paraEsts[3, putthere] <- coef(fstepR)[i]
    }
  }
  c(mean(abs((paraEsts[2,] - paraEsts[1,]) / pmax(abs(paraEsts[2,]), abs(paraEsts[1,]))), na.rm=T),
  mean(abs((paraEsts[3,] - paraEsts[1,]) / pmax(abs(paraEsts[3,]), abs(paraEsts[1,]))), na.rm=T))
})
rowMeans(relBiasSim)

dats <- stack(as.data.frame(t(relBiasSim)))
library(ggplot2)
#pdf("../figures/preselectionbias.pdf", width=5, height=5)
ggplot(dats, aes(y=values, x=ind)) + geom_boxplot(fill="gray", size=0.25) + 
  ylim(0, 1) +
  xlab("type of predictor preselection") +
  ylab("relative bias") + 
  scale_x_discrete(labels=c("correlation", "randomForest")) +
  geom_hline(yintercept=0) +
  annotate(geom="curve", x=2, y=0.2, xend=1.5, yend=0.01, 
    curvature=.3, arrow=arrow(length = unit(2, "mm"))  ) +
  annotate(geom = "text", x=2.01, y=.2, label = "unbiased", hjust = "left") +
  theme_minimal()
#dev.off()

##### real data: single analysis #####
# reference: stepwise (BIC) GLM
RHS1 <- paste0(colnames(BCWD)[-c(1,22:31)], collapse=" + ") # first-order interactions
RHS2 <- paste0(colnames(BCWD)[-c(1,22:31)], collapse="^2) + I(") # quadratic terms
# f <- formula(paste0("Diagnosis ~ (", RHS1, ")^2 + I(", RHS2, "^2)")) # interactions, and quadratic: MUCH TOO SLOW in step!
f <- formula(paste0("Diagnosis ~ ", RHS1, " + I(", RHS2, "^2)")) # no interactions, only quadratic
# f <- formula(paste0("Diagnosis ~ (", RHS1, ")^2")) # interactions, no quadratic: MUCH TOO SLOW in step!
fglm <- glm(f, family=binomial, data=BCWD) # ignore warning: this is only for setting up the full model, not for evaluation
system.time(fstep <- step(fglm, k=log(569), trace=F) ) # takes 12 minutes (no interactions); 10s with 3er variables excluded; 
summary(fstep)

# two-step GLM: always pick k best predictors, then step:
k <- 5
## screen for correlation with response (see Murray & Conner 2009 for justification)
sort(abs(cor(BCWD)[,1]))
varsC <- names(sort(abs(cor(BCWD)[,1]), decr=T))[2:(k+1)]
RHS1C <- paste0(varsC, collapse=" + ") # first-order interactions
RHS2C <- paste0(varsC, collapse="^2) + I(") # quadratic terms
fC <- formula(paste0("Diagnosis ~ (", RHS1C, ")^2 + I(", RHS2C, "^2)")) # as formula for GLM
# fC <- formula(paste0("Diagnosis ~ ", RHS1C, "+ I(", RHS2C, "^2)")) # as formula for GLM
fglmC <- glm(fC, family=binomial, data=BCWD)
system.time(fstepC <- step(fglmC, k=log(569), trace=F) ) # 3s for quadratic; 
summary(fstepC)


## two-step ML:
library(ranger)
fr <- ranger(Diagnosis ~., data=BCWD, importance="impurity")
sort(importance(fr), decr=T)[1:k]
varsR <- names(sort(importance(fr), decr=T)[1:k])

RHS1R <- paste0(varsR, collapse=" + ") # first-order interactions
RHS2R <- paste0(varsR, collapse="^2) + I(") # quadratic terms
#fR <- formula(paste0("Diagnosis ~ (", RHS1R, ")^2 + I(", RHS2R, "^2)")) # as formula for GLM
fR <- formula(paste0("Diagnosis ~ ", RHS1R, "+ I(", RHS2R, "^2)")) # as formula for GLM
fglmR <- glm(fR, family=binomial, data=BCWD)
system.time(fstepR <- step(fglmR, k=log(569), trace=F) )
summary(fstepR)

##### multiple runs #####



#### MC error propagation ####
set.seed(123)
x <- runif(100, 0, 1)
y <- 2 + 2*x + rnorm(100, sd=1)
plot(x, y)
fm <- lm(y ~ x)
summary(fm)

# version 1: analytical error propagation:
preds1 <- predict(fm, newdata=list(x=0.9), se.fit=T)
# version 2: MC sampling (without VCOV)
preds2 <- rnorm(1000, mean=coef(fm)[1], sd=summary(fm)$coefficient[1,2]) + 0.9* rnorm(1000, mean=coef(fm)[2], sd=summary(fm)$coefficient[2,2])
# version 3: MC sampling (with VCOV)
Sigma <- vcov(fm)
cov2cor(Sigma)
library(mvtnorm)
coefs <- rmvnorm(1000, mean=coef(fm), sigma=Sigma)
preds3 <- coefs[,1] + 0.9*coefs[,2]

summary(preds2)
summary(preds3)

plot(density(preds2))
lines(density(preds3))
curve(dnorm(x, preds1$fit, preds1$se.fit), add=T)

library(ggplot2)
dats <- data.frame("pred"=c(preds2, preds3), "type"=rep(c("without VCOV", "with VCOV"), each=1000))
#pdf(file="../figures/MCwithOutVCOV.pdf", width=5, height=5)
ggplot(dats, aes(x=pred, fill=type)) + geom_density(color="white") + 
  scale_fill_manual(values=c("gray90", "gray60")) +
  stat_function(fun=function(x) dnorm(x, mean=preds1$fit, sd=preds1$se.fit), linewidth=1, color="grey10", n=501) +
  xlab("estimate") + ylab("density") +
  annotate("text", 4, 2.3, label="analytical", color="black") +
  theme_minimal() +
  theme(legend.position = c(0.2, 0.9))
#dev.off()



#### workflow study ####
# camels data from here: https://gdex.ucar.edu/dataset/camels/file.html

hy <- read.csv("../Data/CAMEL/camels_hydro.txt", sep=";")
veg <- read.csv("../Data/CAMEL/camels_vege.txt", sep=";")
clim <- read.csv("../Data/CAMEL/camels_clim.txt", sep=";")
soil <- read.csv("../Data/CAMEL/camels_soil.txt", sep=";")
topo <- read.csv("../Data/CAMEL/camels_topo.txt", sep=";")
geol <- read.csv("../Data/CAMEL/camels_geol.txt", sep=";", stringsAsFactors = T)

camels <- Reduce(function(x, y) merge(x, y, all=TRUE), list(hy, veg, clim, soil, topo, geol))  
summary(camels)
sort(colnames(camels))

library(mice)
md.pattern(camels)
camelsI <- complete(mice(camels, method="rf", print=F, m=5))
md.pattern(camelsI) # fine now

#write.csv(camelsI, file="../Data/camelsImputed.csv")


##### Start here ! #####

camelsI <- read.csv("../Data/camelsImputed.csv", stringsAsFactors = T, row.names=1)
# remove categorical, not for fundamental, but for practical reasons (correlations among them an continuous require just that bit too much coding):
camelsI <- camelsI[, -c(21, 31, 34, 52, 54)]
# also remove q5 and q95 as predictors for q_mean: that is too tautological!
camelsI <- camelsI[, -c(1, 7, 8)]
# summary(camelsI)
# colnames(camelsI)

#library(ranger)
#fr <- ranger(q_mean ~ ., data=camelsI, importance="impurity")
#fr <- ranger(q_mean ~ elev_mean + slope_mean + area_gages2 +
#              p_mean + pet_mean + aridity + p_seasonality + frac_snow + low_prec_freq +
#              q5 + # obvious correction, possibly use offset
#              frac_forest + lai_diff + dom_land_cover + root_depth_50 +
#              soil_porosity + sand_frac + organic_frac +
#              geol_1st_class + carbonate_rocks_frac + geol_porostiy,
#              data=camelsI, importance="impurity")
#sort(importance(fr))

# so: preselection substantially reduces fit, from .87 to .67 (r^2 OOB)

# So imagine that mean runoff (q_mean) was affected only by area times precipitation (total water input into a catchment), of which some proportion evaporates (aridity and/or fraction forest) or feeds into groundwater without discharging at the gauge. Then we would have Q = fraction_forest * aridity * area_gages2 * p_mean. Yes, this is a gross simplification, but then this is an illustration, not a scientific research paper.
# Then log10(q_mean) = porosity + area_gages2 + p_mean + error

## get estimates for this model:
summary(glm(log10(q_mean) ~ frac_forest  + aridity + area_gages2 + p_mean, data=camelsI))


# targets: estimates; uncertainty of estimates; prediction to some (random) catchment
# demonstrate that uncertainty is MUCH too low without bootstrapping
# steps: 1. imputation; 2. collinearity pre-selection; 3. variable selection using RF (combine 2 and 3 in one step: choose, for two collinear predictors, the one with higher importance)

# shoot 15% holes into the data set, to create a need for imputation:
camelsHoles <- camelsI
set.seed(123)
for (i in c(1:50)){
  shot <- sample(671, floor(0.15*671))
  camelsHoles[shot, i] <- NA
}
summary(camelsHoles)

library(mice)
library(ranger)
# step 1: imputation
camelsPlugged <- complete(mice(camelsHoles, method="rf", print=F, m=1))
#summary(camelsPlugged)
# step 2: collinearity
# this need not make sense, just present a strategy someone may come up with!
# 1 - identify variables that contribute 90% to total variable importance
# 2 - compute Spearman cor among those
# 3 - of a pair with |r|>0.7, pick the one with higher importance
# 4 - use them in GLM
fr <- ranger(log10(q_mean) ~ ., data=camelsPlugged, importance="impurity")
imps <- sort(importance(fr), decreasing=T)
contribTo100 <- cumsum(imps) / sum(imps)
cutoff <- which.min(abs(contribTo100 - 0.9)) # closest to .9
subCamelsPlugged <- camelsPlugged[, names(contribTo100)[1:cutoff] ]

cors <- cor(subCamelsPlugged, method="sp")

highcor <- which(abs(cors) > 0.7 & cors != 1, arr.ind=T)
dealWith <- which(highcor[,1] < highcor[,2])
removeHighCollinear <- if (length(dealWith) > 0) TRUE else FALSE
while (removeHighCollinear){
  pickThese <- highcor[dealWith[1],]
  # cors[pickThese[1], pickThese[2]] # check
  remove <- which.min(imps[colnames(cors)[pickThese]])
  subCamelsPlugged <- subCamelsPlugged[,-which(colnames(subCamelsPlugged) == names(remove))]
  cors <- cor(subCamelsPlugged, method="sp")
  highcor <- which(abs(cors) > 0.7 & cors != 1, arr.ind=T)
  dealWith <- which(highcor[,1] < highcor[,2])
  removeHighCollinear <- if (length(dealWith) > 0)  TRUE else FALSE
}
# summary(subCamelsPlugged)
RHS <- paste(colnames(subCamelsPlugged), collapse=" + ")
RHS2 <- paste0("I(", paste(colnames(subCamelsPlugged), collapse="^2) + I("), "^2)")
f <- as.formula(paste("log10(q_mean) ~", RHS, "+", RHS2))
fm <- glm(f, data=cbind.data.frame("q_mean"=camelsPlugged$q_mean, subCamelsPlugged))
fmstep <- step(fm, trace=F)
# summary(fmstep)

#out: 
resSingleRun <- rbind(summary(fmstep)$coefficients[,1:2], "gauge1"=unlist(predict(fmstep, newdata=camelsI[1,], se.fit=T))[1:2])


##### show that uncertainty of prediction (and estimates) is substantially underestimated #####
# bootstrap workflow:

workflow <- function(data){
  # step 1: imputation
  camelsPlugged <- complete(mice(data, method="rf", print=F, m=1)) # !! data instead of camelsHoles
  # step 2: quantification of importance
  fr <- ranger(log10(q_mean) ~ ., data=camelsPlugged, importance="impurity")
  imps <- sort(importance(fr), decreasing=T)
  # step 3: dimensional reduction:
  contribTo100 <- cumsum(imps) / sum(imps)
  cutoff <- which.min(abs(contribTo100 - 0.9)) # closest to .9
  subCamelsPlugged <- camelsPlugged[, names(contribTo100)[1:cutoff] ]
  # step 4: collinearity:
  cors <- cor(subCamelsPlugged, method="sp")
  highcor <- which(abs(cors) > 0.7 & cors != 1, arr.ind=T)
  dealWith <- which(highcor[,1] < highcor[,2])
  removeHighCollinear <- if (length(dealWith) > 0) TRUE else FALSE
  while (removeHighCollinear){
    pickThese <- highcor[dealWith[1],]
    # cors[pickThese[1], pickThese[2]] # check
    remove <- which.min(imps[colnames(cors)[pickThese]])
    subCamelsPlugged <- subCamelsPlugged[,-which(colnames(subCamelsPlugged) == names(remove))]
    cors <- cor(subCamelsPlugged, method="sp")
    highcor <- which(abs(cors) > 0.7 & cors != 1, arr.ind=T)
    dealWith <- which(highcor[,1] < highcor[,2])
    removeHighCollinear <- if (length(dealWith) > 0)  TRUE else FALSE
  }
  # step 5: fit model
  RHS <- paste(colnames(subCamelsPlugged), collapse=" + ")
  RHS2 <- paste0("I(", paste(colnames(subCamelsPlugged), collapse="^2) + I("), "^2)")
  f <- as.formula(paste("log10(q_mean) ~", RHS, "+", RHS2))
  fm <- lm(f, data=cbind.data.frame("q_mean"=camelsPlugged$q_mean, subCamelsPlugged))
  # step 6: simplify
  fmstep <- step(fm, trace=F)
  # step 7: save estimates and prediction
  as.matrix(rbind(summary(fmstep)$coefficients[,1:2],  
                  "gauge1"=unlist(predict(fmstep, newdata=camelsI[1,], se.fit=T))[1:2]))
}
system.time(print(fullrun <- workflow(data=camelsHoles))) # 6 s

set.seed(123)
ResultsBootstrapWorkflow <- replicate(100, workflow(data=camelsHoles[sample(671, size=671, replace=T),]))

# hm: every now and then the imputation does not yield a full data set; 
# instead of finding and fixing the issue, we outmanouvre it.
# since this is not the point of the analysis, and would be confusing for the reader, we pretend the above code had worked.
set.seed(123)
for (i in 541:545){
  out <- try(workflow(data=camelsHoles[sample(671, size=671, replace=T),]))
  if (inherits(out, "try-error")) next
  save(out, file=paste0("ResultsBoostrapWorkflow_", i, ".Rdata"))
  rm(out)
}
length(grep("Results", list.files())) # check how many successful runs there are

# put into list:
ResultsBootstrapWorkflow <- list()
for (i in 1:300){
  try(load(paste0("ResultsBoostrapWorkflow_", i, ".Rdata"))) # loads "out"
  if (!exists("out")) out <- NULL 
  ResultsBootstrapWorkflow[[i]] <- out
  rm(out)
}
# remove blanks:
ResultsBootstrapWorkflow <- Filter(length, ResultsBootstrapWorkflow) # https://stackoverflow.com/questions/19023446/remove-empty-elements-from-list-with-character0



# collate estimates (only main effects, not intercept):
res.arr <- array(0, dim=c(length(ResultsBootstrapWorkflow), 50, 2), 
                 dimnames=list(1:length(ResultsBootstrapWorkflow), c(colnames(camelsI)[-1], "gauge1"), 
                               c("estimate", "se"))) 
                # replicates, 49 estimates plus prediction, mean & se
for (i in 1:dim(res.arr)[1]){ # loop over replicates
  for (j in 1:50){ # loop through effects
    findthis <- dimnames(res.arr)[[2]][j]
    puthere <- which(rownames(ResultsBootstrapWorkflow[[i]]) %in% findthis)
    if (length(puthere) > 0){
      res.arr[i, j, ] <- ResultsBootstrapWorkflow[[i]][findthis, ]
    }
  }
}


res.arr 
image(t(res.arr[,,1] != 0))

boxplot(res.arr[,-50,1])
bootstrappedSEs <- apply(res.arr[-50,,1], 2, sd) # sd of estimates
singleRunSEs <- fullrun[-c(1, nrow(fullrun)),2] 
bsOrder <- match(names(singleRunSEs), names(bootstrappedSEs))
plot(bootstrappedSEs[bsOrder], singleRunSEs)
abline(0,1)

resSingleRun[18,]
# c(mean(res.arr[, 50, 1]), mean(res.arr[, 50, 2])) # mean SE of any single run TOO COMPLICATED to explain
c(mean(res.arr[, 50, 1]), sd(res.arr[, 50, 1])) # from bootstraps

# So: 
# image plot: the majority of predictors are not in the simplified model (fine)
# single run vs. bootstrapped runs: mean self-reported SE when bootstrapping the data is 1.5 times that of the single run [LEAVE OUT: TOO COMPLICATED]
# single run vs. bootstrapped SE: bootstrapped SE is > 2 times that of self-reported on single run!
# single run vs. bootstrapped estimates: 1:1 plot with at least 100 bootstraps!

# plot barplot of proportion of times the predictor was in the bootstrapped model; use colour to indicate those predictors that were in the single-run model
dats <- data.frame(incidence=colSums(res.arr[,-50,1] != 0)/dim(res.arr)[1]*100)
dats$type <- factor(rownames(dats), levels=rownames(dats)[order(dats$incidence)])
summary(dats)
library(ggplot2)
#pdf(file="../figures/workflowIncidence.pdf", height=7, width=5)
ggplot(dats, aes(x=incidence, type)) + 
  geom_col() + 
  ylab("") + xlab("predictor occurs in percentage of bootstraps") + 
  annotate("point", x=rep(-1, 9), y=49-c(0, 2, 5, 1, 3, 4, 9, 6, 18), colour="grey70") +
  theme_minimal()
#dev.off()

# comparison of SEs
CVone <- fullrun[-1,2] / abs(fullrun[-1,1]) * 100 # exclude intercept
CVbs <- apply(res.arr[-51, , 1], 2, sd) / abs(colMeans(res.arr[-51, , 1])) * 100
bsOrder <- match(names(CVone), names(CVbs))
cbind(CVone[c(1:11, 22)], CVbs[na.omit(bsOrder)])

dats <- data.frame("CV"=c(CVone[c(1:11, 22)], CVbs[na.omit(bsOrder)]), "approach"=rep(c("full", "bootstrap"), each=12), "coef"=rep(names(CVone[c(1:11, 22)]), 2))
dats$approach <- factor(dats$approach, levels=c("full", "bootstrap"))

pdf(file="../figures/workflowCoefCV.pdf", height=5, width=5)
ggplot(dats, aes(approach, CV, group=coef)) + geom_line() + geom_point() +  
  ylab("coefficient of variation [%]") + xlab("") + 
  scale_y_continuous(trans='log10') +
  scale_x_discrete(expand = expansion(add = c(0.2, .5))) +
  annotate("text", x=rep(2.02, 12), y=dats$CV[13:24]+c(0,0,-5,0,5,-5,0,0,0,5,0,0), label=dats$coef[1:12], hjust = 0, size=2.7) + 
  theme_minimal()
dev.off()


 ##### Bias analysis with simulate data #####
# We like to think that runoff as being affected by rainfall, forest cover, aridity and area (since precipitation was measured at the gauge and may hence be unrepresentative for the entire catchment).
# We get estimates for this model to use for the simulation:
summary(lm(log10(q_mean) ~ p_mean + frac_forest  + aridity + area_gages2, data=camelsI))

y <- with(camelsI, 10^(-0.4 + 0.135*p_mean + 0.366*frac_forest - 0.275*aridity -0.00004*area_gages2 + rnorm(671, 0, .2)))

#plot(y, camelsI$q_mean, log="")
#abline(0,1)
#cor(y, camelsI$q_mean)^2

datos <- camelsHoles
datos$q_mean <- y
fullrun.sim <- workflow(data=datos)
fullrun.sim

set.seed(123)
camelsSims <- replicate(1, {
  y <- with(camelsI, 10^(-0.4 + 0.135*p_mean + 0.366*frac_forest - 0.275*aridity -0.00004*area_gages2 + rnorm(671, 0, .2)))
  datos$q_mean <- y
  fullrun.sim <- workflow(data=datos)
  fullrun.sim
})
#save(camelsSims, file="camelsSims.Rdata")
load(file="camelsSims.Rdata")

sim.res.arr <- array(0, dim=c(length(camelsSims), 50, 2), 
                 dimnames=list(1:length(camelsSims), c(colnames(camelsI)[-1], "gauge1"), 
                               c("estimate", "se"))) 
# replicates, 49 estimates plus prediction, mean & se
for (i in 1:dim(sim.res.arr)[1]){ # loop over replicates
  for (j in 1:50){ # loop through effects
    findthis <- dimnames(sim.res.arr)[[2]][j]
    puthere <- which(rownames(camelsSims[[i]]) %in% findthis)
    if (length(puthere) > 0){
      sim.res.arr[i, j, ] <- camelsSims[[i]][findthis, ]
    }
  }
}
sim.res.arr[,,1]
# dimnames(sim.res.arr) #: p_mean: 19, aridity: 23; area_gages2: 43; frac_forest: 11
colMeans(sim.res.arr[, c(19, 11, 23, 43, 50),1])
# true prediction value:
with(camelsI, (-0.4 + 0.135*p_mean[1] + 0.366*frac_forest[1] -0.275*aridity[1] -0.00004*area_gages2[1])) # 0.0903

apply(camelsI[, c(20, 12, 24, 44)], 2, sd) / sd(datos$q_mean) * colMeans(sim.res.arr[, c(19, 11, 23, 43),1])

# hm: (estimate - truth)/truth *100% = relative bias
boxplot(sim.res.arr[, c(19, 11, 23, 43), 1])

dats <- stack(as.data.frame(sim.res.arr[, c(19, 11, 23, 43, 50), 1]))
dats$truth <- rep(c(0.135, 0.366, -0.275, -0.00004, 0.0903), each=100)
dats$relBias <- (dats$values - dats$truth) / dats$truth
dats$relBias2 <- (dats$values - dats$truth) / dats$values
summary(dats)
#pdf(file="../figures/workflowCoefRelBias.pdf", height=5, width=5)
ggplot(dats, aes(ind, relBias)) + 
  geom_hline(yintercept=0, color="gray") + 
  geom_boxplot(fill="gray") + 
  xlab("coefficient")  + ylab("relative bias of estimate") + theme_minimal()
#dev.off()
