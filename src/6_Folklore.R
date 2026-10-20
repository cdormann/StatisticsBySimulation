# 6: Folklore #


#### Model selection helps interpretation ####
library(faraway)
summary(africa)
# impute missing values:
library(mice)
africa.imp <- complete(mice(africa, m=5, method="pmm"))
summary(africa.imp)
  
plot(africa.imp) # no obvious non-linearities

library(ranger)
fr <- ranger(miltcoup ~., data=africa.imp, importance = "impurity")
sort(ranger::importance(fr), decreasing=T)

fm <- glm(miltcoup ~ (oligarchy + pollib + parties + pctvote + popn + size + numelec)^2, data=africa.imp, family=poisson)
summary(fm)
fmstep <- step(fm, trace=F, k=log(47)) # BIC
summary(fmstep)

library(MuMIn)
options(na.action = "na.fail")
fd <- dredge(fm, m.lim=c(0,5)) # takes a few minutes: 3
subset(fd, delta < 4)


##### simulate a data set #####
library(faraway)
summary(africa)
library(mice)
africa.imp <- complete(mice(africa, m=5, method="pmm"))

attach(africa.imp)
set.seed(123)
coups <- rpois(nrow(africa.imp), lambda=exp(-.7 + .15*oligarchy -.0018*size + 0.03*numelec + .0002*numelec*size))
dats <- cbind.data.frame(coups, africa.imp)
fst <- step(ffull <- glm(coups ~ (oligarchy + pollib + parties + pctvote + popn + size + numelec)^2, data=dats, family=poisson, maxit=100), trace=F)
summary(ffull)
summary(fst)

# count spurious effects, correct identifications and wrong signs ("s-error"), ? also order of magnitude error?
# compare to full model!!

coef(ffull)
coef(fst)

# matrix of estimates:
coef.mat <- matrix(NA, ncol=28, nrow=2)
colnames(coef.mat) <- names(coef(ffull)[-1]) # kick out intercept!
coef.mat[1,] <- coef(ffull)[-1]
coef.mat[2, colnames(coef.mat) %in% names(coef(fst)[-1])] <- coef(fst)[-1]

# matrix of p-values:
p.mat <- matrix(NA, ncol=28, nrow=2)
colnames(p.mat) <- names(coef(ffull)[-1]) # kick out intercept!
p.mat[1,] <- summary(ffull)$coefficients[-1, 4]
p.mat[2, colnames(p.mat) %in% names(coef(fst)[-1])]  <- summary(fst)$coefficients[-1, 4]

# coef.mat analysis:
truth <- coef(ffull)[-1] # without intercept
truth[] <- 0
truth[c("oligarchy", "size", "numelec", "size:numelec")] <- c(.15, -.0018, .003, .0002)

# true positives
realIndex <- which(truth != 0)
sum(p.mat[1, realIndex] < 0.05) 
sum(p.mat[2, realIndex] < 0.05)

(sum(p.mat[1, realIndex] < 0.05))/28  # true positive rate
(sum(p.mat[2, realIndex] < 0.05))/28  # true positive rate

# false positives
sum(p.mat[1, (1:28)[-realIndex]] < 0.05) # later need to subtract 1 for the (trivial) intercept!
sum(p.mat[2, (1:28)[-realIndex]] < 0.05, na.rm=T) # later need to subtract 1 for the (trivial) intercept!

sum(p.mat[1, (1:28)[-realIndex]] < 0.05)/28 # error rate
sum(p.mat[2, (1:28)[-realIndex]] < 0.05, na.rm=T)/28 # error rate

# READ UP HERE!
# https://cran.r-project.org/web/packages/retrodesign/vignettes/Intro_To_retrodesign.html

# signifcant only:
signifs <-  (p.mat < 0.05 ) * matrix(truth != 0, nrow=2, ncol=28, byrow=T) 
# sign error: only for significant findings!!
sum(sign(coef.mat[1, signifs[1,]==1]) != sign(truth[signifs[1,]==1])) 
sum(sign(coef.mat[2, signifs[2,]==1]) != sign(truth[signifs[2,]==1]), na.rm=T) 

# order-of-magnitude error: only for significant findings!!
ratio1 <- truth / coef.mat[1, signifs[1,]==1] 
sum(abs(ratio1) > 10 | abs(ratio1) < .1) # wrong oder of magnitude
ratio2 <- truth[signifs[2,]==1] / coef.mat[2, signifs[2,]==1] # later need to subtract 1 for the (trivial) intercept!
sum(abs(ratio2) > 10 | abs(ratio2) < .1, na.rm=T) # wrong oder of magnitude

eval.full.and.step <- function(p.mat, coef.mat, truth){
  # coef.mat analysis function
  # p.mat, coef.mat and truth as prepared above; 
  # number of columns of both matrices and length of truth vector must be the same
  
  # true positives
  realIndex <- which(truth != 0)
  
  truePosRateFull <- (sum(p.mat[1, realIndex] < 0.05))/28  # true positive rate
  truePosRateStep <- sum(p.mat[2, realIndex] < 0.05, na.rm=T)/28  # true positive rate
  
  # false positives
  falsePosRateFull <- sum(p.mat[1, (1:28)[-realIndex]] < 0.05)/28 # error rate
  falsePosRateStep <- sum(p.mat[2, (1:28)[-realIndex]] < 0.05, na.rm=T)/28 # error rate
  
  # signifcant only:
  signifs <-  (p.mat < 0.05 ) * matrix(truth != 0, nrow=2, ncol=28, byrow=T) 
  
  # sign/magn error: only for significant findings!!
  if (sum(signifs[1,], na.rm=T) != 0){
    SignErrorFull <- sum(sign(coef.mat[1, signifs[1,]==1]) != sign(truth[signifs[1,]==1])) / sum(signifs[1,], na.rm=T)
    ratio1 <- truth[signifs[1,]==1] / coef.mat[1, signifs[1,]==1] 
    MagnErrorFull <- sum(abs(ratio1) > 10 | abs(ratio1) < .1)  /  sum(signifs[1,], na.rm=T) # wrong oder of magnitude
  } else {
    SignErrorFull <- 0
    MagnErrorFull <- 0
  }
  if (sum(signifs[2,], na.rm=T) != 0){
    # S-error
    SignErrorStep <- sum(sign(coef.mat[2, signifs[2,]==1]) != sign(truth[signifs[2,]==1]), na.rm=T) / sum(signifs[2,], na.rm=T)
    # M-error
    ratio2 <- truth[signifs[2,]==1] / coef.mat[2, signifs[2,]==1] 
    MagnErrorStep <- sum(abs(ratio2) > 10 | abs(ratio2) < .1, na.rm=T) / sum(signifs[2,], na.rm=T)
    
  } else {
    SignErrorStep <- 0
    MagnErrorStep <- 0
  }
  
  return(c(truePosRateFull, truePosRateStep, falsePosRateFull, falsePosRateStep, SignErrorFull, SignErrorStep, MagnErrorFull, MagnErrorStep))
  
}
eval.full.and.step(p.mat, coef.mat, truth) # test

##### structural error many runs #####
set.seed(1234) 
# occassionally throws an error during stepping! maybe have to resort to for loop ....
# that may actually be a nice example for a miss ... seed 1234 immidiately throws an error

# coef.mat analysis:
truth <- coef(ffull)[-1] # without intercept
truth[] <- 0
truth[c("oligarchy", "size", "numelec", "size:numelec")] <- c(.15, -.0018, .003, .0002)

sims <- replicate(1, {
  coups <- rpois(nrow(africa.imp), lambda=exp(-.7 + .15*oligarchy -.0018*size + 0.03*numelec + .0002*numelec*size))
  dats <- cbind.data.frame(coups, africa.imp)
  fst <- step(ffull <- glm(coups ~ (oligarchy + pollib + parties + pctvote + popn + size + numelec)^2, data=dats, family=poisson, maxit=100), trace=F)
  
  # matrix of estimates:
  coef.mat <- matrix(NA, ncol=28, nrow=2)
  colnames(coef.mat) <- names(coef(ffull)[-1]) # kick out intercept!
  coef.mat[1,] <- coef(ffull)[-1]
  coef.mat[2, colnames(coef.mat) %in% names(coef(fst)[-1])] <- coef(fst)[-1]
  
  # matrix of p-values:
  p.mat <- matrix(NA, ncol=28, nrow=2)
  colnames(p.mat) <- names(coef(ffull)[-1]) # kick out intercept!
  p.mat[1,] <- summary(ffull)$coefficients[-1, 4]
  p.mat[2, colnames(p.mat) %in% names(coef(fst)[-1])]  <- summary(fst)$coefficients[-1, 4]
  
  eval.full.and.step(p.mat, coef.mat, truth)
})
rowMeans(sims)
# save(sims, file="fullStepSimus.Rdata")
#hm keeps throwing errors; use a for-loop:
sims <- matrix(NA, nrow=1000, ncol=8)
set.seed(123)
for (i in 1:nrow(sims)){
  coups <- rpois(nrow(africa.imp), lambda=exp(-.7 + .15*oligarchy -.0018*size + 0.03*numelec + .0002*numelec*size))
  dats <- cbind.data.frame(coups, africa.imp)
  fst <- try(step(ffull <- glm(coups ~ (oligarchy + pollib + parties + pctvote + popn + size + numelec)^2, data=dats, family=poisson, maxit=100), trace=F))
  if (inherits(fst, "try-error")) next
  # matrix of estimates:
  coef.mat <- matrix(NA, ncol=28, nrow=2)
  colnames(coef.mat) <- names(coef(ffull)[-1]) # kick out intercept!
  coef.mat[1,] <- coef(ffull)[-1]
  coef.mat[2, colnames(coef.mat) %in% names(coef(fst)[-1])] <- coef(fst)[-1]
  
  # matrix of p-values:
  p.mat <- matrix(NA, ncol=28, nrow=2)
  colnames(p.mat) <- names(coef(ffull)[-1]) # kick out intercept!
  p.mat[1,] <- summary(ffull)$coefficients[-1, 4]
  p.mat[2, colnames(p.mat) %in% names(coef(fst)[-1])]  <- summary(fst)$coefficients[-1, 4]
  
  sims[i,] <- eval.full.and.step(p.mat, coef.mat, truth)
  cat(i, " ")
}
# only 2 errors in output!
which(is.na(rowSums(sims)))
i <- which(is.na(rowSums(sims)))[1] # define i and run body of for-loops!
#save(sims, file="fullStepSimus.Rdata")
load(file="fullStepSimus.Rdata")
apply(sims, 2, summary)
colnames(sims) <- c("truePosRateFull", "truePosRateStep", "falsePosRateFull", "falsePosRateStep", "SignErrorFull", "SignErrorStep", "MagnErrorFull", "MagnErrorStep")

dats <- stack(as.data.frame(sims))
dats$fullstep <- substr(as.character(dats$ind), nchar(as.character(dats$ind)) -4+1, 16)
dats$index <- substr(as.character(dats$ind), 1, nchar(as.character(dats$ind)) -4)
dats$index <- factor(dats$index, levels=c("truePosRate", "falsePosRate", "SignError", "MagnError") , labels=c("true positive", "false positive", "S-error", "M-error"), ordered=T)
head(dats); tail(dats)
library(ggplot2)
  pdf(file="../figures/modelSelectionStructureRates.pdf", width=5, height=5)
  ggplot(dats, aes(y=values, x=index, fill=fullstep)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=grey.colors(2)) + 
    scale_y_continuous(name="rate", limits=c(0, .6), trans="sqrt") +
    guides(fill=guide_legend(title="")) +
    theme_minimal() +
    theme(legend.position = c(.8, 0.8), axis.title.x=element_blank())
  dev.off()
  

  

#### Select one of two correlated predictors ####
soil <- read.csv("../Data/soil.csv")
soil <- as.data.frame(scale(soil[,-c(1,2,3,16)]))
head(soil)
round(cor(soil), 2) # not correlated with Ca or pH: Mn, HumicMatter, Zn: use as reference
set.seed(123)
# soil$S <- round(with(soil, 15 + .5*Ca - 1*pH + rnorm(nrow(soil), mean=0, sd=.2)))
soil$S <- round(with(soil, 15 + .5*Ca + rnorm(nrow(soil), mean=0, sd=.2)))

library(ranger)

# with both:
summary(lm(S ~ Ca + Mg + K + Na + P + Cu + Zn + Mn + HumicMatter + Density + pH + ExchAc, data=soil))
sort(importance(ranger(S ~ Ca + Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, importance="impurity")), decreasing=T)
importance(ranger(S ~ Ca + Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, importance="impurity"))["Ca"]

# with Ca only: 
summary(lm(S ~ Ca + K + Na + P + Cu + Zn + Mn + HumicMatter + Density + pH + ExchAc, data=soil))
importance(ranger(S ~ Ca + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, importance="impurity"))

# with Mg only: 
summary(lm(S ~ Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil))
sort(importance(ranger(S ~ Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, importance="impurity")), decreasing=T)

importance(ranger(S ~ Ca + Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, importance="impurity"))[c("Ca", "Mg")]
importance(ranger(S ~ Ca +K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, importance="impurity"))["Ca"]
importance(ranger(S ~ Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, importance="impurity"))["Mg"]

c(1.14, 0.99) /sum(c(1.14, 0.99)) # single models
imp <- importance(ranger(S ~ Ca + Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, importance="impurity"))[c("Ca", "Mg")]
imp/sum(imp)
# replicate to get more stable estimate
# point: mean importance in full model == mean importance in combined model

set.seed(123)
importances.simu <- replicate(1000, {
  soil$S <- round(with(soil, 15 + .5*Ca + rnorm(nrow(soil), mean=0, sd=.2)))
  c(coef(lm(S ~ Ca + Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil))[c("Ca", "Mg")],
  coef(lm(S ~ Ca + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil))["Ca"],
  coef(lm(S ~ Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil))["Mg"],
  importance(ranger(S ~ Ca + Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, importance="impurity"))[c("Ca", "Mg")],
  importance(ranger(S ~ Ca +K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, importance="impurity"))["Ca"],
  importance(ranger(S ~ Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, importance="impurity"))["Mg"]
  )
})
#save(importances.simu, file="importances.simu.all8estimates.Rdata")
summary(t(importances.simu))
boxplot(t(importances.simu))
# LM: Ca > Mg?
sum(importances.simu[1,] > importances.simu[2,]) # 639
# rf: Ca > Mg?
sum(importances.simu[5,] > importances.simu[6,]) # 562
# 0.5 + mean( (importances.simu[5,] - importances.simu[6,]) /  (importances.simu[5,] + importances.simu[6,]))

# average importance in full rf model:
mean(importances.simu[5,] /  (importances.simu[5,] +importances.simu[6,])) 
mean(importances.simu[6,] /  (importances.simu[5,] +importances.simu[6,])) 

# average importance, if comparing single models
mean(importances.simu[7,] /  (importances.simu[7,] +importances.simu[8,])) 
mean(importances.simu[8,] /  (importances.simu[7,] +importances.simu[8,]))


# # reduced version of this simulation:
# set.seed(123)
# importances.simu <- replicate(1000, {
#   soil$S <- round(with(soil, 15 + .5*Ca + rnorm(nrow(soil), mean=0, sd=.2)))
#   c(importance(ranger(S ~ Ca + Mg + K+ Na+ P+ Cu+ Zn+ Mn+ HumicMatter+ Density+ pH+ ExchAc, data=soil, importance="impurity"))[c("Ca", "Mg")],
#     coef(lm(S ~ Ca + Mg + K+ Na+ P+ Cu+ Zn+ Mn+ HumicMatter+ Density+ pH+ ExchAc, data=soil))[c("Ca", "Mg")])
# })
# sum(importances.simu[1,] > importances.simu[2,]) # 57.1
# 
# sum(importances.simu[3,] > importances.simu[4,]) # 62.4
#  

# analyse the actual estimates and their distribution, to demonstrate that with model selection estimates are unduly "pseudo-precise"

load(file="importances.simu.all8estimates.Rdata") # those from the previous simu, not the one only counting
summary(importances.simu[1,])
summary(importances.simu[2,])

# show estimates of LM:
plot(density(importances.simu[1,]), ylim=c(0, 1.6))
lines(density(importances.simu[2,]), col="grey")
lines(density(importances.simu[3,]), col="black", lty=2)
lines(density(importances.simu[4,]), col="grey", lty=2)
abline(v=c(0,1), col=c("grey", "black"))

# 
library(ggplot2)
dats <- as.data.frame(t(importances.simu)[1:100, 1:4])
colnames(dats) <- c("Ca.lm", "Mg.lm", "Ca.solo", "Mg.solo")

pdf(file="../figures/collinearLM.pdf", width=5, height=5)

ggplot(dats) + 
  geom_density(aes(x=Ca.lm), size=.8) +
  geom_density(aes(x=Mg.lm), color="grey", size=.8) +
  geom_density(aes(x=Ca.solo)) +
  geom_density(aes(x=Mg.solo), color="grey") +
  geom_vline(xintercept=c(0, .5), linetype=2, colour=c("grey70", "black")) +
  xlab(label="model estimates") + 
  geom_text(data = data.frame(x=1.1, y=1.4, label = "[Ca]-only:\n[Ca] effect"),
            mapping = aes(x = x, y = y, label = label),
            size = 3, alpha = 1, inherit.aes = FALSE) +
  geom_text(data = data.frame(x=1.25, y=0.9, label = "[Mg]-only:\n[Mg] effect"),
            mapping = aes(x = x, y = y, label = label),
            size = 3, alpha = 1, inherit.aes = FALSE) +
  geom_text(data = data.frame(x=2.15, y=0.305, label = "[Ca] & [Mg]:\n[Ca] effect"),
            mapping = aes(x = x, y = y, label = label),
            size = 3, alpha = 1, inherit.aes = FALSE) +
  geom_text(data = data.frame(x=-1.8, y=0.27, label = "[Ca] & [Mg]:\n[Mg] effect"),
            mapping = aes(x = x, y = y, label = label),
            size = 3, alpha = 1, inherit.aes = FALSE) +
  theme_minimal()

dev.off()



## inter- vs extrapolation
# make simulation slightly more complex:
#soil$S <- round(with(soil, 15 + 5*Ca - 3*Ca*Ca+ rnorm(nrow(soil), mean=0, sd=.5)))
soil <- read.csv("../Data/soil.csv")
plot(Diversity ~ Ca, data=soil)

flm <- glm(Diversity ~ Ca + Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, family=gaussian)
summary(flm)
fr <- ranger(Diversity ~ Ca + Mg + K+Na+P+Cu+Zn+Mn+HumicMatter+Density+pH+ExchAc, data=soil, keep.inbag=T)

interextra <- data.frame("Ca"=seq(0, 1, len=101), t(colMeans(soil)))
head(interextra)

preds.glm <- predict(flm, newdata=interextra, se.fit=T)
preds.rf <- predict(fr, data=interextra, type="se")
#lines(interextra$Ca, exp(preds.glm$fit))
#lines(interextra$Ca, preds.rf$predictions)

pred <- data.frame("Ca"= interextra$Ca,"lmpred"=(preds.glm$fit), "lmlow"=(preds.glm$fit-2*preds.glm$se.fit), "lmhigh"=(preds.glm$fit+2*preds.glm$se.fit), "rfpred"=preds.rf$predictions, "rflow"=preds.rf$predictions - 2*preds.rf$se, "rfhigh"=preds.rf$predictions + 2*preds.rf$se)

library(ggplot2)
ggplot(soil, aes(x=Ca, y=Diversity)) + 
  geom_point() +
  geom_line(data=pred, aes(x=Ca, y=lmpred)) +
  geom_line(data=pred, aes(x=Ca, y=lmlow), linetype=2) +
  geom_line(data=pred, aes(x=Ca, y=lmhigh), linetype=2) +
  geom_line(data=pred, aes(x=Ca, y=rfpred)) +
  geom_line(data=pred, aes(x=Ca, y=rflow), linetype=2) +
  geom_line(data=pred, aes(x=Ca, y=rfhigh), linetype=2) +
  ylim(0, .4) +
  theme_minimal()
## hm; not a good example: too few data; 
## either remove several predictors or use different example





#### Default detrending ####
# start taken from 
# https://www.r-bloggers.com/2023/01/making-non-stationary-data-stationary/

library(aTSA)
plot(AirPassengers)
adf.test(AirPassengers)
# p>0.05 indicates problem
plot(diff(AirPassengers))
adf.test(diff(AirPassengers))
# plot cumulative variance, i.e. whether it levels off or not:
library(TidyDensity)
plot(cvar(diff(AirPassengers)), type = "l")
# no it doesn't
plot(log(AirPassengers))
plot(diff(log(AirPassengers)))
adf.test(diff(log(AirPassengers)))
plot(cvar(diff(log(AirPassengers))), type = "l")
# looks fine (no need for double-diff, as in blog post)

plot(acf(diff(log(AirPassengers)), type="partial"))

# "The detrending is the estimation and removal of the tree's natural biological growth trend." (from dplr::detrend.series)

# Ahrestani et al. (Sci Rep 2013) detrended time series before running a density-dependent population model, i.e. removing the potential signal of increasing density: sounds like a bad idea to me.
#(page 7:) "As detrending a time-series may be necessary to satisfy the assumption of statio- narity of a population time-series, and to maintain consistency across all populations, we therefore detrended all the population time-series and ran the Bayesian state space models on the detrended time series."
# Garcia-Carreras et al. (2011) detrend all time series before analysing the effect of climate on the spectral composition of GPDD time series: "All time series were linearly detrended" (p. 1044) This DOES make sense, as the focus is on fluctuations, not the trend itself.

# idea: simulate for 100 time series their variability and trend without/with density dependence; analyse with/out detrending and compare error rates
# --> should we detrend time series to investigate density dependence?



#### Detrending and changepoint ####

#Ideal detrending and differencing example (AME) fig:idealdetrend

t <- seq(1901,2000,1)

#linear trend
a<-0
b<-2
mu <- a + b*t

set.seed(12345)
y <- mu + rnorm(100,0,10)
y1 <- y-mu

par(mfrow=c(2,3))
plot(y~t, xlab="Year", type="l", ylab="Y")
abline(lm(y~t))
plot(y1~t, xlab="Year", ylab="Residuals", type="l", pch=19)
plot(diff(y)~t[-1], xlab="Year", ylab="First differences", type="l", pch=19)
pacf(y, main="", ci.col="gray")
pacf(y1, main="", ci.col="gray")
pacf(diff(y), ci.col="gray", main="")


library(ggplot2)
library(ggfortify)
library(gridExtra)
pdf(file="../figures/6_idealdetrend.pdf", width=12, height=8)
g1 <- ggplot(cbind.data.frame("Year"=t, "Y"=y), aes(Year, Y)) + 
  geom_line() + 
  theme_minimal()
g2 <- ggplot(cbind.data.frame("Year"=t, "residuals"=y1), aes(Year, residuals)) + 
  geom_line() + 
  theme_minimal()
g3 <- ggplot(cbind.data.frame("Year"=t[-1], "diff"=diff(y)), aes(Year, diff)) + 
  geom_line() + 
  theme_minimal()
g4 <- autoplot(pacf(y, plot=F), conf.int.colour="grey", ylab="partial ACF") + theme_minimal()
g5 <- autoplot(pacf(y1, plot=F), conf.int.colour="grey", ylab="partial ACF") + theme_minimal()
g6 <- autoplot(pacf(diff(y), plot=F), conf.int.colour="grey", ylab="partial ACF") + theme_minimal()
grid.arrange(g1, g2, g3, g4, g5, g6, ncol = 3)
dev.off()

###NTL data for example graph of turbidity
###Read in directly from EDI

# Package ID: knb-lter-ntl.31.32 Cataloging System:https://pasta.edirepository.org.
# Data set title: North Temperate Lakes LTER: Secchi Disk Depth; Other Auxiliary Base Crew Sample Data 1981 - current.
# Data set creator:  John Magnuson - University of Wisconsin-Madison 
# Data set creator:  Stephen Carpenter - University of Wisconsin-Madison 
# Data set creator:  Emily Stanley - University of Wisconsin-Madison 
# Metadata Provider:  NTL Information Manager - University of Wisconsin-Madison 
# Contact:    -  NTL LTER  - ntl.infomgr@gmail.com
# Stylesheet v2.11 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu 

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-ntl/31/32/d01c782e0601d2217b94dd614444bd33" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")

                   
 NTL <-read.csv(infile1,header=F 
          ,skip=1
            ,sep=","  
                ,quot='"' 
        , col.names=c(
                    "lakeid",     
                    "year4",     
                    "daynum",     
                    "sampledate",     
                    "sta",     
                    "secview",     
                    "secnview",     
                    "timeon",     
                    "timeoff",     
                    "airtemp",     
                    "windir",     
                    "windspd",     
                    "waveht",     
                    "cloud",     
                    "ice"    ), check.names=TRUE)
               
unlink(infile1)
		    
# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings
                
if (class(NTL$lakeid)!="factor") NTL$lakeid<- as.factor(NTL$lakeid)
if (class(NTL$year4)=="factor") NTL$year4 <-as.numeric(levels(NTL$year4))[as.integer(NTL$year4) ]               
if (class(NTL$year4)=="character") NTL$year4 <-as.numeric(NTL$year4)
if (class(NTL$daynum)=="factor") NTL$daynum <-as.numeric(levels(NTL$daynum))[as.integer(NTL$daynum) ]               
if (class(NTL$daynum)=="character") NTL$daynum <-as.numeric(NTL$daynum)                                   
# attempting to convert NTL$sampledate dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%d"
tmp1sampledate<-as.Date(NTL$sampledate,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1sampledate) == length(tmp1sampledate[!is.na(tmp1sampledate)])){NTL$sampledate <- tmp1sampledate } else {print("Date conversion failed for NTL$sampledate. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1sampledate) 
if (class(NTL$sta)!="factor") NTL$sta<- as.factor(NTL$sta)
if (class(NTL$secview)=="factor") NTL$secview <-as.numeric(levels(NTL$secview))[as.integer(NTL$secview) ]               
if (class(NTL$secview)=="character") NTL$secview <-as.numeric(NTL$secview)
if (class(NTL$secnview)=="factor") NTL$secnview <-as.numeric(levels(NTL$secnview))[as.integer(NTL$secnview) ]               
if (class(NTL$secnview)=="character") NTL$secnview <-as.numeric(NTL$secnview)
if (class(NTL$timeon)=="factor") NTL$timeon <-as.numeric(levels(NTL$timeon))[as.integer(NTL$timeon) ]               
if (class(NTL$timeon)=="character") NTL$timeon <-as.numeric(NTL$timeon)
if (class(NTL$timeoff)=="factor") NTL$timeoff <-as.numeric(levels(NTL$timeoff))[as.integer(NTL$timeoff) ]               
if (class(NTL$timeoff)=="character") NTL$timeoff <-as.numeric(NTL$timeoff)
if (class(NTL$airtemp)=="factor") NTL$airtemp <-as.numeric(levels(NTL$airtemp))[as.integer(NTL$airtemp) ]               
if (class(NTL$airtemp)=="character") NTL$airtemp <-as.numeric(NTL$airtemp)
if (class(NTL$windir)!="factor") NTL$windir<- as.factor(NTL$windir)
if (class(NTL$windspd)=="factor") NTL$windspd <-as.numeric(levels(NTL$windspd))[as.integer(NTL$windspd) ]               
if (class(NTL$windspd)=="character") NTL$windspd <-as.numeric(NTL$windspd)
if (class(NTL$waveht)=="factor") NTL$waveht <-as.numeric(levels(NTL$waveht))[as.integer(NTL$waveht) ]               
if (class(NTL$waveht)=="character") NTL$waveht <-as.numeric(NTL$waveht)
if (class(NTL$cloud)=="factor") NTL$cloud <-as.numeric(levels(NTL$cloud))[as.integer(NTL$cloud) ]               
if (class(NTL$cloud)=="character") NTL$cloud <-as.numeric(NTL$cloud)
if (class(NTL$ice)!="factor") NTL$ice<- as.factor(NTL$ice)


NTL.TR <- NTL[NTL$lakeid == "TR",c(1:6)]

NTL.TR.annual <- aggregate(secview ~ year4, data=NTL.TR, mean, na.rm=TRUE)
plot(secview~year4, data=NTL.TR.annual, xlab="Year", 
	ylab="Secchi Disk visibility (m)", type="o", pch=19)

library(ggplot2)
pdf(file="../figures/6_NTLsecchi.pdf", width=8, height=5)
ggplot(NTL.TR.annual, aes(year4, secview)) + 
  geom_line() + geom_point() + xlab(label="Year") + ylab(label="Secchi disk visibility [m]") + 
  theme_minimal()
dev.off()


# think: clear to turbid lake
library(changepoint)
library(dplR)
library(ggplot2)


time <- 1970:2021
rangit <- function(x) (x-min(x)) / (max(x) - min(x))
set.seed(1234)
jump <- .5 # jump size
lightextinction <- rangit(sapply((time-1970)/100, \(x) -1.2*x^5 + .5*x^4 -3*x^3 + 2.75*x^2 -1*x + .21) + rnorm(length(time), 0, .1) )
lightextinction <- lightextinction + rep(c(-jump, jump), each=26) + jump

#(-1.2*(1/100 * time-1.970)^5 + .5*(1/100 *time-1.970)^4 -3*(1/100 *time-1.970)^3 + 2.75*(1/100 *time-1.970)^2 -1*(1/100 *time-1.970))
#curve(-1.2*x^5 + .5*x^4 -3*x^3 + 2.75*x^2 -1*x, from=0, to=.52)
plot(time, lightextinction, type="s")

pdf(file="../figures/6_detrendingLakeExample.pdf", width=8, height=5)
ggplot(data.frame(time, lightextinction), aes(x=time, y=lightextinction)) + 
  geom_line() + xlab(label="Year") + ylab(label="Secchi depth visibility [m]") +
  theme_minimal()
dev.off()

ts.raw <- lightextinction 
ts.det <- detrend.series(ts.raw)[,-4] # remove AR, as it causes NAs
chpts <- apply(cbind(ts.raw, ts.det), 2, cpt.mean)
# str(chpts[[7]])
sapply(chpts, function(x) x@cpts[x@cpts != 52]) # get changepoints

# do this 100 times:
set.seed(123)
jump <- 1 # % jump
chpts.simu <- replicate(100, {
  lightextinction <- rangit(sapply((time-1970)/100, \(x) -1.2*x^5 + .5*x^4 -3*x^3 + 2.75*x^2 -1*x + .21) + rnorm(length(time), 0, .1) )
  lightextinction <- lightextinction + rep(c(-jump, jump), each=26) + jump
  ts.raw <- lightextinction 
  ts.det <- detrend.series(ts.raw)[,-4] # remove AR, as it causes NAs
  chpts <- apply(cbind(ts.raw, ts.det), 2, cpt.mean)
  sapply(chpts, function(x) x@cpts[x@cpts != 52])
})
chpts.simu
# compute number of times a changepoint was detected
apply(chpts.simu, 1, \(x) length(unlist(x)))

# do this for several values of jump:
jumps <- seq(0, 1, by=0.05)
manyjumps <- matrix(NA, ncol=length(jumps), nrow=7)
set.seed(123)
for (i in seq_along(jumps)){
  jump <- jumps[i]
  chpts.simu <- replicate(100, {
    lightextinction <- rangit(sapply((time-1970)/100, \(x) -1.2*x^5 + .5*x^4 -3*x^3 + 2.75*x^2 -1*x + .21) + rnorm(length(time), 0, .1) )
    lightextinction <- lightextinction + rep(c(-jump, jump), each=26) + jump
    ts.raw <- lightextinction 
    ts.det <- detrend.series(ts.raw)[,-4] # remove AR, as it causes NAs
    chpts <- apply(cbind(ts.raw, ts.det), 2, cpt.mean)
    sapply(chpts, function(x) x@cpts[x@cpts != 52])
  })
  manyjumps[,i] <- apply(chpts.simu, 1, \(x) length(unlist(x)))
  cat(i, " ")
}
manyjumps
save(manyjumps, file="manyjumps.RData")
matplot(t(manyjumps), type="o")

library(ggplot2)
dats <- data.frame(jumps, t(manyjumps))
colnames(dats) <- c("size", rownames(chpts.simu))

pdf(file="../figures/changepointjumps.pdf", width=5, height=5)
ggplot(dats, aes(x=size, y=ts.raw)) + 
  geom_line(aes(y=ModNegExp), color="gray40", linewidth=.8) + 
  geom_line(aes(y=AgeDepSpline), color="gray60", linewidth=.8) + 
  geom_line(linewidth=.8) +
  xlab(label="jump size") + ylab(label="percentage detected") +
  annotate("text", x=c(.3, .58, .8), y=c(95, 70, 50), label=c("age-dep. spline", "raw data", "neg. exp."), size=5, color=c("grey60", "black", "grey40")) +
  theme_minimal()
dev.off()  
  



#### AIRLINE analysis ####
air <- read.csv("../Data/AirPassengers/BTSairpassengers.csv", header=TRUE, skip=1)

#create the time series; #domestic US passengers only, in millions

air.ts <- ts(air[,3]/1000000, start=c(2002,10), frequency=12)
plot(air.ts, xlab="Year", ylab="Millions of passengers")

pdf(file="../figures/6_AirPassengers-timeseries.pdf", width=10, height=5)
autoplot(air.ts, xlab="Year", ylab="Millions of passengers") + theme_minimal()
dev.off()


pdf(file="../figures/6_AirPassengers-linear-stl.pdf", height=6)
#ggfortify:::autoplot.ts(stl(air.ts, s.window="per", t.window=1000), 
#     labels=c("Passengers\n(millions)","Seasonal","Trend","Residuals"),
#     set.pars=list(mar = c(0, 6, 0, 6), oma = c(6, 0, 4, 0),
#                   tck = -0.01, mfrow = c(4, 1), las=0)) + 
#  theme_minimal()
forecast::autoplot(stl(air.ts, s.window="per", t.window=1000)) + 
  xlab(label="Year") +
  theme_minimal()
dev.off()


#linear trend


#loess smooth for trend 
pdf(file="../figures/6_AirPassengers-robust-stl.pdf", height=6)
#autoplot(stl(air.ts, s.window="per", robust=TRUE),
#     labels=c("Passengers\n(millions)","Seasonal","Trend","Residuals"),
#     set.pars=list(mar = c(0, 6, 0, 6), oma = c(6, 0, 4, 0),
#                   tck = -0.01, mfrow = c(4, 1), las=0)) +
#  theme_minimal()
forecast::autoplot(stl(air.ts, s.window="per")) + 
  xlab(label="Year") +
  theme_minimal()
dev.off()


#compare linear trend to loess trend; information in the range bars
#from https://stats.stackexchange.com/questions/7876/interpreting-range-bars-in-rs-plot-stl
#
#
#linear: linear trend has week explanatory power here (range bar is very large
#relative to data range bar); residuals are contributing about the same information
#as the data; seasonal pattern is next on the list.
#
#loess: magnitude of trend and residual range bars are about the same as
#those of data bar, suggesting that the trend and residuals are about as important
#as the data themselves, and the contribution of the seasonal component to the 
#time series is relatively small.

plot(decompose(air.ts, type="multiplicative"))


  




#### MC time series: van Marle et al. ####
library(readxl)
dats <- read_excel("../Data/vanMarle2023/Marle_et_al_Nature_AirborneFraction_Datasheet.xlsx", sheet=6, skip=6)
colnames(dats) <- c("year", "AFnew", "AFnew.error", "AFcorrnew", "AFcorrnew.error", "AFHN", "AFHN.error", "AFHNcorr", "AFHNcorr.error", "AFGCP", "AFGCP.error", "AFGCPcorr", "AFGCPcorr.error")
head(dats)
# the python code is a bit murky: line 595 calls a function "data_selector" that is not defined; it produces a tuple with information on
# dCO2, FF and lu, from which AF is computed as dCO2/(FF+lu) (line 600)
# In the re-analysis, Bennedsen et al. use the data imported above (their "main_Figure1.m", line 34)
# that is what we do here as well, implicitly assuming that the slightly different procedure of van Marle et al. is equivalent

# plot data:
library(ggplot2)
pdf(file="../figures/vanMarleRawData.pdf", width=8, height=5)
ggplot(dats, aes(x=year, y=AFnew)) + 
  geom_line(aes(x=year, y=AFHN), color="grey30", size=.5) + 
  geom_line(aes(x=year, y=AFGCP), color="grey10", size=.5) + 
  geom_segment(aes(xend=year, yend=AFnew + 2*AFnew.error, y=AFnew - 2*AFnew.error), color="grey70") +
  ylab(label = "airborne fraction") +
  geom_line(color="grey30", size=.8) + 
  geom_point(color="grey30") +
  theme_minimal()
dev.off()

# analyse for linear trend:
summary(fm <- lm(AFnew ~ year, data=dats))
acf(residuals(fm), type="partial") # no significant TA; plot not shown

library(Kendall)
MannKendall(dats$AFnew)

curve(dnorm(x, mean=coef(fm)[2], sd=summary(fm)$coefficients[2,2]), from=-0.01, 0.01)
abline(v=0)
pnorm(0, mean=-0.0006233, sd=summary(fm)$coefficients[2,2])

# weighted regression
summary(lm(AFnew ~ year, data=dats, weights = 1/AFnew.error^2))
pnorm(0, mean=coef(fm)[2], sd=0.0009900) # 0.793

# let's say we don't trust the normality of the estimates and want to simulate them:

# simulate a new data set, according to van Marle et al.:
dats$AFprime <- dats$AFnew + rnorm(nrow(dats), mean=0, sd=dats$AFnew.error)
# reanalyse:
summary(lm(AFprime ~ year, data=dats))

# do this 10000 times and store estimate of slope:
set.seed(123)
slopes <- replicate(10000, {
  dats$AFprime <- dats$AFnew + rnorm(nrow(dats), mean=0, sd=dats$AFnew.error)
  summary(lm(AFprime ~ year, data=dats))$coefficients[2,1]
})
mean(slopes) # -0.0008
sd(slopes)
mean(slopes< 0) # 0.823
mean(slopes <= -0.0008096) # 49% have more negative slopes

plot(density(slopes))
curve(dnorm(x, mean=-0.0008096, sd=0.0009126), add=T, col="green")
abline(v=0)
# looks very similar: normality approximation seems to hold


# analysis as implied by the criticism of Bennedsen et al.:
# two ways to remove possible trend:
# (a) use residuals of original regression as Y
# (b) randomise years before analysis
# we'll use (a), as that yields only one element for simulation, not two:
dats$AFresids <- residuals(lm(AFprime ~ year, data=dats))
#plot(AFnew ~ year, dats, type="l")  
#plot(AFresids ~ year, dats, type="l")  
# add noise:
dats$AFresidsprime <- dats$AFresids + rnorm(nrow(dats), mean=0, sd=dats$AFnew.error)
# reanalyse:
summary(lm(AFresidsprime ~ year, data=dats))


# second point of Bennedsen et al.: compare test statistics, not p-values:

# van Marle report a "probability of negative trend", based on 10000 noisified time series with a trend;
# Alternative: how likely is a negative trend of at least that size for noisy but trend-less data?

set.seed(123)
resid.slopes <- replicate(10000, {
  dats$AFresidsprime <- dats$AFresids + rnorm(nrow(dats), mean=0, sd=dats$AFnew.error)
  summary(lm(AFresidsprime ~ year, data=dats))$coefficients[2,1]
})
plot(density(resid.slopes)) 
abline(v=-0.0008096)
mean(resid.slopes > -0.0008096) # 0.174

# 17% of no-trend data would have this slope of steeper; meaning 83% would have less negative (or positive) slope

plot(density(slopes)) # 80% negative
lines(density(resid.slopes), col="red") # 50% negative
curve(dnorm(x, mean=-0.0008096, sd=0.0009126), add=T, col="green") # 80% negative
abline(v=c(0, -0.0008), lty=1:2)

toplot <- cbind.data.frame(slopes, resid.slopes)
pdf(file="../figures/trendAnalysis.pdf", width=5, height=5)
ggplot(toplot) + 
  geom_density(aes(x=slopes), colour="grey70", size=1) + 
  geom_density(aes(x=resid.slopes), size=1) + 
  geom_vline(xintercept=c(0, -0.0008), linetype=1, colour=c("grey70", "black"), size=.8) +
  geom_function(fun = dnorm, args = list(mean = -0.0008096, sd = .0009126), colour = "grey20", linetype=2) +
  annotate(geom="text", x=-0.0015, y=400, label="noisy trend\n simulations", hjust="right", size=5, color="grey50") +
  annotate(geom="text", x=0.0006, y=400, label="no trend\n simulations", hjust="left", size=5) +
  annotate(geom="curve", curvature=0, x=-0.0025, y=250, xend=-.0019, yend=220, size=.5, color="grey20", linetype=2) +
  annotate(geom="text", x=-0.0025, y=250, label="traditional", hjust="right", size=4, color="grey20") +
  theme_minimal()
dev.off()

# Bayesian approach:
(0.001 + 0.0002)/2
-0.0008+0.0006
sqrt(2*(0.0005/2)^2 ) # sd of prior  
sqrt(0.0009^2+0.001^2) # sd of convolution
pnorm(0, mean=-0.0002, sd=0.001)

# using the simulations:
dnorm(-0.002, mean=0.0006, sd=0.001)
dnorm(slopes[1:2], mean=0.0006, sd=0.001)
dnorm(slopes[1:2], mean=0.0006, sd=0.001) * slopes[1:2]
weights <- dnorm(slopes, mean=0.0006, sd=0.001)
sum(weights * slopes) / sum(weights)# -0.00019
sum(weights * slopes) / sum(weights) *10000# -0.00019


10000/sum(weights)

library(ggplot2)
toplot <- data.frame(slopes, weights=dnorm(slopes, mean=0.0006, sd=0.001), ones=1)

pdf(file="../figures/trendBayes.pdf", width=5, height=5)
ggplot(toplot, aes(slopes)) + #, after_stat(density)
  geom_histogram(fill="grey30", alpha = 0.3) +  # slopes
  geom_histogram(aes(slopes, weight=weights*0.0058), fill="grey30", alpha = 0.5) + # weighted slopes
  #geom_function(fun = dnorm, args = list(mean = 0.0006, sd = .001), colour = "black", linetype=1) +
  geom_function(fun = \(x) 1.946*dnorm(x, mean=0.0006, sd = .001), colour = "black", linetype=1) + # convert density into counts!
  geom_vline(xintercept=0, linetype=1, colour="black", size=.8) +
  ylab(label="number of simulations") +
  theme_minimal()
dev.off()

# Bennedsen et al. argue that "p noisy under trend" is not a measure of "p of trend", as it is unclear what "p noisy without trend" would look like.
# That is a LRT: P(trend) / P(no trend) = P(this slope or more negative | trend)  /  P(this slope or more negative | no trend) =
# K = 0.493 / 0.174 = 2.833, i.e. bits= ln(k)=log(2.833)=1.04 = "barely worth mentioning" (https://en.wikipedia.org/wiki/Bayes_factor)
# 
# dchisq(2.833, df=1) # 0.057 BUT: chisq is the assumption we want to avoid using our simulations!


# Strikes me that these two analyses quantify something different:
# van Marle: probability that observed data exhibit a negative trend (84%)
# (the simu assumes that y.true + noise = y.obs; then use y.obs - noise to simulate y.true; since noise is symmetric around 0: y.obs - noise = y.obs + noise; fine!)
# Bennedsen: probability of observed negative trend being spurious (17%)

# Final argument: previous studies claimed a positive trend, i.e. the prior is shifted to the right
# Bayes:
# P(theta | data) = P(data | theta) * P(theta) / P(data)
# Use GCP and HN as priors!

simulate priors; multiply each by 1000 MCs above; divide by sum of all products (to sum to 1): gives probability of a negative slope, given data and previous data

library(brms)
prior <- prior("normal(.1, 1)", class=b)
fm <- brm(AFnew*1000 ~ year, data=dats, prior=prior)
summary(fm)


numerator <- function(x) dnorm(x, mean = -0.0008, sd = 0.0009) * dnorm(x, mean=0.0006, sd=0.001)
denominator <- integrate(numerator, -Inf, Inf) # 172
post <- function(x) numerator(x)/denominator$value
curve(post(x), -.002, .002)
#curve(dnorm(x, mean=-0.0008+0.0006, sd=sqrt(0.0009^2 + 0.0006^2)), add=T, col="red")
integrate(post, -Inf, 0)
integrate(post, 0, Inf)
# find mean:
optimize(post, interval=c(-.1, .1), maximum=T)






#### I need machine learning for that ####



#### Categorisation's loss of information ####
# colour categorisation: https://en.wikipedia.org/wiki/Visible_spectrum
# effect of colour on impulsivity: https://theconversation.com/does-colour-really-affect-our-mind-and-body-a-professor-of-colour-science-explains-84382
# on percieved effect of drugs: https://www.bmj.com/content/313/7072/1624?ref=fugo-blog

#--> blue, green: tranquilising; yellow, red, orange: stimulant

# put aside the problem of lightness typically coming with colours
#categorical five vs continuous wavelength: linear effect of lambda?
#vpredict to brown, white, pink, cyan: requires RGB, i.e. as many parameters as three categories!
  

library(effects)
data(TitanicSurvival)
summary(TitanicSurvival)
library(mice)
Ti <- complete(mice(TitanicSurvival, method="pmm", print=F))
Ti$survived <- ifelse(Ti$survived == "yes", 1, 0)

library(ggplot2)
pdf(file="../figures/TitanicSexAge.pdf", width=5, height=5)
ggplot(Ti, aes(age, survived, group = sex)) + 
  geom_jitter(aes(colour=sex), width = 1, height = 0.05, size=.5) +
  geom_smooth(aes(colour=sex), method="gam",
              method.args=list(family="binomial")) +
  scale_colour_grey() +
  theme_minimal() +
  theme(legend.position = c(0.1, 0.75), axis.title.x=element_blank())
dev.off()

library(Hmisc)
ageCat <- cut2(Ti$age, cuts=c(0, 18, 100))
table(Ti$survived, ageCat)

table(Ti$survived, ageCat) / matrix(tapply(Ti$survived, ageCat, length), nrow=2, ncol=2, byrow=T)

tapply(Ti$survived, list(ageCat, Ti$sex), \(x) sum(x)/length(x))
# lower survival for girls than women
# slightly higher survival for boys than men

summary(glm(survived ~ sex*age, family=binomial, data=Ti))
summary(glm(survived ~ sex*ageCat, family=binomial, data=Ti))


###### simulate categorisation ######
# simulate data
set.seed(1234)
isMale <- ifelse(Ti$sex == "male", 1, 0)
theAge <- Ti$age
theSex <- Ti$sex
survSim <- rbinom(1309, prob=plogis(0 - 1*isMale + 0.01*theAge) , size=1) #-0.05*isMale*theAge

# cut at several ages
ageCuts <- 1:21
p.vec <- NA # vector to collect p-values
ellCat.vec <- madCat.vec <- NA # ditto
for (i in 1:21){
  ageCat <- cut2(theAge, cuts=c(0, ageCuts[i], 100))
  fcat <- glm(survSim ~ ageCat, family=binomial)
  p.vec[i] <- summary(fcat)$coefficients[2,4]
  fitCat <- predict(fcat, type="response")
  ellCat.vec[i] <- sum(dbinom(survSim, prob=fitCat, size=1, log=T))
  madCat.vec[i] <- mean(abs(survSim - fitCat))
  rm(fcat, ageCat)
}
fcont <- glm(survSim ~ theAge, family=binomial)
pCont <- summary(fcont)$coefficients[2,4]
fitCont <- predict(fcont, type="response")
ellCont <- sum(dbinom(survSim, prob=fitCat, size=1, log=T))
madCont <- mean(abs(survSim - fitCont))

plot(1:21, p.vec, type="l", lwd=2)
abline(h=pCont)

pdf(file="../figures/TitanicAgeCutoffSimu1Run.pdf", width=9, height=5)
ggplot(as.data.frame(p.vec), aes(x=1:21, y=p.vec)) + 
  geom_line() + 
  labs(x="age threshold for 'young'") +
  labs(y="p-value of age effect") +
  scale_x_continuous(breaks=1:21, labels=as.character(1:21)) +
  geom_hline(yintercept=pCont, col="grey") + 
  theme_minimal()
dev.off()

#plot(1:21, ellCat.vec, type="l", lwd=2, ylim=c(-899, -894))
#abline(h=ellCont)

# compute error rates
(isSignificant <- ifelse(c(p.vec, pCont) < 0.05, 1, 0))
# MAD
(MAD <- c(madCat.vec, madCont))


##### replicate many times #####
isMale <- ifelse(Ti$sex == "male", 1, 0)
theAge <- Ti$age
theSex <- Ti$sex
ageCuts <- 1:21

set.seed(123)
sims <- replicate(100, {
  # simulate data
  survSim <- rbinom(1309, prob=plogis(0 - 1*isMale + 0.02*theAge) , size=1) #-0.05*isMale*theAge
  
  # cut at several ages
  p.vec <- NA
  ellCat.vec <- madCat.vec <- NA
  for (i in 1:21){
    ageCat <- cut2(theAge, cuts=c(0, ageCuts[i], 100))
    fcat <- glm(survSim ~ ageCat, family=binomial)
    p.vec[i] <- summary(fcat)$coefficients[2,4]
    fitCat <- predict(fcat, type="response")
    ellCat.vec[i] <- sum(dbinom(survSim, prob=fitCat, size=1, log=T))
    madCat.vec[i] <- mean(abs(survSim - fitCat))
    rm(fcat, ageCat, fitCat)
  }
  fcont <- glm(survSim ~ theAge, family=binomial)
  pCont <- summary(fcont)$coefficients[2,4]
  fitCont <- predict(fcont, type="response")
  ellCont <- sum(dbinom(survSim, prob=fitCat, size=1, log=T))
  madCont <- mean(abs(survSim - fitCont))
  
  isSignificant <- ifelse(c(p.vec, pCont) < 0.05, 1, 0)
  MAD <- c(madCat.vec, madCont)
  rbind(isSignificant, MAD)
})
rowMeans(sims[1,,])
rowMeans(sims[2,,])


# plot MAD:
dats <- stack(as.data.frame(t(sims[2,,])))
library(ggplot2)
pdf(file="../figures/TitanicAgeCutoffSimu.pdf", width=9, height=5)
ggplot(dats, aes(x=ind, y=values)) + 
  geom_boxplot(fill="grey80")  + 
  labs(x="age threshold for 'young'") +
  labs(y="error (mean absolute difference)") +
  scale_x_discrete(labels=c(as.character(1:21), "continuous")) +
  theme_minimal() 
dev.off()


dats <- data.frame("vals"=rowMeans(sims[1,,]), "thresh"=c(as.character(1:21), "continuous"))
dats$thresh <- factor(dats$thresh, levels=dats$thresh)
pdf(file="../figures/TitanicAgeCutoffSimuPvalues.pdf", width=9, height=5)
ggplot(dats, aes(x=thresh, y=vals)) + 
  geom_col(fill="grey") +
  labs(x="age threshold for 'young'") +
  labs(y="proportion significant") +
  theme_minimal() 
dev.off()






#### unused code ####
# same as Titanic, but with a uniform predictor and normal response
# compares regression, cuts into equidistant bins, and cuts into equal-sized bins

set.seed(123)
x <- round(runif(34), 2)
y <- 2 + 2*x + rnorm(34, 0, 1)
hist(x)
# make g groups
library(Hmisc)
ind <- cut2(x, g=3)
xequal <- factor(ind, labels=LETTERS[1:3])
table(xequal)
# make g equal-distant:
ind2 <- cut2(x, cuts=seq(0, 1, len=3+1))
table(ind2)
xequi <- factor(ind2, labels=LETTERS[1:3])
table(xequi)

summary(fm1 <- lm(y ~ x))
summary(fm2 <- aov(y ~ xequal))
summary(fm3 <- aov(y ~ xequi))

xnew <- seq(0,1,len=21)

predict(fm1, newdata=data.frame(x=xnew))
predict(fm2, newdata=data.frame(xequal="B"))
predict(fm3, newdata=data.frame(xequi="B"))


