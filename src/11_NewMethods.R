## New methods



#### GDF ####
set.seed(123)
x1 <- 0:100
x2 <- sample(0:100)
#x3 <- sample(0:100)
x3 <- x2 + rnorm(101, sd=20) #highly correlated with x2 --> affects GDF?
cor(x2, x3)
x4 <- sample(0:100)
x5 <- sample(0:100)
x6 <- sample(0:100)
y <- rnorm(101, mean=2 + 2*x1 -0.5*x2 + 1*x3, sd=30)
dats <- cbind.data.frame(y, x1, x2, x3, x4, x5, x6)
#plot(x, y)
summary(fm <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6, data=dats))

##### gdf.lm #####

gdf.lm <- function(object, N=200, noise=0.1){
  # computes GDF for a linear model, following Elder 2003
  # CFD 27-Oct-2022
  # add some noise to the y and re-run original model, N times, storing both perturbed y (=yprime) and fitted yhat:
  yextracted <- model.frame(object)[,1]
  n <- length(yextracted)
  preds.prime <- replicate(N, {
    LHS <- all.vars(object$call)[1]
    RHS <- paste(names(coef(object)[-1]), collapse = ' + ') 
    yprime <- yextracted + rnorm(n, sd=noise*sd(yextracted))
    newfm <- update(object, yprime ~ ., data=cbind.data.frame("yprime"=yprime, object$model)) 
    cbind(yprime, "yhat"=fitted.values(newfm))
  }
  ) # an array of n x 2 x N
  elasticities <- apply(preds.prime, 1, function(x) coef(lm(x[2,] ~ x[1,])))
  c("GDF"=unname(rowSums(elasticities)[2])) # GDF!
}
gdf.lm(fm)

library(mgcv)
fgam1 <- gam(y ~ s(x1) + s(x2) + s(x3) + s(x4) + s(x5) + s(x6))
fgam2 <- gam(y ~ s(x1, bs="cs") + s(x2, bs="cs") + s(x3, bs="cs") + s(x4, bs="cs") + s(x5, bs="cs") + s(x6, bs="cs"))
summary(fgam1) # 7 (6 + intercept)
summary(fgam2) # 7 (6 + intercept)
gdf.lm(fgam1)
gdf.lm(fgam2)

library(randomForest)
library(ranger)
fra <- ranger(y ~ x1 + x2 + x3 + x4 + x5 + x6, data=dats) # no update for ranger
frf <- randomForest(y ~ x1 + x2 + x3 + x4 + x5 + x6, data=dats) # no update for ranger

gdf.rf <- function(object, N=50, noise=0.1, data, type="response"){
  # computes GDF for a randomForest or ranger model, following Elder 2003
  # since neither stores the data in the object (by default), data have to be handed over as argument;
  # note that ranger is typically substantially faster than randomForest
  # CFD 12-Nov-2022
  # add some noise to the y and re-run original model, N times, storing both perturbed y (=yprime) and fitted yhat:
  if (!any(class(object) %in% c("ranger", "randomForest"))) stop("This function currently only handles ranger and randomForest objects!")
  
  if (is(object, "randomForest")){
    # for randomForest objects
    yextracted <- model.frame(object)[,1]
    RHS <- strsplit(as.character(object$call$formula), "-", fixed=T)[[3]]
  }
  if (is(object, "ranger")){
    # for ranger objects
    split1 <- strsplit(as.character(object$call), "-", fixed=T)
    split2 <- strsplit(split1[[2]], " ~ ", fixed=TRUE)
    RHS <- split2[[1]][2]
    yextracted <- data[, split2[[1]][1]]
  }
  n <- length(yextracted)
  if (is(object, "randomForest")){
    preds.prime <- replicate(N, {
      yprime <- yextracted + rnorm(n, sd=noise*sd(yextracted))
      dats <- cbind.data.frame("yprime"=yprime, data)
      newfrf <- randomForest(as.formula(paste("yprime ~ ", RHS)), data=dats)
      fitted.values.rf <- unname(predict(newfrf, newdata=dats, type=type))
      cbind(yprime, "yhat"=fitted.values.rf)
    }) # an array of n x 2 x N
  }
  if (is(object, "ranger")){
    preds.prime <- replicate(N, {
      yprime <- yextracted + rnorm(n, sd=noise*sd(yextracted))
      dats <- cbind.data.frame("yprime"=yprime, data)
      newfrf <- ranger(as.formula(paste("yprime ~ ", RHS)), data=dats)
      fitted.values.rf <- predict(newfrf, data=dats, type=type)
      cbind(yprime, "yhat"=fitted.values.rf$predictions)
    }) # an array of n x 2 x N
  }
  #lm(preds.prime[1,2,] ~ preds.prime[1,1,]) # Yehat vs Ye for each data point
  elasticities <- apply(preds.prime, 1, function(x) coef(lm(x[2,] ~ x[1,])))
  #elasticities <- 1
  #for (i in 1:n){ # loop over data points
  #   elasticities[i] <- coef(lm(preds.prime[i,2,] ~ preds.prime[i,1,]))[2] # Yehat vs Ye for 
  #}
  c("GDF"=unname(rowSums(elasticities)[2])) # GDF!
}

gdf.rf(fra, data=dats)
gdf.rf(frf, data=dats)

library(nnet)
set.seed(123)
fnn1 <- nnet(y ~ x1 + x2 + x3 + x4 + x5 + x6, data=dats, linout=TRUE, size=7, decay=.1, maxit=1E5, trace=F)
fnn2 <- nnet(y ~ x1 + x2 + x3 + x4 + x5 + x6, data=dats, linout=TRUE, size=30, decay=.1, maxit=1E5, trace=F)

gdf.nnet <- function(object, N=200, noise=0.1, data){
  # computes GDF for a nnet neural network, following Elder 2003
  # CFD 27-Oct-2022
  # add some noise to the y and re-run original model, N times, storing both perturbed y (=yprime) and fitted yhat:
  
  # some sanity checks:
  if (!any(class(object) %in% "nnet")) stop("This function currently only handles nnet objects!")
  if (!object$call$linout) stop("This approach to GDF only makes sense for regressions. Since linout=F, this is not the case here.")
  if (object$convergence != 0) stop("Seems like this neural net did not converge!")
  
  yextracted <- model.frame(object)[,1]
  n <- length(yextracted)
  preds.prime <- replicate(N, {
    LHS <- all.vars(object$call)[1]
    RHS <- strsplit(as.character(object$call$formula), "-", fixed=T)[[3]]
    yprime <- yextracted + rnorm(n, sd=noise*sd(yextracted))
    dats <- cbind.data.frame("yprime"=yprime, data)
    newfnn <- update(object, yprime ~ ., data=dats, trace=F) 
    cbind(yprime, "yhat"=fitted.values(newfnn))
  }
  ) # an array of n x 2 x N
  elasticities <- apply(preds.prime, 1, function(x) coef(lm(x[2,] ~ x[1,])))
  c("GDF"=unname(rowSums(elasticities)[2])) # GDF!
}

(gdf.fnn1 <- gdf.nnet(fnn1, data=dats))
(gdf.fnn2 <- gdf.nnet(fnn2, data=dats))


# compute likelihood fit:
ell <- function(obs, preds){
  sum(dnorm(obs, mean=preds, sd=sd(obs-preds), log=T))
}
-2*ell(dats$y, predict(fm)) + 2*gdf.lm(fm)
-2*ell(dats$y, predict(fgam1)) +2*gdf.lm(fgam1)
-2*ell(dats$y, predict(fgam2)) +2*gdf.lm(fgam2)
-2*ell(dats$y, predict(fra, data=dats)$prediction) +2*gdf.rf(fra, data=dats)
-2*ell(dats$y, predict(fnn1)) +2*gdf.fnn1
-2*ell(dats$y, predict(fnn2)) +2*gdf.fnn2



#### propagation of uncertainty: averaging ####
set.seed(123)
x1 <- rnorm(1000, mean=-3, sd=3)
x2 <- rnorm(1000, mean=3, sd=5)

# mean of two normals:
mean.sim <- 0.5 * sample(x1, 1000, replace=T) + 0.5 * sample(x2, 1000, replace=T)
mean(mean.sim)
sd(mean.sim)

# one might think that the mean-distribution may also have a "mean" standard deviation, i.e. (5+5)/2 = 5
# that is not the case; variances sum, not standard deviations;
# and for the weighted sum we are looking at a specific case of ...
# propagation of uncertainty (https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulae):
# var(aA + bB) = a^2*var(A) + b^2*var(B) + 2abcov(A, B)
# since cor(x1, x2) =0, only the first two terms are relevant:
(sdmean <- sqrt(0.5^2 * var(x1) + 0.5^2 * var(x2)))

library(ggplot2)
dats <- as.data.frame(mean.sim)
pdf(file="../figures/errorpropAverage.pdf", width=5, height=5)
ggplot(dats, aes(x=mean.sim)) + geom_density(linewidth=1) + 
  geom_function(fun = dnorm, args = list(mean = -3, sd = 3), colour = "grey20", linetype=2) +
  geom_function(fun = dnorm, args = list(mean = 3, sd = 5), colour = "grey20", linetype=2) +
  geom_function(fun = dnorm, args = list(mean = 0, sd = sdmean), colour = "grey50", linetype=1) +
  theme_minimal()
dev.off()



#### Shannon ####
library(vegan)
c1 <- c(1000, 10, 10, 10, 1)
c2 <- c(10, 10, 10, 1) # abundant species lost
d1 <- c(10, 10, 10, 10, 1)
d2 <- c(10, 10, 10, 10, 1, 1000) #abundant species gained 

round(sapply(list(c1, c2, d1, d2), diversity), 2)

exp(sapply(list(c1, c2, d1, d2), diversity)) # no, exp does not help: try harder


##### Bolger analysis and DHARMa ####
bolger <- read.delim("../Data/bolger.txt") 
fm <- glm(RODENTSP ~  PERSHRUB, family=binomial, data=bolger)
summary(fm)
# plot(predict(fm), residuals(fm), las=1, pch=16, xlab="link-scale prediction", ylab="")
pdf(file="../figures/binomDiagnosticPlot.pdf", width=5, height=5.5)
plot(predict(fm), rstudent(fm), las=1, pch=16, xlab="link-scale prediction", ylab="studentised residuals")
abline(h=0, col="grey")
dev.off()

library(DHARMa)
par(xpd=T)
sim <- simulateResiduals(fm)
pdf(file="../figures/binomDHARMaPlot.pdf", width=10, height=5.5)
plot(sim)
dev.off()




#### code not used ####

#### new idea tested: spiking calibration curve ####
# 5-parameter asymmetric sigmoidal:
# {y}_j=a+\frac{d-a}{{\left(1+{\left(\frac{x_j}{c}\right)}^b\right)}^g} # https://link.springer.com/article/10.1208/s12248-017-0159-4
##### generated data: #####
a <- 100 # upper asymptote
b <- 2   # growth factor
c <- 40  # concentration at inflection point = half-distance d to a
d <- 10  # lower asymptote
g <- 1   # asymmetry factor

cal <- function(par=c("a"=100, "b"=3, "c"=40, "d"=10, "g"=1.5), x){
  a <- par["a"] # upper asymptote
  b <- par["b"] # growth factor
  c <- par["c"] # concentration at inflection point = half-distance d to a
  d <- par["d"] # lower asymptote
  g <- par["g"] # asymmetry factor
  
  a + (d-a) / (1 + (x/c)^b)^g
}
set.seed(123)
x <- seq(0, 100, len=100)
y <- cal(x=x) + rnorm(length(x), mean=0, sd=1)
plot(x, y, type="p")

##### spike with errors #####
p <- 0.1 # proportion of errors
# replace random points with values of 100:
y[sample(seq_along(x), round(p*length(x)))] <- 100
plot(x, y, type="l")

##### idea: analyse as mixture model of normal and binomial #####
# show that optim finds correct values without spiking:
optfun <- function(par, x, y){
  yhat <- cal(par, x)
  -sum(dnorm(y, mean=yhat, sd=par["sigma"], log=T))
}
optim(par=c("a"=100, "b"=3, "c"=40, "d"=10, "g"=1.5, "sigma"=1), fn=optfun, x=x, y=y, method="BFGS") # fine

optfunSpike <- function(par, x, y){
  yhat <- cal(par, x)
  p <- plogis(par["p"])
  sigma <- exp(par["sigma"])
  
  ell <- ifelse(y==100, {p*dbinom(1, size=1, prob=p) + (1-p)*dnorm(100, mean=yhat, sd=sigma)}, {dnorm(y, mean=yhat, sd=sigma)})
  -sum(ell)
}
optfunSpike(par=c("a"=100, "b"=3, "c"=40, "d"=10, "g"=1.5, "p"=.1, "sigma"=1), x=x, y=y)

(op <- optim(par=c("a"=100, "b"=3, "c"=40, "d"=10, "g"=1.5, "p"=-2.2, "sigma"=0), fn=optfunSpike, x=x, y=y, method="SANN"))
plot(x,y)
curve(cal(op$par, x), add=T, col="orange")

# So: works, in principle, but is extremely sensitive to starting values and hence not usable (not even with SANN)!



##### product of two distributions (e.g. in Bayes) #####
# to illustrate a failed example of understanding something (not didactically very useful unless a generic problem!)
# Idea: Convolution gives product of two distributions as required in Bayes (spoiler: it does not)
# (This idea may well be returned by an internet search!)
# in this case, we weight x1 by x2, which is equivalent to a convolution:
# N()+N()=N((mu1+mu2)/2, sd=sqrt(sum(sigma^2)))
# we use x2 as prior, and x1 as likelihood
set.seed(123)
x1 <- rnorm(1000, mean=-3, sd=3)
w <- dnorm(x1, mean=3, sd=5) # probability density of x1-observations under x2 (=N(3,5)) prior
# product of the two functions needs to be standardised to have an area of 1; so we need to compute what the area of their product is:
# following the example of https://stats.stackexchange.com/questions/86472/posterior-very-different-to-prior-and-likelihood:
numerator <- function(x) dnorm(x, mean = -3, sd = 3)*dnorm(x, mean=3, sd=5)
denominator <- integrate(numerator, -Inf, Inf) # 0.040
post <- function(x) numerator(x)/denominator$value

# so, for each observation we now also have a weight!
dats <- data.frame(x1, w)
ggplot(dats) + geom_density(aes(x=x1, weight=w), linewidth=1, colour="black") + 
  geom_function(fun = dnorm, args = list(mean = -3, sd = 3), colour = "grey20", linetype=2) +
  geom_function(fun = dnorm, args = list(mean = 3, sd = 5), colour = "grey20", linetype=2) +
  # geom_function(fun = dnorm, args = list(mean = -2.9, sd = .09), colour = "green", linetype=2) +
  geom_function(fun = post, colour = "grey50") +
  theme_minimal()

# see conjugate prior computation for known sigma^2: https://en.wikipedia.org/wiki/Conjugate_prior
# prior: N(3, 5); likelihood: N(-3, 3), 1000 data points:
# mean:
(1/5^2 + 1000/3^2)^(-1) * (3/5^2 + sum(x1)/3^2) # -2.7
# sd:
sqrt((1/5^2 + 1000/3^2)^(-1)) # 0.095 
# no, can't be!
# error in thinking: "3" is not the standard deviation of the likelihood! 
# The sigma of the underlying distribution from which the mean is estimated is VERY wide, so that the resulting standard error of 3 emerges with 1000 samples, i.e. 3*sqrt(1000) = 95!!



