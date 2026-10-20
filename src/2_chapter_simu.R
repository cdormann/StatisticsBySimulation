## R code for chapter introducing simulations


#### Polya's urn ####
set.seed(123)
urn <- c("r", "b") # one red, one blue
N <- 100
for (i in 1:N) { urn <- c(urn, sample(urn, size=1)) }
propBlue <- cumsum(ifelse(urn[-c(1:2)]=="b", 1, 0))/c(1:100)
#plot(propBlue, type="o", pch=16, las=1)
library(ggplot2)
ggplot(data.frame(blue=propBlue, trials=1:100), aes(x=trials, y=blue)) + 
  geom_point() + geom_line() + labs(y="proportion blue") + theme_minimal()
table(urn[-c(1,2)])
table(urn)["b"] / N

set.seed(321)
urn2 <- c("r", "b") # one red, one blue
N <- 100
for (i in 1:N) { urn2 <- c(urn2, sample(urn2, size=1)) }
table(urn2[-c(1,2)])
propBlue2 <- cumsum(ifelse(urn2[-c(1:2)]=="b", 1, 0))/c(1:100)
runs <- data.frame(blue=c(propBlue, propBlue2), trials=rep(1:100, 2), run=rep(1:2, each=100))
pdf(file="../figures/blues.pdf")
ggplot(runs, aes(x=trials, y=blue, group=run)) + geom_point() + 
  geom_line() + labs(y="proportion blue") + theme_minimal()
dev.off()

blueProportion <- 1
for (j in 1:5000){
  set.seed(2+j)
  urn <- c("r", "b") # one red, one blue
  for (i in 1:N) { urn <- c(urn, sample(urn, size=1)) }
  blueProportion[j] <- table(urn[-c(1:2)])["b"]/ N
}
sum(is.na(blueProportion)) / 5000 # 1% without blue!
sum(blueProportion == 1, na.rm=T) / 5000 # 1% only blue!
#hist(blueProportion)
pdf(file="../figures/bluesHisto.pdf")
ggplot(data.frame(blue=blueProportion), aes(x=blue)) + geom_histogram(bins=20) + theme_minimal()
dev.off()




#### Confidence interval ####
mean.true <- 10
sd <- 2

obs <- rnorm(10, mean=mean.true, sd=sd)
fact <- qt(0.975, df=9) # Find the appropriate quantile of the t distribution for computing the 95% CI
fact

c(mean(obs) - fact*sd(obs)/sqrt(10), mean(obs), mean(obs) + fact*sd(obs)/sqrt(10))

t.test(obs)$conf.int # Confirms the previous computation

set.seed(12345)
simu <- replicate(100, { # 100 Simulations
  obs <- rnorm(10, mean=mean.true, sd=sd)
  c("lower"=mean(obs) - fact*sd(obs)/sqrt(10), "mean"=mean(obs), "upper"=mean(obs) + fact*sd(obs)/sqrt(10))
})
simu <- cbind.data.frame("N"=1:100, t(simu))
length(which(simu$lower > 10 | simu$upper < 10))

# not show in code box:
# find non-overlapping samples
out <- which(simu$lower > 10 | simu$upper < 10)
g1 <- simu[out, ]

library(ggplot2)
pdf(file="../figures/SEplot.pdf")
ggplot(simu, aes(mean, N))+        # ggplot2 plot with confidence intervals
  geom_errorbar(aes(xmin = lower, xmax = upper), color="grey70") +
  geom_point(color="grey40") +
  geom_vline(xintercept=10, color="grey50") +
  geom_point(data=g1, colour="black") +
  geom_errorbar(data=g1, aes(xmin = lower, xmax = upper), colour="black") +
  theme_minimal() # theme_classic()
dev.off()



#### KS and Lilliefors ####

library(nortest)
# once to see that the code works:
A <- rnorm(10) # we could use any other mean and sd, of course
ks.test(A, "pnorm", mean=mean(A), sd=sd(A))

lillie.test(A) 

# now repeat 1000 times:
set.seed(234)
res <- replicate(1000, {
  A <- rnorm(10)
  c(ks.test(A, "pnorm", mean=mean(A), sd=sd(A))$p.value, lillie.test(A)$p.value)
})
rowSums(res < 0.05)



#### logistic regression uncertainty #### 
set.seed(345)
## this if we want to simulate data and fit:
sims <- replicate(100, {
  X <- runif(123, 0 ,1)
  Y <- rbinom(123, size=1, prob=plogis(1 - 4*X))
  coef(glm(Y ~ X, family=binomial))
})

## this if we want to start with a fitted model:
set.seed(345)
X <- runif(123, 0 ,1)
Y <- rbinom(123, size=1, prob=plogis(1 - 2*X))
fm <- glm(Y ~ X, family=binomial)
vcov(fm)
library(mvtnorm)
sims <- rmvnorm(n=100, mean=coef(fm), sigma=vcov(fm))

# par(las=1, pch=16, mar=c(4,4,1,1))
#curve(plogis(sims[1,1] + sims[2,1]*x), add=F, lwd=2, col=rgb(.2, .2, .2, .3), ylim=c(0,1), ylab="y")
#sapply(2:100, \(k) curve(plogis(sims[1,k] + sims[2,k]*x), add=T, lwd=2, col=rgb(.2, .2, .2, .3)))

library(purrr) # for pmap and pipe
library(dplyr) # for tibble
library(ggplot2)

fun <- function(x, a, b) plogis(a + b*x)
coefs <- as.data.frame(sims) 
colnames(coefs) <- c("intercept", "slope")

curves <- coefs %>% 
  pmap(function(...){
    dots <- tibble(...)
    stat_function(data = dots, aes(0),
                  fun = function(x) fun(x, dots$intercept, dots$slope), 
                  xlim = c(0, 1),
                  color=rgb(.2, .2, .2, .3))
  })

Xnew <- seq(0, 1, len=100)
preds <- predict(fm, newdata=data.frame("X"=Xnew), se.fit=T)
CIs <- cbind.data.frame("X"=Xnew, "upper"=plogis(preds$fit+2*preds$se.fit), "lower"=plogis(preds$fit-2*preds$se.fit), "mean"=plogis(preds$fit))

pdf(file="../figures/bernoulliMultiple.pdf", width=6, height=5)
ggplot() + 
  curves + scale_y_continuous(expand = expansion(), limits = c(0, 1))  +  theme_minimal()
dev.off()

# with underlying CI ribbon:
ggplot(CIs, aes(x=X)) + 
  geom_ribbon(aes(ymax=upper, ymin=lower, xmin=0, xmax=1), fill="grey80") + geom_line(aes(x=X, y=mean)) +
  curves + scale_y_continuous(expand = expansion(), limits = c(0, 1))  +  theme_minimal()

# only classical CI:
pdf(file="../figures/bernoulliCI.pdf", width=6, height=5)
ggplot(CIs, aes(x=X, y=mean)) + 
  geom_ribbon(aes(ymax=upper, ymin=lower, xmin=0, xmax=1), fill="grey80") + 
  geom_line() +  
  theme_minimal()
dev.off()

## this for the bootstrap:
set.seed(345)
X <- runif(123, 0 ,1)
Y <- rbinom(123, size=1, prob=plogis(1 - 2*X))
dats <- cbind.data.frame(X, Y)
fm <- glm(Y ~ X, family=binomial, data=dats)
#single bootstrap run:
coef(glm(Y ~ X, family=binomial, data=dats[sample(nrow(dats), replace=T),]))
# 100 replicates:
boots <- replicate(100, coef(glm(Y ~ X, family=binomial, data=dats[sample(nrow(dats), replace=T),])))

## plot for bootstraps:
coefs.boot <- as.data.frame(t(boots))
colnames(coefs.boot) <- c("intercept", "slope")

curves <- coefs.boot %>% 
  pmap(function(...){
    dots <- tibble(...)
    stat_function(data = dots, aes(0),
                  fun = function(x) fun(x, dots$intercept, dots$slope), 
                  xlim = c(0, 1),
                  color=rgb(.2, .2, .2, .3))
  })

pdf(file="../figures/bernoulliBoots.pdf", width=6, height=5)
ggplot() + curves + scale_y_continuous(expand = expansion(), limits = c(0, 1)) + theme_minimal()
dev.off()



#### t-test and normality ####
#hist(rlnorm(100, 2.3, .5))
set.seed(123)
obs <- data.frame("group"=rep(c("A", "B"), each=10), "surface"=round(c(rlnorm(10, 2, 0.8), rlnorm(10, 2.9, 1) ), 1 ) )
write.csv(obs, file="tTestData.csv", row.names = F)

obs <- read.csv("../Data/tTestData.csv")
with(obs, t.test(surface ~ group)) 
#boxplot(surface ~ group, data=obs, las=1)
#hist(residuals(lm(surface ~ group, data=obs)))


mean(obs$surface)
sd(obs$surface)

set.seed(12345)
# how often is a difference described as significant?
simu <- replicate(1000, {
  sim <- data.frame("group"=rep(c("A", "B"), each=10), "surface"=rnorm(20, mean(obs$surface), sd(obs$surface))) 
  with(sim, t.test(surface ~ group))$p.value
})
sum(simu < 0.05)/1000 # 0.05, as it should be

# get parameters for lognormal:
library(MASS)
fit.lognormal <- fitdistr(obs$surface, "log-normal")

# same now for lognormally generate data:
simu.lognormal <- replicate(1000, {
  sim <- data.frame("group"=rep(c("A", "B"), each=10), "surface"=rlnorm(20, fit.lognormal$estimate[1], fit.lognormal$estimate[2])) 
  with(sim, t.test(surface ~ group))$p.value
})
sum(simu.lognormal < 0.05)/1000 # 0.027, while p-value of test was 0.17!


# now make a density plot of the actual t-statistic for both sampling schemes:
simu <- replicate(1000, {
  sim <- data.frame("group"=rep(c("A", "B"), each=10), "surface"=rnorm(20, mean(obs$surface), sd(obs$surface))) 
  with(sim, t.test(surface ~ group))$statistic
})

plot(density(simu), las=1)
curve(dt(x, df=18), add=T, col="green")
# fits perfectly


simu.lognormal <- replicate(1000, {
  sim <- data.frame("group"=rep(c("A", "B"), each=10), "surface"=rlnorm(20, fit.lognormal$estimate[1], fit.lognormal$estimate[2])) 
  with(sim, t.test(surface ~ group))$statistic
})

lines(density(simu.lognormal), las=1, col="red")

d <- data.frame("distr"=rep(c("normal", "log-normal"), each=1000), "t"=c(simu, simu.lognormal))
library(ggplot2)
pdf(file="../figures/tTestPlot.pdf", width=5, height=5)
ggplot(d, aes(x=t, color=distr)) +
  stat_function(fun=\(x) dt(x, df=18), show.legend=F, color="grey90", linetype=1, size=1.25) +
  geom_density(show.legend=F, size=1.25) +
  scale_color_grey() +
  geom_vline(xintercept=-2.2404, color="grey", size=1.25) +
  theme_classic()
dev.off()

# # now visualise 1000 sample-means for normal and lognormal in a boxplot
# diff(tapply(rnorm(20, mean(obs$surface), sd(obs$surface)), obs$group, mean))
# simu <- replicate(1000, diff(tapply(rnorm(20, mean(obs$surface), sd(obs$surface)), obs$group, mean)))
# simu.lognormal <- replicate(1000, diff(tapply(rlnorm(20, fit.lognormal$estimate[1], fit.lognormal$estimate[2]), obs$group, mean)))
# boxplot(cbind(simu, simu.lognormal))
# apply(cbind(simu, simu.lognormal), 2, sd)




#### temporal autocorrelation ####
set.seed(12345)
rho = 0.7; sigma = 2; ts.mean.target <- 0
a <- ts.mean.target * (1-rho)
y <- rnorm(1, ts.mean.target, sigma) 
for (i in 2:100) y[i] = a + rho* y[i-1] + rnorm(1, 0, sigma)

plot(y, type="l") 
abline(h=ts.mean.target)
acf(y)
acf(y, type="partial")


library(ggfortify)
p1 <- autoplot(acf(y, plot = FALSE), conf.int.fill = '#0000FF', conf.int.value = 0.8, conf.int.type = 'ma') 
print(p1) 

pdf(file="../figures/tempAutocorPACF.pdf")
forecast::ggtsdisplay(y, lag.max=10, theme=theme_minimal(), ylab="y")
dev.off()


##### wrongly assuming independence #####
summary(lm(y ~ 1))

##### correctly accounting for temporal autocorrelation #####
library(nlme)
summary(gls(y ~ 1, correlation=corAR1(form = ~1)))






# or, using the multivariate approach:
sigma = 2
rho = 0.7
R = diag(rep(1, 100))
diag(R[-nrow(R), -1]) <- sqrt(rho)
diag(R[-1, -nrow(R)]) <- sqrt(rho)
R[1:5, 1:5]
Sigma <- diag(rep(sigma, 100)) %*% R %*%  diag(rep(sigma, 100))
Sigma[1:5, 1:5]
# this thing is typically NOT positive definite: grrrr

library(mvtnorm)
#set.seed(345)
y <- as.vector(rmvnorm(n=1, mean=rep(0, 100), sigma=Sigma))
acf(y)
acf(y, type="partial")

