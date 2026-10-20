
sapply(c("effectsize", "ggplot2", "mgcv", "mvtnorm", "tidyr"), library, character.only = TRUE)

#### prospective sample size for a t-test ####
# this is to demonstrate that simulation will yield similar values as power.t.test

# decide on: delta, sigma, alpha and power
# I want to be able to detect a difference between 2 groups of 5 units, when the standard deviation (or noise) is around 4. Alpha and power stay at their typical values of 0.05 and 0.8, respectively.
# In other words, I want to know the sample size that gives me a 0.05-significance in 80% of the cases that a true difference of at least delta (with some sd) exists.

delta=5; sigma=4; alpha=0.05; power=0.8
# initial run:
N <- 10 
y <-  rep(c(0, delta), each=N) + rnorm(2*N, mean=delta, sd=sigma)
x <- rep(c("A", "B"), each=N)
t.test(y ~ x, var.equal=T)

# many runs:
Ps <- replicate(1000, {
  y <-  rep(c(0, delta), each=N) + rnorm(2*N, mean=delta, sd=sigma)
  x <- rep(c("A", "B"), each=N)
  t.test(y ~ x, var.equal=T)$p.value
})
sum(Ps < 0.05)

# many runs for different values of N:
set.seed(123)
prop.less.5percent <- 1
for (i in 3:20){
  N <- i
  Ps <- replicate(1000, {
    y <-  rep(c(0, delta), each=N) + rnorm(2*N, mean=delta, sd=sigma)
    x <- rep(c("A", "B"), each=N)
    t.test(y ~ x, var.equal=T)$p.value
  })
  prop.less.5percent[i] <- sum(Ps < 0.05)
}
plot(3:20, prop.less.5percent[3:20], type="o", pch=16, las=1, ylab="proportion of P<0.05", xlab="sample size")
abline(h=c(0.5, 0.4, 0.3, 0.2, 0.1)*1000, col="blue")
text(20, 550, "power:")
text(20, c(0.5, 0.4, 0.3, 0.2, 0.1)*1000, c("50%", "40%", "30%", "20%", "10%"))
# with 11 or so samples (per group), 80% of the simulations return P-values lower than 0.05.
# with 18 or more samples per group we reach the 0.95 threshold for power

# compare to power.t.test with low values for power
power.t.test(delta=5, sd=4, sig.level=0.05, power=0.8)
# --> n=11, so very similar

## --> simulation works fine
library(ggplot2)
dats <- data.frame(x=3:20, y=prop.less.5percent[3:20]/1000)
pdf(file="../figures/prospective_ttest_simu.pdf", width=5, height=5)
ggplot(dats, aes(x, y)) + geom_point(size=2) + geom_line() + 
  xlab("sample size") + ylab("proportion of tests significant at P<0.05") +
  geom_hline(yintercept=c(.5, .6, .7, .8, .9), col="grey20", linewidth=.75) + 
  annotate("text", x=20, y=0.475, label= "Power  ") + 
  annotate("text", x=20, y=c(.5, .6,.7, .8, .9) + .015, label= c("50%", "60%", "70%", "80%", "90%")) + 
  theme_bw()
dev.off()



#### effect size: back and forth ####
bp.control <- 140
s <- 10
es <- 0.5
bp.treat <- bp.control - s*es
y <- rnorm(1000, mean=rep(c(bp.control, bp.treat), each=500), sd=s)
t.test(y[1:500], y[501:1000]) 
#!! effect size != t-value !!
a <- y[1:500]
b <- y[501:1000]

d.cohen <- function(x1, x2){
  # compute Cohen's d from two samples (https://en.wikipedia.org/wiki/Effect_size#Cohen's_d)
  s <- sqrt(((length(x1)-1)*var(x1) + (length(x2)-1) * var(x2)) / (length(c(x1,x2)-2)))
  (mean(x1) - mean(x2))/s
}
d.cohen(a,b)


# approximately this can be computed from the t-statistic:
#    t = delta/s = delta/ s * sqrt(N) 
#<-> d = delta/s = t / sqrt(N)
# but this yields too low values:
t.test(a, b)$statistic / sqrt(500)
# because means are estimated with error, and t is underestimated

ts <- replicate(1000, {
  N=100
  y <- rnorm(2*N, mean=rep(c(bp.control, bp.treat), each=N), sd=s)
  t.test(y[1:N], y[(N+1):(2*N)])$statistic / sqrt(N/2)
})
summary(ts) # i.e. biased independent of N!

# using N/2 rather than N yields correct estimates (whyever: neither single-tailed nor var.equal makes a difference): this is also the approach of effectsize::t_to_d
library(effectsize)
t_to_d(t.test(a, b)$statistic, t.test(a, b)$parameter, ci=NULL)

# Hedges correction (does not make up the sqrt(2) factor):
1 - (3/(4*(2*N) - 9))



#### How many samples do I need for a regression slope to become significantly different from 0? ####
# situation: we want to go out and sample N locations/subjects for their values of x and y; x will be a value between 0 and 100, and y will be normally distributed but dependent on x; given that we want to be able to detect a slope of value "slope", and there is noise of value "sd", how many samples do we need to have a power of 80% to detect a slope significantly different from 0 at p=0.05?
# similar to t-test, decide on:
# slope, p, power, sd to estimate n
sd=1
slope=1
N <- 20

set.seed(123)
x <- runif(N, 0, 1)
y <- slope * x + rnorm(N, mean=0, sd=sd)

X <- as.matrix(cbind(1, x))
solve(t(X) %*% X) %*% t(X) %*% y # estimates of intercept and slope
# same as:
coef(lm(y ~ x))

# run R simulations for sample size N:
R <- 100
N <- 20
set.seed(123)
beta <- Ps <- NA
for (r in 1:R){
  x <- runif(N, 0, 1)
  y <- slope * x + rnorm(N, mean=0, sd=sd)
  fm <- lm(y ~ x)
  beta[r] <- coef(fm)[2]
  Ps[r]   <- summary(fm)$coefficients[2, 4]
  rm(x, y, fm)
}
quantile(beta, c(0.025, 0.5, 0.975))
sum(Ps < 0.05)/R


# run simulations for R runs of sample size N:
R <- 500
N <- seq(10, 200, by=10)
set.seed(123)
beta.mat <- Ps.mat <- matrix(NA, nrow=max(N), ncol=R)
for (n in N){
  for (r in 1:R){
    x <- runif(n, 0, 1)
    y <- slope * x + rnorm(n, mean=0, sd=sd)
    fm <- lm(y ~ x)
    beta.mat[n, r] <- coef(fm)[2]
    Ps.mat[n, r]   <- summary(fm)$coefficients[2, 4]
    rm(x, y, fm)
  }
}
quants <- apply(beta.mat, 1, quantile, c(0.025, 0.5, 0.975), na.rm=T)
plot(N, quants[2, N], ylim=c(-3, 4), pch=16)
abline(h=1, col="grey")
arrows(N, quants[1, N], N, quants[3, N], angle=90, code=3, length = .02)
# where I expected the 95% CI, but where it isn't:
#curve(1 + qt(0.975, df=x)*1/sqrt(x), add=T, col="cornflowerblue")
#curve(1 - qt(0.975, df=x)*1/sqrt(x), add=T, col="cornflowerblue")


(prop.sign <- apply(Ps.mat < 0.05, 1, sum, na.rm=T)/R  )
plot(N, prop.sign[N], las=1, xlab="sample size", ylab="proportion of significant slopes", type="o", pch=16)
abline(h=0.8, col="grey")

# This analysis can be repeated for different slope or noise levels.


#### simulate detectable effect size for a given sample size  #### 


#### link effect size (d, g, t, Glass' delta, ...) to regression R2, correlation r and standardised regression coefficient ####
# stand. regr. coefs. won't wash, since not applicable to two groups !?


#### Correlation simulation ####
# reimplements Schönbrodt & Perugini 2013

# 1. Simulate a bivariate Gaussian distribution with 1,00,000 cases and a specified correlation q (the ‘‘population’’).
# 2. Draw B = 100,000 bootstrap samples1 with nmax = 1000 cases from this population.
# 3. For each bootstrap sample calculate the correlation for every sample size n, starting from nmin = 20 up to nmax = 1000, with a step size of 1 (‘‘trajectory of the correlation’’).
# 4. Calculate the POS for each bootstrapped trajectory. For that purpose trace back the trajectory from nmax (99.98 % of all trajecto- ries terminated within the COS at n = 1000) until it breaks the COS for the first time. The sample size of this break is recorded as the POS for this trajectory.

phi = 0.3
B = 1e4 # Schönbrodt: 1e5 takes quite a while, 1e4 is already a lot
N = 20:1000
w=0.1
# first we do everything once only:
# ad 1.:
library(mvtnorm)
population <- rmvnorm(1e5, mean=c(0,0), sigma=matrix(c(1, phi, phi, 1), 2, 2))
# ad 2.:
boots <- population[sample(1000, replace=T), ]
# ad 3.: 
cors <- NA
for (i in N) cors[i] <- cor(boots[1:i,1], boots[1:i,2])
#cors <- sapply(boots, \(x) cor(x[,1], x[,2]))
# ad 4.:
(POS <- max(which((cors > (0.3+w) | cors < (0.3-w)))))

plot(N, cors[N], type="l", las=1, xlim=c(20, 200))
abline(h=0.3, col="grey")
abline(h=c(0.3+w, 0.3-w), lty=2)
abline(v=POS)

## now do this repeatedly (B times):
set.seed(123)
boots <- replicate(B, population[sample(1000, replace=T), ]) # array: 1000 x 2 x B

cors.mat <- matrix(NA, nrow=max(N), ncol=B)
for (i in 1:B){
  cors <- NA
  for (j in N) cors[j] <- cor(boots[1:j, 1, i], boots[1:j, 2, i])
  cors.mat[,i] <- cors
  rm(cors)
}
cors.mat <- cors.mat[-c(1:19),] # remove N=1:19
POSvec <- apply(cors.mat, 2, \(x) max(which((x > (0.3+w) | x < (0.3-w)))) ) + 19
POSvec[which(is.infinite(POSvec))] <- 20 # runs always within COS
summary(POSvec)
CIs <- apply(cors.mat, 1, quantile, c(0.025, 0.975))
#save(cors.mat, POSvec, w, CIs, file="schonbrodt.Rdata")
load(file="schonbrodt.Rdata")

pdf(file="Schonbrodtfig.pdf", width=10, height=7)
par(mar=c(4,4,1,1))
plot(N, cors.mat[,1], type="n", lwd=2, las=1, xlim=c(20, 300), ylab="correlation", ylim=c(-.1, .7))
for (i in 1:40) lines(N, cors.mat[,i], type="l", lwd=0.5, col="grey")
lines(N, cors.mat[,11], type="l", lwd=2, col="black")
lines(N, CIs[1,], type="l", lwd=1, lty=3, col="black")
lines(N, CIs[2,], type="l", lwd=1, lty=3, col="black")
abline(h=0.3, col="black")
abline(h=c(0.3+w, 0.3-w), lty=2)
abline(v=POSvec[11], col="grey")
grid()
arrows(220, 0.525, 210, 0.41, length=.1)
text(220, 0.55, "point of stability (POS):\n n=209", adj=0, cex=1)
dev.off()


### variation on this approach ###
## now do this repeatedly (B times):
set.seed(123)
B=1e4
cors.mat2 <- matrix(NA, nrow=max(N), ncol=B)
for (i in N){
  cors.mat2[i,] <- replicate(B, {
    X <- rmvnorm(i, mean=c(0,0), sigma=matrix(c(1, phi, phi, 1), 2, 2))  
    cor(X)[1,2]
  })
}
cors.mat2[1:21, 1:21]
cors.mat2 <- cors.mat2[-c(1:19),] # remove N=1:19
POSvec2 <- apply(cors.mat2, 2, \(x) max(which((x > (0.3+w) | x < (0.3-w)))) ) + 19
summary(POSvec2)
CIs2 <- apply(cors.mat2, 1, quantile, c(0.025, 0.975))
sd(cors.mat2[, 1000])
# save(cors.mat2, POSvec2, w, CIs2, file="schonbrodt2.Rdata")


pdf(file="Schonbrodtfig2.pdf", width=14, height=7)
par(mar=c(4,4,1,1))
plot(N, cors.mat2[,1], type="n", lwd=2, las=1, xlim=c(20, 600), ylab="correlation", ylim=c(-.1, .7))
for (i in 1:40) lines(N, cors.mat2[,i], type="l", lwd=0.5, col="grey")
lines(N, cors.mat2[,11], type="l", lwd=2, col="black")
lines(N, CIs2[1,], type="l", lwd=1, lty=3, col="black")
lines(N, CIs2[2,], type="l", lwd=1, lty=3, col="black")
abline(h=0.3, col="black")
abline(h=c(0.3+w, 0.3-w), lty=2)
abline(v=POSvec2[11], col="grey")
grid()
arrows(590, 0.525, 608, 0.41, length=.1)
text(520, 0.535, "point of stability (POS):\n n=609", adj=0, cex=1)
dev.off()




#### Recruit until significant? ####
# effect of washing hands on patient death
# patients randomly assigned to +/- hand washing (https://www.sciencedirect.com/topics/medicine-and-dentistry/sequential-analysis)
# alternative case study: Fisher's tea milk testing (https://en.wikipedia.org/wiki/Fisher's_exact_test)

set.seed(123)
pC <- 0.18 # probability of death in controls
pT <- 0.01 # probability of death in treatment (washing hands)
control <- sample(c(0,1), size=1000, prob=c(1-pC, pC), replace=T)
treat <- sample(c(0,1), size=1000, prob=c(1-pT, pT), replace=T)

Co <- table(control[1:50])
Tr <- table(treat[1:50])

chisq.test(cbind(Co, Tr))

# run one sequential simulation and look at rates when it becomes significant:
pC <- 0.18 # probability of death in controls
pT <- 0.01 # probability of death in treatment (washing hands)

set.seed(1234)
Ns <- seq(10, 1000, by=50)
simu1res <- matrix(NA, nrow=length(Ns), ncol=4)
colnames(simu1res) <- c("N", "p-value", "control rate", "treatment rate")
for (i in seq_along(Ns)){
  control <- sample(c(0,1), size=Ns[i], prob=c(1-pC, pC), replace=T)
  treat <- sample(c(0,1), size=Ns[i], prob=c(1-pT, pT), replace=T)
  Co <- table(control)
  Tr <- table(treat)
  simu1res[i,] <- c("N"=Ns[i], "p-value"=suppressWarnings(chisq.test(cbind(Co, Tr))$p.value), "control rate"=Co[2]/sum(Co), "treatment rate"=Tr[2]/sum(Tr))
}
simu1res[is.na(simu1res)] <- 0
# now identify N when first time p<0.05:
firstSignif <- min(which(simu1res[, 2] < 0.05))
simu1res[firstSignif, ]
# significant at around N=110
# mean control rate estimate = 0.16 (true: 0.18)
# mean treat rate estimate = 0.0018 (true: 0.01)

# plot estimates over time: note that never are they biased, unless we use the stopping rule!
colnames(simu1res)[c(3,4)] <- c("controlRate", "treatmentRate")
library(ggplot2)
pdf(file="../figures/stoppingRuleSingleRunEstimate.pdf", width=5, height=5)
ggplot(as.data.frame(simu1res), aes(x=N, y=controlRate)) + 
  geom_smooth(color="white") +
  geom_point() +
  geom_smooth(aes(y=treatmentRate), color="white") +
  geom_point(aes(y=treatmentRate), color="grey") +
  geom_segment(aes(x=10, y=0.01, xend=970, yend=0.01), color="grey60") +   
  geom_segment(aes(x=10, y=0.18, xend=970, yend=0.18), color="black") +   
  xlab(label="sample size") + 
  ylab(label="infection rate estimate") + 
  theme_minimal()
dev.off()

# now do this many times:
set.seed(1234)
Ns <- seq(20, 1000, by=20)
#colnames(simu1res) <- c("N", "p-value", "control rate", "treatment rate")
simus <- replicate(100, {
  simu1res <- matrix(NA, nrow=length(Ns), ncol=4)
  for (i in seq_along(Ns)){
    control <- sample(c(0,1), size=Ns[i], prob=c(1-pC, pC), replace=T)
    treat <- sample(c(0,1), size=Ns[i], prob=c(1-pT, pT), replace=T)
    Co <- table(control)
    Tr <- table(treat)
    simu1res[i,] <- c("N"=Ns[i], "p-value"=suppressWarnings(chisq.test(cbind(Co, Tr))$p.value), "control rate"=Co[2]/sum(Co), "treatment rate"=Tr[2]/sum(Tr))
  }
  simu1res[is.na(simu1res)] <- 0
  # now identify N when first time p<0.05:
  firstSignif <- min(which(simu1res[, 2] < 0.05))
  c(simu1res[firstSignif, ], simu1res[nrow(simu1res),])
})  
rowMeans(simus)
# significant at around N=37
# mean control rate estimate = 0.157 (true: 0.18): unbiased!
# mean treat rate estimate = 0.0032 (true: 0.01)
# full length of trial:
# mean control rate estimate = 0.182 (true: 0.18)
# mean treat rate estimate = 0.010 (true: 0.01)
library(ggplot2)
library(tidyr)
dats <- pivot_longer(as.data.frame(t(simus)), cols = c(3,4,7,8), names_to="rates", values_to="estimate") 
dats$rates <- factor(dats$rates, levels = c("V3", "V7", "V4", "V8"), labels=c("C stop" ,"C full", "T stop", "T full"), ordered = TRUE)
#head(t(simus))
# give names to variables for plotting!
#summary(dats)
pdf(file="../figures/stoppingRule.pdf", width=5, height=5)
ggplot(dats, aes(y=estimate, x=rates)) + 
  scale_y_continuous(trans='sqrt') +
  geom_jitter(size=0.7) +
  geom_segment(aes(x = 2.5, y = pT, xend = 4.5, yend = pT), color="grey", linewidth=1) +
  geom_segment(aes(x = .5, y = pC, xend = 2.5, yend = pC), color="grey", linewidth=1) +
  geom_boxplot(outlier.size = -1, fill="transparent") + 
  theme_minimal()
dev.off()



## now do this for varying values of pT
# show that bias depends on prevalence: the closer to 0.5, the less bias
pC <- 0.18
pT <- seq(0.01, 0.18, by=0.01)
set.seed(123)
Ns <- seq(10, 1000, by=20)

manysimus <- matrix(NA, nrow=length(pT), ncol=8)
for (j in seq_along(pT)){
  simus <- replicate(100, {
    simu1res <- matrix(NA, nrow=length(Ns), ncol=4)
    for (i in seq_along(Ns)){
      control <- sample(c(0,1), size=Ns[i], prob=c(1-pC, pC), replace=T)
      treat <- sample(c(0,1), size=Ns[i], prob=c(1-pT[j], pT[j]), replace=T)
      Co <- table(control)
      Tr <- table(treat)
      simu1res[i,] <- c("N"=Ns[i], "p-value"=suppressWarnings(chisq.test(cbind(Co, Tr))$p.value), "control rate"=Co[2]/sum(Co), "treatment rate"=Tr[2]/sum(Tr))
    }
    simu1res[is.na(simu1res)] <- 0
    # now identify N when first time p<0.05:
    firstSignif <- min(which(simu1res[, 2] < 0.05))
    c(simu1res[firstSignif, ], simu1res[nrow(simu1res), ])
  })  
  manysimus[j,] <- rowMeans(simus, na.rm=T)
}
manysimus
colnames(manysimus) <- c("Nstop", "pstop", "Cstop", "Tstop", "Nfull", "pfull", "Cfull", "Tfull")
manysimus <- cbind.data.frame("TtrueRate"=pT, manysimus)
plot(manysimus[, c(4, 8)])
abline(0,1)
# bias ~ prevalence:
plot(pT, manysimus[,4] - pT)
abline(h=0)

out <- cbind(manysimus[, c(1, 2, 4, 5)], "bias stopping"=manysimus[,5]-manysimus[,1], manysimus[, c(8,9)], "bias full"=manysimus[,9]-manysimus[,1])
xtable::xtable(out, digits=c(0, 2, 0, 3, 3, 3, 3, 3, 4), type = "latex") # possibly use print(..., file="...") to export to file


pdf(file="../figures/mortalityEstimatesStopping.pdf", width=5, height=5)
ggplot(as.data.frame(manysimus), aes(x=TtrueRate, y=Tstop-pT)) + 
  geom_hline(yintercept=0, color="grey") + 
  geom_point() +
  geom_point(aes(x=TtrueRate, y=Tfull-pT), color="grey") +
  xlab("true mortality rate with hand washing") + 
  ylab("bias in rate estimate for stopping rule data") + 
  annotate("text", x=c(0.01, 0.01), y=c(-0.021, -0.002), label=c("stopping rule", "full data set"), size=5, color=c("black", "grey"), adj=0) +
  theme_minimal()
dev.off()






#### time series: how long? ####
# specify question precisely!
# After how many months, starting in month 100 (because of missing data before), will the increasing trend in CO2 be detectable?
# problems: * also non-linear trend? * will correcting for seasonality improve (linear) trend detection? * Is switch between methods allowed (i.e. start with lm, then move to gam)?

lm(co2 ~ months, data=data.frame("months"=1:length(co2), "co2"=co2)) #: 0.11/month!

# First, we observe a strong seasonal pattern that we need to estimate in order to simulate it.
dats <- data.frame("months"=1:24, "month"=rep(1:12, times=2), CO2=co2[1:24])
fnls <- nls(CO2 ~ a*sin(2*pi/12*month+b) + c, start=list(a=1, b=0, c=316), data=dats)
summary(fnls)
plot(CO2 ~ months, data=dats, type="l")
curve(2.29*sin(2*pi/12*x -0.77) + 316, add=T)
#abline(v=c(10, 22))

library(ggplot2)
pdf(file="../figures/co2_2yearsWithSin.pdf", width=5, height=5)
ggplot(dats, aes(x=months, y=CO2)) + geom_point() + 
  geom_line() +
  stat_function(fun=\(x) 2.29*sin(2*pi/12*x -0.77) + 316, color="grey") + 
  ylab(label="atmospheric CO2 concentration [ppm]") +
  xlab(label="months since January 1959") +
  theme_minimal()
dev.off()

noise <- sd(residuals(fnls)) * 1.5 # it is unlikely that the data contain the possible underlying noise, therefore multiply by 50% safety margin ...
noise
acf(residuals(fnls), plot=F)
acf(residuals(fnls))
rho <- acf(residuals(fnls), plot=F)$acf[2]

# simulate 8 years without trend:
ms <- 1:96 # months since start
m <- rep(1:12, len=length(ms)) # month of the year
set.seed(123)
noises <- rnorm(1, 0, noise)
for (i in 2:length(ms)) noises[i] <- rho*noises[i-1] + rnorm(1, 0, noise)
acf(noises, plot=F)[0:8] # fine
co2hat <- 2.29*sin(2*pi/12*ms -0.77) + 316 + noises

plot(co2hat, type="l")

# now add linear trend of varying degree (e.g. between 0 and 0.2°C/month)

noises <- rnorm(1, 0, noise)
for (i in 2:length(ms)) noises[i] <- rho*noises[i-1] + rnorm(1, mean=0.01*ms[i], noise)
co2hat <- 2.29*sin(2*pi/12*ms -0.77) + 316 + noises
#plot(co2hat, type="l")

dats10 <- rbind.data.frame(dats, cbind("months"=25:(24+length(ms)), "month"=m, "CO2"=co2hat))

library(mgcv)
fgam0 <- gam(CO2 ~ months, data=dats10)
summary(fgam0)
# note that accounting for seasonality increases sensitivity:
fgam <- gam(CO2 ~ s(month, k=12) + months, data=dats10)
summary(fgam)
acf(residuals(fgam), plot=T) # not shown

fgamm <- gamm(CO2 ~ s(month, k=12) + months, correlation=corAR1(form = ~ months), data=dats10)
summary(fgamm$gam)
summary(fgamm$lme)
#plot simulated data and fitted GAM
# add bar to symbolise the residual SE in bottom right corner
# also: make x and y labels same as in previous
predgamm <- predict(fgamm$gam, se.fit=T)
library(ggplot2)
pdf(file="../figures/co2_10yearsWithGAMM.pdf", width=10, height=5)
ggplot(dats10, aes(x=months, y=CO2)) + 
  geom_ribbon(aes(ymin = predgamm$fit - 2*predgamm$se.fit, 
                    ymax = predgamm$fit + 2*predgamm$se.fit), fill="gray80") +
  geom_line(color="grey40") + 
  geom_point(color="grey40") +
  geom_line(data=dats[1:24,]) + 
  geom_point(data=dats[1:24,]) + 
  ylab(label="atmospheric CO2 concentration [ppm]") +
  xlab(label="months since January 1959") +
  theme_minimal()
dev.off()

# now real simu:
N <- 1 # number of additional years of monitoring
trend <- 0.01
months  <- 1:(12*(2+N)) # months
moty <- rep(1:12, len=length(months)) # month of the year
set.seed(123)
CO2 <- round(2.29*sin(2*pi/12*months -0.77) + 316 + rnorm(length(months), 0.01*months, noise), 2)
dats <- data.frame(months, moty, CO2)
fgamm <- gamm(CO2 ~ s(moty, k=12) + months, correlation=corAR1(form = ~ months), data=dats)
summary(fgamm$gam)$p.table[2,4]

# now many times:
set.seed(123)
simu.one <- replicate(100, {
  CO2 <- round(2.29*sin(2*pi/12*months -0.77) + 316 + rnorm(length(months), 0.01*months, noise), 2)
  dats <- data.frame(months, moty, CO2)
  fgamm <- gamm(CO2 ~ s(moty, k=12) + months, correlation=corAR1(form = ~ months), data=dats)
  summary(fgamm$gam)$p.table[2,4]
  })
#simu.one
mean(simu.one < 0.05)


# now for many N values:
set.seed(123)
Ns <- 1:10 # number of years
simu.res01 <- NA
for (i in seq_along(Ns)){
  simu.one <- replicate(100, {
    months  <- 1:(12*(2+Ns[i])) # months
    moty <- rep(1:12, len=length(months)) # month of the year
    CO2 <- round(2.29*sin(2*pi/12*months -0.77) + 316 + rnorm(length(months), 0.01*months, noise), 2)
    dats <- data.frame(months, moty, CO2)
    fgamm <- gamm(CO2 ~ s(moty, k=12) + months, correlation=corAR1(form = ~ months), data=dats)
    summary(fgamm$gam)$p.table[2,4]
  })
  simu.res01[i] <- mean(simu.one)
}
round(simu.res01, 2)
plot(Ns, simu.res01, type="o")
abline(h=0.05)


# at this point is is worth digressing on how to parallalise a for loop!!

# now for many trends and Ns values:
set.seed(123)
trends <- seq(0, 0.2, len=21) # possible trends
Ns <- 1:10 # number of years
simu.mat <- matrix(NA, nrow=length(Ns), ncol=length(trends))
for (j in seq_along(trends)){ # this takes quite a few minutes!
  for (i in seq_along(Ns)){
    simu.one <- replicate(1000, {
      months  <- 1:(12*(2+Ns[i])) # months
      moty <- rep(1:12, len=length(months)) # month of the year
      CO2 <- round(2.29*sin(2*pi/12*months -0.77) + 316 + rnorm(length(months), trends[j]*months, noise), 2)
      dats <- data.frame(months, moty, CO2)
      fgamm <- gamm(CO2 ~ s(moty, k=12) + months, correlation=corAR1(form = ~ months), data=dats)
      summary(fgamm$gam)$p.table[2,4]
    })
    simu.mat[i,j] <- mean(simu.one) # across the 1000 simulations per trend and year
  }
}
#save(simu.mat, file="CO2simuMat.Rdata")
load(file="CO2simuMat.Rdata")
#plot(Ns, simu.mat[,5], type="o")

# convert: proportion of significant (i.e. 1-simu); stack for ggplot
library(tidyr)
Ns <- 1:10 # number of years
dats <- pivot_longer(as.data.frame(simu.mat), cols = 1:21, names_to="trend", values_to="signif") 
pdf(file="../figures/co2_trendN.pdf", width=5, height=5)
ggplot(cbind(Ns=rep(Ns, each=21), dats), aes(x=Ns, y=signif, color=trend)) + 
  expand_limits(y=.5) +
  annotate('rect', xmin=1, xmax=2, ymin=0, ymax=.52, alpha=.2, fill='grey') +
  annotate('text', label=c(0, 0.01, 0.02, 0.03, 0.04, 0.05), x=rep(1.8, 6), y=1-c(0.532, 0.58, 0.72, 0.84, .92, .96), size=2) +
  geom_line(show.legend = FALSE) +
  scale_colour_grey(start=0, end=.75, na.value = "green") +
  geom_hline(yintercept=0.025, color="grey") + 
  scale_x_continuous(trans='sqrt', breaks=c(2,3,4,5,6,7,8,9,10)) +
  xlab(label="number of years monitored") +
  ylab(label="mean p-value for trend") +
  theme_minimal()
dev.off()



#### prior-posterior-W-boson ####
library(ggplot2)
x <- seq(80.32, 80.48, len=3)
y <- dnorm(x, mean=80, sd=.01)
# posterior: normalised product of the two:
numerator <- function(x) dnorm(x, mean=80.379, sd=0.012) * dnorm(x, mean=80.433, sd=0.009)
denominator <- integrate(numerator, 80, 81) # 0.040
post <- function(x) numerator(x)/denominator$value

dats <- data.frame("mass"=x, "density"=y)
pdf(file="../figures/Wboson.pdf", width=5, height=5)
ggplot(data=dats, aes(x=mass, y=density)) + 
  geom_vline(xintercept=80.36, color="grey", linewidth=1.2) +
  stat_function(fun=\(x) dnorm(x, mean=80.379, sd=0.012), color="grey40", n=1001, linewidth=1.2) +   
  stat_function(fun=\(x) dnorm(x, mean=80.433, sd=0.009), color="grey60", n=1001, linewidth=1.2) +   
  stat_function(fun=\(x) post(x), color="black", n=1001, linewidth=1.2) +   
  stat_function(fun=\(x) dnorm(x, mean=80.360, sd=0.016), color="grey60", n=1001, linewidth=0.5, linetype=2) +   
  annotate('text', label="2021 prior", x=80.38, y=35, size=4) +
  annotate('text', label="2022 TEVATRON", x=80.46, y=22, size=4) +
  annotate('text', label="2022 posterior", x=80.44, y=52, size=4) +
  annotate('text', label="2023 ATLAS", x=80.33, y=12, size=2, color="grey60") +
  theme_minimal()
dev.off()
