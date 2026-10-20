setwd("/Users/carsten/Library/CloudStorage/Dropbox/FooledByData_SEE_dropbox/Rcode")

sapply(c("raster", "sp", "ggplot2", "faraway", "rgbif", "maps"), library, character.only = TRUE)

#### Normality tests ####
# simulate ceramic shard thickness following Daszkiewicz et al. 2019
# range of thickness: 2.6 - 12.7 mm
# cultures: Early Pre-Roman Iron Age: -700 : 1000; prerom: -500 : -1; early emperor: 1 : 375; most in 
# see: https://de.wikipedia.org/wiki/Vorlage:Bronzezeit%2BEisenzeit
N <- 113
set.seed(123)
years <- runif(N, -300, 400) # Fig. 2
# thickness: 8mm --> 5 mm
#curve(8 - 3*300/700 - 3/700*x, -300, 400)
thickness <-  6.71  - 0.0043 * years + rnorm(N, 0, 2)
library(ggplot2)
dats <- data.frame("thickness"=thickness, "year"=years)
ggplot(dats, aes(x=years, y=thickness)) + geom_point() + theme_minimal()

N <- 113
set.seed(123)
years <- runif(N, -300, 400) # Fig. 2
# thickness: 8mm --> 5 mm
#curve(8 - 3*300/700 - 3/700*x, -300, 400)
thickness <-  6.71  - 0.0043 * years + runif(N, -3, 3)

fm <- lm(thickness ~ years, data=dats)
shapiro.test(residuals(fm))

rm(fm, thickness, years)
set.seed(123)
N=113
Ps <- replicate(1000, {
  years <- runif(N, -300, 400)
  thickness <-  6.71  - 0.0043*years + rnorm(N, 0, 2)
  dats <- data.frame("thickness"=thickness, "year"=years)
  fm <- lm(thickness ~ years, data=dats)
  shapiro.test(residuals(fm))$p.value
})
sum(Ps < 0.05)

set.seed(123)
N=1000
Ps <- replicate(1000, {
  years <- runif(N, -300, 400)
  thickness <-  6.71  - 0.0043*years + runif(N, -3, 3)
  dats <- data.frame("thickness"=thickness, "year"=years)
  fm <- lm(thickness ~ years, data=dats)
  shapiro.test(residuals(fm))$p.value
})
sum(Ps < 0.05)


#### Heteroscedasticity  #### 
library(faraway)
pdf(file="../figures/butterfat.pdf", width=5, height=5)
ggplot(butterfat, aes(y=Butterfat, x=Breed)) + geom_boxplot(fill="grey") + labs(y="Butterfat [% mass]") + theme_minimal()
dev.off()
fm <- aov(Butterfat ~ Breed, data=butterfat)
summary(fm)
fligner.test(Butterfat ~ Breed, data=butterfat) # problem!
car::leveneTest(fm) # problem!

TukeyHSD(fm)
with(butterfat, tapply(Butterfat, Breed, var))
# Ad 1.:
# Assume same (average) variance for all:
mean(with(butterfat, tapply(Butterfat, Breed, sd))) # sd=0.394

set.seed(123)
butterfat$Sim <- rnorm(100, mean=rep(with(butterfat, tapply(Butterfat, Breed, mean)), each=20), sd=0.394)
fmSim <- aov(Sim ~ Breed, data=butterfat)
#summary(fmSim)[[1]]$"Pr(>F)"[1]
summary(fmSim)[[1]]$F[1]

# ANOVA not assuming equal variances:
oneway.test(Butterfat ~ Breed, data=butterfat)

# many times:
set.seed(234)
fmList <- list()
for (i in 1:1000){
  butterfat$Sim <- rnorm(100, mean=rep(with(butterfat, tapply(Butterfat, Breed, mean)), each=20), sd=0.394)
  fmList[[i]] <- aov(Sim ~ Breed, data=butterfat)
}
Ps <- sapply(fmList, \(x) summary(x)[[1]]$"Pr(>F)"[1])
Fs <- sapply(fmList, \(x) summary(x)[[1]]$F[1])
mean(Fs)
pdf(file="../figures/heteroscedasticityANOVA.pdf", width=5, height=5)
ggplot(data.frame("x"=Fs), aes(x=x)) + geom_density(colour="grey", size=1.5) + labs(x=expression(paste(italic(F),"-value of ANOVA"))) + geom_vline(xintercept=summary(fm)[[1]]$F[1], linewidth=1.5) + geom_vline(xintercept=qf(0.05, 4, 95), linetype=2) + theme_minimal() 
dev.off()

pf(57.3, 4, 95, lower.tail=F)
pf(49.8, 4, 95, lower.tail=F)


# Does that affect a specific pairwise comparison, in particular one close to significance (Canadian-Ayrshire)?
TukeySim <- sapply(fmList, \(x) TukeyHSD(x)$Breed[1, 4])
ggplot(data.frame("x"=TukeySim), aes(x=x)) + geom_density(colour="grey", size=1.5) + labs(x=expression(paste(italic(P),"-value of Tukey's HSD Canadian vs. Ayreshire"))) + geom_vline(xintercept=0.038, size=1.5) + geom_vline(xintercept=0.05, linetype=2) + theme_minimal() 
summary(TukeySim)
sum(TukeySim >0.05) / length(TukeySim)

#### Multiple Testing, False Discovery Rate and Family-Level error rates ####
# builds on Quantitude Episode S4E16: https://quantitudepod.org/s4e16-false-discovery/

# For $k=10$, this value is $p=1-(1-0.05)^10 = 0.41$.

tenner <- sample(size=10, c("different","same"), prob=c(0.05, 0.95), replace=T)
any(tenner == "different")

# Let us now repeat this many times, to see how often we discover a (non-existent) difference:

set.seed(123)
simu10er <- replicate(10000, {
  tenner <- sample(size=10, c("different","same"), prob=c(0.05, 0.95), replace=T)
  any(tenner == "different")
})
sum(simu10er)/length(simu10er)

#For 10 comparison the critical $p$-value would thus be $1- 0.95^{1/10} = 0.0051$.

set.seed(1234)
FWER <- NA
j=1
for (i in 10^seq(-4, -1, by=0.1)) {# log-scaled sequence
  simu10er <- replicate(1000, {
    tenner <- sample(size=10, c("different","same"), prob=c(i, 1-i), replace=T)
    any(tenner == "different")
  })
  FWER[j] <- sum(simu10er)/length(simu10er)
  j <- j + 1 # count up for each loop
}
plot(10^seq(-4, -1, by=0.1), FWER, type="o", pch=16, log="x")
abline(h=0.05)
abline(v=0.0051)

library(ggplot2)
library(latex2exp)
dats <- cbind.data.frame(FWER, "p"=10^seq(-4, -1, by=0.1))
#pdf(file="../figures/FWERplot.pdf", 6, 6)
pdf(file="FWERplot.pdf", 6, 6)
ggplot(dats, aes(y=FWER, x=p)) + 
  scale_x_continuous(trans='log10', breaks=as.vector(c(1,2,3,4,5,6,7,8,9) %o% 10^c(-c(4:2)))) +
  ylim(0, 0.38) + 
  xlab(TeX(r'(Critical value $\alpha\prime$ for single test)'))+
  ylab("Family-wise Error Rate") +
  geom_line() + 
  geom_vline(xintercept=0.0051, color="grey") + 
  geom_hline(yintercept=0.05, color="grey") + 
  theme_minimal() + 
  theme(axis.text.x  = element_text(angle=90))
dev.off()



# breaks=c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1), 
#          labels=c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)
#Now, an important assumption here is that the different tests are uncorrelated. But if we take 10 different measures on the same subject and analyse each for its relation to the predictor, then often these 10 measures are correlated: taller people have larger hands, are heavier and wear darker clothes (say). So then we cannot use the above computation: if one test is significant, it is more likely that the others are, too (and similarly for not significant). 
#
#We still can solve this problem analytically: what is the overall probability of a false positive for $k$ correlated tests?
#
#
#[go through the statsexchange solution here!]
#
#But how do we simulate this? Actually, the simulations are a bit involved, and we borrow them from the R-package bindata. The key idea is to simulate multivariate normal data that are then converted into 0s and 1s. To do so, we first have to understand how the correlation of the binary is related to that of the normal variates.
# \href{https://stats.stackexchange.com/questions/284996/generating-correlated-binomial-random-variables#285008}{solution posted at stackexchange}.
#
#The first thing to realise is that correlation is about the joint distribution of two random variables. If they were uncorrelated, $P(X, Y) = 0$. So for correlated random Bernoulli variables, $X$ and $Y$, the probability of observing a specific combination is non-zero. Let $Pr(X=1) = p$, $Pr(Y=1)=q$ and their correlation $=\rho$.
#
#Remember that the variance of a Bernoulli variable is $var(X) = p(1-p)$. The equation for Pearson's correlation scales the difference between joint and marginal expectation by their variances:
# \begin{equation}
# \rho_{X, Y} = \frac{E(X,Y) - E(X)E(Y)}{\sqrt{(E(X^2) - E(X)^2)(E(Y^2) - E(Y)^2)}}
# \end{equation}
# We can now insert our Bernoulli-notation, with $a$ being the joint probability:
# \begin{equation}
# \rho_{X, Y} = \frac{a - pq}{\sqrt{(p(1-p)q(1-q)}}
# \end{equation}  
# So if we want to specify the correlation, $\rho$, then we need to first compute the joint probability $Pr((X,Y)=(0,0))=a$ from $\rho, p$ and $q$. Rearranging the last equation yields:
# \begin{equation}
# a = (1-p)(1-q) + \rho_{X, Y} \sqrt{(p(1-p)q(1-q)}
# \end{equation}  
# From this we can now compute the probability for the Bernoulli-realisation of each of the four possible combinations of $X$ and $Y$:
#   \begin{array}
#   Pr((X,Y)=(0,0)) = a =  (1-p)(1-q) + \rho_{X, Y} \sqrt{(p(1-p)q(1-q)}
#   Pr((X,Y)=(1,0)) = ?? = 1 - q - a 
#   Pr((X,Y)=(0,1)) = ?? = 1 - p - a 
#   Pr((X,Y)=(1,1)) = ?? = p + q + a - 1
#   \end{array}
# 
# 
# %The probability of observing $(X,Y) = (0,0)$ would then be $Pr((X,Y) = (0,0)) = (1-p)(1-q) - \rho \sqrt{pq(1-p)(1-q)}$, i.e. the individual probability for $X=0$ and $Y=0$ minus the joint probability $Pr( (X,Y) = (0,0) ) = \rho \sqrt{pq(1-p)(1-q)}$. 
# %Let us start with observing the event $X=0$ and $Y=0$: $Pr( (X,Y) = (0,0) ) = a$. 
# %If $Pr(X=1) = p$ and $Pr(Y=1)=q$, then $Pr( (X,Y) = (1,0) ) = p(1-q)-a = 1- q- a$.
# 
# # bincorr2commonprob converts the matrix of pairwise correlations into joint probabilities:
# bincor <- matrix(0.9, 5, 5)
# diag(bincor) <- 1
# bindata::bincorr2commonprob(margprob=rep(0.9, 5), bincorr=bincor)  
# rmvbin(n=10, margprob=rep(0.5, 5), bincorr=bincor)  


# If we do multiple comparisons (e.g. because we have many factor levels in a post-hoc comparison; or because we have measured many characteristics of our target), the probability of false positives builds up:
#### FDR a) simu p of false pos  ####
# 10 comparisons move p from 0.05 to 0.4
1 - (1-0.05)^10 # = 0.40; for independent tests
# now simulate this:
tenner <- sample(size=10, c(T,F), prob=c(0.05, 0.95), replace=T)
any(tenner == T)

simu10er <- replicate(1000, {
  tenner <- sample(size=10, c(T,F), prob=c(0.05, 0.95), replace=T)
  any(tenner == T)
})
sum(simu10er)/length(simu10er)

# now include correlation between samples:
corr <- 0.8 # has to be larger than 0!
(first <- sample(size=1, c(T,F), prob=c(0.05, 0.95), replace=T))
rest <- if (first){
  sample(size=9, c(T, F), prob=c(corr, 1-corr), replace=T)
} else {
  sample(size=9, c(T, F), prob=c(1-corr, corr), replace=T)
}
all <- c(first, rest)
any(all == T)

# now replicate

# alternative: use the approach of whuber on stackexchange:
# https://stats.stackexchange.com/questions/284996/generating-correlated-binomial-random-variables#285008
# (see also: https://stats.stackexchange.com/questions/324158/how-to-generate-correlated-bernoulli-variables)



#### FDR b) mixing differences and non-differences ####
# show that distributions of p-values differ for nulls and working hypothesis being true:
## (1) t-test for null vs t-test for different
set.seed(1234)
N <- 20
dats0 <- data.frame("y"=rnorm(2*N), "treat"=rep(LETTERS[1:2], each=N))
t.test(y ~ treat, data=dats0)$p.value

dats1 <- data.frame("y"=rnorm(2*N, mean=rep(c(1, 1.2), each=N)), "treat"=rep(LETTERS[1:2], each=N))
t.test(y ~ treat, data=dats1)$p.value


p0s <- replicate(1000, {
  dats0 <- data.frame("y"=rnorm(2*N), "treat"=rep(LETTERS[1:2], each=N))
  t.test(y ~ treat, data=dats0)$p.value
})
hist(p0s) # uniform

p1s <- replicate(1000, {
  dats1 <- data.frame("y"=rnorm(2*N, mean=rep(c(1, 1.3), each=N)), "treat"=rep(LETTERS[1:2], each=N))
  t.test(y ~ treat, data=dats1)$p.value
})
hist(p1s) # skewed

# mix true nulls and false nulls:
set.seed(12)
prop0true <- 0.7 # in 70% of the cases the null is true, i.e. no difference
pmix <- c(
  replicate(prop0true*10000, {
    datsmix <- data.frame("y"=rnorm(2*N), "treat"=rep(LETTERS[1:2], each=N))
    t.test(y ~ treat, data=datsmix)$p.value
  }),
  replicate((1-prop0true)*10000, {
    datsmix <- data.frame("y"=rnorm(2*N, mean=rep(c(1, 1.2), each=N)), "treat"=rep(LETTERS[1:2], each=N))
    t.test(y ~ treat, data=datsmix)$p.value
  })
)
hist(pmix)

pdf("../figures/mixingPvalues.pdf", width=7, height=7)
par(mar=c(4,4,1,1))
hist(pmix, freq=T, las=1, main="", xlab="p-value", col="grey90", border="white", ylab="number (out of 10000 replicates)", breaks=seq(0, 1, by=0.05))
hist(pmix[1:7000], freq=T, las=1, main="", xlab="p-value", add=T, col="grey50", border="white", breaks=seq(0, 1, by=0.05))
hist(pmix[7001:10000], freq=T, las=1, main="", xlab="p-value", add=T, col="grey30", border="white", breaks=seq(0, 1, by=0.05))
legend("topright", col=c("grey90", "grey50", "grey30"), pch=15, bty="n", legend=c("total", "H0 true", "H1 true"), cex=1.5)
dev.off()


#lines(density(pmix[1:700], from=0, to=1)$x, 0.7*density(pmix[1:700], from=0, to=1)$y, col="steelblue1", lwd=2)
#lines(density(pmix[701:1000], from=0, to=1)$x, 0.3*density(pmix[701:1000], from=0, to=1)$y, col="steelblue3", lwd=2)
#lines(density(pmix, from=0, to=1), col="black", lwd=2)

library(ggplot2)
df <- data.frame(group=c(rep("H0 true", prop0true*1000), rep("H1 true", (1-prop0true)*1000)), value=pmix)
ggplot(data=df, aes(x=value, fill=group)) + 
  geom_histogram() +
  theme_minimal()
# not correct: counts don't add to true histogram! see "base" graphs above

##### other stuff #####

obs <- c(1, 2, 3, 4)
obs <- c(1, 1, 1, 1)
P <- c(.25, .25, .25, .25)
P <- c(.5, .3, .1, .1)
P <- c(.1, .2, .3, .4)
dmultinom(obs, size=10, prob=P, log=T)

fun <- function(theta){
  P <- seq(0.1, 2, len=4)^theta * c(.1, .2, .3, .4)
  -dmultinom(obs, size=4, prob=P, log=T)
}
optimize(f=fun, interval=c(-5, 5))
fun(0.00107)



# case study: Rastelli et al. 2021: 22 indices (in 6 groups) of happiness used in clinical study; employ FDR (in some, undescribed way); no open data
# no control: pre-post comparison (and 1 vs 2 follow up)
N <- 29 # plus 9 dropouts
#! Simulate to explain the idea of mixing null and real effect; with mixture estimated from p-values.
#So: 22 p-values from different mixtures, and then proportion significant vs FDR-corrected proportion significant
# https://www.graphpad.com/support/faq/using-simulations-to-compute-the-false-discovery-rate/
prop.null <- 0.9
sim <- matrix(rnorm(2*N), ncol=2)
t.test(sim)$p.value

simus <- replicate(1000, {
  sim <- matrix(rnorm(2*N), ncol=2)
  # for proportion of non-nulls add some clear true effect:
  if (sample(1000, 1) < (1-prop.null) * 1000) sim[,2] <- sim[,2] + 2
  t.test(sim)$p.value
})
sum(simus<0.05)
(sum(simus<0.05) - ((1-prop.null) * 1000) ) / 1000 # (FP / N)
# roughly 0.05, i.e. for a single test we have nominal FDR

simus.multiple <- replicate(1000, {
  sim <- matrix(rnorm(2*N), ncol=2)
  # for proportion of non-nulls add some clear true effect:
  if (sample(1000, 1) < (1-prop.null) * 1000) sim[,2] <- sim[,2] + 2
  t.test(sim)$p.value
})


# now multiple tests:
N <- 29
k <- 22
onesim <- replicate(k, {
  sim <- matrix(rnorm(2*N), ncol=2)
  t.test(sim)$p.value
})
# count if ANY p<0.05 (family-wise significance)
any (onesim < 0.05)

# repeat this many times:
set.seed(321)
simus22 <- replicate(1000, {
  onesim <- replicate(k, {
    sim <- matrix(rnorm(2*N), ncol=2)
    t.test(sim)$p.value
  })
  any (onesim < 0.05)
})
mean(simus22) # 0.678 instead of 0.05: same as in analytical solution

# compute FDR = false positives / total positives; since all positives are false: FDR = sum(simu22) / sum(simus22) = 1


# now add proportion of null tests to this code:
prop.null <- 0.1
set.seed(321)
simus22 <- replicate(1000, {
  effect <- FALSE
  if (sample(1000, 1) < (1-prop.null) * 1000) effect <- TRUE # add true effect?
  if (effect) { # real effect
      onesim <- replicate(k, {
        sim <- matrix(rnorm(2*N), ncol=2)
        sim[,2] <- sim[,2] + 2
        t.test(sim)$p.value 
      })
  } else{ # no effect
      onesim <- replicate(k, {
        sim <- matrix(rnorm(2*N), ncol=2)
        t.test(sim)$p.value 
      })
  }
  any (onesim < 0.05)
})
mean(simus22) # 0.851

# FDR:
(sum(simus22) - (1- prop.null) * 1000) / sum(simus22) # 0.41

# now let's control the FDR using the Benjamini-Hochberg procedure:
onesim <- replicate(k, {
  sim <- matrix(rnorm(2*N), ncol=2)
  t.test(sim)$p.value
})
# 1. sort in ascending order
sort(onesim)
plot(sort(onesim))
abline(0, 0.05/22)
# 2. find largest k such that P_k < k / m * alpha (m=22, alpha=0.05)
ks <- sort(onesim) < (1:22 / 22 * 0.05)
signifs <- max(0, which(ks == TRUE)) # largest k
signifs > 0

# repeat many times:
k <- 22
set.seed(321)
simusFDR <- replicate(1000, {
  onesim <- replicate(k, {
    sim <- matrix(rnorm(2*N), ncol=2)
    t.test(sim)$p.value
  })
  ks <- sort(onesim) < (1:k / k * 0.05)
  signifs <- max(0, which(ks == TRUE)) # largest of ks
  signifs > 0
})
mean(simusFDR) # 0.041: nice!

# finally: add some real effects:
prop.true <- 0.5
k = 22
set.seed(321)
simusFDRtrue <- replicate(1000, {
  effect <- FALSE
  if (runif(1) < prop.true) effect <- TRUE # add true effect?
  if (effect) { # real effect
    onesim <- replicate(k, {
      sim <- matrix(rnorm(2*N), ncol=2)
      sim[,2] <- sim[,2] + 2
      t.test(sim)$p.value 
    })
  } else{ # no effect
    onesim <- replicate(k, {
      sim <- matrix(rnorm(2*N), ncol=2)
      t.test(sim)$p.value 
    })
  }
  ks <- sort(onesim) < (1:k / k * 0.05)
  signifs <- max(0, which(ks == TRUE)) # largest k
  signifs > 0
})
mean(simusFDRtrue)
# FDR:
(sum(simusFDRtrue) - prop.true * 1000) / sum(simusFDRtrue)



#### Regression dilution ####
set.seed(321)
x <- runif(100)
y <- 2 + 2 * x + rnorm(100)

plot(x, y)
abline(lm(y~x))

for (i in 1:100){
  x2 <- x + rnorm(100, sd=0.1)
  abline(lm(y ~ x2), col="grey")
}

for (i in 1:100){
  x2 <- x + rnorm(100, sd=0.5)
  abline(lm(y ~ x2), col="grey")
}

dats <- data.frame(X=x, Y=y)

library(ggplot2)
ggplot(dats, aes(x=X, y=Y)) + 
  geom_point(size=1.5) + 
  theme_classic()


# possible example:
## GDP has high variability in global south as informal sector is not quantified: bias and heterogeneity!
## blood pressure in patients measured every day --> used to predict stroke risk
## ground water distance effect on pumping cost (Mieno & Brozovic 2016)
## impact of immigration on wages (Aydemir & Borjas 2011): wage ~ immigration share per county
set.seed(43210)
immi <- runif(500, 0.2, 0.4) # proportion immigrant in Canada's 5360 "labour cells" (Table 4)
slope <- -0.1 # (no intercept; impact on log10 weekly wage)
lgwage <- .1 * mean(immi) + slope * immi + rnorm(length(immi), sd=0.005) # intercept only to have values centred on 0
plot(immi, lgwage)
abline(h=0)
summary(lm(lgwage ~ immi)) # Aydemir's model has an R2 of 0.982, but the partial contribution of immi is not stated
coef(lm(lgwage ~ immi))
dats <- data.frame(immi = immi, lgwage = lgwage)

x <- immi + rnorm(length(immi), sd=.05)
plot(x, lgwage)
abline(h=0)
lm(lgwage ~ x)

dilut05 <- as.data.frame(t(replicate(20, coef(lm(lgwage ~ I(immi + rnorm(length(immi), sd=.05)))))))
dilut1 <- as.data.frame(t(replicate(20, coef(lm(lgwage ~ I(immi + rnorm(length(immi), sd=.1)))))))
colnames(dilut1) <- colnames(dilut05) <- c("Intercept", "slope")


library(ggplot2)

pdf(file="../figures/dilution.pdf", width=5.1, height=5)
ggplot(dats, aes(x=immi, y=lgwage)) +
  geom_abline(data=dilut05, aes(intercept=Intercept, slope=slope), col="grey") + 
  geom_abline(data=dilut1, aes(intercept=Intercept, slope=slope), col="grey") + 
  geom_text(x=.403, y=-0.0015, label="  0.1", size=4, check_overlap = TRUE) +
  geom_text(x=.403, y=-0.0075, label="0.05", size=4, check_overlap = TRUE) +
  labs(x="proportion of immigrants", y="centred log wage") + 
  geom_point(size = 1.5) + 
  geom_smooth(method = "lm", col = "black") +
  theme_classic()
dev.off()



# risk of stroke as a function of dyastolic blood pressure: from MacMahon et al. 1990, Fig. 1 left
# just to get some baseline for data simulation
x <- c(76, 84, 91, 98, 105)
y <- c(0.3, 0.51, 1, 1.8, 3.8)
plot(x, log10(y))
abline(lm(log10(y)~ x))
lm(log10(y )~ x)
# simulation would require data and analysis at patient level: too involved as example!

# references:
## MacMahon sets the tone, using stroke as case study
## Frost & Thompson compare different methods to correct for dilution bias, find MacMahon's approach to work best, but not be extendable to multivariate regression
## Hausman 2001 for econometrics
## https://en.wikipedia.org/wiki/Regression_dilution for various other references and extensions
## https://en.wikipedia.org/wiki/Errors-in-variables_models is the more general webpage for this topic











#### loss of gradient: using intended or realised range? ####

##### real data: intended vs realised #####
jena <- read.csv("../Data/JenaExperiment/487_3_data.csv", stringsAsFactors = T)
summary(jena)
jena22 <- jena[jena$year==2022 & jena$season=="Total",]
table(jena22$sowndiv)
plot(jena22[,c("sowndiv", "realizedSR")])
abline(0,1)
cor(jena22[, "sowndiv"], jena22[,"realizedSR"])

Sint <- jena22$sowndiv
Sreal <- jena22$realizedSR

with(jena[jena$year==2022 & jena$season=="Total",], table(sowndiv, realizedSR))
# simu based on these values: 
Sint <- rep(c(1,2,4,8,16,60), times=c(14, 16, 16, 16, 14, 4) ) # from jena487
Sreal <- c(sample(c(1,2), 14, prob=c(0.9, 0.1), replace=T), #1
           sample(c(2,3,4), 16, prob=c(.85, .1, .05), replace=T), #2
           sample(c(1,2,3,4), 16, prob=c(1/16, 5/16, 1/2, 1/8), replace=T), #4
           sample(2:8, 16, prob=c(1,1,2,1,6,2,3)/16, replace=T), #8
           sample(6:12, 14, prob=c(1,0,0,2,6,3,2)/14, replace=T), #16
           sample(c(29,30,33,38), size=4, prob=rep(1/4, 4), replace=T)) #60

##### simulated loss of richness in ignorance of the real data #####
Sint <- rep(c(1,2,4,8,16,60), times=c(16, 16, 16, 16, 14, 4) ) # from Weisser et al. 2017
# 46/60 established (page 17), i.e. 14/60=23% did not
set.seed(123)
Sreal <- c(sample(c(1,2,3,4), 16, prob=c(0.9, 0.1), replace=T), #1
           sample(c(2,3), 16, prob=c(0.9, 0.1), replace=T), #2
           sample(c(3,4,5), 16, prob=c(0.2, 0.7, 0.1), replace=T), #4
           sample(c(5,6,7), 16, prob=c(0.25, 0.5, 0.25), replace=T), #8
           sample(11:14, 14, prob=c(.2, .4, .2, .2), replace=T), #16
           sample(25:34, size=4, prob=rep(0.1, 10), replace=T)) #60
plot(Sint, Sreal)
abline(0,1)
cor(Sint, Sreal) # "closely related"

##### modelling the effect of BD on EF #####

K=1 # +/- 0.5 95%CI 
Ymax=2 # /- 0.5 95%CI   # Cardinale et al. 2006
effect.real <- Ymax*Sreal/(K+Sreal) + rnorm(length(Sreal), 0, .15*Ymax) # SD easily 15%: Fig. 18 in Weisser
effect.int <- Ymax*Sint/(K+Sint) + rnorm(length(Sreal), 0, .15*Ymax)

#plot(effect.int, effect.real)
#abline(0,1)
#cor(effect.int, effect.real)

plot(effect.real ~ Sreal, log="x", xlim=c(1, 60), xlab="plant species richness", ylab="ecosystem function", las=1)
points(effect.int ~ Sint, pch=16, col="blue")

library(ggplot2)
dats <- cbind.data.frame(Sreal, Sint, effect.real, effect.int)
pdf(file="../figures/BDEFJena1.pdf", width=5, height=5)
ggplot(dats) + 
  geom_segment(aes(x=ifelse(Sint==Sreal, NA, Sint), y=effect.int, xend=Sreal, yend=effect.real), arrow = arrow(length=unit(.15, "inches"), angle=15, type="closed"), colour="grey80", linejoin="mitre") +
  geom_point(aes(y=effect.int, x=Sint), colour="grey50") + 
  geom_point(aes(y=effect.real, x=Sreal)) + 
  scale_x_continuous(trans='log2') +
  labs(x = "species richness", y = "ecosystem function") +
  theme_minimal()
dev.off()

#nls(effect.int ~ Ymax*Sint/(K+Sint), start=list("Ymax"=2, "K"=1)) # ideal
summary(fmJ <- nls(effect.real ~ Ymax*Sint/(K+Sint), start=list("Ymax"=2, "K"=1))) # what they did
summary(fmC <- nls(effect.real ~ Ymax*Sreal/(K+Sreal), start=list("Ymax"=2, "K"=1))) # what they should have done

pdf(file="../figures/BDEFJena2.pdf", width=5, height=5)
ggplot(dats) + 
  geom_point(aes(y=effect.int, x=Sint), colour="grey50") + 
  geom_point(aes(y=effect.real, x=Sreal)) + 
  geom_function(fun = \(x) 2*x/(1+x), linetype=1) +
  geom_function(fun = \(x) coef(fmJ)[1]*x/(coef(fmJ)[2]+x), colour = "grey70", linetype=2, linewidth=1) +
  geom_function(fun = \(x) coef(fmC)[1]*x/(coef(fmC)[2]+x), colour = "black", linetype=2, linewidth=1) +
  scale_x_continuous(trans='log2') +
  labs(x = "species richness", y = "ecosystem function") +
  annotate("point", x=12, y=c(.6, .7), col=c("grey70", "black")) +
  annotate("text", x=13, y=c(.7, .6), label=c("realised richness", "intended richness"), hjust=0, col=c("black", "grey60")) +
  theme_minimal()
dev.off()

# summary: no big deal, although biased in the asymptote; effect in either case weaker than real in the way they analysed the data
# check for effect of error on x, i.e. rerun with variable x-values;
# this builds on the regression dilution section

summary(fmJ)$coefficients[1:2, 1:2]
summary(fmC)$coefficients[1:2, 1:2] # repeat, use weighted averaging to compute final values

K=1 # +/- 0.5 95%CI     # Cardinale et al. 2006
Ymax=2 # /- 0.5 95%CI   
N = 1000
simu.results <- data.frame("Ymax.real"=rep(NA, N), "K.real"=NA, "Ymax.int"=NA, "K.int"=NA)
Sint <- rep(c(1,2,4,8,16,60), times=c(14, 16, 16, 16, 14, 4) ) # from jena487
set.seed(123)
for (i in 1:N){
  # generate Sreal-values:
  Sreal.simu <-  c(sample(c(1,2), 14, prob=c(0.9, 0.1), replace=T), #1
                    sample(c(2,3,4), 16, prob=c(.85, .1, .05), replace=T), #2
                    sample(c(1,2,3,4), 16, prob=c(1/16, 5/16, 1/2, 1/8), replace=T), #4
                    sample(2:8, 16, prob=c(1,1,2,1,6,2,3)/16, replace=T), #8
                    sample(6:12, 14, prob=c(1,0,0,2,6,3,2)/14, replace=T), #16
                    sample(c(29,30,33,38), size=4, prob=rep(1/4, 4), replace=T)) #60
  # Sreal.simu <- Sreal # using observed realised S 
  # compute real effect:
  effect.real.simu <- Ymax*Sreal.simu/(K+Sreal.simu) + rnorm(length(Sreal.simu), 0, .15*Ymax)
  # fit MM-kinetics, once on real, once on intended S:
  fmC <- nls(effect.real.simu ~ Ymax * Sreal.simu / (K + Sreal.simu), start=list("Ymax"=2, "K"=1))
  fmJ <- nls(effect.real.simu ~ Ymax * Sint/(K + Sint), start=list("Ymax"=2, "K"=1))
  simu.results[i,] <- c(coef(fmC), coef(fmJ))
  rm(Sreal.simu, fmC, fmJ)
}
summary(simu.results)

dats <- stack(simu.results)
dats$coef <- rep(rep(c("Ymax", "K"), each=N), 2)
dats$treat <- rep(c("real", "intended"), each=2*N)

pdf(file="../figures/BDEFJena3.pdf", width=5, height=5)
ggplot(dats, aes(y=values, x=coef, fill=treat)) +
  geom_boxplot(color= c("black", "grey70", "black", "grey70")) + 
  scale_fill_manual(values=c("grey70", "black")) +
  labs(x = "", y = "coefficient") +
  scale_x_discrete(breaks = NULL) +
  guides(fill=guide_legend(title="S value used\nin analysis")) +
  annotate("text", x=c(1, 2), y=.2, label=c("K* = 1", "Ymax* = 2"), size=5) +
  theme_minimal() + 
  theme(legend.position = c(.8, 0.5))
dev.off()
# slight underestimation of BD effect!



#### Mixed-effect model non-normal hyperparameter #### 
set.seed(12345)
k=20
classes <- rep(LETTERS[1:15], each=k)
# draw random means from log-normal:
curve(dlnorm(x, meanlog = 1, sdlog=2), from=0, to=15)
randomEffectClassNorm <- rep(rnorm(k, mean = 0, sd=10), each=15)
randomEffectClassNorm <- randomEffectClassNorm - mean(randomEffectClassNorm) # centre random effect
randomEffectClass <- rep(rlnorm(k, meanlog = 1, sdlog=2), each=15)
randomEffectClass <- randomEffectClass - mean(randomEffectClass) # centre random effect
plot(randomEffectClass)
x <- seq(0, 100, len=15*k)
y <- round(0 + randomEffectClass + 1*x + rnorm(k*15, mean=0, sd=10), 1)
yn <- round(0 + randomEffectClassNorm + 1*x + rnorm(k*15, mean=0, sd=10), 1)
y2 <- round(0 + randomEffectClass + 0.1*randomEffectClass*x + rnorm(k*15, mean=0, sd=10), 1)
dats <- data.frame("Y"=y, "Yn"=yn, "X"=x, "group"=randomEffectClass, "groupID"=classes)

library(ggplot2)
ggplot(dats, aes(X, Y, colour=classes)) + geom_point() + theme_minimal()
ggplot(dats, aes(X, Yn, colour=classes)) + geom_point() + theme_minimal()
#ggplot(dats, aes(X, Y2, colour=classes)) + geom_point() + theme_minimal()

library(lme4)
summary(flme <- lmer(Y ~ X + (1|classes), data=dats))
summary(flmen <- lmer(Yn ~ X + (1|classes), data=dats))
#summary(flme2 <- lmer(Y2 ~ X + (X|classes), data=dats))
library(lattice)
dotplot(ranef(flme, condVar=TRUE))

## repeat 1000 times to evaluate how often the grand intercept is inside 95% CI for lognormal and normal, respectively; should be widely out for lognormal



# pupils, classes, schools -> does non-normal class random effect bias estimates of school ranking?
# We use a simplified version of PISA simulations: students in classes in schools in provinces.

# https://rdrr.io/cran/simstudy/f/vignettes/clustered.Rmd
library(simstudy)
# start at highest level, then move down from there!
gen.school <- defData(varname="s0", dist = "normal", formula = 0, variance = 3, id = "idSchool") # the names of the school
gen.school <- defData(gen.school, varname = "nClasses", dist = "noZeroPoisson", formula = 3) # their attribute: number of classes

set.seed(123)
dtSchool <- genData(8, gen.school)
dtSchool

gen.class <- defDataAdd(varname = "c0", dist = "normal", formula = 0, variance = 2) # the name of the class
gen.class <- defDataAdd(gen.class, varname = "nStudents", dist = "noZeroPoisson", formula = 20) # the attribute: number of students

dtClass <- genCluster(dtSchool, "idSchool", numIndsVar = "nClasses",level1ID = "idClass") # ...
dtClass <- addColumns(gen.class, dtClass)
dtClass

gen.student <- defDataAdd(varname="Male", dist="binary", formula=0.5)
gen.student <- defDataAdd(gen.student, varname="age", dist = "uniform", formula="14; 16.5")
gen.student <- defDataAdd(gen.student, varname="test", dist = "normal",
                          formula = "50 - 5*Male + s0 + c0", variance = 2)
dtStudent <- genCluster(dtClass,cLevelVar="idClass", numIndsVar = "nStudents",
                        level1ID = "idChild")

dtStudent <- addColumns(gen.student, dtStudent)

library(lmerTest)
summary(lmer(test ~ Male + age + (1|idSchool/idClass), data=dtStudent))


## province - school - class - student, with effect of school, i.e. their formula != 0
library(simstudy)
library(lmerTest)
# start at highest level, then move down from there!
# the names of the provinces and their effect
gen.province <- defData(varname="Province", dist = "normal", formula = 0, variance = 4, id = "idProvince")
#exaProvince <- genData(16, gen.province)
#exaProvince
gen.province <- defData(gen.province, varname = "nSchools", dist = "noZeroPoisson", formula = 20) # their attribute: number of schools
exaProvince <- genData(16, gen.province)
exaProvince
# effect of school, and ID:
gen.school <- defDataAdd(varname = "School", dist = "normal", formula = "Province", variance = 6) 
gen.school <- defData(gen.school, varname = "affluence", dist = "normal", formula = 2, variance=1) # attribute: GDP
gen.school <- defDataAdd(gen.school, varname = "nClasses", dist = "noZeroPoisson", formula = 4) # attribute: number of classes
# make nested data, and assign ID to schools
exaSchool <- genCluster(exaProvince, "idProvince", numIndsVar = "nSchools", level1ID = "idSchool") # ...
exaSchool <- addColumns(gen.school, exaSchool)
exaSchool
# effect of class, and ID:
gen.class <- defDataAdd(varname = "Class", dist = "normal", formula = "affluence*1 - 1 + School", variance = 8) # HERE EFFECT OF AFFLUENCE!!
gen.class <- defDataAdd(gen.class, varname = "nStudents", dist = "noZeroPoisson", formula = 22) # attribute: number of students
# make nested data, and assign ID to schools
exaClass <- genCluster(exaSchool, "idSchool", numIndsVar = "nClasses", level1ID = "idClass") 
exaClass <- addColumns(gen.class, exaClass)
exaClass
# effect of student, and ID:
gen.student <- defDataAdd(varname="test", dist = "normal", formula = "0 + Class", variance = 10)
simStudent <- genCluster(exaClass,cLevelVar="idClass", numIndsVar = "nStudents", level1ID = "idStudent")
simStudent <- addColumns(gen.student, simStudent)
simStudent
summary(simStudent)

summary(fm <- lmer(test ~ affluence + (1|idProvince/idSchool/idClass), data=simStudent))
#extract variances:
as.data.frame(VarCorr(fm))[-c(2,3,5)]
# extract estimates for intercept and effect of affluence:
fm@beta
# combine them in single data.frame:
res <- rbind(as.data.frame(VarCorr(fm))[-c(2,3,5)], data.frame("grp"=c("intercept", "affluence"), "vcov"=fm@beta))
colnames(res) <- c("term", "estimate")
res


#summary(lmer(test ~ affluence + (1|idProvince/idSchool/idClass), data=simStudent))$coefficient[2,1]
# expands to:
#summary(lmer(test ~ GDPprovince + (1|idProvince) + (1|idProvince : idSchool) + (1|idProvince : idSchool : idClass) , data=simStudent))

#run another simulation by repeating the gen... and later addColumns steps:
sims <- replicate(100, {
  exaProvince <- genData(16, gen.province)
  exaSchool <- genCluster(exaProvince, "idProvince", numIndsVar = "nSchools", level1ID = "idSchool") # ...
  exaSchool <- addColumns(gen.school, exaSchool)
  exaClass <- genCluster(exaSchool, "idSchool", numIndsVar = "nClasses", level1ID = "idClass") 
  exaClass <- addColumns(gen.class, exaClass)
  simStudent <- genCluster(exaClass, cLevelVar="idClass", numIndsVar = "nStudents", level1ID = "idStudent")
  simStudent <- addColumns(gen.student, simStudent)
  print("one run almost done")
  fm <- lmer(test ~ affluence + (1|idProvince/idSchool/idClass), data=simStudent)
  c(as.data.frame(VarCorr(fm))[,4], fm@beta)
})
rownames(sims) <- c(as.data.frame(VarCorr(fm))[,1], "intercept", "affluence")
summary(t(sims))
plot(density(sims[4,]))



# now non-normal for hyperparameter on class:
# use mixture of two normals (https://www.rdatagen.net/post/adding-mixture-distributions-to-simstudy/):
plot(density(c(rnorm(0.8*1E6, mean=-1, sd=sqrt(8)), rnorm(0.2*1E6, mean=7, sd=sqrt(8)))))
dats <- data.frame(x=c(rnorm(0.8*1E6, mean=-1, sd=sqrt(8)), rnorm(0.2*1E6, mean=7, sd=sqrt(8))))
library(ggplot2)

# compute mean and variance of mixture: https://en.wikipedia.org/wiki/Mixture_distribution
## mu:
(mu <- 0.8*-1 + 0.2*7)
##var:
0.8*(8 + 1) + 0.2*(8+49) - mu^2

pdf(file="../figures/MEnormmix.pdf", width=5.1, height=5)
ggplot(dats, aes(x=x)) + geom_density() + theme_minimal()
dev.off()

gen.class <- defDataAdd(varname = "Class1", dist = "normal", formula = "affluence*1 - 1 + School", variance = 8) # HERE EFFECT OF AFFLUENCE!!
gen.class <- defDataAdd(gen.class, varname = "Class2", dist="normal", formula = "affluence*1 + 7 + School", variance = 8)
gen.class <- defDataAdd(gen.class, varname="Class", formula="Class1 | 0.8 + Class2 | 0.2", dist="mixture")
gen.class <- defDataAdd(gen.class, varname = "nStudents", dist = "noZeroPoisson", formula = 22) # attribute: number of students

# analysis: seems effect is nicely absorbed by variance estimate, not affluence estimate!

sims2 <- replicate(100, {
  exaProvince <- genData(16, gen.province)
  exaSchool <- genCluster(exaProvince, "idProvince", numIndsVar = "nSchools", level1ID = "idSchool") # ...
  exaSchool <- addColumns(gen.school, exaSchool)
  exaClass <- genCluster(exaSchool, "idSchool", numIndsVar = "nClasses", level1ID = "idClass") 
  exaClass <- addColumns(gen.class, exaClass)
  simStudent <- genCluster(exaClass, cLevelVar="idClass", numIndsVar = "nStudents", level1ID = "idStudent")
  simStudent <- addColumns(gen.student, simStudent)
  print("one run almost done")
  fm <- lmer(test ~ affluence + (1|idProvince/idSchool/idClass), data=simStudent)
  c(as.data.frame(VarCorr(fm))[,4], fm@beta)
})
rownames(sims2) <- c(as.data.frame(VarCorr(fm))[,1], "intercept", "affluence")
#save(sims, sims2, file="schooltest100Simulations.RData")
summary(t(sims2))

plot(density(sims2[6,]))


# diagnostics:
lattice::dotplot(ranef(fm, condVar=TRUE, whichel="idClass:(idSchool:idProvince)"))
library(DHARMa)
simres <- simulateResiduals(fm)
testUniformity(simres, plot=F)
plot(simres)


library(ggplot2)
dats <- rbind.data.frame(t(sims), t(sims2))
dats
dats1 <- stack(dats)
dats1$distr <- rep(rep(c("normal", "non-normal"), each=100), 6)

pdf(file="../figures/5_MEestimates.pdf", width=8, height=5)
plot1 <- ggplot(dats1, aes(y=values, x=ind, fill=distr)) + geom_boxplot() + scale_fill_grey(start=0.5) + 
  scale_x_discrete(labels=c("class", "school", "province", "student", "Intercept", "affluence"), name="") +
  scale_y_continuous(name="estimate") +
  theme_minimal()
plot1 + 
  theme(strip.text.x = element_blank(), legend.title=element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.9,.85)  )
dev.off()



#### Extrapolation with changing correlation structure ####
# Red imported fire ant: use S.Am. distribution to fit, everywhere else to predict
# Centaurea maculosa (Brönnimann et al. 2007)
#download worldclim
#a) demonstrate "shift":
#* get data points for native range
#* get data points for invaded range(s)
#* do PCA of all occurrences ("is there a shift?") 
#b) draw random points and do PCA: same "shift"?
#c) place niche, i.e. draw points as function of predictor, e.g. plogis(y) = temp*1 - temp^2*0.5; plot PCA and show "niche shift"
# no need to fit models: point to Dormann et al. 2013 for that

# Download worldclim (https://gis.stackexchange.com/questions/227585/using-r-to-extract-data-from-worldclim)
library(raster)
r <- getData("worldclim", var="bio", res=10) # note: only down to -60°S
# in the future:
#library(geodata)
#r <- worldclim_global(var="bio", res=10, path=getwd()) # note: different column names!!
library(sp)

#r <- r[[c(1,12)]]
#names(r) <- c("Temp","Prec")
#points <- spsample(as(r@extent, 'SpatialPolygons'),n=100, type="random") 
#extract(r, )
# Get fire ant native range polygon: (https://www.r-bloggers.com/2021/03/downloading-and-cleaning-gbif-data-with-r/; https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html)
library(rgbif)
#gbif_fireant <- occ_search(scientificName = "Solenopsis invicta", hasCoordinate = TRUE, limit=200000) # requires no registration, but is discouraged ...
#save(gbif_fireant, file="gbif_fireant.Rdata")

load("gbif_fireant.Rdata")
fireant_coords <- gbif_fireant$data[ , c("decimalLongitude", "decimalLatitude", "occurrenceStatus", "continent", "coordinateUncertaintyInMeters", "institutionCode", "references", "occurrenceRemarks")]
summary(as.factor(fireant_coords$continent))
library(maps)
map("world", xlim = range(fireant_coords$decimalLongitude), ylim = range(fireant_coords$decimalLatitude))  # if the map doesn't appear right at first, run this command again
points(fireant_coords[ , c("decimalLongitude", "decimalLatitude")], pch = ".", col="red", cex=4)

climAllFireAntLocations <- extract(r, fireant_coords[ , c("decimalLongitude", "decimalLatitude")])
summary(climAllFireAntLocations)
# which(rowSums(is.na(climAllFireAntLocations)) > 0) # same as following: 
#[1]  1287  1915  5363 13538 13854 13871 13883
which(is.na(climAllFireAntLocations), arr.ind=T) # quite a few: remove them!
climAllFireAntLocations[-unique(which(is.na(climAllFireAntLocations), arr.ind=T)[,1]), ]
missing.index <- unique(which(is.na(climAllFireAntLocations), arr.ind=T)[,1])
climAllFireAntLocations <- climAllFireAntLocations[-missing.index, ]
summary(climAllFireAntLocations)

fireAntPCA <- prcomp(climAllFireAntLocations, scale=T)
indexSA <- which(gbif_fireant$data[-missing.index,]$continent=="SOUTH_AMERICA")
indexNA <- which(gbif_fireant$data[-missing.index,]$continent=="NORTH_AMERICA")
indexAS <- which(gbif_fireant$data[-missing.index,]$continent=="ASIA")

pdf(file="../figures/CLIMPCA.pdf", width=6, height=6)
par(mar=c(1,1,1,1))
biplot(fireAntPCA, xlabs=rep(".", nrow(fireAntPCA$x)), col="grey")
points(fireAntPCA$x[indexNA, 1], fireAntPCA$x[indexNA, 2], col="grey20", pch=15, cex=.2)
points(fireAntPCA$x[indexSA, 1], fireAntPCA$x[indexSA, 2], col="grey80", pch=15, cex=.2)
#points(fireAntPCA$x[indexAS, 1], fireAntPCA$x[indexAS, 2], col="green", pch=15)
dev.off()

pdf(file="../figures/5_CLIMPCA.pdf", width=6, height=6)
library(ggbiplot)
ggbiplot(fireAntPCA, alpha=0, varname.size=0.00001, varname.color="white") + 
  geom_point(data=as.data.frame(fireAntPCA$x[indexNA, ]), aes(x=PC1, y=PC2), colour="grey") + 
  geom_point(data=as.data.frame(fireAntPCA$x[indexSA, ]), aes(x=PC1, y=PC2), colour="black") + 
  theme_minimal()
dev.off()

# position on PC1 is rather different, that on PC2 not:
plot(density(fireAntPCA$x[indexNA, 1]), col="blue")
lines(density(fireAntPCA$x[indexSA, 1]), col="red")
lines(density(fireAntPCA$x[indexAS, 1]), col="green")

pdf(file="../figures/CLIMcorimages.pdf", width=8, height=4.25)
par(mfrow=c(1,2), mar=c(1,1,1,1))
image(cor(climAllFireAntLocations[indexNA, c(1,2,5,6, 12:14,16,17)], use="pairwise.complete.obs", method="sp"), col=grey.colors(50), axes=F)
image(cor(climAllFireAntLocations[indexSA, c(1,2,5,6, 12:14,16,17)], use="pairwise.complete.obs", method="sp"), col=grey.colors(50), axes=F)
dev.off()

# are these niche shifts?
mapa <- map("world", regions = c("Brazil", "Uruguay", "Argentina", "Suriname", "French Guiana", "Colombia", "Venezuela", "Bolivia", "Ecuador", "Chile", "Paraguay", "Peru", "Guyana"), fill=T, ylim=c(-60,20))
mapb <- maptools::map2SpatialPolygons(mapa, IDs=mapa$names, proj4string=CRS("+proj=longlat"))
mapc <- map("state", fill=T) # only 48!
mapd <- maptools::map2SpatialPolygons(mapc, IDs=mapc$names, proj4string=CRS("+proj=longlat"))
randomNativepoints <- spsample(mapb, n=324, type="random") 
randomUSpoints <- spsample(mapd, n=13553, type="random") 
#randomPoints <- rbind(mapb, mapd)
# plot(randomPoints) # turns spatialPoints into spatialPolygons (somehow)
# use manual approach to combine the two data sets:
randomPoints <- randomNativepoints
randomPoints@coords <- rbind(randomNativepoints@coords, randomUSpoints@coords)
randomPoints@bbox <- apply(cbind(randomNativepoints@bbox, randomUSpoints@bbox), 1, range)
plot(randomPoints) # works
climAllRandomLocations <- extract(r, randomPoints) # for terra::extract: vect(randomPoints)
summary(climAllRandomLocations)
randomPCA <- prcomp(climAllRandomLocations, scale=T)
randomPCA$rotation
biplot(randomPCA, xlim=c(-0.004, 0.006), ylim=c(-.004, 0.004))
points(randomPCA$x[325:13877, 1], randomPCA$x[325:13877, 2], col="blue", pch=15)
points(randomPCA$x[1:324, 1], randomPCA$x[1:324, 2], col="red", pch=15)
# yields Fig. 5.9
plot(density(randomPCA$x[325:13877, 1]), col="blue")
lines(density(randomPCA$x[1:324, 1]), col="red")



## repeatedly:
set.seed(123)
PC1s <- replicate(25, {
  randomNativepoints <- spsample(mapb, n=324, type="random") 
  randomUSpoints <- spsample(mapd, n=13553, type="random") 
  randomPoints <- randomNativepoints
  randomPoints@coords <- rbind(randomNativepoints@coords, randomUSpoints@coords)
  randomPoints@bbox <- apply(cbind(randomNativepoints@bbox, randomUSpoints@bbox), 1, range)
  climAllRandomLocations <- extract(r, randomPoints)
#  if (nrow(climAllRandomLocations) != 324+13553) stop()
  randomPCA <- prcomp(na.omit(climAllRandomLocations), scale=T)
  randomPCA$x[,1]
})
# this could be a list, if some data points are NA and are hence omitted by na.omit

# make Fig. 5.11:
pdf(file="../figures/CLIMdensitiesPC1.pdf", width=5, height=5)
par(mar=c(3,4,1,1))
plot(density(PC1s[325:13877, 1]), col="grey20", xlim=c(-20, 7), las=1, main="")
lines(density(PC1s[1:324, 1]), col="grey80")
for (i in 2:25){
  lines(density(PC1s[325:13877, i]), col="grey20")
  lines(density(PC1s[1:324, i]), col="grey80")
}
legend("topleft", bty="n", legend=c("USA48", "South America"), lty=1, col=c("grey20", "grey80"), cex=1)
dev.off()



##### fit GLM to PA fire ant #####
# we gloss over a few things: random points are not the same as absences; no interactions in the GLM, no standardisation of predictors, no removal of very highly correlated predictors, calibration of randomForest, ...

library(rgbif)
#gbif_fireant <- occ_search(scientificName = "Solenopsis invicta", hasCoordinate = TRUE, limit=200000) # requires no registration, but is discouraged ...
load("gbif_fireant.Rdata")
fireant_coords <- gbif_fireant$data[ , c("decimalLongitude", "decimalLatitude", "occurrenceStatus", "continent", "coordinateUncertaintyInMeters", "institutionCode", "references", "occurrenceRemarks")]
climAllFireAntLocations <- extract(r, fireant_coords[ , c("decimalLongitude", "decimalLatitude")])
missing.index <- unique(which(is.na(climAllFireAntLocations), arr.ind=T)[,1])
indexSA <- which(gbif_fireant$data[-missing.index,]$continent=="SOUTH_AMERICA")
indexUSA <- which(gbif_fireant$data[-missing.index,]$continent=="NORTH_AMERICA")

set.seed(123)
randomNativePoints <- spsample(mapb, n=length(indexSA), type="random") 
systematicUSPoints <- spsample(mapd, n=length(indexUSA), type="regular") 
climAllSARandomLocations <- extract(r, randomNativePoints)
climAllUSASystematicLocations <- extract(r, systematicUSPoints)
datsSA <- na.omit(data.frame("PAfireant"=rep(c(0,1), each=length(indexSA)), rbind(climAllSARandomLocations, climAllFireAntLocations[indexSA,])))
datsUSA <- na.omit(data.frame("PAfireant"=rep(c(0,1), times=c(nrow(climAllUSASystematicLocations), length(indexUSA))), rbind(climAllUSASystematicLocations, climAllFireAntLocations[indexUSA,])))
# fit the model to south American data, only to get an idea which X are important:
f <- as.formula(paste("PAfireant ~ ", paste0("bio", 1:19, collapse=" + "), "+", paste0("I(bio", 1:19, "^2)", collapse=" + ")))
fm <- glm(f, data=datsSA, family=binomial)
fmstep <- step(fm, trace=F)
summary(fmstep)


##### simulate similar species ####
set.seed(123)
ySA <-  plogis(with(datsSA[,-1], -150.2 + .6* bio1 + 1.2*bio10 + -.003 * bio10^2 - 0.7*bio11 + 0.001 * bio11^2))
mean(ySA)
hist(ySA)
##### fit randomForest and predict to USA #####
yUSA <-  plogis(with(data.frame(datsUSA[,-1]), -150.2 + .6* bio1 + 1.2*bio10 + -.003 * bio10^2 - 0.7*bio11 + 0.001 * bio11^2))
mean(yUSA)

library(ranger)
fr <- ranger(as.factor(ySA>0.5) ~., data=datsSA)
preds2USA <- predict(fr, data=datsUSA)$predictions # I'd rather use probabilities ...
summary(preds2USA)

pdf(file="../figures/CLIMfakeantmaps.pdf", width=15, height=6)
par(mfrow=c(1,2))
map("state", fill=F)
points(systematicUSPoints, col=ifelse(preds2USA==TRUE, "white", "black"), pch=15, cex=.5)
mtext("prediction from SA", side=3, line=0, adj=0, cex=1.5)
map("state", fill=F, add=T)

map("state", fill=F)
points(systematicUSPoints, col=ifelse(yUSA < mean(yUSA), "white", "black"), pch=15) # use prevalence as threshold!
mtext("simulated truth", side=3, line=0, adj=1, cex=1.5)
map("state", fill=F, add=T)
dev.off()