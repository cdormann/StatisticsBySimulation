library(pwr)
library(retrodesign)
library(ggplot2)
library(Hmisc)

#1. Reproduce simulation from Lane and Dunlap (1978)

#2 groups
#3 effect sizes (D = 4, 8, 16, corresponding to 0.25, .5, 1 sd units)
#4 sample sizes (n = 5, 10, 15 ,20 per group; I did 5,10,50,100)
#Data drawn from normal distributions
#Analysis using one-way ANOVA; parameter is omega^2
#5000 replicates

library(effectsize)
library(parameters)

set.seed(123)
output <- data.frame("n"=numeric(), "effect"=numeric(), "pval"=numeric(), "omega2"=numeric()) 
n <- c(5, 10, 25, 50, 75, 100) #number of individuals per group
eff <- c(4, 8, 16) #actual difference between groups
nreps <- 5000 #number of simulations
m <- 1 #counter
for (i in 1:length(n)) {
  for (j in 1:length(eff)) { 
    for (k in 1:nreps) {
      group1 <- rnorm(n[i], 0, 16)
      group2 <- rnorm(n[i], eff[j], 16)
      temp <- stack(data.frame(group1, group2))
      names(temp) <- c("IQ", "group")
      a <- model_parameters(aov(IQ ~ group, data=temp),
                            effectsize_type = "omega")
      output[m,1] <- n[i]
      output[m,2] <- eff[j]
      output[m,3] <- a$p[1]
      output[m,4] <- a$Omega2[1]
      m <- m + 1
    }}}

# Stuff for plotting
# get mean values for each p value threshold
pvals <- c(.001, .01, .05, .1, 1.0) #cut points for p values
output2 <- NULL
for (i in 1:length(pvals)) {
  output2a <- aggregate(omega2~n+effect, data=output[output$pval<pvals[i],], mean)
  output2a$pval <- pvals[i]
  output2 <- rbind(output2,output2a)
}
output2$pval <- as.factor(output2$pval)
output2$effect <- as.factor(output2$effect)

#create labels to annotate plot

omegas <- round(output2[output2$n==5 & output2$effect==4, ]$omega2, 2)
ann_text <- data.frame(n=c(-9, -7, -6, -2, -2), omega2=omegas + 0.02, 
                     effect=factor(rep(4, length(omegas)), levels=c("4", "8", "16")), 
                     pval=as.factor(pvals), label=c("0.001", "0.01", "0.05", "0.1", "1.0"))
supp.labs <- c("Effect = 0.25 SD", "Effect = 0.5 SD", "Effect = 1.0 SD")
names(supp.labs) <- c("4","8","16")
reflines <- data.frame(effect=factor(c(4,8,16)), yint=c(.0154, .0588, .2))

#plot it
pdf("../figures/9_LaneDunlap.pdf", width=14, height=5)
ggplot(output2, aes(x=n, y=omega2, group=pval)) +
  geom_line(aes(linetype=pval), show.legend = FALSE) +
  geom_hline(data=reflines, aes(yintercept=yint), color="grey") +
  geom_point() +
  scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dashed")) +
  facet_grid(. ~ effect, scales="free_x", labeller=labeller(effect=supp.labs)) +
  geom_text(data = ann_text, mapping = aes(x=n, y=omega2, label=label), 
            size=3, hjust="left") +
  labs(x="Group (sample) size", y=expression(omega^2)) +
  theme_minimal() +	
  theme(axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=16),
        axis.text = element_text(size=11),
        strip.text = element_text(size=11)) 
dev.off()



#### MDD simulation (Mair et al. 2020)
# This is the code found in the zenodo/github repository for the Mair et al. paper!
source("MDD.R", local = T)

##### Simulation 1: Can the MDD discriminate between true and false negative tests? #####
# Define number of experiments to be simulated:
size <- 2000

# Define experimental/ real world parameters:
experiments <- data.frame(
  effect = c(rep("no",size/2), rep("yes",size/2)),
  effectsize = c(rep(0,size/2), rep(0.5,size/2)), #1/2 no effect, 1/2 effect of 0.5 
  theta = runif(size, max = 50, min=1),
  sampleSize = sample(3:10, size, replace = T),
  abund = sample(10:100, size, replace=T)
)
experiments$variance <- experiments$theta * experiments$abund

# Prepare empty output dataframe:
out <- data.frame(
  p = rep(NA,size),
  pMDD = rep(NA,size),
  pMDE = rep(NA,size),
  pCI = rep(NA,size)
)

# Simulate and analyze experiments:
for (i in 1:size)
{
  N <- experiments$sampleSize[i]
  theta <- experiments$theta[i]
  abund <- experiments$abund[i]
  effectsize <- experiments$effectsize[i]
  
  # Simulate experiments:
  control <- rnbinom(n = N, mu = abund, size = abund/(theta-1))
  tmean <- abund * (1 - effectsize)
  treatment <- rnbinom(n = N, mu = tmean, size = tmean/(theta-1))
  
  # log-transform data prior to analysis:
  transc <- log((2*control)+1)
  transt <- log((2*treatment)+1)
  
  # Perform tests and calculations on transformed data:
  test <- t.test(transc, transt, alternative = "greater", var.equal = T)
  mdd.ln <- MDD(N1 = N, N2 = N, variance1 = var(c(transc-mean(transc),transt-mean(transt))),
                alpha=0.05, two.sided = F, var.equal = T)
  postpower <- power.t.test(n = N, delta = NULL, sd = sd(c(transc-mean(transc),transt-mean(transt))), 
                            sig.level = 0.05, alternative="one.sided", type = "two.sample", 
                            power = 0.8)
  upperCI <- mean(transc)-mean(transt) + 
    qt(0.05/1, df = 2*N-2, lower.tail = FALSE) * sqrt(var(c(transc-mean(transc),transt-mean(transt)))) * sqrt(2/N)
  
  # Back-transformation:
  mdd.abu <- (exp(mean(transc))-exp(mean(transc) - mdd.ln$mdd))/2
  mde.abu <- (exp(mean(transc))-exp(mean(transc) - postpower$delta))/2
  upperCI.abu <- (exp(mean(transc))-exp(mean(transc) - upperCI))/2
  
  # Relation to control mean:
  pMDD <- 100 * mdd.abu / ((exp(mean(transc))-1)/2)
  pMDE <- 100 * mde.abu / ((exp(mean(transc))-1)/2)
  pCI <- 100 * upperCI.abu / ((exp(mean(transc))-1)/2)
  
  # Store results in output dataframe:
  out[i,1] <- test$p.value
  out[i,2] <- pMDD
  out[i,3] <- pMDE
  out[i,4] <- pCI
}

# Combine experimental and output dataframe:
all <- cbind(experiments, out)

# Subset dataframe: include only non-significant results at alpha = 0.05
nonsig <- all[all$p >= 0.05,]

##### now turn to code in Figure4.R #####
### Note: Run code in 'Simulation1.R' first to obtain dataframe 'nonsig'.

## Calculation of error rates ----------------

# Define threshold levels:
thresholds = seq(from = 30, to = 100, by = 10)

# Prepare output dataframe:
ErrorRates <- data.frame(
  threshold = thresholds,
  mdd.ftrust = rep(NA,length(thresholds)),
  mdd.fmist = rep(NA,length(thresholds)),
  post.ftrust = rep(NA,length(thresholds)),
  post.fmist = rep(NA,length(thresholds)),
  CI.ftrust = rep(NA,length(thresholds)),
  CI.fmist = rep(NA,length(thresholds))
)

# Calculate false trust and false mistrust rates:
for (i in 1:length(thresholds))
{
  mdd.ftrust <- length(which(nonsig$pMDD <= ErrorRates$threshold[i] & 
                               nonsig$effect == "yes"))/length(which(nonsig$effect == "yes"))
  mdd.fmist <- length(which(nonsig$pMDD > ErrorRates$threshold[i] & 
                              nonsig$effect == "no"))/length(which(nonsig$effect == "no"))
  post.ftrust <- length(which(nonsig$pMDE <= ErrorRates$threshold[i] & 
                                nonsig$effect == "yes"))/length(which(nonsig$effect == "yes"))
  post.fmist <- length(which(nonsig$pMDE > ErrorRates$threshold[i] & 
                               nonsig$effect == "no"))/length(which(nonsig$effect == "no")) 
  CI.ftrust <- length(which(nonsig$pCI <= ErrorRates$threshold[i] & 
                              nonsig$effect == "yes"))/length(which(nonsig$effect == "yes"))
  CI.fmist <- length(which(nonsig$pCI > ErrorRates$threshold[i] & 
                             nonsig$effect == "no"))/length(which(nonsig$effect == "no")) 
  ErrorRates[i,2] <- mdd.ftrust
  ErrorRates[i,3] <- mdd.fmist
  ErrorRates[i,4] <- post.ftrust
  ErrorRates[i,5] <- post.fmist
  ErrorRates[i,6] <- CI.ftrust
  ErrorRates[i,7] <- CI.fmist
}


# pMDD:
plot(ErrorRates$mdd.ftrust, ErrorRates$mdd.fmist, xlim = c(-0.1, 1), ylim = c(-0.05,1), 
     type = "l", xlab = "False trust rate", ylab = "", bty = "l", las = 1, 
     cex.lab = 1.5, cex.axis = 1.5, col = "darkred", lwd = 2)
mtext("False mistrust rate", side = 2, line = 3.5, cex = 1.5)
text(ErrorRates$mdd.ftrust, ErrorRates$mdd.fmist, labels = ErrorRates$threshold, 
     cex= 1.5, col = "darkred", font = 4, pos = c(4, 4, 4, rep(3, 5)))
points(ErrorRates$mdd.ftrust, ErrorRates$mdd.fmist, pch = 20, col = "darkred")

# pMDE:
lines(ErrorRates$post.ftrust, ErrorRates$post.fmist, col = "darkgreen", lty = 3, lwd = 2)
text(ErrorRates$post.ftrust, ErrorRates$post.fmist, labels=ErrorRates$threshold, 
     cex= 1.5, col = "darkgreen", pos = c(4, 2, 4, 4, 3, 3, 3, 3))
points(ErrorRates$post.ftrust, ErrorRates$post.fmist, pch = 20, col = "darkgreen")

# pCI:
lines(ErrorRates$CI.ftrust, ErrorRates$CI.fmist, col = "darkblue", lty = 2, lwd = 2)
text(ErrorRates$CI.ftrust, ErrorRates$CI.fmist, labels=ErrorRates$threshold, 
     cex= 1.5, col = "darkblue", pos = c(rep(2, 4), rep(1, 4)))
points(ErrorRates$CI.ftrust, ErrorRates$CI.fmist, pch = 20, col = "darkblue")

# Add legend:
arrows(0.4 , 0.8, 0.5, 0.8, angle = 0, col = "darkblue", lty = 2, lwd = 2, length = 0)
arrows(0.4 , 0.7, 0.5, 0.7, angle = 0, col = "darkgreen", lty = 3, lwd = 2, length = 0)
arrows(0.4 , 0.6, 0.5, 0.6, angle = 0, col = "darkred", lty = 1, lwd = 2, length = 0)
text(0.52 , 0.8, "pCI", col = "darkblue", pos=4, offset= 0, cex = 1.5)
text(0.52 , 0.7, "pMDE", col = "darkgreen", pos=4, offset= 0, cex = 1.5)
text(0.52 , 0.6, "pMDD", col = "darkred", pos=4, offset= 0, cex = 1.5)

par(op)


library(ggplot2)
head(ErrorRates)
pdf(file="../figures/9_CIMDDMDE.pdf")
ggplot(ErrorRates, aes(mdd.ftrust, mdd.fmist)) + geom_line(size=1) +
  geom_point() +
  annotate("text", x=ErrorRates$mdd.ftrust-.02, y=ErrorRates$mdd.fmist, 
           label=ErrorRates$threshold) +
  geom_line(aes(post.ftrust, post.fmist), colour="grey", size=1) + 
  geom_point(aes(post.ftrust, post.fmist), colour="grey") +
  annotate("text", x=ErrorRates$post.ftrust+.03, y=ErrorRates$post.fmist, label=ErrorRates$threshold, color="grey") +
  geom_line(aes(CI.ftrust, CI.fmist), colour="grey30", linetype=2, size=1) +
  geom_point(aes(CI.ftrust, CI.fmist), colour="grey30") +
  annotate("text", x=ErrorRates$CI.ftrust-.02, y=ErrorRates$CI.fmist, 
           label=ErrorRates$threshold, color="grey30") +
  xlab("False trust rate") + ylab("False mistrust rate") +
  # legend:
  annotate("text", x=rep(0.76, 3), y=c(.85, .81, .77), label=c("pMDD", "pMDE", "pCI"), col=c("black", "grey", "grey30"), size=5, adj=0) + 
  annotate("segment", x=rep(0.65, 3), xend=rep(0.75, 3), y=c(.85, .81, .77), yend=c(.85, .81, .77), linetype=c(1, 1, 2), size=1, colour=c("black", "grey", "grey30")) + 
  theme_minimal()
dev.off()




#### Retrospective power Hoenig & Heisey 2001 ####
#2. Reproduce Fig from Hoenig and Heisey (2001) with inset from Nuzzo (2020)
#Relationship of retrospective power to p-value for post-hoc power for alpha = 0.05

#pvalue is a random variable P on [0,1]
#cdf of the pvalue is Pr(P <= p) = G_d(p), d is parameter value
#For the one-sample z test of H0: mu <=0 vs mu >0, data are rnorm(mu, sigma)
#G_d(p) = pnorm(Zp-d), where Zp = 100(1-p)th percentile of the standard normal distribution

#let percentiles be 0 -> 100 (so probabilities are 0 to 1)

qnorm(seq(0,1,.1))
#Zscores
zscores <- qnorm(seq(0,1,.005))
#pvalues:
pvals <- 1 - pnorm(zscores)
#pvals <- 1-pnorm(qnorm(seq(0,1,.005)))

#power
#power <- 1-(pnorm(qnorm(.95)-qnorm(seq(0,1,.005))))
power <- 1 - (pnorm(qnorm(.95) - zscores))

1 - (pnorm(qnorm(.95)-qnorm(.95)))

pdf(file="../figures/9_pvals_vs_power.pdf", width=5.5, height=6)
plot(pvals, power, type="l", xlab=expression(paste(italic("p"), "-value")), 
        ylab="Retrospective power", las=1)
   lines(c(.05,.05), c(-.1,.5), lty=2)
   lines(c(-.1,0.05), c(.5,.5), lty=2)
   temp <- subplot(
     plot(pvals, power, type="l", xlim=c(0,0.1),xlab=expression(paste(italic("p"), "-value")), 
          ylab="Retrospective power", las=1),
     x=c(.5,1), y=c(.5,1), pars=list(cex=.75)
   )
   
   op <- par(no.readonly=TRUE)
   par(temp)
   lines(c(.005,.005), c(-.025,.8), lty=2)
   lines(c(-.0025,0.005), c(.8,.8), lty=2)
   par(op)
dev.off()
   
library(ggplot2)
library(cowplot)
# repeat from above:
zscores <- qnorm(seq(0,1,.005))
pvals <- 1 - pnorm(zscores)
power <- 1 - (pnorm(qnorm(.95) - zscores))
dats <- data.frame(retro=power, pvals=pvals)
main.plot <- ggplot(dats, aes(pvals, power)) +
  geom_line() + 
  xlab(expression(paste(italic(p), "-value"))) +
  ylab("Retrospective power") +
  annotate("segment", x=0.05, xend=0.05, y=0, yend=0.5, linetype=2) +
  annotate("segment", x=0, xend=0.05, y=0.5, yend=0.5, linetype=2) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0), limits = c(0, 1)) +
  scale_x_continuous(expand=c(0,0.01), limits = c(0, 1))

inset.plot <- ggplot(dats, aes(pvals, power)) +
  geom_line() + 
  xlab(expression(paste(italic(p), "-value"))) +
  ylab("Retrospective power") +
  annotate("segment", x=0.0065, xend=0.0065, y=0.8, yend=0, linetype=2) +
  annotate("segment", x=0, xend=0.0065, y=0.8, yend=0.8, linetype=2) +
  theme_classic() +
  scale_y_continuous(expand=c(0,0), limits = c(0, 1)) +
  scale_x_continuous(expand=c(0,0), limits = c(0, .1))

pdf(file="../figures/9_pvals_vs_power.pdf", width=5.5, height=6)
ggdraw() +
  draw_plot(main.plot) +
  draw_plot(inset.plot, x = 0.44, y = .45, width = .54, height = .53)
dev.off()


#retrodesign figure

D_range <- c(seq(0,1,.01),seq(1,10,.1),100)
n <- length(D_range)
power <- rep(NA, n)
typeS <- rep(NA, n)
exaggeration <- rep(NA, n)

for (i in 1:n){
  a <<- retrodesign(D_range[i], 1)
  power[i] <- a$power
  typeS[i] <- a$type_s
  exaggeration[i] <- a$type_m
}



pdf("SMpower.pdf", height=2.5, width=3)
par(mfrow=c(1,2), mar=c(4,4,4,2), mgp=c(1.7,.5,0), tck=-.01)
plot(power, typeS, type="l", xlim=c(0,1.05), ylim=c(0,0.54), xaxs="i", yaxs="i",
     xlab="Power", ylab="Type-S error rate", bty="l", cex.axis=.9, cex.lab=.9, las=1)
#segments(.07, .01, .07, .1, col="gray")
plot(power, exaggeration, type="l", xlim=c(0,1.05), ylim=c(.5,11.5), xaxs="i", yaxs="i",
     xlab="Power", ylab="Type-M error", bty="l", yaxt="n", cex.axis=.9, cex.lab=.9)
axis(2, seq(1,11,2), las=1)
#segments(.05, 1.25, .55, 1.25, col="gray")
dev.off()

set.seed(123)
retrodesign(.15, .055)
simplot<-sim_plot(.15,.055, alpha = 0.05, df = Inf, n.sims = 5000, gg = TRUE)

simplot +
  ylim(-.05,.35) +
  theme_minimal()

   
   
   
   
#### post hoc prediction uncertainty ####

#Idea: 
#  fit to a data set; 
#  predict to new data (interpolation & extrapolation); 
#  get SE of model prediction in different ways;
#  compare coverage of model prediction relative to truth.


library(ranger)

set.seed(1234)
X <- runif(100, 0 ,1)
Y <- rnorm(100, mean=1 - 4*X, sd=1)
plot(X, Y)

## fit models:
fglm <- lm(Y ~ X)
franger <- ranger(Y ~ X, data=cbind.data.frame(X=X, Y=Y), keep.inbag=T)
fqrf <- ranger(Y ~ X, data=cbind.data.frame(X=X, Y=Y), keep.inbag=T, quantreg=T) # quantile regression
fbootrf <- replicate(25, ranger(Y ~ X, data=cbind.data.frame(X=X, Y=Y)[sample(100, replace=T),]), simplify=F) # bootstrap of ranger


par(mfrow=c(1,2))
plot(Y, predict(fglm))
  abline(0,1)
plot(Y, predict(franger, data=data.frame(X=X))$prediction)
  abline(0,1)

newX <- c(0.5, 1.1) # seq(0, 1.2, len=101)
truth <- 1 - 4*newX

# predict:
predsGLM <- predict(fglm, newdata=data.frame(X=newX), se.fit=T)
predsranger <- predict(franger, data=data.frame("X"=newX), type="se")
predsqrf <- predict(fqrf, data=data.frame("X"=newX), type="quantiles", quantiles=c(0.025, 0.5, 0.975))
predsbootrf <- sapply(fbootrf, \(x) predict(x, data=data.frame("X"=newX))$prediction)

# residuals:
predsGLM$fit - truth
predsranger$predictions - truth  # jack
predsqrf$predictions[,2] - truth # quantile
rowMeans(predsbootrf) - truth # boot

# RMSE
sqrt(mean( (predsGLM$fit - truth)^2 ))
sqrt(mean( (predsranger$predictions - truth)^2 ))
sqrt(mean( (predsqrf$predictions[,2] - truth)^2 )) # median prediction
sqrt(mean( (rowMeans(predsbootrf) - truth)^2 )) # boot



plot(newX, truth, type="p", xlim=c(0.25, 1.35), ylim=c(-4, -.5))
curve(1-4*x, add=T, col="grey")
points(newX, predsGLM$fit, col="darkgreen", pch=16)
points(newX, predsranger$predictions, col="red", pch=16)
points(newX, predsqrf$predictions[,2], col="blue", pch=16)
points(newX, rowMeans(predsbootrf), col="pink", pch=16) # same as ranger (as it should be)

# ProbIntTrans
pnorm(truth, mean=predsGLM$fit, sd=predsGLM$se.fit)
pnorm(truth, mean=predsranger$predictions, sd=predsranger$se)
ifelse(truth < predsqrf$predictions[,1] | truth > predsqrf$predictions[,3], "out", "in") # truth in quantiles?
pnorm(truth, mean=rowMeans(predsbootrf), sd=apply(predsbootrf, 1, sd)) # boot




# Y <- rbinom(100, size=1, prob=plogis(1 - 4*X))
# ## fit models:
# fglm <- glm(Y ~ X, family=binomial)
# franger <- ranger(Y ~ X, data=cbind.data.frame(X=X, Y=Y), keep.inbag=T, classification=T, probability=T)
# plot(X, Y)
# plot(X, franger$predictions[,1])
# 
# 
# newX <- c(0.5, 1.1)
# newX <- seq(0, 1.2, len=101)
# truth <- 1 - 4*newX
# # predict:
# predsGLM <- predict(fglm, newdata=data.frame(X=newX), se.fit=T)
# predsranger <- predict(franger, data=data.frame("X"=newX))
# points(newX, predsranger$predictions, col="red")
# 
# # relative error
# (predsGLM$fit - truth)/truth
# (qlogis(predsranger$predictions[,1]) - truth)/truth
# 
# # ProbIntTrans
# pnorm(truth, mean=predsGLM$fit, sd=predsGLM$se.fit)
# pnorm(truth, mean=predsGLM$fit, sd=predsGLM$se.fit)





#### Power approach paradox ####

#Example 1

effect1 <- 0.01
effect2 <- 1
var <- 1.1
effect3 <- 0.25
varlarge <- 2
varsmall <- 0.3


plot(effect1, 3.5, pch=19,ylim=c(-1,5), xlim=c(-3,3),axes=F, xlab="Effect size", ylab="")
points(effect2, 4, pch=19)
abline(v=0, lty=2, color="gray")
lines(c(effect1+var, effect1-var), c(3.5,3.5))
lines(c(effect2+var, effect2-var), c(4,4))

points(effect3, 2, pch=15)
points(effect3, 1.5, pch=15)
lines(c(effect3+varlarge, effect3-varlarge), c(2,2))
lines(c(effect3+varsmall, effect3-varsmall), c(1.5,1.5))

points(0,0, pch=18, cex=1.25)
points(0,-.5, pch=18, cex=1.25)
lines(c(varsmall, -1*varsmall), c(0,0))
lines(c(var, -1*var), c(-.5,-.5))

text(-2.5, 3.75, "Case I")
text(-2.5, 1.75, "Case II")
text(-2.5, -.25, "Case III")

text(2.5,4, "Exp. 1", cex=.75)
text(2.5,3.5, "Exp. 2", cex=.75)
text(2.5,2, "Exp. 1", cex=.75)
text(2.5,1.5, "Exp. 2", cex=.75)
text(2.5,0, "Exp. 1", cex=.75)
text(2.5,-.5, "Exp. 2", cex=.75)

axis(1, at=c(-2, -1, 0, 1, 2), lwd.ticks=0, lwd=0)


