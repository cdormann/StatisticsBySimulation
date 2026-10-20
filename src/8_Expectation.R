#### Plotting ####

sapply(c("brms", "cowplot", "DHARMa", "faraway", "future.apply", "ggfortify", "ggplot2", "MASS", "mgcv", "mvtnorm", "ranger"), library, character.only = TRUE)

# simu and make a complicated model:

# https://archive.ics.uci.edu/dataset/186/wine+quality
wine <- read.csv("../Data/wine+quality/winequality-red.csv", sep=";")
head(wine)
summary(wine)

#curve(3*x - 0.15*x^2, from=4.6, to=15.9)
set.seed(123)
wine$preference <- with(wine, -10 + 3*fixed.acidity - 0.15*fixed.acidity^2 + 0.4*fixed.acidity*residual.sugar  + 1*residual.sugar - 0.01*total.sulfur.dioxide + 1*alcohol + rnorm(nrow(wine), 0, 1))
range(wine$preference)

summary(fm <- lm(preference ~ fixed.acidity + I(fixed.acidity^2) + fixed.acidity:residual.sugar + residual.sugar + total.sulfur.dioxide + alcohol + chlorides + sulphates, data=wine))


##### conditional #####
newFA <- data.frame("fixed.acidity"=seq(min(wine$fixed.acidity), max(wine$fixed.acidity), len=100), t(apply(wine[,-1], 2, median)))
newFA$predsFA <- predict(fm, newdata=newFA)

#plot(preference ~ fixed.acidity, data=wine)
#lines(newFA$fixed.acidity, predsFA)

library(ggplot2)
ggplot(wine, aes(x=fixed.acidity, y=preference)) + geom_point() + 
  geom_line(data=newFA, aes(x=fixed.acidity, y=predsFA), linewidth=2, col="grey") +
  xlab(label="fixed acidity") + ylab("wine preferences") +
  ylim(10, 60) +
  theme_minimal()

##### marginal #####
# produce a data set, where the 100 evaluations are performed on every combination; so, for the combination of values in the first row, evaluate the 100 levels of "fixed.acidity", then repeat for each other row in the data set;
# predictions are then the average across all of these!
seqFA <- seq(min(wine$fixed.acidity), max(wine$fixed.acidity), len=100)
newFAmarg <- data.frame("fixed.acidity"=rep(seqFA, times=nrow(wine)), wine[rep(1:nrow(wine), each=100), 2:13] )
#head(newFAmarg)
preds <- predict(fm, newdata=newFAmarg)
predsMeans <- tapply(preds, newFAmarg$fixed.acidity, mean)
dats <- cbind.data.frame("fixed.acidity"=seqFA, "preference"=predsMeans)

pdf(file="../figures/plottingCondMarg.pdf", width=5, height=5)
ggplot(wine, aes(x=fixed.acidity, y=preference)) + geom_point(color="grey40") + 
  geom_line(data=dats, aes(x=fixed.acidity, y=preference), linewidth=1.2, col="grey5") + # marginal
  geom_line(data=newFA, aes(x=fixed.acidity, y=predsFA), linewidth=1.2, col="grey") + # conditional
  ylim(10, 60) +
  xlab(label="fixed acidity") + ylab("wine preferences") +
  annotate("text", 14.5, 25, label="conditional", col="grey50") +
  annotate("text", 15, 31, label="marginal", col="grey5") +
  theme_minimal()
dev.off()


 ##### partial dependence #####
# predict to data, but keep focal variable at its mean; thereby account for all others; plot observed - this prediction against original focal variable
wineFAmean <- wine
wineFAmean$fixed.acidity <- mean(wine$fixed.acidity)
wine$FApredsPartial <- predict(fm, newdata=wineFAmean) # add to orginal data!
wine$FAResidualsPartial <- wine$preference - wine$FApredsPartial

newFApdp <- data.frame("fixed.acidity"=seq(min(wine$fixed.acidity), max(wine$fixed.acidity), len=100), t(colMeans(wine[,-1])))
newFApdp$predsFA <- predict(fm, newdata=newFApdp) - mean(wine$preference)


plot(newDF$preference - predsPartial ~ fixed.acidity, data=wine) # what mgcv does, I think
lines(newFApdp$fixed.acidity, newFApdp$predsFA)

pdf(file="../figures/plottingPDP.pdf", width=5, height=5)
ggplot(wine, aes(x=fixed.acidity, y=FAResidualsPartial)) + geom_point(color="grey40") +
  xlab(label="fixed acidity") + ylab("partial effect on wine preferences") +
  geom_line(data=newFApdp, aes(x=fixed.acidity, y=predsFA), linewidth=1.2, col="grey") +
  ylim(-8, 8)+
  theme_minimal()
dev.off()

##### partial residual plot #####
wineFAprp <- data.frame("fixed.acidity"=wine$fixed.acidity, t(colMeans(wine[,-1])))
wineFAprp$partialRes <- predict(fm, newdata=wineFAprp) + residuals(fm)
newFAprp <- data.frame("fixed.acidity"=seq(min(wine$fixed.acidity), max(wine$fixed.acidity), len=100), t(colMeans(wine[,-1])))
newFAprp$predictedPreference <- predict(fm, newdata=newFAprp)
  
plot(wine$fixed.acidity, partialRes)
lines(newFA$fixed.acidity, predict(fm, newdata=newFA))

pdf(file="../figures/plottingPRP.pdf", width=5, height=5)
ggplot(wineFAprp, aes(x=fixed.acidity, y=partialRes)) + geom_point(color="grey40") + 
  geom_line(data=newFAprp, aes(x=fixed.acidity, y=predictedPreference), linewidth=1.2, col="grey") +
  xlab(label="fixed acidity") + ylab("partial effect on wine preferences") +
  theme_minimal()
dev.off()





#### Model diagnostics ####
# uses the "wine" model from the previous section"
library(ggfortify)
pdf(file="../figures/8_wineDiagnosticsLM.pdf")
autoplot(fm, smooth.colour="grey") + theme_minimal()
dev.off()

library(DHARMa)
pdf(file="../figures/8_wineDiagnosticsDHARMa.pdf", width=9, height=5.5)
plot(simulateResiduals(fm))
dev.off()


## or, using simple invented data:
set.seed(123)
x <- runif(100, 0, 1)
y <- -x^2 + x + rnorm(100, 0 , .05)
plot(x,y, las=1)

fm <- lm(y ~ x)
#summary(fm)

fm.fitted <- predict(fm)
fm.residuals <- y - fm.fitted
# library(DHARMa)
# plot(simulateResiduals(fm))

library(ggplot2)
pdf(file="../figures/diagnosticsLM.pdf", width=5, height=5)
ggplot(data.frame("fitted"=fm.fitted, "residuals"=fm.residuals), aes(fitted, residuals)) + geom_point() + 
  geom_hline(yintercept=0, col="grey") +
  theme_minimal()
dev.off()

##### simulate from model #####
sim1 <- predict(fm) + rnorm(100, mean=0, sd=sd(residuals(fm)))
#plot(sim1, y)
#abline(0,1)
sims <- replicate(200, {predict(fm) + rnorm(100, mean=0, sd=sd(residuals(fm)))})

##### compute quantile of observation #####
ecdf(sims[1,])(y[1])
plot(density(sims[1,]))
abline(v=y[1])

pdf(file="../figures/diagnosticsDens.pdf", width=5, height=5)
ggplot(data.frame("dens"=sims[1, ]), aes(x=dens)) + geom_density() + 
  geom_vline(xintercept=y[1], color="grey") + 
  xlab(label="response value") +
  theme_minimal()
dev.off()

pdf(file="../figures/diagnosticsECDF.pdf", width=5, height=5)
ggplot(data.frame("dens"=sims[1, ]), aes(x=dens)) + 
#  geom_density() + 
  stat_ecdf(geom="line") + 
  geom_vline(xintercept=y[1], color="grey") + 
  xlab(label="response value") + ylab("cumulative density") +
  theme_minimal()
dev.off()


# for all 100 data points:
qs <- 0
for (i in 1:100) qs[i] <- ecdf(sims[i,])(y[i])
summary(qs)
plot(ecdf(qs))
curve(punif(x, 0, 1), add=T)

#hist(qs)
#plot(fm.fitted, qs)
pdf(file="../figures/diagnosticsLMsim.pdf", width=5, height=5)
ggplot(data.frame("fitted"=fm.fitted, "residuals"=qs), aes(fitted, residuals)) + geom_point() + 
  geom_hline(yintercept=0.5, col="grey") +
  theme_minimal()
dev.off()

cor(fm.residuals, qs)
cor(fm.residuals, qs, method="sp")
plot(fm.residuals, qs)
curve(pnorm(x, mean=0, sd=sd(residuals(fm))), add=T)

pdf(file="../figures/diagnosticsResidsOLSvssim.pdf", width=5, height=5)
ggplot(data.frame("OLSresid"=fm.residuals, "quantresid"=qs), aes(OLSresid, quantresid)) + 
  geom_point() + xlab("OLS residuals") + ylab("simulated residuals") +
  stat_function(fun=pnorm, geom="line", args=list(mean=0, sd=sd(fm.residuals)), col="grey") + 
  theme_minimal()
dev.off()

plot((1:100)/(100 + 1), sort(qs))
dats <- data.frame("expected"=(1:100)/(100+1), "observed"=sort(qs))
pdf(file="../figures/diagnosticsQQ.pdf", width=5, height=5)
ggplot(dats, aes(observed, expected)) + 
  stat_function(fun=punif, geom="line", args=list(0, 1), col="grey") + 
  geom_point() +  theme_minimal()
dev.off()

##### diagnostics a la DHARMa #####

ks.test(qs, "punif", 0, 1) # rather uniform

#var(fm.residuals) / mean(apply(sims, 2, var))
summary(apply(sims, 2, function(x) var(fm.residuals) / var(x)))

range(qs)
sum(qs==0) + sum(qs==1)


plot(fm.fitted, qs)

pdf(file="../figures/diagnosticsQuantResids.pdf", width=5, height=5)
ggplot(data.frame("fitted"=fm.fitted, "quantresid"=qs), aes(fitted, quantresid)) + 
  xlab("rank-transformed fitted values") + ylab("simulated residuals") +
  geom_smooth(method="gam", color="grey40") + geom_point() + 
  theme_minimal()
dev.off()


##### random forest simulations #####
BCWD <- read.csv(file="../Data/breast+cancer+wisconsin+diagnostic/BCWD.csv", stringsAsFactors = T, row.names = 1)
summary(BCWD)
BCWD$Diagnosis <- ifelse(BCWD$Diagnosis == "B", 0, 1) # recode as 0 (benign) and 1 (malign)
library(ranger)
fr <- ranger(factor(Diagnosis) ~., data=BCWD, probability=T)
fr.fitted <- predict(fr, data=BCWD)$prediction[,2] # prob of y=1
n <- nrow(BCWD)

set.seed(123)
sim1 <- rbinom(nrow(BCWD), p=fr.fitted, size=1)
sims <- replicate(200, rbinom(n, p=fr.fitted, size=1))

# add normal noise to simulations and observations when computing quantile residuals:
qr <- 0
for (i in 1:n) qr[i] <- ecdf(sims[i,] + runif(200, -0.5, 0.5))(BCWD$Diagnosis[i] + runif(1, -0.5, 0.5))
summary(qr)


dats <- data.frame("expected"=(1:n)/(n+1), "observed"=sort(qr))

pdf(file="../figures/diagnosticsQQbinary.pdf", width=5, height=5)
ggplot(dats, aes(observed, expected)) + 
  stat_function(fun=punif, geom="line", args=list(0, 1), col="grey") + 
  geom_point() +  theme_minimal()
dev.off()

ks.test(qr, "punif", 0, 1) # rather uniform

range(qr)
sum(qr==0) + sum(qr==1)

pdf(file="../figures/diagnosticsQuantResidsbinary.pdf", width=5, height=5)
ggplot(data.frame("fitted"=rank(fr.fitted), "quantresid"=qr), aes(fitted, quantresid)) + 
  xlab("fitted values") + ylab("simulated residuals") +
  geom_smooth(method="gam", color="grey40") + geom_point() + 
  theme_minimal()
dev.off()


# rewrite plotQQunif  to use grey colours:
plotQQunif <- function (simulationOutput, testUniformity = T, testOutliers = T, 
          testDispersion = T, ...){
  a <- list(...)
  a$pch = DHARMa:::checkDots("pch", 2, ...)
  a$bty = DHARMa:::checkDots("bty", "n", ...)
  a$logscale = DHARMa:::checkDots("logscale", F, ...)
  a$col = DHARMa:::checkDots("col", "black", ...)
  a$lcol="grey"
  a$main = DHARMa:::checkDots("main", "QQ plot residuals", ...)
  a$cex.main = DHARMa:::checkDots("cex.main", 1, ...)
  a$xlim = DHARMa:::checkDots("xlim", c(0, 1), ...)
  a$ylim = DHARMa:::checkDots("ylim", c(0, 1), ...)
  simulationOutput = DHARMa:::ensureDHARMa(simulationOutput, convert = "Model")
  do.call(gap::qqunif, append(list(simulationOutput$scaledResiduals), 
                              a))
  if (testUniformity == TRUE) {
    temp = DHARMa:::testUniformity(simulationOutput, plot = F)
    legend("topleft", c(paste("KS test: p=", round(temp$p.value, 
                     digits = 5)), paste("Deviation ", ifelse(temp$p.value < 0.05, "significant", "n.s."))), 
           text.col = ifelse(temp$p.value < 0.05, "grey", "black"), bty = "n")
  }
  if (testOutliers == TRUE) {
    temp = DHARMa:::testOutliers(simulationOutput, plot = F)
    legend("bottomright", c(paste("Outlier test: p=", round(temp$p.value, 
                                                            digits = 5)), paste("Deviation ", ifelse(temp$p.value < 0.05, "significant", "n.s."))), text.col = ifelse(temp$p.value < 0.05, "grey", "black"), bty = "n")
  }
  if (testDispersion == TRUE) {
    temp = DHARMa:::testDispersion(simulationOutput, plot = F)
    legend("center", c(paste("Dispersion test: p=", round(temp$p.value, 
                                                          digits = 5)), 
                       paste("Deviation ", ifelse(temp$p.value < 0.05, "significant", "n.s."))), text.col = ifelse(temp$p.value < 0.05, "grey", "black"), bty = "n")
  }
}

# rewrite plotResiduals to use grey colours:
plotResiduals <- function (simulationOutput, form = NULL, quantreg = NULL, rank = T, 
          asFactor = NULL, smoothScatter = NULL, quantiles = c(0.25, 0.5, 0.75), ...){
  a <- list(...)
  a$ylab = DHARMa:::checkDots("ylab", "DHARMa residual", ...)
  a$xlab = DHARMa:::checkDots("xlab", ifelse(is.null(form), "Model predictions", 
                                    gsub(".*[$]", "", deparse(substitute(form)))), ...)
  if (rank == T)  a$xlab = paste(a$xlab, "(rank transformed)")
  simulationOutput = DHARMa:::ensureDHARMa(simulationOutput, convert = T)
  res = simulationOutput$scaledResiduals
  if (inherits(form, "DHARMa")) stop("DHARMa::plotResiduals > argument form cannot be of class DHARMa. Note that the syntax of plotResiduals has changed since DHARMa 0.3.0. See ?plotResiduals.")
  pred = DHARMa:::ensurePredictor(simulationOutput, form)
  if (!is.factor(pred)) {
    if (rank == T) {
      pred = rank(pred, ties.method = "average")
      pred = pred/max(pred)
      a$xlim = DHARMa:::checkDots("xlim", c(0, 1), ...)
    }
    nuniq = length(unique(pred))
    ndata = length(pred)
    if (is.null(asFactor)) asFactor = (nuniq == 1) | (nuniq < 10 & ndata/nuniq > 10)
    if (asFactor) pred = factor(pred)
  }
  if (is.null(quantreg)) 
    if (length(res) > 2000) quantreg = FALSE
  else quantreg = TRUE
  switchScatter = 10000
  if (is.null(smoothScatter)) 
    if (length(res) > switchScatter)  smoothScatter = TRUE
  else smoothScatter = FALSE
  blackcol = rgb(0, 0, 0, alpha = max(0.1, 1 - 3 * length(res)/switchScatter))
  if (is.factor(pred)) {
    DHARMa:::testCategorical(simulationOutput = simulationOutput, 
                    catPred = pred, quantiles = quantiles)
  } else if (smoothScatter == TRUE) {
    defaultCol = ifelse(res == 0 | res == 1, 2, blackcol)
    do.call(graphics::smoothScatter, append(list(x = pred, 
                                                 y = res, ylim = c(0, 1), axes = FALSE, 
                                                 colramp = colorRampPalette(c("white", "darkgrey"))), a))
    points(pred[defaultCol == 2], res[defaultCol == 2], col = "grey", cex = 0.5)
    axis(1)
    axis(2, at = c(0, quantiles, 1))
  } else {
    defaultCol = ifelse(res == 0 | res == 1, 2, blackcol)
    defaultPch = ifelse(res == 0 | res == 1, 8, 1)
    a$col = DHARMa:::checkDots("col", defaultCol, ...)
    a$pch = DHARMa:::checkDots("pch", defaultPch, ...)
    do.call(plot, append(list(res ~ pred, ylim = c(0, 1), axes = FALSE), a))
    axis(1)
    axis(2, at = c(0, quantiles, 1))
  }
  main = DHARMa:::checkDots("main", ifelse(is.null(form), "Residual vs. predicted", 
                                  "Residual vs. predictor"), ...)
  out = NULL
  if (is.numeric(pred)) {
    if (quantreg == F) {
      title(main = main, cex.main = 1)
      abline(h = quantiles, col = "black", lwd = 0.5, lty = 2)
      try({
        lines(smooth.spline(pred, res, df = 10), lty = 2, 
              lwd = 2, col = "grey")
        abline(h = 0.5, col = "grey", lwd = 2)
      }, silent = T)
    }
    else {
      out = DHARMa:::testQuantiles(simulationOutput, pred, quantiles = quantiles, 
                          plot = F)
      if (any(out$pvals < 0.05, na.rm = TRUE)) {
        main = paste(main, "Quantile deviations detected (grey curves)", 
                     sep = "\n")
        if (out$p.value <= 0.05) {
          main = paste(main, "Combined adjusted quantile test significant", 
                       sep = "\n")
        }
        else {
          main = paste(main, "Combined adjusted quantile test n.s.", 
                       sep = "\n")
        }
        maincol = "grey"
      }
      else {
        main = paste(main, "No significant problems detected", 
                     sep = "\n")
        maincol = "black"
      }
      title(main = main, cex.main = 0.8, col.main = maincol)
      for (i in 1:length(quantiles)) {
        lineCol = ifelse(out$pvals[i] <= 0.05 & !(is.na(out$pvals[i])), 
                         "grey", "black")
        filCol = ifelse(out$pvals[i] <= 0.05 & !(is.na(out$pvals[i])), 
                        "grey", "#00000020")
        abline(h = quantiles[i], col = lineCol, lwd = 0.5, 
               lty = 2)
        polygon(c(out$predictions$pred, rev(out$predictions$pred)), 
                c(out$predictions[, 2 * i] - out$predictions[, 
                                                             2 * i + 1], rev(out$predictions[, 2 * i] + 
                                                                               out$predictions[, 2 * i + 1])), col = "#00000020", 
                border = F)
        lines(out$predictions$pred, out$predictions[, 
                                                    2 * i], col = lineCol, lwd = 2)
      }
    }
  }
  invisible(out)
}

pdf(file="../figures/nhtempDharmaResiduals.pdf", width=8.5, height=5)
par(xpd=F, mfrow=c(1,2))
plotQQunif(createDHARMa(simulatedResponse=sims, observedResponse = BCWD$Diagnosis, fittedPredictedResponse = fr.fitted, integerResponse = T))
plotResiduals(createDHARMa(simulatedResponse=sims, observedResponse = BCWD$Diagnosis, fittedPredictedResponse = fr.fitted, integerResponse = T))
dev.off()


#### Confidence and Prediction Intervals #### 
# # data from here: https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/ 
# temp <- read.table("../Data/NOAAweatherSpokane/CRNH0203-2023-WA_Spokane_17_SSW.txt", na.strings = "-9999.0")
# summary(temp)
# temp <- temp[, c(4, 5, 10)] # 4: local day, 5: local time, 10: mean hourly temperature
# colnames(temp) <- c("DOY", "TOD", "T") 
# temp$TOD <- sprintf("%04d", temp$TOD)
# temp$TODhours <- as.numeric(temp$TOD)/100
# temp$datetime <- as.POSIXct(paste(temp$DOY, temp$TOD), format="%Y%m%d %H%M")
# #head(temp)
# #which(temp$DOY==20230630)
# spring <- temp[2160:4339,]
# tail(spring)
# plot(T ~ datetime, data=spring, type="l")
# 
# library(mgcv)
# fm <- gam(T ~ s(DOY) + s(TODhours, bs="cc", k=24), data=spring)
# # plot(fm) # looks good
# predict(fm, newdata=data.frame("DOY"=20230630, "TODhours"=12), se.fit=T)

data(trees)  
fm <- lm(Volume ~ Girth, data=trees)
summary(fm)

##### confidence interval #####
# predict to 20 inch tree:
predsLM <- predict(fm, newdata=list("Girth"=20), interval="confidence")

# MC from parameters:
library(mvtnorm)
Sig <- vcov(fm)
b <- rmvnorm(1000, mean=coef(fm), sigma=Sig)
predsMC <- c(1, 20) %*% t(b)
mean(predsMC)
quantile(predsMC, c(0.025, 0.975))

# Analytically: CI from equation 8.1
sigma <- summary(fm)$sigma
X <- cbind(1, trees$Girth) # include intercept!
Z <- c(1, 20)

# mean plus/minus this:
(mean.analytical <- c(1, 20) %*% coef(fm))
SECI.analytical <- qt(0.025, df=fm$df.residual) * as.numeric(sqrt(sigma^2 * t(Z) %*% solve(t(X) %*% X) %*% Z))

rbind("LM"=as.vector(predsLM), 
      "Analytical"=mean.analytical + c(0, 1, -1) * SECI.analytical,
      "MC"=c(mean(predsMC), quantile(predsMC, c(0.025, 0.975)))
      ) 


##### prediction interval #####
predsLM.PI <- predict(fm, newdata=list("Girth"=20), interval="prediction")
predsMC.PI <- rnorm(1000, mean=predsMC, sd(residuals(fm)))

##### analytical: PI from eqn. 8.2 #####
SEPI.analytical <- qt(0.025, df=fm$df.residual) * as.numeric(sqrt(sigma^2 * (1 + t(Z) %*% solve(t(X) %*% X) %*% Z)))

rbind("LM"=as.vector(predsLM.PI), 
      "Analytical"= c(mean.analytical + c(0, 1, -1) * SEPI.analytical),
      "MC"=c(mean(predsMC.PI), quantile(predsMC.PI, c(0.025, 0.975)))
)
 


library(faraway)
data(africa)
library(MASS)
fmnb <- glm.nb(miltcoup ~ oligarchy, data=africa)
summary(fmnb)
pred.fmnb <- predict(fmnb, newdata=list("oligarchy"=12), se.fit=T)
# confidence interval:
exp(c(pred.fmnb$fit-2*pred.fmnb$se.fit, pred.fmnb$fit+2*pred.fmnb$se.fit))
# or, by foot:
Sig.nb <- vcov(fmnb)  
set.seed(123)
b.nb <- rmvnorm(1000, mean=coef(fmnb), sigma=Sig.nb)
predsMC.nb <- c(1, 12) %*% t(b.nb) # link-scale MC predictions
# exp(predsMC.nb) # response-scale
quantile(exp(predsMC.nb), c(0.025, 0.975))

# prediction interval:
realizations.nb <- rnbinom(1000, mu=exp(predsMC.nb), size=fmnb$theta)
mean(realizations.nb); quantile(realizations.nb, c(0.025, 0.975))
table(realizations.nb)/1000
# compared to distribution without parameter uncertainty:
round(dnbinom(0:14, mu=exp(pred.fmnb$fit), size=fmnb$theta), 3)

library(ggplot2)
dats <- data.frame("count"=c(as.numeric(names(table(realizations.nb))), 0:14), 
                   "type"=rep(c("included", "excluded"), times=c(13, 15)), 
                   "probability"=c(table(realizations.nb)/1000, dnbinom(0:14, mu=exp(pred.fmnb$fit), size=fmnb$theta))) 
pdf(file="../figures/predictionIntervalnegbin.pdf", width=5, height=5)
ggplot(dats, aes(x=count, y=probability, fill=type)) + 
  geom_col(position=position_dodge()) + 
  scale_fill_manual("parameter uncertainty", values=c("grey20", "grey70")) + 
  theme_minimal() + theme(legend.position = c(0.8, 0.9))
dev.off()



# Bernoulli example
library(faraway)
data(pima)
library(mgcv)
fgam <- gam(test ~ s(insulin), data=pima, family=binomial)
summary(fgam)
#plot(fgam)

(pred.gam.link.scale <- predict(fgam, newdata=list(insulin=500), se.fit=T))
# confidence:
plogis(c(pred.gam.link.scale$fit-2*pred.gam.link.scale$se.fit, pred.gam.link.scale$fit+2*pred.gam.link.scale$se.fit))
set.seed(123)
pred.gam.response.scale <- plogis(rnorm(1000, pred.gam.link.scale$fit, sd=pred.gam.link.scale$se.fit))
# prediction:
table(rbinom(1000, prob=pred.gam.response.scale, size=1))


##### Confidence and prediction for machine learning ####

library(faraway)
data(wcgs)
#summary(wcgs)
wcgs <- as.data.frame(scale(na.omit(wcgs[, -c(7, 9:13)])))
#summary(wcgs)
# 7 continuous, 5 categorical
MM <- model.matrix(~(age+height+weight+sdp+dbp+chol+cigs)^2, data=wcgs)
b <- rep(c(1, -1), len=29)
truth <- MM %*% b 
set.seed(1234)
wcgs$y <- as.vector(truth + rnorm(nrow(MM), 0, 2))
summary(fm <- lm(y ~ (age+height+weight+sdp+dbp+chol+cigs)^2, data=wcgs))
predict(fm, newdata=wcgs[1,], se.fit=T)

library(ranger)
fr <- ranger(y ~ ., data=wcgs, keep.inbag=T, quantreg=T)
fr
str(predict(fr, data=wcgs[1,], type="se", se.method="infjack")) # default
str(predict(fr, data=wcgs[1,], type="se", se.method="jack"))

predsse <- predict(fr, data=wcgs[1,], type="se", se.method="infjack")
c(predsse$prediction-1.96*predsse$se, predsse$prediction+1.96*predsse$se)

system.time(frbs <- replicate(500, {
  ranger(y ~ ., data=wcgs[sample(nrow(wcgs), replace=T),])
  }, simplify=F)) # 400s for 500 runs

library(future.apply) # parallelize replicate!
# length(availableWorkers())
plan(multisession, workers = 8)
system.time(frbs <- future_replicate(100, {
  list(ranger(y ~ ., data=wcgs[sample(nrow(wcgs), replace=T),]))
})) # 160s
# need to use list() because future_replicate automatically simplifies!
#frbs[[1]]

predsbs <- sapply(frbs, \(x) predict(x, data=wcgs[1,])$prediction)
mean(predsbs); sd(predsbs)
quantile(predsbs, c(0.025, 0.975))

# quantile prediction with ranger
predict(fr, data=wcgs[1,], type="quantiles", quantiles=c(0.025, 0.975))$prediction
str(predsq <- predict(fr, data=wcgs[1,], type="quantiles", quantiles=seq(0.01, 0.99, len=100)))

# visualise for one data point prediction:
(MM %*% b)[1,] # mean, truth
curve(dnorm(x, mean=1.933, sd=1.987), from=-3, to=8, lwd=3, ylim=c(0, .65)) # fit LM
curve(dnorm(x, mean=1.77, sd=2), lwd=2, ylim=c(0, .65), lty=2, add=T) # truth
abline(v= (MM %*% b)[1,], lty=2)
curve(dnorm(x, mean=1.32, sd=0.833), add=T, col="grey", lwd=3) # RFse
lines(density(predsbs), col="green", lwd=3, xpd=T) # RFbs
lines(predsq$predictions, pmin(seq(0.01, 0.99, len=100), 1-seq(0.01, 0.99, len=100)), lwd=3, col="red") # RFq
dats2 <- as.data.frame(density(predsbs)[1:2])
#plot(1:1000, pmin(seq(0.025, 0.975, len=1000), 1-seq(0.025, 0.975, len=1000)))
#plot(1:1000, predsq$predictions)

dats <- data.frame("X"=as.vector(predsq$predictions), "Y"=pmin(seq(0.01, 0.99, len=100), 1-seq(0.01, 0.99, len=100)), "boots"=predsbs)
library(ggplot2)
pdf(file="../figures/CImachlearn_example.pdf", width=6, height=5)
ggplot(dats, aes(x=X, y=Y)) + xlim(-3, 8) + 
  geom_function(fun = dnorm, args=list(mean=1.77, sd=2), linewidth=1.5, linetype=1) + 
  annotate("text", 0.3, 0.18, label="truth") +
#  geom_vline(xintercept=1.77) +
  geom_function(fun = dnorm, args=list(mean=1.933, sd=1.987), col="gray55", linewidth=1.5) + 
  annotate("text", .8, 0.11, label="linear model") +
  geom_function(fun = dnorm, args=list(mean=0.95, sd=0.723), col="gray65", linewidth=1.5) + 
  annotate("text", 2.8, .45, label="infinite jackknife") +
  geom_line(linewidth=1.5, col="gray75") +
  annotate("text", 4.5, 0.25, label="quantile regression") +
  geom_density(data=dats2, aes(x, y), stat="identity", linewidth=1.5, col="gray85") +
  annotate("text", 4.2, 0.02, label="bootstrapping") +
  ylab(label="density") + xlab(label="predicted value") +
  theme_minimal()
dev.off()


# count times the truth is in the 95% CI:
inCV <- data.frame("LM"=NA, "RFse"=NA, "RFbs"=NA, "RFqu"=NA)
# LM:
predsLM <- predict(fm, newdata=wcgs, se.fit=T)
CI_lm <- cbind(predsLM$fit -1.96*predsLM$se.fit, predsLM$fit + 1.96*predsLM$se.fit)
inCV[, "LM"] <- sum(truth > CI_lm[,1] & truth < CI_lm[,2])
# RFse:
predsRFse <- predict(fr, data=wcgs, type="se")
CI_RFse <- cbind(predsRFse$predictions - 1.96*predsRFse$se, predsRFse$predictions + 1.96*predsRFse$se)
inCV[, "RFse"] <- sum(truth > CI_RFse[,1] & truth < CI_RFse[,2])
# RFbs:
predsbs <- sapply(frbs, \(x) predict(x, data=wcgs)$prediction)
CI_RFbs <- t(apply(predsbs, 1, quantile,  c(0.025, 0.975)))
inCV[, "RFbs"] <- sum(truth > CI_RFbs[,1] & truth < CI_RFbs[,2])
# RFqu:
predsRFqu <- predict(fr, data=wcgs, type="quantiles", quantiles=c(0.025, 0.975))
CI_RFqu <- predsRFqu$predictions
inCV[, "RFqu"] <- sum(truth > CI_RFqu[,1] & truth < CI_RFqu[,2])
# get average coverage:
inCV/nrow(wcgs)
#        LM      RFse      RFbs      RFqu
# 0.9506684 0.9042011 0.6791852 0.9977721


##### repeat simulation #####
# takes a while due to bootstrapping; parallel over replicates
# 2.5 hours on my (fast) laptop

# library(future.apply) # parallelize replicate!
# plan(multisession, workers = 8)

reps <- 100
set.seed(123)
before <- Sys.time()
inCVs <- future_replicate(reps, { 
  inCV <- c("LM"=1, "RFse"=NA, "RFbs"=NA, "RFqu"=NA)
  b <- runif(29, -2, 2) * rep(c(1, -1), len=29)
  truth <- MM %*% b 
  wcgs$y <- as.vector(truth + rnorm(nrow(MM), 0, 2))
  fm <- lm(y ~ (age+height+weight+sdp+dbp+chol+cigs)^2, data=wcgs)
  fr <- ranger(y ~ ., data=wcgs, keep.inbag=T, quantreg=T)
  frbs <- replicate(200, {
    ranger(y ~ ., data=wcgs[sample(nrow(wcgs), replace=T),])  
  }, simplify=FALSE)
  # extract CIs:
  # LM:
  predsLM <- predict(fm, newdata=wcgs, se.fit=T)
  CI_lm <- cbind(predsLM$fit -1.96*predsLM$se.fit, predsLM$fit + 1.96*predsLM$se.fit)
  inCV["LM"] <- sum(truth > CI_lm[,1] & truth < CI_lm[,2])
  # RFse:
  predsRFse <- predict(fr, data=wcgs, type="se")
  CI_RFse <- cbind(predsRFse$predictions - 1.96*predsRFse$se, predsRFse$predictions + 1.96*predsRFse$se)
  inCV["RFse"] <- sum(truth > CI_RFse[,1] & truth < CI_RFse[,2])
  # RFbs:
  predsbs <- sapply(frbs, \(x) predict(x, data=wcgs)$prediction)
  CI_RFbs <- t(apply(predsbs, 1, quantile,  c(0.025, 0.975)))
  inCV["RFbs"] <- sum(truth > CI_RFbs[,1] & truth < CI_RFbs[,2])
  # RFqu:
  predsRFqu <- predict(fr, data=wcgs, type="quantiles", quantiles=c(0.025, 0.975))
  CI_RFqu <- predsRFqu$predictions
  inCV["RFqu"] <- sum(truth > CI_RFqu[,1] & truth < CI_RFqu[,2])
  inCV
})
after <- Sys.time()
after-before
# save(inCVs, file="inCV.Rdata")
load(file="inCV.Rdata")

apply(inCVs/nrow(wcgs), 1, summary)

dats <- stack(as.data.frame(t(inCVs)/nrow(wcgs)))
library(ggplot2)
pdf(file="../figures/CImachlearn_coverage.pdf", width=6, height=5)
ggplot(dats, aes(ind, values)) + 
  geom_hline(yintercept=0.95, linewidth=1.5, col="grey") + 
  geom_boxplot(fill="grey") + 
  ylab(label="confidence interval coverage (nominal: 0.95)") + 
  xlab(label="") +
  theme_minimal()
dev.off()



#### Alice's Restaurant ####

restaurant <- data.frame("Name"=c("Blue Lagoon", "Joe's Finger Food", "Beanache", "Real Saigon", "Tête a Tête"),
                         "Formality"=c(5, 3, 7, 7, 9),
                         "Enjoyment"=c(3, 6, 4, 6, 8))

plot(Enjoyment ~ Formality, data=restaurant, pch=16)

# library(brms)
# fm <- brm(Enjoyment ~ Formality, data=restaurant)

fm <- lm(Enjoyment ~ Formality, data=restaurant)
summary(fm)
cov2cor(vcov(fm))

library(mvtnorm)
set.seed(1234)
priorSimus <- rmvnorm(1000, mean=coef(fm), sigma=vcov(fm))

plot(Enjoyment ~ Formality, data=restaurant, pch=16, las=1, xlim=c(0,10), ylim=c(0,10))
for (i in 1:20){ curve(priorSimus[i,1] + priorSimus[i,2] * x, add=T, col="grey") }

# new impressions:
newRestaurant <- data.frame("Name"=c("Paradise", "Backstreet"),
                            "Formality"=c(2, 5),
                            "Enjoyment"=c(2, 7)) # yes, Alice likes it rough

priorExpectation <- sapply(1:1000, \(x) priorSimus[x,1] + priorSimus[x,2] * newRestaurant$Formality)

par(mfrow=c(1,2))
plot(density(priorExpectation[1,]))
abline(v=newRestaurant$Enjoyment[1], lty=2)
plot(density(priorExpectation[2,]))
abline(v=newRestaurant$Enjoyment[2], lty=2)

# now compute posterior!
library(brms)
fmbrms <- brm(Enjoyment ~ Formality, data=restaurant)
fixef(fmbrms)
coef(fm)
# from stanvar example:
bprior <- prior(multi_normal(M, V), class = "b")
stanvars <- stanvar(unname(coef(fm)), "M", scode = "  vector[2] M;") +
  stanvar(unname(vcov(fm)), "V", scode = "  matrix[2, 2] V;")
make_stancode(Enjoyment ~ Formality, data=newRestaurant, prior=bprior, stanvars = stanvars)
fit <- brm(Enjoyment ~ 0 + Intercept + Formality, data=newRestaurant, prior=bprior, stanvars = stanvars) # hack; see ?prior -> Details 1.
fit

# visualise prior and posterior; re-fit one data point at a time to show which data point let to most learning?
posts <- as.data.frame(posterior_predict(fit))
colnames(posts) <- c("Paradise", "Backstreet")
posts.stack <- stack(posts)
colnames(posts.stack) <- c("draw", "Restaurant")
# simulate from prior:
priorPred <- predict(fit, sample_prior="only")



library(ggplot2)
library(cowplot) # for two-panel plot

#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colours for colour-blinds

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") # colours for colour-blinds



pdf(file="../figures/restaurantsPriorPosterior.pdf", width=10, height=5)
p1 <- ggplot(as.data.frame(posts)[,1, drop=F], aes(x=Paradise)) + 
  xlim(0, 10) +
  geom_density(linewidth=1, color="black") +
  geom_vline(xintercept=newRestaurant$Enjoyment[1], color="grey80", linetype="dashed", size=1) +
  stat_function(fun=function(x) dnorm(x, mean=priorPred[1,1], sd=priorPred[1,2]), color="grey60", size=1) +
  ylab("enjoyment") +
  theme_minimal()

p2 <- ggplot(as.data.frame(posts)[,2, drop=F], aes(x=Backstreet)) + 
  xlim(0, 10) +
  geom_density(linewidth=1, color="black") +
  geom_vline(xintercept=newRestaurant$Enjoyment[2], color="grey", linetype="dashed", size=1) +
  stat_function(fun=function(x) dnorm(x, mean=priorPred[2,1], sd=priorPred[2,2]), color="grey60", size=1) +
  ylab("enjoyment") +
  theme_minimal()
plot_grid(p1, p2, align = "h", ncol = 2, rel_heights = c(.5, .5))
dev.off()

# plot regression before and after:
restaurants <- rbind(restaurant, newRestaurant)
restaurants$time <- rep(c("before", "after"), times=c(5, 2))

pdf(file="../figures/restaurant.pdf", width=5, height=5)
ggplot(restaurants, aes(x=Formality, y=Enjoyment, color=time)) +
#  scale_colour_manual(values=cbPalette[-1]) +
  scale_colour_manual(values=c("black", "grey60")) +
  geom_smooth(data=restaurant, method = "lm", fullrange = F, color="grey60", fill="grey80") +
  geom_smooth(method = "lm", fullrange = F, color="black", fill="grey80") +
  scale_x_continuous(limits = c(1, 10.5), expand = c(0,0)) +
  annotate("text", x=8.2, y=2.75, label= "before", color = "grey60", size=5) +
  annotate("text", x=2.4, y=3, label= "after", color = "black", size=5) +
  geom_point(size = 2) +
  ylab("enjoyment") +
  xlab("formality") +
  theme_minimal() +
  theme(legend.position = "none") # no legend, needs to go after theme_!
dev.off()

#And create an example residual plot
fm <- lm(Enjoyment ~ Formality, data=restaurants)
plot(fm, panel = function(x, y) panel.smooth(x, y, col.smooth="black"))




#### fore & hind casting ####
library(faraway)
data(globwarm)

library(mgcv)
globwarm$year <- as.numeric(rownames(globwarm))
train <- globwarm[globwarm$year > 1880 & globwarm$year < 1975, ]
before <- globwarm[globwarm$year <= 1880 & globwarm$year > 1855, ]
after <- globwarm[globwarm$year >= 1975, ]
fm <- gam(nhtemp ~ s(year), data=train)
summary(fm)

## model diagnostics (for later):
acf(residuals(fm))
library(DHARMa)
sim <- simulateResiduals(fm)


B <- 100
set.seed(123)
modellist <- replicate(B, simplify=F, # or use for-loop
                       {
                         fmbs <- gam(nhtemp ~ s(year), 
                                     data=train[sample(nrow(train), replace=T), ])
                       })

predsHind <- sapply(modellist,  predict, newdata=before, simplify=F)
predsFore <- sapply(modellist,  predict, newdata=after, simplify=F)

RMSEs <- sapply(predsHind, function(x) 
  sqrt(mean((before$nhtemp - x)^2)) )
mean(RMSEs)

pdf(file="../figures/nhtempHindFore.pdf", width=6, height=4)
par(mar=c(4,4,1,1))
plot(nhtemp ~ year, data=globwarm, type="l", las=1,
     ylab="annual mean temperature anomaly [°C]", lwd=2, ylim=c(-0.8, 0.8), xlim=c(1850, 2000))
lines(1881:1974, predict(fm), lwd=3, col="grey")
abline(v=c(1880, 1975))
abline(h=0, col="grey")
for (i in 1:B){
  lines(1856:1880, predsHind[[i]], col=rgb(.5, .5, .5, .5)) 
}
for (i in 1:B){
  lines(1975:2000, predsFore[[i]], col=rgb(.5, .5, .5, .5)) 
}
arrows(1875, 0.2, 1875, mean(RMSEs)+0.2, code=3, angle=90, length=.1)
dev.off()


## now re-analyse 10 years later and compare:
train2 <- globwarm[globwarm$year > 1880 & globwarm$year < 1985, ]
after2 <- globwarm[globwarm$year >= 1985, ]
fm2 <- gam(nhtemp ~ s(year), data=train2)

set.seed(123)
modellist2 <- replicate(B, simplify=F, # or use for-loop
                       {
                         fmbs <- gam(nhtemp ~ s(year), 
                                     data=train2[sample(nrow(train2), replace=T), ])
                       })
predsHind2 <- sapply(modellist2,  predict, newdata=before, simplify=F)
predsFore2 <- sapply(modellist2,  predict, newdata=after2, simplify=F)


pdf(file="../figures/nhtempFore10yearsLater.pdf", width=6, height=4)
par(mar=c(4,4,1,1))
plot(nhtemp ~ year, data=globwarm, type="l", las=1,
     ylab="annual mean temperature anomaly [°C]", lwd=2, ylim=c(-0.8, 0.8), xlim=c(1850, 2000))
lines(1881:1974, predict(fm), lwd=3, col="grey", lty=3)
lines(1881:1984, predict(fm2), lwd=3, col="grey")
abline(v=c(1880, 1975, 1985))
abline(h=0, col="grey")
# old
for (i in 1:B){
  lines(1856:1880, predsHind[[i]], col=rgb(.5, .5, .5, .5), lty=3) 
}
for (i in 1:B){
  lines(1975:2000, predsFore[[i]], col=rgb(.5, .5, .5, .5), lty=3) 
}

# new
for (i in 1:B){
  lines(1856:1880, predsHind2[[i]], col=rgb(.5, .5, .5, .5)) 
}
for (i in 1:B){
  lines(1985:2000, predsFore2[[i]], col=rgb(.5, .5, .5, .5)) 
}
#arrows(1875, 0.2, 1875, mean(RMSEs)+0.2, code=3, angle=90, length=.1)
dev.off()
