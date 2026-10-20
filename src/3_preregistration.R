#########################################################################################################
#Chapter 3
#########################################################################################################

#########################################################################################################
#load required libraries
library(ggplot2)
library(reshape2)
library(pwr)
#########################################################################################################

#########################################################################################################
#initial figure (blank figure, labelled axes only)

PerTruth <-seq(0,1,.01)
ITE <- PerTruth-PerTruth^2

df <- data.frame(PerTruth, ITE)

pdf("../figures/3_blankplot.pdf", width=5, height=5)

ggplot(df, aes(x=PerTruth, y=ITE)) +
geom_blank() +
coord_cartesian(xlim = c(0,1), ylim=c(0,.5))+
labs(x="Initial plausibility", y = "Change in plausibility following repetition\n(Illusory Truth Effect)") +
theme_bw()+
theme(axis.title=element_text(size=8))
 
dev.off()
#########################################################################################################

#########################################################################################################
#Reproduce simulation in Fazio et al. (2019)
#
#Simulations converted to R from Matlab code in preregistration: fazio_pennycoock_rand_code.m
#As preregistered on https://osf.io/w4k2c/
#
#Simulate four different possibilities: 
#	1. repetition effect constant and independent of baseline plausibility
#	2. repetition effect increases linearly with plausibility
#	3. repetition effect decreases linearly with plausibility 
#	4. repetition effect peaks at the midpoint of plausibility
#
#For three different degrees of "fluency" *f_i* - the subject's perceived accuracy of a novel statment
#
###Matlab note:
#Most Matlab code has similar syntax to R
#But one notable difference is in drawing random numbers.
#
#For "novel" subjects - those seeing statement for the first time (half of total subjects N) - 
#the (normal) random draw of perceived accuracy scores in Fazio et al. is a 
#N/2 rows x 1 column matrix with mean = plausibility *p_i* and sd fixed at 0.5
#Matlab code is (randn(N/2,1)*sd+i)
#corresponding R code is: rnorm(N/2, plausibility[i], sd)
#
#
#Code conversion, with embedded comments paralleling Fazio et al. 2019
#
#specify different relationships between plausibility and repetition (possibilities 1-4, above)
#  0:  constant repetition effect size
# -1: repetition effect size decreases with plausibility
# +1: repetition effect size increases with plausibility
# +2: repetition effect size is maximal at midpoint of plausibility

skewes <- c(0, -1, 1, 2) 	#vector of repetition effects

#Different sizes of the basic effect of repetition ("fluencies") to simulate

fluencies <- c(.1, .2, .4)	#vector of fluencies

#Define baseline plausibility ratings per Fazio et al.

plausibility <- seq(-1,2,.01)	#vector of plausibilities 

#set constants
N <- 100000 #number of subjects (Fazio et al. use 100,000). For this sample size, actual run time is < 30 seconds.
sd <- 0.5 #sd of normal distribution from which accuracy scores are sampled (Fazio et al call this "mu")

#define blank variables for storing output
scoresR <-  array(NA, c(length(plausibility), length(fluencies)+1, length(skewes)))
	dimnames(scoresR)[[3]] <- paste("RepEffect.",skewes, sep="")
	dimnames(scoresR)[[2]] <-c("Plausibility", paste("Fluencies.",fluencies,sep=""))
	scoresR[,1,]<- plausibility #first column is plausibility vector 
scoresN <-  array(NA, c(length(plausibility), length(fluencies)+1, length(skewes)))
	dimnames(scoresN)[[3]] <- paste("RepEffect.",skewes, sep="")
	dimnames(scoresN)[[2]] <-c("Plausibility", paste("Fluencies.",fluencies,sep=""))
	scoresN[,1,]<- plausibility #first column is plausibility vector
NFSample <- NULL 	#novel participants (N = 50,000 of them)
FSample <-NULL	#participants (N = 50,000) who see a statement a 2nd time

#Loop it over skewes, fluencies, and participants (this could be streamlined, for sure)

for (i in 1:length(skewes)) {
	for (j in 1:length(fluencies)) {
		for (k in 1:length(plausibility)) {
			if(skewes[i]==0) {	#this first ifelse series applies the functions
							#relating repetition effect size to plausibility
				fluencyEffect <- fluencies[j]
			}else if(skewes[i]==-1) {
				fluencyEffect <- fluencies[j] * (1 - plausibility[k])
			}else if(skewes[i]==1) {
				fluencyEffect <- fluencies[j] * plausibility[k]
			}else if(skewes[i]==2) {
				fluencyEffect <- fluencies[j] * 2 * (.5 - abs(.5 - plausibility[k]))
			}
			#Simulate 50,000 novel subjects:
				NFsample <- rnorm(N/2, plausibility[k], sd)
				NFsample <- NFsample > 0.5
				scoresN[k,j+1,i] <- mean(NFsample)
			#Simulate 50,000 repeated subjects:
				Fsample <- rnorm(N/2, (plausibility[k] + fluencyEffect), sd)
				Fsample <- Fsample > 0.5
				scoresR[k,j+1,i] <- mean(Fsample)
		}
	}
}

#Organize data to prepare for plotting

#For each of the four repetition effects (here "skewes", array index i) Fazio et al. plot illusory truth effect (scoresR[k]-scoresN[k]) vs. 
#perceived truth ((scoresR[k]+scoresN[k])/2) for each level of fluency (here each of the three columns after plausibility of array index j)

#So, we want a 3-D array with matrices grouped by the repetition effect, rows for each random replicate,
#and 8 coluns: repetition effect size [1], plausibility [2], 
#perceived truths [3-5] and illusory truth effects [6-8]

plotScores <- array(NA, c(length(plausibility), 2*length(fluencies)+2, length(skewes)))	#empty array
dimnames(plotScores)[[3]] <- paste("RepEffect.",skewes, sep="")
dimnames(plotScores)[[2]] <-c("Skew","Plausibility", paste("PerTruth.",fluencies,sep=""), paste("ITE.",fluencies,sep=""))
for (i in 1:length(skewes)) plotScores[,1,i]<-rep(skewes[i],length(plausibility)) #first column is repetition effect size
plotScores[,2,]<- plausibility #second column is plausibility vector but we don't use this in plotting, just for record-keeping and debugging. 

#Computations of perceived truth (x-axis: (R+N)/2) and illusory truth effect (y-axis: R-N):

for (i in 3:(length(fluencies)+2)) {
	plotScores[,i,] <- (scoresR[,i-1,]+scoresN[,i-1,])/2	#Perceived Truth: (R+N)/2
	plotScores[,(i+3),] <- scoresR[,i-1,]-scoresN[,i-1,]	#Illusory Truth Effect: R-N
}

#Now, we are going to want to plot, for each skewness (array index i), each ITE against its corresponding PerTruth
#ggplot doesn't work well with arrays (or I don't know how), so
#First, we convert the array into a long-form dataframe so that ggplot2 can deal with it effectively
#(Note: this could be indexed on array index i so that we could easily adjust for different numbers of repetition effect scenarios)

plotScoresLong <- data.frame(rbind(plotScores[,,1],plotScores[,,2],plotScores[,,3],plotScores[,,4]))

#A little more data wrangling to get this well-organized for ggplot

plotScoresLong$Skew <- as.factor(plotScoresLong$Skew)
PerTruths <- melt(plotScoresLong[,c(1,3:5)], id.vars=1, measure.vars=c(2:4))
PerTruths$variable <- as.factor(substr(PerTruths$variable,10,12))
names(PerTruths) <- c("Skew", "Fluency", "PerTruth")
ITE <- melt(plotScoresLong[,c(1,6:8)], id.vars=1, measure.vars=c(2:4))
ITE$variable <- as.factor(substr(ITE$variable,5,7))
names(ITE) <- c("Skew", "Fluency", "ITE")

PlotSims <-cbind(PerTruths, ITE)		#I used cbind() here because merge() and left_join() multiplied rows for no obvious reason
PlotSims<- PlotSims[,c(1:3,6)]		#Keep only the columns we need for plotting

PlotSims$Skew <- factor(PlotSims$Skew, levels=c("0", "1", "-1", "2"))	#Re-order factor levels to match conceptual order of repetition effects

#plots, finally
pdf("../figures/3_simplots.pdf", width=8, height=8)

ggplot(PlotSims, aes(x=PerTruth, y=ITE, color=Fluency, linetype=Fluency)) +
scale_colour_manual(values=c("black", "black", "black")) +
scale_linetype_manual(values=c(3,2,1)) +
stat_smooth(span=.1, se=FALSE) +
facet_wrap(vars(Skew), ncol=2,
		labeller = as_labeller(c("0"="Constant repetition effect",
                             "1"="Positive repetition effect",
                             "-1"="Negative repetition effect",
                             "2"="Midpoint repetition effect"))) +
labs(x="Perceived truth ([R+N]/2)", y="Illusory Truth Effect (R-N)") +
theme_bw()

dev.off()

########################################################################################################
########################################################################################################
#Sample size calculation for Fazio et al 2019
#value of d (Hedges' g) from Dechene et al. 2010 is 0.39, 95% ci = [.32, .47]
#for within-items, homogeneous context (as used by Fazio et al. 2019)
#assume alpha = 0.05, power = 0.08
#Fazio et al. use t-test (or equivalent)

pwr.t.test(n = NULL, d = .49, sig.level = 0.05, power = .8, type="paired")

pwr.t.test(n = 500, d= NULL, sig.level = 0.05, power = .8, type="paired")


#suggests n = 66 per group; Fazio used n = 500 participants, but evaluated
#relatively small number of statements (40 "true", 40 "false") per individual

###power for general linear model
## polynomial regression a + bx + cx^2 has 2 independent variables (x, x^2); u = 2
## for effects size f^2 "medium" (= .15) to correspond with medium e.s. of illusory truth effect,
## what is sample size needed for power = .8, alpha = .5? This is v + 2 + 1, where v is error df (N = u + v +1)

pwr.f2.test(u = 2, v = NULL, f2 = .15, sig.level = 0.05, power = .8)

## returns 64.3, so n = 80 is approximately ok.

#finally, we might ask what difference (effect size) could be detected from 0.5 (location of peak) for n = 80
#power = .8, alpha = .5

pwr.t.test(n=80, d=NULL, sig.level=.5, power=.8)

#returns 3 = 23%

#### code for figure 3.3 ####
library(pwr)
library(ggplot2)

## Two independent samples (power)
## Exercise 2.1 p. 40 from Cohen (1988)
d <- .3
a <- pwr.t.test(d=d,n=10,sig.level=c(.001,.005,.01,.05,.1),type="two.sample",alternative="two.sided")

## Two independent samples (sample size)
## Exercise 2.10 p. 59
pwr.t.test(d=0.3, power=0.75, sig.level=0.05, type="two.sample", alternative="greater")

df <- data.frame(matrix(ncol = 4, nrow = 0))
for (n in c(10, seq(100,1000,100))){
  a<-pwr.t.test(d=.2,n=n,sig.level=seq(.001,.1,.001),type="two.sample",alternative="two.sided")
  b<-cbind(a$d, a$n, a$sig.level, a$power)
  #print(b)
  df <- rbind(df,b)
}
names(df) <- c("d", "n", "alpha", "power")

pdf("../figures/3_analytical_power.pdf", width=5.5, height=5)

ggplot(df, aes(x=alpha, y=power, group=n)) + 
  geom_smooth(colour="black", se=FALSE) +
  theme_bw() +
  geom_hline(yintercept=0.8, linetype="dashed") +
  geom_vline(xintercept=0.05, linetype="dashed") +
  geom_text(label="N=10", x=.0975, y=.09) +
  geom_text(label="N=100", x=.0975, y=.36) +
  geom_text(label="N=400", x=.0975, y=.84) + 
  theme_minimal()
dev.off() 







#Results coded as proportion of statements rated as "true" across different "bins"
#of plausibility
#
#Data are in \Data\FazioExperiment\StatementEffects.csv 
#and at https://osf.io/9zjre (not a direct link) 

statement.effects <- read.csv("../Data/FazioExperiment/StatementEffects.csv")

pdf("../figures/3_dataplot.pdf", width=5, height=5)
ggplot(statement.effects, aes(x=perceived_truth, y=illusory_truth_effect)) +
  geom_point() +
  geom_smooth(color="black", lty=2, fill="black") +
  geom_smooth(color="black", method=lm, formula=y ~poly(x,2), fill="white") +
  geom_smooth(color="black", lty=2, se=FALSE) +
  labs(x="Perceived truth ([R+N]/2)", y="Illusory Truth Effect (R-N)") +
  theme_minimal()
dev.off()




##########################################
#identify peaks in models

#Loess estimate from stack overflow 
#https://stackoverflow.com/questions/12183137/calculate-min-max-slope-of-loess-fitted-curve-with-r


    #loess model
	loess.out <- loess(illusory_truth_effect~perceived_truth, data=statement.effects, degree=2)
    # Then the first difference
    	first_diff <- diff(loess.out$fitted)
    # Then the corresponding x and y values for the minima and maxima (only interested in the minimum, here)
   	 res <- cbind(loess.out$x[c(which.min(first_diff), which.max(first_diff))], 
            loess.out$fitted[c(which.min(first_diff), which.max(first_diff))])
    	colnames(res) <- c("x", "y")
    	rownames(res) <- c("min", "max")
    	res
	max_x <-res[1,1]



	#polynomial model
	 lm.out <- lm(illusory_truth_effect~poly(perceived_truth, 2, raw=T), data=statement.effects)
	 summary(lm.out)
	
	cf=coef(lm.out)
	# model fits
	perceived_truth = seq(min(statement.effects$perceived_truth), max(statement.effects$perceived_truth), length=1000)
	pred = data.frame(perceived_truth, illusory_truth_effect=predict(lm.out, newdata=data.frame(perceived_truth)))

	max_y <- max(pred[2])
	max_x <- pred[pred[2]==max_y,1]

	max_x
par(mfrow=c(2,2))
plot(lm.out)

plot(resid(loess.out) ~ fitted(loess.out))
qqnorm(resid(loess.out))
qqline(resid(loess.out))
hist(resid(loess.out))
#####################################
##Compute CIs on estimated peaks. Fazio et al. use bootstrap with 5K replicates, so we do too.

library(boot)

peak_function.lm <- function(formula, data, indices) {
	d <- data[indices,]
	lm.out <- lm(formula, data=d)
	perceived_truth = seq(min(d$perceived_truth), max(d$perceived_truth), length=1000)
	pred = data.frame(perceived_truth, illusory_truth_effect=predict(lm.out, newdata=data.frame(perceived_truth)))
	max_y <- max(pred[2])
	max_x <- pred[pred[2]==max_y,1]
	return(max_x)
}

peak_function.loess <- function(formula, degree, span, data,indices) {
	d <- data[indices,]
	loess.out <- loess(formula, data=d, degree=degree, span=span)
    # Then the first difference
    	first_diff <- diff(loess.out$fitted)
    # Then the corresponding x value for the minimum first_diff
     	max_x <- loess.out$x[which.min(first_diff)]
	return(max_x)
}	

set.seed(1234)

statement.effects <- read.csv("../Data/FazioExperiment/StatementEffects.csv")

reps.lm <- boot(data=statement.effects, statistic=peak_function.lm, R=5000, 
	formula=illusory_truth_effect~poly(perceived_truth, 2, raw=T))

reps.lm

plot(reps.lm)
lm.ci <-boot.ci(reps.lm, type="bca")

set.seed(1234)
reps.loess <- boot(data=statement.effects, statistic=peak_function.loess, R=5000, 
	formula=illusory_truth_effect~perceived_truth, degree=2, span=.75)

reps.loess

plot(reps.loess)
loess.ci<-boot.ci(reps.loess, type="bca")

matrix(rbind(unlist(lm.ci)[c(1,2,7,8)], unlist(loess.ci)[c(1,2,7,8)]), nrow=2,dimnames = list(c("lm", "loess"), c("Reps", "Estimate", "Lower", "Upper"))) 


all.data <-read.csv("../Data/FazioExperiment/full.range.csv")
t.test(all.data$prp_true_rep, all.data$prp_true_new, paired=TRUE)

x1 <- mean(all.data$prp_true_rep)
x2 <- mean(all.data$prp_true_new)
s1 <- sd(all.data$prp_true_rep)
s2 <- sd(all.data$prp_true_new)
n <-503
swithin <- sqrt(((n-1)*(s1^2 + s2^2))/(n*2-2))
Cohens.d <- (x2-x1)/swithin
Hedges.g <- Cohens.d*(1-(3/(4*502-1)))

pwr.t.test(n=80, d=NULL, sig.level=.05, power=.8)

lm.ci68 <- boot.ci(reps.lm, conf=.68, type="bca")
loess.ci68 <- boot.ci(reps.loess, type="perc")
ci68s <- matrix(rbind(unlist(lm.ci68)[c(1,2,7,8)], unlist(loess.ci68)[c(1,2,7,8)]), nrow=2,dimnames = list(c("lm", "loess"), c("Reps", "Estimate", "Lower", "Upper"))) 
x1 <- as.numeric(ci68s[1,2])
x2 <- as.numeric(ci68s[2,2])
s1 <- as.numeric(ci68s[1,4])-as.numeric(ci68s[1,3])
s2 <- as.numeric(ci68s[2,4])-as.numeric(ci68s[2,3])
n <- 5000
swithin <- sqrt(((n-1)*(s1^2 + s2^2))/(n*2-2))
Cohens.d <- abs(x2-x1)/swithin
Hedges.g <- Cohens.d*(1-(3/(4*n-1)))

CohensD <- function(x1,x2, s1, s2, n) {
	swithin <- sqrt(((n-1)*(s1^2 + s2^2))/(n*2-2))
	d <- abs(x2-x1)/swithin
	return(d)
}





