# Import data 
setwd("U:/Desktop/S_AdultCuriosity/DataAnalysis")
data = read.csv("datalp.csv",head=T)

# For some trials it is not possible to estimate Learning Progress (LP), e.g., in the first trial
data$LP[is.nan(data$LP)]=NA
# Remove data points for all variables to make AIC comparable
data$Time[is.nan(data$LP)]=NA 
data$novelty[is.nan(data$LP)]=NA
data$expLP[is.nan(data$expLP)]=NA

# Import packages
require(R.basic);
require(mgcv); # for additive models and for binomial model
require('qpcR'); # for plotting
require(mgcViz);

################################################################################
### FIT Additive Model (see Figure 2) ##########################################
data$subj <- factor(data$subj) # subjects as factor to add it as random effect

models=data.frame() # set up matrix where all models are saved
modname=data.frame()
b1 <- bam(switch ~ 1 + s(subj,bs="re"), family=binomial, data=data, discrete=TRUE, nthreads=2)
sumb=summary(b1); models[1,1]=sumb$r.sq; models[1,2]=sumb$dev.expl; models[1,3]=sumb$residual.df; models[1,4]=b1[["aic"]]
modname[1,1]="Null Model"

b2 <- bam(switch ~ s(LP) + s(subj,bs="re"), family=binomial, data=data, discrete=TRUE, nthreads=2)
sumb=summary(b2); models[2,1]=sumb$r.sq; models[2,2]=sumb$dev.expl; models[2,3]=sumb$residual.df; models[2,4]=b2[["aic"]]
modname[2,1]="Learning Progress"

b3 <- bam(switch ~ s(novelty)+ s(subj,bs="re"), family=binomial, data=data, discrete=TRUE, nthreads=2)
sumb=summary(b3); models[3,1]=sumb$r.sq; models[3,2]=sumb$dev.expl; models[3,3]=sumb$residual.df; models[3,4]=b3[["aic"]]
modname[3,1]="Novelty"

b4 <- bam(switch ~ s(Time) + s(subj,bs="re"), family=binomial, data=data, discrete=TRUE, nthreads=2)
sumb=summary(b4); models[4,1]=sumb$r.sq; models[4,2]=sumb$dev.expl; models[4,3]=sumb$residual.df; models[4,4]=b4[["aic"]]
modname[4,1]="Random Search"

bicmod=BIC(b1,b2,b3,b4) # interprer BIC with care because it does not keep into account that params are additive (AIC does)

aicstats=akaike.weights(models[,4]) #  get AIC weights
modelfinal=cbind(modname,models,aicstats[["deltaAIC"]],aicstats[["rel.LL"]],aicstats[["weights"]])
modelfinal=cbind(modelfinal,bicmod)
columns=c("Model","Adjusted R^2","Deviance Explained","residual df","AIC","deltaAIC","Relative Likelihood","AIC weights", "df", "BIC")
colnames(modelfinal)=columns
modelfinal=modelfinal[order(modelfinal$AIC),] #order from best to worst

# plot best model
b2=getViz(b2)
pdf("LP.pdf", width=1.75, height=2.15,colormodel="cmyk")
plot( sm(b2, 1) )+
  l_ciPoly(fill = "#d1e6fa") +
  l_fitLine(colour = "black", size=0.4) +
  l_fitLine() + geom_hline(yintercept = 0, linetype = 5,size=0.4,colour='black') +
  #geom_smooth(methsummary(b2od="loess",span = 0.3,colour = "black", size=1.3) +
  ylab("Hazard (a.u.)") +
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(breaks=seq(-1.4,1.41,.4), limits=c(-1,1)) +
  xlab("Learning Progress")+
  theme_classic() +
  theme(
    axis.line = element_line(colour = 'black', size = 0.4),
    axis.ticks = element_line(colour = "black", size = 0.4),
    axis.text = element_text(color = "black",size=6),
    axis.text.x = element_text(size = 7),
    axis.title.y = element_text(size = 7),
    axis.title.x = element_text(size = 7))
dev.off()

################################################################################
### Fit binomial model (Figure 4) ##############################################

# Import data
EIGdata = read.csv("EIGdata.csv",head=T)
EIGdata$subject <- factor(EIGdata$subject) #subjects as factor

b1 <- bam(choice ~ DeltaNovelty + deltaEIG+s(subject,bs="re"), family=binomial, data=EIGdata, discrete=TRUE, nthreads=2)
summary(b1)

# plot results
df=data.frame(name=c('Novelty','Learning Progress'), coef=c(0.09667, 0.08707), se=c(0.01912, 0.03120))
df$name <- factor(df$name, levels = c('Novelty','Learning Progress'))
limits <- aes(ymax = coef + se, ymin=coef - se)
p <- ggplot(df, aes(fill=name, y=coef, x=name))
p + geom_bar(colour="black", stat="identity", position=position_dodge(), size=0.8, width=0.8) +
  guides(fill=FALSE) +
  ylab("Beta coefficients (a.u.)") + # Set axis labels
  ggtitle("Saccadic Latencies") +     # Set title
  scale_x_discrete(labels=c('Novelty','Expected Learning Progress'))+
  scale_y_continuous(breaks=seq(0,0.15,0.05), limits=c(0, 0.15))+
  theme_bw()+
  scale_fill_manual(values=c("#d1e6fa", "#d1e6fa")) +
  geom_errorbar(limits, size=1, width=0.2)+
  theme(
    axis.text = element_text(color = "black",size=15),
    axis.text.x = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black",size=1))

################################################################################
### Subjective preference analyses (Figure 5) ##################################
data = read.csv("datalp13.csv",head=T)
# get overall time spent on each sequence (totalTime), and average prediction error (avgPE)
freqmat=data.frame() #cols= subj, char, tottrials, avgKL
for (i in 1:60) {
  thisdata <- subset(data, data$subj==i)
  for (n in 1:3) {
    freqmat[(i-1)*3+n,1]=i #subj
    freqmat[(i-1)*3+n,2]=n #char
    freqmat[(i-1)*3+n,3]=sum(thisdata$char==n) #totalTime
    freqmat[(i-1)*3+n,4]=mean(subset(thisdata,thisdata$char==n)$PE, na.rm = T)#avgPE
  }
}
colnames(freqmat)=c("subj","char","totalTime","avgPE")
freqmat$totalTime[freqmat$totalTime==0]=NA #doesn't make sense to check preference if they have never experienced a character
freqmat$char <- factor(freqmat$char) 
freqmat$subj <- factor(freqmat$subj) 

# used for plots but NOT for analysis
ggstatsplot::ggbetweenstats(
  data = freqmat,
  x = char,
  y = avgPE,
  xlab = "Conditions",
  ylab = "Average Prediction Error",
  sample.size.label = FALSE,
  nboot = 1000,
  ggtheme = ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(), axis.line = element_line(colour = "black"),
                           axis.title.y = element_text(face="plain")),
  messages = FALSE
)+ggplot2::scale_color_manual(values = c("#4c72b0", "#c44e52", "#55a868"))+
  ggplot2::scale_x_discrete(labels= c('Intermediate', 'High Volatility','High Noise'))

# same
ggstatsplot::ggbetweenstats(
  data = freqmat,
  x = char,
  y = totalTime,
  xlab = "Conditions",
  ylab = "Total Time",
  sample.size.label = FALSE,
  nboot = 1000,
  ggtheme = ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(), axis.line = element_line(colour = "black"),
                           axis.title.y = element_text(face="plain")),
  messages = FALSE
)+ggplot2::scale_color_manual(values = c("#4c72b0", "#c44e52", "#55a868"))+
  ggplot2::scale_x_discrete(labels= c('Intermediate', 'High Volatility','High Noise'))

freqmat$char <- factor(freqmat$char) 
freqmat$subj <- factor(freqmat$subj) 
# GLMs to predict effect of character type over avgPE and totalTime with subjects as random factor
library(lmerTest)
library(modelbased)
fit <- lmer(totalTime ~ char + (1|subj), data=freqmat)
anova(fit)
estimate_contrasts(fit)

fit <- lmer(avgPE ~ char + (1|subj), data=freqmat)
anova(fit)
estimate_contrasts(fit)

# get the bayes factor for such models
library(BayesFactor)
timedata = na.omit(freqmat)
full_BF = lmBF(totalTime ~ char, data = timedata, whichRandom = 'subj')
null_BF = lmBF(totalTime ~ ., data = timedata, whichRandom = 'subj')
full_BF / null_BF# The Bayes factor in favor of the full model

full_BF = lmBF(avgPE ~ char, data = timedata, whichRandom = 'subj')
null_BF = lmBF(avgPE ~ ., data = timedata, whichRandom = 'subj')
null_BF / full_BF# The Bayes factor in favor of the NULL model

#At last, analyze polldata (see plot in python
dat <- data.frame(Most  = c(17,12,9), 
                  Least  = c(5,12,18), 
                  row.names = c("Intermediate", "Volatile", "Noisy"))
chi.t <- chisq.test(dat)

#plot
mosaicplot(chi.t$observed, cex.axis =1 , main = "Observed counts")
