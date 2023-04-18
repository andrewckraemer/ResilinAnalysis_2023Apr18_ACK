rm(list=ls())
library(RRPP)
library(lme4)
library(car)

data_main<-read.csv('ddPCR_Keime_Main.csv')
data_main$Age<-factor(data_main$Age, levels = c("Pupae", "Hatchling", "Nurse", "Forager", "Drone")) # reorder levels to correspond with age stage sequence

#Resilin by age class. Does not consider effect of hive, which is complicated by the fact that not every age class is sampled from every hive (result: hive is ignored in all analyses)
shapiro.test(data_main$Resilin)
leveneTest(data_main$Resilin~data_main$Age)
lm1<-lm.rrpp(data_main$Resilin~data_main$Age) #create linear model - we use the randomized residual permutation procedure because the pupae class is much more variable than the other groups
summary(lm1) #calculate results of statistical analysis for main model
summary(pairwise(lm1,groups=data_main$Age[-c(71,73,91,92)])) #post-hoc pairwise analysis

par( family="serif")
boxplot(data_main$Resilin~data_main$Age,col=c('cornflowerblue','orange','forestgreen','firebrick','mediumorchid1'),ylim=c(0,6.5), ylab="Gene expression (pro-resilin/eIF3-S8)", xlab="Age Class")
segments(1,6.3,3.5,6.3,lwd=2)
segments(1,6.3,1,6.1,lwd=2)
segments(3.5,6.3,3.5,6.1,lwd=2)
segments(2,6.1,5,6.1,lwd=2)
text(2.25,6.5,'***',cex=1.7)

# CTCF Data
rm(list=ls())
data_CTCF<-read.csv('ddPCR_Keime_CTCF_end2022.csv')
data_CTCF$Age<-factor(data_CTCF$Age, levels = c("Hatchling", "Nest","Forager","Drone")) # reorder levels to correspond with age stage sequence
data_CTCF$Resilin
shapiro.test(data_CTCF$CTCF) #not normally distributed with a natural log transformation

# CTCF values are best predicted by a model with the main effects of Age and Joint Type as well as an interaction term of age by joint type
lm.full<-lm.rrpp(data_CTCF$CTCF~data_CTCF$Age*data_CTCF$Joint.Type)
lm.reduced<-lm.rrpp(data_CTCF$CTCF~data_CTCF$Age+data_CTCF$Joint.Type)
summary(anova(lm.reduced, lm.full))
summary(lm.full)

# Joint Type 1m-cu LFV
lmculfv<-subset(data_CTCF,data_CTCF$Joint.Type=='1m-cu LFV')
#boxplot(lmculfv$CTCF~ lmculfv$Age)
summary(lm.rrpp(lmculfv$CTCF~ lmculfv$Age))
summary(pairwise(lm.rrpp(lmculfv$CTCF~ lmculfv$Age),groups=lmculfv$Age))$summary.table

# Joint Type 1m-cu RFV
lmcurfv<-subset(data_CTCF,data_CTCF$Joint.Type=='1m-cu RFV')
#boxplot(lmcurfv$CTCF~ lmcurfv$Age)
summary(lm.rrpp(lmcurfv$CTCF~ lmcurfv$Age))
summary(pairwise(lm.rrpp(lmcurfv$CTCF~ lmcurfv$Age),groups= lmcurfv$Age))$summary.table

# Joint Type Cu-V LFV
cuvlfv<-subset(data_CTCF,data_CTCF$Joint.Type=='Cu-V LFV')
#boxplot(log(cuvlfv$CTCF)~ cuvlfv$Age)
summary(aov(log(cuvlfv$CTCF)~ cuvlfv$Age))
summary(pairwise(lm.rrpp(cuvlfv$CTCF~ cuvlfv$Age),groups= cuvlfv$Age))$summary.table

# Joint Type Cu-V RFV
cuvrfv<-subset(data_CTCF,data_CTCF$Joint.Type=='Cu-V RFV')
#boxplot(log(cuvrfv$CTCF)~ cuvrfv$Age)
summary(aov(log(cuvrfv$CTCF)~ cuvrfv$Age))
summary(pairwise(lm.rrpp(cuvrfv$CTCF~ cuvrfv$Age),groups= cuvrfv$Age))$summary.table

# Joint Type Cu-V LHV
cuvlhv<-subset(data_CTCF,data_CTCF$Joint.Type=='Cu-V LHV')
#boxplot(log(cuvlhv$CTCF)~ cuvlhv$Age)
summary(aov(log(cuvlhv$CTCF)~ cuvlhv$Age))
summary(pairwise(lm.rrpp(cuvlhv$CTCF~ cuvlhv$Age),groups= cuvlhv$Age))$summary.table

# Joint Type Cu-V RHV
cuvrhv<-subset(data_CTCF,data_CTCF$Joint.Type=='Cu-V RHV')
boxplot(cuvrhv$CTCF~ cuvrhv$Age)
summary(aov(log(cuvrhv$CTCF)~ cuvrhv$Age))
summary(pairwise(lm.rrpp(cuvrhv$CTCF~ cuvrhv$Age),groups= cuvrhv$Age))$summary.table

par(mfrow=c(3,2), tcl=-0.5, family="serif", omi=c(0.2,0.6,0,0))

par(mai=c(0.3,0.4,0.4,0.1))
boxplot(lmculfv$CTCF~ lmculfv$Age,main='1m-cu LFV',col=c('orange','forestgreen','firebrick','darkolivegreen2','dodgerblue3'),ylim=c(-31000,430000),xlab="",ylab="",)
segments(1,390000,3,390000,lwd=2)
segments(1,390000,1,370000,lwd=2)
segments(3,390000,3,370000,lwd=2)
segments(2,370000,4,370000,lwd=2)
text(2,420000,'***',cex=1.7)

par(mai=c(0.3,0.4,0.4,0.1))
boxplot(lmcurfv$CTCF~ lmcurfv $Age,main='1m-cu RFV',col=c('orange','forestgreen','firebrick','darkolivegreen2','dodgerblue3'),ylim=c(-31000,430000),xlab="",ylab="")
segments(1,390000,3,390000,lwd=2)
segments(1,390000,1,370000,lwd=2)
segments(3,390000,3,370000,lwd=2)
segments(2,370000,4,370000,lwd=2)
text(2,420000,'***',cex=1.7)
segments(3,280000,4,280000,lwd=2)
segments(3,280000,3,260000,lwd=2)
segments(4,280000,4,260000,lwd=2)
text(3.5,310000,'***',cex=1.7)

par(mai=c(0.3,0.4,0.4,0.1))
boxplot(cuvlfv $CTCF~ cuvlfv $Age,main='Cu-V LFV',col=c('orange','forestgreen','firebrick','darkolivegreen2','dodgerblue3'),ylim=c(-31000,430000),xlab="",ylab="")
segments(1,390000,3,390000,lwd=2)
segments(1,390000,1,370000,lwd=2)
segments(3,390000,3,370000,lwd=2)
segments(2,370000,4,370000,lwd=2)
text(2,420000,'***',cex=1.7)

par(mai=c(0.3,0.4,0.4,0.1))
boxplot(cuvrfv $CTCF~ cuvrfv $Age,main='Cu-V RFV',col=c('orange','forestgreen','firebrick','darkolivegreen2','dodgerblue3'),ylim=c(-31000,430000),xlab="",ylab="")
segments(1,390000,3,390000,lwd=2)
segments(1,390000,1,370000,lwd=2)
segments(3,390000,3,370000,lwd=2)
segments(2,370000,4,370000,lwd=2)
text(2,420000,'***',cex=1.7)

par(mai=c(0.3,0.4,0.4,0.1))
boxplot(cuvlhv $CTCF~ cuvlhv $Age,main='Cu-V LHV',col=c('orange','forestgreen','firebrick','darkolivegreen2','dodgerblue3'),ylim=c(-31000,430000),xlab="",ylab="")
segments(1,390000,3,390000,lwd=2)
segments(1,390000,1,370000,lwd=2)
segments(3,390000,3,370000,lwd=2)
segments(2,370000,4,370000,lwd=2)
text(2,420000,'***',cex=1.7)

par(mai=c(0.3,0.4,0.4,0.1))
boxplot(cuvrhv $CTCF~ cuvrhv $Age,main='Cu-V RHV',col=c('orange','forestgreen','firebrick','darkolivegreen2','dodgerblue3'),ylim=c(-31000,430000),xlab="",ylab="")
segments(1,390000,3,390000,lwd=2)
segments(1,390000,1,370000,lwd=2)
segments(3,390000,3,370000,lwd=2)
segments(2,370000,4,370000,lwd=2)
text(2,420000,'***',cex=1.7)

par(xpd=NA)
text(x=c(-5.7),y=900000,"Corrected Total Cell Fluorescence (CTCF)",cex=1.5,font=2,srt=90)
