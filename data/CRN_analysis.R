##### Analysis Overview Zettersten & Lupyan (2018) ####


#load libraries
library(datasets) #version 3.3.1
library(graphics) #version 3.3.1
library(grDevices) #version 3.3.1
library(methods) #version 3.3.1
library(stats) #version 3.3.1
library(utils) #version 3.3.1

library(Matrix) #version 1.2-10
library(cowplot) #version 0.8.0
library(lme4) #version 1.1-14
library(plyr) #version 1.8.4
source('summarizeData.R') #helper functions for summarizing data

#read in data
d=read.csv("CRN_data.txt")

####Experiment 1A####

#summarizing subject performance
subj_1A=ddply(subset(d,experiment=="1A"),.(subject,condition),summarize,
              accuracy=mean(isRight),
              rt=mean(RT[RT<=5000]))

#overall accuracy
#high condition
mean(subj_1A$accuracy[subj_1A$condition=="high"])
confint(lm(accuracy~1,data=subset(subj_1A,condition=="high")))
#low condition
mean(subj_1A$accuracy[subj_1A$condition=="low"])
confint(lm(accuracy~1,data=subset(subj_1A,condition=="low")))

#overall RTs
#high condition
mean(subj_1A$rt[subj_1A$condition=="high"])
confint(lm(rt~1,data=subset(subj_1A,condition=="high")))
#low condition
mean(subj_1A$rt[subj_1A$condition=="low"])
confint(lm(rt~1,data=subset(subj_1A,condition=="low")))

#overall model fit
mLearn_1A=glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment=="1A"),family=binomial)
summary(mLearn_1A)
confint(mLearn_1A,method="Wald")[2:3,]

#interaction with trial number
mTrialInteraction_1A=glmer(isRight~conditionC*totalTrialNumC+ (1+totalTrialNumC|subject),data=subset(d,experiment=="1A"),family=binomial)
summary(mTrialInteraction_1A)
confint(mTrialInteraction_1A,method="Wald")[4:7,]


####Experiment 1B####

#summarizing subject performance
subj_1B=ddply(subset(d,experiment=="1B"),.(subject,condition),summarize,
              accuracy=mean(isRight),
              rt=mean(RT[RT<=5000]))

#overall accuracy
#high condition
mean(subj_1B$accuracy[subj_1B$condition=="high"])
confint(lm(accuracy~1,data=subset(subj_1B,condition=="high")))
#low condition
mean(subj_1B$accuracy[subj_1B$condition=="low"])
confint(lm(accuracy~1,data=subset(subj_1B,condition=="low")))

#overall RTs
#high condition
mean(subj_1B$rt[subj_1B$condition=="high"])
confint(lm(rt~1,data=subset(subj_1B,condition=="high")))
#low condition
mean(subj_1B$rt[subj_1B$condition=="low"])
confint(lm(rt~1,data=subset(subj_1B,condition=="low")))

#overall model fit
mLearn_1B=glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment=="1B"),family=binomial)
summary(mLearn_1B)
confint(mLearn_1B,method="Wald")[2:3,]

#interaction with trial number
mTrialInteraction_1B=glmer(isRight~conditionC*totalTrialNumC+ (1+totalTrialNumC|subject),data=subset(d,experiment=="1B"),family=binomial)
summary(mTrialInteraction_1B)
confint(mTrialInteraction_1B,method="Wald")[4:7,]

####Experiment 2A####

#summarizing subject performance
subj_2A=ddply(subset(d,experiment=="2A"),.(subject,condition),summarize,
              accuracy=mean(isRight),
              rt=mean(RT[RT<=5000]))

#overall accuracy
#high condition
mean(subj_2A$accuracy[subj_2A$condition=="high"])
confint(lm(accuracy~1,data=subset(subj_2A,condition=="high")))
#low condition
mean(subj_2A$accuracy[subj_2A$condition=="low"])
confint(lm(accuracy~1,data=subset(subj_2A,condition=="low")))

#overall RTs
#high condition
mean(subj_2A$rt[subj_2A$condition=="high"])
confint(lm(rt~1,data=subset(subj_2A,condition=="high")))
#low condition
mean(subj_2A$rt[subj_2A$condition=="low"])
confint(lm(rt~1,data=subset(subj_2A,condition=="low")))

#overall model fit
mLearn_2A=glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment=="2A"),family=binomial)
summary(mLearn_2A)
confint(mLearn_2A,method="Wald")[2:3,]

#interaction with trial number
mTrialInteraction_2A=glmer(isRight~conditionC*totalTrialNumC+ (1+totalTrialNumC|subject),data=subset(d,experiment=="2A"),family=binomial)
summary(mTrialInteraction_2A)
confint(mTrialInteraction_2A,method="Wald")[4:7,]

####Experiment 2B####

#summarizing subject performance
subj_2B=ddply(subset(d,experiment=="2B"),.(subject,condition),summarize,
              accuracy=mean(isRight),
              rt=mean(RT[RT<=5000]))

#overall accuracy
#high condition
mean(subj_2B$accuracy[subj_2B$condition=="high"])
confint(lm(accuracy~1,data=subset(subj_2B,condition=="high")))
#low condition
mean(subj_2B$accuracy[subj_2B$condition=="low"])
confint(lm(accuracy~1,data=subset(subj_2B,condition=="low")))

#overall RTs
#high condition
mean(subj_2B$rt[subj_2B$condition=="high"])
confint(lm(rt~1,data=subset(subj_2B,condition=="high")))
#low condition
mean(subj_2B$rt[subj_2B$condition=="low"])
confint(lm(rt~1,data=subset(subj_2B,condition=="low")))

#overall model fit
mLearn_2B=glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment=="2B"),family=binomial)
summary(mLearn_2B)
confint(mLearn_2B,method="Wald")[2:3,]

#interaction with trial number
mTrialInteraction_2B=glmer(isRight~conditionC*totalTrialNumC+ (1+totalTrialNumC|subject),data=subset(d,experiment=="2B"),family=binomial)
summary(mTrialInteraction_2B)
confint(mTrialInteraction_2B,method="Wald")[4:7,]

####Experiment 3####

#summarizing subject performance
subj_3=ddply(subset(d,experiment=="3"),.(subject,condition),summarize,
              accuracy=mean(isRight),
              rt=mean(RT[RT<=5000]))

#overall accuracy
#high condition
mean(subj_3$accuracy[subj_3$condition=="high"])
confint(lm(accuracy~1,data=subset(subj_3,condition=="high")))
#low condition
mean(subj_3$accuracy[subj_3$condition=="low"])
confint(lm(accuracy~1,data=subset(subj_3,condition=="low")))

#overall RTs
#high condition
mean(subj_3$rt[subj_3$condition=="high"])
confint(lm(rt~1,data=subset(subj_3,condition=="high")))
#low condition
mean(subj_3$rt[subj_3$condition=="low"])
confint(lm(rt~1,data=subset(subj_3,condition=="low")))

#overall model fit
mLearn_3=glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment=="3"),family=binomial)
summary(mLearn_3)
confint(mLearn_3,method="Wald")[2:3,]

#interaction with trial number
#model with random slope for totalTrialNumC does not converge
mTrialInteraction_3=glmer(isRight~conditionC*totalTrialNumC+ (1|subject),data=subset(d,experiment=="3"),family=binomial)
summary(mTrialInteraction_3)
confint(mTrialInteraction_3,method="Wald")[2:5,]


#### GRAPHING ####

#Experiment 1
subj_blocks=ddply(d,.(experiment,subject,condition,blockNum),summarize,
                           accuracy=mean(isRight))

#### Experiment 1 ####
subjOverall_byBlock_exp1A = summarySEwithin(subset(subj_blocks,experiment=="1A"),"accuracy",betweenvars=c("condition"),withinvars=c("blockNum"),idvar="subject")
subjOverall_byBlock_exp1B= summarySEwithin(subset(subj_blocks,experiment=="1B"),"accuracy",betweenvars=c("condition"),withinvars=c("blockNum"),idvar="subject")

## 1A
#plot by block
p1A = ggplot(subjOverall_byBlock_exp1A, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block\n\nExperiment 1A")+
  ylab("Training Accuracy")+
  ylim(c(0.48,1))+
  geom_hline(yintercept=0.5,linetype="dotted",size=1.1)+
  theme(legend.position=c(.7, .2),legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))

#1B
p1B = ggplot(subjOverall_byBlock_exp1B, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block\n\nExperiment 1B")+
  ylab("Training Accuracy")+
  ylim(c(0.48,1))+
  geom_hline(yintercept=0.5,linetype="dotted",size=1.1)+
  theme(legend.position=c(.7, .2), legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))

#plot together
p_exp1=plot_grid(p1A,p1B,ncol=2,labels=c("A","B"))
p_exp1

#### Experiment 2 ####
subjOverall_byBlock_exp2A = summarySEwithin(subset(subj_blocks,experiment=="2A"),"accuracy",betweenvars=c("condition"),withinvars=c("blockNum"),idvar="subject")
subjOverall_byBlock_exp2B= summarySEwithin(subset(subj_blocks,experiment=="2B"),"accuracy",betweenvars=c("condition"),withinvars=c("blockNum"),idvar="subject")

## 1A
#plot by block
p2A = ggplot(subjOverall_byBlock_exp2A, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block\n\nExperiment 2A")+
  ylab("Training Accuracy")+
  ylim(c(0.48,1))+
  geom_hline(yintercept=0.5,linetype="dotted",size=1.1)+
  theme(legend.position=c(.7, .2),legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))

#2B
p2B = ggplot(subjOverall_byBlock_exp2B, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block\n\nExperiment 2B")+
  ylab("Training Accuracy")+
  ylim(c(0.48,1))+
  geom_hline(yintercept=0.5,linetype="dotted",size=1.1)+
  theme(legend.position=c(.7, .2), legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))

#plot together
p_exp2=plot_grid(p2A,p2B,ncol=2,labels=c("A","B"))
p_exp2

#### Experiment 3 ####
subjOverall_byBlock_exp3 = summarySEwithin(subset(subj_blocks,experiment=="3"),"accuracy",betweenvars=c("condition"),withinvars=c("blockNum"),idvar="subject")

p3 = ggplot(subjOverall_byBlock_exp3, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block\n\nExperiment 3")+
  ylab("Training Accuracy")+
  ylim(c(0.45,0.775))+
  geom_hline(yintercept=0.5,linetype="dotted",size=1.1)+
  theme(legend.position=c(.2, .8), legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))
p3
