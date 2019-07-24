##### Analysis Overview ####

#load libraries
library(datasets) #version 3.5.2
library(graphics) #version 3.5.2
library(grDevices) #version 3.5.2
library(methods) #version 3.5.2
library(stats) #version 3.5.2
library(utils) #version 3.5.2
library(car) #version 3.0-2
library(Matrix) #version 1.2-15
library(cowplot) #version 0.9.3
library(lme4) #version 1.1-18-1
library(plyr) #version 1.8.4
library(tidyverse) #version 1.2.1
library(effsize) #version 0.7.4
source('summarizeData.R') #helper functions for summarizing data
source('vif.mer.R') #function for checking for multicollinearity

#### DATA ####

#read in data experimental data
d <- read.csv("../data/CRN_data.txt")
#read in verbal strategy coding data
verbal_strategy <- read.csv("../data/CRN_verbal_strategies.csv")
#read in trial-by-trial data from color discriminability norming study
color_rt <- read.csv("../data/color_discriminability_rt_data.csv")
#read in data on individual color properties
color_properties <- read.csv("../data/color_properties.csv")
#read in data on color pair discriminability
color_properties_discrim <- read.csv("../data/color_properties_discriminability.csv")
#read in data on individual shape properties
shape_properties <- read.csv("../data/shape_properties.csv")
#read in data on shape pair discriminability
shape_properties_discrim <- read.csv("../data/shape_properties_discriminability.csv")
#read in trial-by-trial data from shape discriminability norming study
shape_rt <- read.csv("../data/shape_discriminability_rt_data.csv")
#read in properties of tangram (Exp 4) shapes
tangram_properties <- read.csv("../data/tangram_properties.csv")

####Experiment 1A####

#summarizing subject performance
subj_1A <- filter(d,experiment=="1A") %>%
  group_by(subject,condition,Age,Gender,L1,totalTimeMinutes) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000]),
    num_excludedTrials_rt=sum(RT>5000),
    total_trials = sum(!is.na(RT)))

#overall descriptives
exp1A_descriptives <- subj_1A %>%
  group_by() %>%
  summarize(
    N=sum(!is.na(unique(subject))), 
    gender_f=sum(Gender=="Female"), 
    avg_age = round(mean(Age,na.rm=T),2), 
    min_age=round(min(Age,na.rm=T),2),
    max_age=round(max(Age,na.rm=T),2), 
    native_lang=sum(grepl("English",L1)),
    num_high_condition=sum(condition=="high"),
    mean_time=mean(totalTimeMinutes),
    sd_time=sd(totalTimeMinutes))
exp1A_descriptives

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
mLearn_1A <- glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment=="1A"),family=binomial)
summary(mLearn_1A)
confint(mLearn_1A,method="Wald")[2:3,]

#model fit controlling for item type and random item effects
mLearn_1A_item_non_conv <- glmer(isRight~conditionC+imageTypeC+ (1+imageTypeC|subject)+(1|imageName),data=subset(d,experiment=="1A"),family=binomial)
#model fails to converge, so simplify by removing random slope for imageTypeC
mLearn_1A_item <- glmer(isRight~conditionC+imageTypeC+ (1|subject)+(1|imageName),data=subset(d,experiment=="1A"),family=binomial)
summary(mLearn_1A_item)
confint(mLearn_1A_item,method="Wald")[3:5,]

#interaction with trial number
mTrialInteraction_1A <- glmer(isRight~conditionC*totalTrialNumC+ (1+totalTrialNumC|subject),data=subset(d,experiment=="1A"),family=binomial)
summary(mTrialInteraction_1A)
confint(mTrialInteraction_1A,method="Wald")[4:7,]

#effect size
effsize::cohen.d(subj_1A$accuracy,subj_1A$condition)

####Experiment 1B####

#summarizing subject performance
subj_1B <- filter(d,experiment=="1B") %>%
  group_by(subject,condition,Age,Gender,L1,totalTimeMinutes) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000]),
    num_excludedTrials_rt=sum(RT>5000),
    total_trials = sum(!is.na(RT)))

#overall descriptives
exp1B_descriptives <- subj_1B %>%
  group_by() %>%
  summarize(
    N=sum(!is.na(unique(subject))), 
    gender_f=sum(Gender=="Female"), 
    avg_age = round(mean(Age,na.rm=T),2), 
    min_age=round(min(Age,na.rm=T),2),
    max_age=round(max(Age,na.rm=T),2), 
    native_lang=sum(grepl("English",L1)),
    num_high_condition=sum(condition=="high"),
    mean_time=mean(totalTimeMinutes),
    sd_time=sd(totalTimeMinutes))
exp1B_descriptives

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
mLearn_1B <- glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment=="1B"),family=binomial)
summary(mLearn_1B)
confint(mLearn_1B,method="Wald")[2:3,]

#model fit controlling for item type and random item effects
mLearn_1B_item <- glmer(isRight~conditionC+imageTypeC+ (1+imageTypeC|subject)+(1|imageName),data=subset(d,experiment=="1B"),family=binomial)
summary(mLearn_1B_item) 
confint(mLearn_1B_item,method="Wald")[5:7,]

#interaction with trial number
mTrialInteraction_1B <- glmer(isRight~conditionC*totalTrialNumC+ (1+totalTrialNumC|subject),data=subset(d,experiment=="1B"),family=binomial)
summary(mTrialInteraction_1B)
confint(mTrialInteraction_1B,method="Wald")[4:7,]

#test for an interaction between experiment 1A and 1B
#center experiment so condition effect is interpretable
d$experiment1AB_centered <- ifelse(d$experiment=="1A",-0.5,
                                   ifelse(d$experiment=="1B",0.5,NA))
mExperimentInteraction_1AB <- glmer(isRight~conditionC*experiment1AB_centered+ (1|subject),data=subset(d,experiment %in% c("1A","1B")),family=binomial)
summary(mExperimentInteraction_1AB)
confint(mExperimentInteraction_1AB,method="Wald")[2:5,]

#effect size
effsize::cohen.d(subj_1B$accuracy,subj_1B$condition)

#Low-performing participants across experiments 1A and 1B
#summarize subject accuracy
subj_1 <- ddply(subset(d,experiment=="1A"|experiment=="1B"),.(subject,experiment,condition),summarize,
                   accuracy=mean(isRight))
#summarizing subject performance across blocks
subj_1_block <- ddply(subset(d,experiment=="1A"|experiment=="1B"),.(subject,experiment,condition,blockNum),summarize,
                    accuracy=mean(isRight))
#participants with lower than 75% accuracy on final block classified as "non learners"
non_learners <- unique(subj_1_block$subject[subj_1_block$accuracy<0.75&subj_1_block$blockNum==3])
subj_1$non_learner <- ifelse(subj_1$subject %in% non_learners,1,0)
#what conditions are non-learners in?
table(subj_1$condition[subj_1$non_learner==1])
#2 high, 16 low nameability condition
binom.test(16,18,p=0.5)
#model fit with non-learner participants removed
mLearn_1 <- glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment %in% c("1A","1B")&!(subject %in% non_learners)),family=binomial)
summary(mLearn_1)
confint(mLearn_1,method="Wald")[2:3,]

#### Experiments 1A-1B - Graph ####

#create dataframe aggregating subject accuracy across each block for each experiment
subj_blocks <- d %>% 
  group_by(experiment,subject,condition,blockNum) %>%
  summarize(accuracy=mean(isRight))

#summarize across subjects for Experiments 1A and 1B
subjOverall_byBlock_exp1A <- summarySEwithin(subset(subj_blocks,experiment=="1A"),"accuracy",betweenvars=c("condition"),withinvars=c("blockNum"),idvar="subject")
subjOverall_byBlock_exp1B <- summarySEwithin(subset(subj_blocks,experiment=="1B"),"accuracy",betweenvars=c("condition"),withinvars=c("blockNum"),idvar="subject")

## 1A
#plot by block
p1A <- ggplot(subjOverall_byBlock_exp1A, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block\nExperiment 1A")+
  ylab("Training Accuracy")+
  ylim(c(0.48,1))+
  geom_hline(yintercept=0.5,linetype="dotted",size=1.1)+
  theme(legend.position=c(.7, .2),legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))

#1B
p1B <-  ggplot(subjOverall_byBlock_exp1B, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block\nExperiment 1B")+
  ylab("Training Accuracy")+
  ylim(c(0.48,1))+
  geom_hline(yintercept=0.5,linetype="dotted",size=1.1)+
  theme(legend.position=c(.7, .2), legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))

#plot together
p_exp1 <- plot_grid(p1A,p1B,ncol=2,labels=c("A","B"))
p_exp1

####Experiment 2A####

#summarizing subject performance
subj_2A <- filter(d,experiment=="2A") %>%
  group_by(subject,condition,Age,Gender,L1) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000&RT>=200]),
    num_excludedTrials_rt=sum(RT>5000|RT<200),
    total_trials = sum(!is.na(RT)))

subj_2A_overall <- filter(d,experiment=="2A") %>%
  group_by(subject,Age,Gender,L1) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000&RT>=200]),
    num_excludedTrials_rt=sum(RT>5000|RT<200),
    total_trials = sum(!is.na(RT)))

#overall descriptives
exp2A_descriptives <- subj_2A_overall %>%
  group_by() %>%
  summarize(N=sum(!is.na(unique(subject))), gender_f=sum(Gender=="Female"), avg_age = round(mean(Age,na.rm=T),2), min_age=round(min(Age,na.rm=T),2),max_age=round(max(Age,na.rm=T),2), native_lang=sum(grepl("English",L1)))
exp2A_descriptives

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
t.test(subj_2A$rt[subj_2A$condition=="low"],subj_2A$rt[subj_2A$condition=="high"],var.equal=T)

#checking for a speed-accuracy tradeoff
m <- lmer(accuracy~rt+(1|subject),data=subj_2A)
summary(m)
ggplot(subj_2A,aes(rt,accuracy,color=condition))+
  geom_point()+
  geom_smooth(method="lm")

#main model fit
mLearn_2A <- glmer(isRight~conditionC+(1+conditionC|subject),data=filter(d, experiment=="2A"),family=binomial)
summary(mLearn_2A)  
confint(mLearn_2A,method="Wald")[4:5,]

#model fit controlling for item type and random item effects
mLearn_2A_item <- glmer(isRight~conditionC+ imageTypeC+(1+conditionC+imageTypeC|subject)+ (1|imageName),data=filter(d, experiment=="2A"),family=binomial,control=glmerControl(optimizer="bobyqa"))
summary(mLearn_2A_item)  
confint(mLearn_2A_item,method="Wald")[8:10,]

#main model fit - excluding very short (<200ms) and very long (>5000ms) RTs
mLearn_2A_rtExclusions <- glmer(isRight~conditionC+(1+conditionC|subject),data=filter(d, experiment=="2A"&(RT>=200&RT<=5000)),family=binomial)
summary(mLearn_2A_rtExclusions)  
confint(mLearn_2A,method="Wald")[4:5,]

#model fit controlling for item type and random item effects - excluding very short (<200ms) and very long (>5000ms) RTs
mLearn_2A_item_rtExclusions <- glmer(isRight~conditionC+ imageTypeC+(1+conditionC+imageTypeC|subject)+ (1|imageName),data=filter(d, experiment=="2A"&(RT>=200&RT<=5000)),family=binomial,control=glmerControl(optimizer="bobyqa"))
summary(mLearn_2A_item_rtExclusions)  
confint(mLearn_2A_item_rtExclusions,method="Wald")[8:10,]

#interaction with trial number
mTrialInteraction_2A <- glmer(isRight~conditionC*totalTrialNumC+(1+conditionC+totalTrialNumC+conditionC:totalTrialNumC|subject),data=filter(d, experiment=="2A"),family=binomial,control=glmerControl(optimizer="bobyqa"))
#model does not converge - simplify by successively removing lower-order random effect slopes
mTrialInteraction_2A <- glmer(isRight~conditionC*totalTrialNumC+(1+conditionC+conditionC:totalTrialNumC|subject),data=filter(d, experiment=="2A"),family=binomial,control=glmerControl(optimizer="bobyqa"))
summary(mTrialInteraction_2A) 
confint(mTrialInteraction_2A,method="Wald")

#effect size
subj_2A_differences <- d %>%
  filter(experiment=="2A") %>%
  group_by(subject) %>%
  summarize(accuracy_difference=mean(isRight[condition=="high"])-mean(isRight[condition=="low"]))
d_z <- mean(subj_2A_differences$accuracy_difference)/sd(subj_2A_differences$accuracy_difference)
psych::d.ci(d_z,n1=39)

####Experiment 2B####

#summarizing subject performance
subj_2B <- filter(d,experiment=="2B") %>%
  group_by(subject,condition,Age,Gender,L1) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000&RT>=200]),
    num_excludedTrials_rt=sum(RT>5000|RT<200),
    total_trials = sum(!is.na(RT)))

subj_2B_overall <- filter(d,experiment=="2B") %>%
  group_by(subject,Age,Gender,L1) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000&RT>=200]),
    num_excludedTrials_rt=sum(RT>5000|RT<200),
    total_trials = sum(!is.na(RT)))

#overall descriptives
exp2B_descriptives <- subj_2B_overall %>%
  group_by() %>%
  summarize(N=sum(!is.na(unique(subject))), gender_f=sum(Gender=="Female"), avg_age = round(mean(Age,na.rm=T),2), min_age=round(min(Age,na.rm=T),2),max_age=round(max(Age,na.rm=T),2), native_lang=sum(grepl("English",L1)))
exp2B_descriptives

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
t.test(subj_2B$rt[subj_2B$condition=="low"],subj_2B$rt[subj_2B$condition=="high"],var.equal=T)
#checking for a speed-accuracy tradeoff
m <- lmer(accuracy~rt+(1|subject),data=subj_2B)
summary(m)
ggplot(subj_2B,aes(rt,accuracy,color=condition))+
  geom_point()+
  geom_smooth(method="lm")

#main model fit
mLearn_2B <- glmer(isRight~conditionC+(1+conditionC|subject),data=filter(d, experiment=="2B"),family=binomial)
summary(mLearn_2B)  
confint(mLearn_2B,method="Wald")[4:5,]

#model fit controlling for item type and random item effects
mLearn_2B_item <- glmer(isRight~conditionC+ imageTypeC+(1+conditionC+imageTypeC|subject)+ (1|imageName),data=filter(d, experiment=="2B"),family=binomial)
summary(mLearn_2B_item)  
confint(mLearn_2A_item,method="Wald")[8:10,]

#main model fit - excluding very short (<200ms) and very long (>5000ms) RTs
mLearn_2B_rtExclusions <- glmer(isRight~conditionC+(1+conditionC|subject),data=filter(d, experiment=="2B"&(RT>=200&RT<=5000)),family=binomial)
summary(mLearn_2B_rtExclusions)  
confint(mLearn_2B,method="Wald")[4:5,]

#model fit controlling for item type and random item effects - excluding very short (<200ms) and very long (>5000ms) RTs
mLearn_2B_item_rtExclusions <- glmer(isRight~conditionC+ imageTypeC+(1+conditionC+imageTypeC|subject)+ (1|imageName),data=filter(d, experiment=="2B"&(RT>=200&RT<=5000)),family=binomial)
summary(mLearn_2B_item_rtExclusions)  
confint(mLearn_2B_item_rtExclusions,method="Wald")[8:10,]

#interaction with trial number
mTrialInteraction_2B <- glmer(isRight~conditionC*totalTrialNumC+(1+conditionC+totalTrialNumC+conditionC:totalTrialNumC|subject),data=filter(d, experiment=="2B"),family=binomial,control=glmerControl(optimizer="bobyqa"))
summary(mTrialInteraction_2B) 
confint(mTrialInteraction_2B,method="Wald")

#effect size
subj_2B_differences <- d %>%
  filter(experiment=="2B") %>%
  group_by(subject) %>%
  summarize(accuracy_difference=mean(isRight[condition=="high"])-mean(isRight[condition=="low"]))
d_z <- mean(subj_2B_differences$accuracy_difference)/sd(subj_2B_differences$accuracy_difference)
psych::d.ci(d_z,n1=39)

#Low-performing participants across experiments 2A and 2B
#summarize subject accuracy
subj_2 <- filter(d,experiment %in% c("2A","2B")) %>%
  group_by(subject,experiment,condition) %>%
  summarize(accuracy=mean(isRight))
#summarizing subject performance
subj_2_block <- filter(d,experiment %in% c("2A","2B")) %>%
  group_by(subject,experiment,condition,blockNum) %>%
  summarize(accuracy=mean(isRight))
#participants with lower than 75% accuracy on final block
non_learners_high <- unique(subj_2_block$subject[subj_2_block$accuracy<0.75&subj_2_block$blockNum==6&subj_2_block$condition=="high"])
non_learners_low <- unique(subj_2_block$subject[subj_2_block$accuracy<0.75&subj_2_block$blockNum==6&subj_2_block$condition=="low"])
subj_2$non_learner_high <- ifelse(subj_2$subject %in% non_learners_high,1,0)
subj_2$non_learner_low <- ifelse(subj_2$subject %in% non_learners_low,1,0)
#classifying non-learners/ low learners
subj_2$non_learner_classification <- ifelse(subj_2$non_learner_high!=subj_2$non_learner_low & subj_2$non_learner_high==1,"non_learner_high_only",
                                            ifelse(subj_2$non_learner_high!=subj_2$non_learner_low & subj_2$non_learner_low==1,"non_learner_low_only",
                                                   ifelse(subj_2$non_learner_high==1&subj_2$non_learner_low==1,"non_learner_both","learner")))
#spread dataframe
subj_2_wide <- subj_2 %>%
  spread(condition,accuracy)
#distribution of non-learners/ low learners
table(subj_2_wide$non_learner_classification)
non_learner_dist <- c(9,3,15)
prob <- c(1/3,1/3,1/3)
#3 high condition, 18 low condition, 9 both
library(EMT)
multinomial.test(non_learner_dist,prob, useChisq = T)
#model fit with these participants removed
m <- glmer(isRight~conditionC+(1+conditionC|subject),data=filter(d, (experiment %in% c("2A","2B"))&!(subject %in% non_learners_low) &!(subject %in% non_learners_high)),family=binomial)
summary(m)  
confint(m,method="Wald")

#predict accuracy from color feature characteristics
m <- glmer(isRight~avg_simpson_withinSubjCentered + imageTypeC+avg_dE2000_withinSubjCentered+avg_saturation_withinSubjCentered+avg_discr_rt_z+(1+avg_simpson_withinSubjCentered+imageTypeC+avg_discr_rt_z+avg_dE2000_withinSubjCentered+avg_saturation_withinSubjCentered|subject)+ (1|imageName),
        data=filter(d, experiment %in% c("2A","2B")),
        family=binomial,
        control=glmerControl(optimizer="bobyqa"))
#does not converge; simplify random effects by removing a less important random slope (imageTypeC)
m <- glmer(isRight~avg_simpson_withinSubjCentered + imageTypeC+avg_dE2000_withinSubjCentered+avg_saturation_withinSubjCentered+avg_discr_rt_z+(1+avg_simpson_withinSubjCentered+avg_discr_rt_z+avg_dE2000_withinSubjCentered+avg_saturation_withinSubjCentered|subject)+ (1|imageName),
        data=filter(d, experiment %in% c("2A","2B")),
        family=binomial,
        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=50000)))
vif.mer(m) #check for multicollinearity
summary(m)
confint(m,method="Wald")[17:22,]

####Experiments 2A-2B - Graph####

#summarize across subjects for Experiments 2A and 2B
subjOverall_byBlock_exp2A <- summarySEwithin(subset(subj_blocks,experiment=="2A"),"accuracy",betweenvars=c(),withinvars=c("condition","blockNum"),idvar="subject")
subjOverall_byBlock_exp2B <- summarySEwithin(subset(subj_blocks,experiment=="2B"),"accuracy",betweenvars=c(),withinvars=c("condition","blockNum"),idvar="subject")

#plot by block
p2A <- ggplot(subjOverall_byBlock_exp2A, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block\nExperiment 2A")+
  ylab("Training Accuracy")+
  scale_y_continuous(breaks=seq(0.4,1,0.1), limits=c(0.38,1))+
  theme(legend.position=c(.7, .3),legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))

#2B
p2B <- ggplot(subjOverall_byBlock_exp2B, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block\nExperiment 2B")+
  ylab("Training Accuracy")+
  scale_y_continuous(breaks=seq(0.4,1,0.1), limits=c(0.38,1))+
  theme(legend.position=c(.7, .3), legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))

#plot together
p_exp2=plot_grid(p2A,p2B,ncol=2,labels=c("A","B"))
p_exp2

####Experiments 1A, 1B, 2A, 2B - Color Properties####

#DE2000

#Color Set 1 - within-prototype dE2000 discriminability
#high nameability
high_prototype_pairs <- c("bluered","bluepurple","purplered","orangeyellow","brownorange","brownyellow")
mean(filter(color_properties_discrim,colorSet=="colorset1"&colorPair %in% high_prototype_pairs)$dE2000)
sd(filter(color_properties_discrim,colorSet=="colorset1"&colorPair %in% high_prototype_pairs)$dE2000)
#low nameability
low_prototype_pairs <- c("pinkturquoise","darkgreenbluepink","darkgreenblueturquoise","mustardneonyellow","lightredmustard","lightredneonyellow")
mean(filter(color_properties_discrim,colorSet=="colorset1"&colorPair %in% low_prototype_pairs)$dE2000)
sd(filter(color_properties_discrim,colorSet=="colorset1"&colorPair %in% low_prototype_pairs)$dE2000)
#compare high and low nameability
t.test(filter(color_properties_discrim,colorSet=="colorset1"&colorPair %in% high_prototype_pairs)$dE2000,filter(color_properties_discrim,colorSet=="colorset1"&colorPair %in% low_prototype_pairs)$dE2000,var.equal=T)

#Color Set 2 - within-prototype dE2000 discriminability
#high nameability
high_prototype_pairs_2 <- c("bluered","blueorange","orangered","brownpurple","browngrey","greypurple")
mean(filter(color_properties_discrim,colorSet=="colorset2"&colorPair %in% high_prototype_pairs_2)$dE2000)
sd(filter(color_properties_discrim,colorSet=="colorset2"&colorPair %in% high_prototype_pairs_2)$dE2000)
#low nameability
low_prototype_pairs_2 <- c("darkblueneonyellow","darkbluepinkpurple","neonyellowpinkpurple","lightgreenpink","lightgreenyellowbrown","pinkyellowbrown")
mean(filter(color_properties_discrim,colorSet=="colorset2"&colorPair %in% low_prototype_pairs_2)$dE2000)
sd(filter(color_properties_discrim,colorSet=="colorset2"&colorPair %in% low_prototype_pairs_2)$dE2000)
#compare high and low nameability
t.test(filter(color_properties_discrim,colorSet=="colorset2"&colorPair %in% high_prototype_pairs_2)$dE2000,filter(color_properties_discrim,colorSet=="colorset2"&colorPair %in% low_prototype_pairs_2)$dE2000,var.equal=T)

#Figure S1
ggplot(color_properties_discrim,aes(nameability, dE2000,color=nameability))+
  geom_violin()+
  geom_boxplot(alpha=0)+
  geom_dotplot(aes(fill=nameability),binaxis="y",stackdir="center")+
  theme(legend.position="none")+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  facet_wrap(~colorSet)

#RT BEHAVIORAL DISCRIMINABILITY

#participant accuracy - overall
color_rt %>%
  group_by(subjCode,colorSet) %>%
  summarize(average_accuracy=mean(isRight)) %>%
  summarySE("average_accuracy",groupvars=c("colorSet"))

#participant accuracy - split by nameability
color_rt %>%
  group_by(subjCode,colorSet,firstStimCategory) %>%
  summarize(average_accuracy=mean(isRight)) %>%
  summarySEwithin( "average_accuracy",betweenvars=c("colorSet"),withinvars=c("firstStimCategory"),idvar="subjCode")

#average reaction times - individual colors
#Figure 5
#convert into long format for same/different trials
color_properties_rt_sameDiff <- color_properties %>%
  unite("different trials",c("rt_different","rt_different_ci")) %>%
  unite("same trials", c("rt_same", "rt_same_ci")) %>%
  gather("isSame","rt_ci",c("different trials","same trials")) %>%
  separate(rt_ci,c("rt","ci"),sep="_") %>%
  mutate(rt=as.numeric(as.character(rt)),ci=as.numeric(as.character(ci)))
#color set 1 plot
p1 <- ggplot(subset(color_properties_rt_sameDiff,colorSet=="colorset1"), aes(colorName,rt, color=nameability))+
  geom_bar(stat="identity",size=1.5,fill="white", position=position_dodge(.95))+
  geom_errorbar(aes(ymin=rt-ci,ymax=rt+ci), width=0.3,size=1.2, position=position_dodge(.95))+
  coord_cartesian(ylim=c(500,700))+
  ggtitle("Color Set 1 (Experiment 2A)")+
  ylab("Reaction Time (in ms)")+
  xlab("Color")+
  scale_x_discrete(limits=c(
    "red",
    "brown",
    "blue",
    "yellow",
    "purple",
    "orange",
    "pink",
    "mustard",
    "turquoise",
    "neonyellow",
    "darkgreenblue",
    "lightred"),
    labels=c("(220, 20, 0) - red",
             "(120, 80, 40) - brown",
             "(30, 90, 210) - blue",
             "(250, 240, 0) - yellow",
             "(130, 30, 180) - purple",
             "(250, 120, 30) - orange",
             "(200, 170, 170) - grey",
             "(170, 160, 40) - mustard",
             "(150, 200, 180) - green",
             "(220, 240, 150) - pale green",
             "(70, 100, 90) - grey green",
             "(200, 100, 70) - brown"))+
  scale_color_brewer(palette="Set1",name="Nameability")+
  theme(legend.position=c(0.2,0.85)) + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))+
  facet_wrap(~isSame)
#color set 2 plot
p2 <- ggplot(subset(color_properties_rt_sameDiff,colorSet=="colorset2"), aes(colorName,rt, color=nameability))+
  geom_bar(stat="identity",size=1.5,fill="white", position=position_dodge(.95))+
  geom_errorbar(aes(ymin=rt-ci,ymax=rt+ci), width=0.3,size=1.2, position=position_dodge(.95))+
  coord_cartesian(ylim=c(500,700))+
  ggtitle("Color Set 2 (Experiment 2B)")+
  ylab("Reaction Time (in ms)")+
  xlab("Color")+
  scale_x_discrete(limits=c(
    "blue",
    "brown",
    "red",
    "purple",
    "grey",
    "orange",
    "neonyellow",
    "pink",
    "pinkpurple",
    "lightgreen",
    "darkblue",
    "yellowbrown"),
    labels=c("(10, 90, 210) - blue",
             "(100, 60, 20) - brown",
             "(240, 0, 10) - red",
             "(130, 30, 180) - purple",
             "(130, 130, 130) - grey",
             "(250, 120, 30) - orange",
             "(160, 180, 20) - green",
             "(200, 160, 180) - light purple",
             "(160, 0, 80) - magenta",
             "(190, 220, 210) - light blue",
             "(50, 80, 100) - blue",
             "(110, 100, 10) - brown"))+
  scale_color_brewer(palette="Set1",name="Nameability")+
  theme(legend.position=c(0.2,0.85)) + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))+
  facet_wrap(~isSame)

#final figure (Fig 5)
plot_grid(p1, p2, labels=c("A","B"), label_size=20, align="h")

#summarize across colors
sumNameabilityRTs <- summarySE(color_properties_rt_sameDiff, measurevar="rt",groupvars=c("colorSet","isSame","nameability")) %>%
  mutate(low_ci=rt-ci,high_ci=rt+ci)
#color set 1
#same trials
t.test(
  color_properties_rt_sameDiff$rt[color_properties_rt_sameDiff$colorSet=="colorset1"&color_properties_rt_sameDiff$isSame=="same trials"&color_properties_rt_sameDiff$nameability=="high"],
  color_properties_rt_sameDiff$rt[color_properties_rt_sameDiff$colorSet=="colorset1"&color_properties_rt_sameDiff$isSame=="same trials"&color_properties_rt_sameDiff$nameability=="low"],
  var.equal=T)
#different trials
t.test(
  color_properties_rt_sameDiff$rt[color_properties_rt_sameDiff$colorSet=="colorset1"&color_properties_rt_sameDiff$isSame=="different trials"&color_properties_rt_sameDiff$nameability=="high"],
  color_properties_rt_sameDiff$rt[color_properties_rt_sameDiff$colorSet=="colorset1"&color_properties_rt_sameDiff$isSame=="different trials"&color_properties_rt_sameDiff$nameability=="low"],
  var.equal=T)

#color set 2
#same trials
t.test(
  color_properties_rt_sameDiff$rt[color_properties_rt_sameDiff$colorSet=="colorset2"&color_properties_rt_sameDiff$isSame=="same trials"&color_properties_rt_sameDiff$nameability=="high"],
  color_properties_rt_sameDiff$rt[color_properties_rt_sameDiff$colorSet=="colorset2"&color_properties_rt_sameDiff$isSame=="same trials"&color_properties_rt_sameDiff$nameability=="low"],
  var.equal=T)
#different trials
t.test(
  color_properties_rt_sameDiff$rt[color_properties_rt_sameDiff$colorSet=="colorset2"&color_properties_rt_sameDiff$isSame=="different trials"&color_properties_rt_sameDiff$nameability=="high"],
  color_properties_rt_sameDiff$rt[color_properties_rt_sameDiff$colorSet=="colorset2"&color_properties_rt_sameDiff$isSame=="different trials"&color_properties_rt_sameDiff$nameability=="low"],
  var.equal=T)

#Color Pair Discriminability RTs
#Figure S2
ggplot(color_properties_discrim,aes(nameability, average_rt,color=nameability))+
  geom_violin()+
  geom_boxplot(alpha=0)+
  geom_dotplot(aes(fill=nameability),binaxis="y",stackdir="center")+
  theme(legend.position="none")+
  ylab("Average Reaction Time (in ms)")+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  facet_wrap(~colorSet)

#SATURATION

#Figure S3
ggplot(color_properties,aes(nameability,saturation,color=nameability))+
  geom_violin()+
  geom_boxplot()+
  geom_dotplot(aes(fill=nameability),binaxis="y",stackdir="center")+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  facet_wrap(~colorSet)

#average saturation for each color set (Exps 1A, 1B, 2A: colorset1; Exp 2B: colorset2)
saturationAverage <- summarySE(color_properties,"saturation",groupvars=c("nameability","colorSet"))
saturationAverage <- saturationAverage %>%
  mutate(low_ci=saturation-ci,high_ci=saturation+ci)
saturationAverage

#comparison between saturation for high and low nameability colors for each colorset
t.test(color_properties$saturation[color_properties$nameability=="high"& color_properties$colorSet=="colorset1"],color_properties$saturation[color_properties$nameability=="low"& color_properties$colorSet=="colorset1"],var.equal=T)
t.test(color_properties$saturation[color_properties$nameability=="high"& color_properties$colorSet=="colorset2"],color_properties$saturation[color_properties$nameability=="low"& color_properties$colorSet=="colorset2"],var.equal=T)

####Experiment 3A####

#summarizing subject performance
subj_3A <- filter(d,experiment=="3A") %>%
  group_by(subject,condition,Age,Gender,L1,totalTimeMinutes) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000]),
    num_excludedTrials_rt=sum(RT>5000),
    total_trials = sum(!is.na(RT)))

#overall descriptives
exp3A_descriptives <- subj_3A %>%
  group_by() %>%
  summarize(
    N=sum(!is.na(unique(subject))), 
    gender_f=sum(Gender=="Female"), 
    avg_age = round(mean(Age,na.rm=T),2), 
    min_age=round(min(Age,na.rm=T),2),
    max_age=round(max(Age,na.rm=T),2), 
    native_lang=sum(grepl("English",L1)),
    num_high_condition=sum(condition=="high"),
    mean_time=mean(totalTimeMinutes),
    sd_time=sd(totalTimeMinutes))
exp3A_descriptives

#overall accuracy
#high condition
mean(subj_3A$accuracy[subj_3A$condition=="high"])
confint(lm(accuracy~1,data=subset(subj_3A,condition=="high")))
#low condition
mean(subj_3A$accuracy[subj_3A$condition=="low"])
confint(lm(accuracy~1,data=subset(subj_3A,condition=="low")))

#overall RTs
#high condition
mean(subj_3A$rt[subj_3A$condition=="high"])
confint(lm(rt~1,data=subset(subj_3A,condition=="high")))
#low condition
mean(subj_3A$rt[subj_3A$condition=="low"])
confint(lm(rt~1,data=subset(subj_3A,condition=="low")))

#overall model fit
mLearn_3A <- glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment=="3A"),family=binomial)
summary(mLearn_3A)
confint(mLearn_3A,method="Wald")[2:3,]

#model fit controlling for item type and random item effects
mLearn_3A_item <- glmer(isRight~conditionC+imageTypeC+ (1+imageTypeC|subject)+(1|imageName),data=subset(d,experiment=="3A"),family=binomial)
summary(mLearn_3A_item)
confint(mLearn_3A_item,method="Wald")[5:7,]

#interaction with trial number
mTrialInteraction_3A <- glmer(isRight~conditionC*totalTrialNumC+ (1+totalTrialNumC|subject),data=subset(d,experiment=="3A"),family=binomial)
summary(mTrialInteraction_3A)
confint(mTrialInteraction_3A,method="Wald")[4:7,]

#effect size
effsize::cohen.d(subj_3A$accuracy,subj_3A$condition)

####Experiment 3B####

#summarizing subject performance
subj_3B <- filter(d,experiment=="3B") %>%
  group_by(subject,condition,Age,Gender,L1,totalTimeMinutes) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000]),
    num_excludedTrials_rt=sum(RT>5000),
    total_trials = sum(!is.na(RT)))

#overall descriptives
exp3B_descriptives <- subj_3B %>%
  group_by() %>%
  summarize(
    N=sum(!is.na(unique(subject))), 
    gender_f=sum(Gender=="Female"), 
    avg_age = round(mean(Age,na.rm=T),2), 
    min_age=round(min(Age,na.rm=T),2),
    max_age=round(max(Age,na.rm=T),2), 
    native_lang=sum(grepl("English",L1)),
    num_high_condition=sum(condition=="high"),
    mean_time=mean(totalTimeMinutes),
    sd_time=sd(totalTimeMinutes))
exp3B_descriptives

#overall accuracy
#high condition
mean(subj_3B$accuracy[subj_3B$condition=="high"])
confint(lm(accuracy~1,data=subset(subj_3B,condition=="high")))
#low condition
mean(subj_3B$accuracy[subj_3B$condition=="low"])
confint(lm(accuracy~1,data=subset(subj_3B,condition=="low")))

#overall RTs
#high condition
mean(subj_3B$rt[subj_3B$condition=="high"])
confint(lm(rt~1,data=subset(subj_3B,condition=="high")))
#low condition
mean(subj_3B$rt[subj_3B$condition=="low"])
confint(lm(rt~1,data=subset(subj_3B,condition=="low")))

#overall model fit
mLearn_3B <- glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment=="3B"),family=binomial)
summary(mLearn_3B)
confint(mLearn_3B,method="Wald")[2:3,]

#model fit controlling for stimulus (no difference in frequency of item type in experiment 3B to control for)
mLearn_3B_item <- glmer(isRight~conditionC+(1|subject)+(1|imageName),data=subset(d,experiment=="3B"),family=binomial)
summary(mLearn_3B_item) 
confint(mLearn_3B_item,method="Wald")[3:4,]

#interaction with trial number
mTrialInteraction_3B <- glmer(isRight~conditionC*totalTrialNumC+ (1+totalTrialNumC|subject),data=subset(d,experiment=="3B"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(mTrialInteraction_3B)
confint(mTrialInteraction_3B,method="Wald")[4:7,]

#test for an interaction between experiment 1A and 1B
#center experiment so condition effect is interpretable
d$experiment3AB_centered <- ifelse(d$experiment=="3A",-0.5,
                                   ifelse(d$experiment=="3B",0.5,NA))
mExperimentInteraction_3AB <- glmer(isRight~conditionC*experiment3AB_centered+ (1|subject),data=subset(d,experiment %in% c("3A","3B")),family=binomial)
summary(mExperimentInteraction_3AB)
confint(mExperimentInteraction_3AB,method="Wald")[2:5,]

#effect size
effsize::cohen.d(subj_3B$accuracy,subj_3B$condition)

#Low-performing participants across experiments 3A and 3B
#summarize subject accuracy
subj_3=ddply(subset(d,experiment=="3A"|experiment=="3B"),.(subject,experiment,condition),summarize,
             accuracy=mean(isRight))
#summarizing subject performance
subj_3_block <- ddply(subset(d,experiment=="3A"|experiment=="3B"),.(subject,experiment,condition,blockNum),summarize,
                   accuracy=mean(isRight))
#participants with lower than 75% accuracy on final block
non_learners <- unique(subj_3_block$subject[subj_3_block$accuracy<0.75&subj_3_block$blockNum==3])
subj_3$non_learner <- ifelse(subj_3$subject %in% non_learners,1,0)
#what conditions did non-learners belong to?
table(subj_3$condition[subj_3$non_learner==1])
#11 high, 36 low nameability condition
binom.test(36,47,p=0.5)
#model fit with these participants removed
mLearn_3 <- glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment %in% c("3A","3B")&!(subject %in% non_learners)),family=binomial)
summary(mLearn_3)
confint(mLearn_3,method="Wald")[2:3,]

#predict accuracy from shape feature characteristics
m <- glmer(isRight~avg_simpson + avg_discr_rt_z+avg_complexity_z+(1+avg_simpson + avg_discr_rt_z+avg_complexity_z|subject)+(1|imageName),
           data=filter(d, experiment %in% c("3A","3B")),
           family=binomial,
           control=glmerControl(optimizer="bobyqa"))
vif.mer(m) #check for multicollinearity
summary(m)
confint(m,method="Wald")[12:15,]

#### Experiments 3A-3B - Graph####
subjOverall_byBlock_exp3A <- summarySEwithin(subset(subj_blocks,experiment=="3A"),"accuracy",betweenvars=c("condition"),withinvars=c("blockNum"),idvar="subject")
subjOverall_byBlock_exp3B <-  summarySEwithin(subset(subj_blocks,experiment=="3B"),"accuracy",betweenvars=c("condition"),withinvars=c("blockNum"),idvar="subject")

#plot by block
p3A <-  ggplot(subjOverall_byBlock_exp3A, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block\nExperiment 3A")+
  ylab("Training Accuracy")+
  ylim(c(0.48,1))+
  geom_hline(yintercept=0.5,linetype="dotted",size=1.1)+
  theme(legend.position=c(.7, .2),legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))

#3B
p3B <- ggplot(subjOverall_byBlock_exp3B, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block\nExperiment 3B")+
  ylab("Training Accuracy")+
  ylim(c(0.48,1))+
  geom_hline(yintercept=0.5,linetype="dotted",size=1.1)+
  theme(legend.position=c(.7, .2), legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))

#plot together
p_exp3 <- plot_grid(p3A,p3B,ncol=2,labels=c("A","B"))
p_exp3

####Experiments 3A, 3B - Shape Properties####

####Nameability, Association Value, and Description Length

#all shapes
shape_properties_vg_shapes<- shape_properties %>%
  group_by(nameability) %>%
  summarize(
    simpson_diversity_mean=mean(simpson_diversity),
    modal_agreement_mean=mean(modal_response_agreement),
    description_length_mean=mean(description_length)
  )

#high vs. low Simpson's diversity index
t.test(shape_properties$simpson_diversity[shape_properties$nameability=="high"],
       shape_properties$simpson_diversity[shape_properties$nameability=="low"],var.equal=T)
#high vs. low modal agreement
t.test(shape_properties$modal_response_agreement[shape_properties$nameability=="high"],
       shape_properties$modal_response_agreement[shape_properties$nameability=="low"],var.equal=T)
#correlation between association value and Simpson's diversity index
cor.test(shape_properties$simpson_diversity,shape_properties$association_value)
#correlation between association value and Simpson's diversity index
cor.test(shape_properties$modal_response_agreement,shape_properties$association_value)

#Exp 3A shapes
shape_properties_exp3A<- shape_properties %>%
  filter(exp3A=="yes") %>%
  group_by(nameability) %>%
  summarize(
    simpson_diversity_mean=mean(simpson_diversity),
    modal_agreement_mean=mean(modal_response_agreement),
    description_length_mean=mean(description_length),
    description_length_ci_high=t.test(description_length)$conf.int[2],
    description_length_ci_low=t.test(description_length)$conf.int[1]
  )
#compare shape description length by nameability
t.test(filter(shape_properties,exp3A=="yes"&nameability=="high")$description_length,filter(shape_properties,exp3A=="yes"&nameability=="low")$description_length,var.equal=T)

#Exp 3B shapes
shape_properties_exp3B<- shape_properties %>%
  filter(exp3B=="yes") %>%
  group_by(nameability) %>%
  summarize(
    simpson_diversity_mean=mean(simpson_diversity),
    modal_agreement_mean=mean(modal_response_agreement),
    description_length_mean=mean(description_length),
    description_length_ci_high=t.test(description_length)$conf.int[2],
    description_length_ci_low=t.test(description_length)$conf.int[1]
  )
#compare shape description length by nameability
t.test(filter(shape_properties,exp3B=="yes"&nameability=="high")$description_length,
       filter(shape_properties,exp3B=="yes"&nameability=="low")$description_length,var.equal=T)

#RT BEHAVIORAL DISCRIMINABILITY

#participant accuracy - split by nameability
shape_rt %>%
  group_by(subjCode,firstStimCategory) %>%
  summarize(average_accuracy=mean(isRight)) %>%
  summarySEwithin( "average_accuracy",withinvars=c("firstStimCategory"),idvar="subjCode")

#reaction times - averaged across all individual shapes
#different trials
shape_properties %>%
  summarySE( "rt_different",groupvars=c("nameability")) %>%
  mutate(low_ci =rt_different-ci,high_ci=rt_different+ci)
#same trials
shape_properties %>%
  summarySE( "rt_same",groupvars=c("nameability")) %>%
  mutate(low_ci =rt_same-ci,high_ci=rt_same+ci)

#reaction times, pairwise comparisons - averaged across experiment 3A
#different trials
shape_properties_discrim %>%
  filter(exp3B=="yes"&isSame=="false") %>%
  summarySE( "average_rt",groupvars=c("nameability")) %>%
  mutate(low_ci = average_rt-ci,high_ci=average_rt+ci)
#range of values
range(filter(shape_properties_discrim, exp3B=="yes"&isSame=="false"&nameability=="high")$average_rt)
range(filter(shape_properties_discrim, exp3B=="yes"&isSame=="false"&nameability=="low")$average_rt)
#t-test
t.test(
  filter(shape_properties_discrim, exp3B=="yes"&isSame=="false"&nameability=="high")$average_rt,
  filter(shape_properties_discrim, exp3B=="yes"&isSame=="false"&nameability=="low")$average_rt,var.equal=T)

####Experiment 4####

#summarizing subject performance
subj_4 <- filter(d,experiment=="4") %>%
  group_by(subject,condition,Age,Gender,L1,totalTimeMinutes) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000]),
    num_excludedTrials_rt=sum(RT>5000),
    total_trials = sum(!is.na(RT)))

#overall descriptives
exp4_descriptives <- subj_4 %>%
  group_by() %>%
  summarize(
    N=sum(!is.na(unique(subject))), 
    gender_f=sum(Gender=="Female"), 
    avg_age = round(mean(Age,na.rm=T),2), 
    min_age=round(min(Age,na.rm=T),2),
    max_age=round(max(Age,na.rm=T),2), 
    native_lang=sum(grepl("English",L1)),
    num_high_condition=sum(condition=="high"),
    mean_time=mean(totalTimeMinutes),
    sd_time=sd(totalTimeMinutes))
exp4_descriptives

#overall accuracy
#high condition
mean(subj_4$accuracy[subj_4$condition=="high"])
confint(lm(accuracy~1,data=subset(subj_4,condition=="high")))
#low condition
mean(subj_4$accuracy[subj_4$condition=="low"])
confint(lm(accuracy~1,data=subset(subj_4,condition=="low")))

#overall RTs
#high condition
mean(subj_4$rt[subj_4$condition=="high"])
confint(lm(rt~1,data=subset(subj_4,condition=="high")))
#low condition
mean(subj_4$rt[subj_4$condition=="low"])
confint(lm(rt~1,data=subset(subj_4,condition=="low")))

#overall model fit
mLearn_4=glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment=="4"),family=binomial)
summary(mLearn_4)
confint(mLearn_4,method="Wald")[2:3,]

#model fit controlling for item type and random item effects
mLearn_4_item <- glmer(isRight~conditionC+imageTypeC+ (1+imageTypeC|subject)+(1|imageName),data=subset(d,experiment=="4"),family=binomial)
summary(mLearn_4_item)
confint(mLearn_4_item,method="Wald")[5:7,]

#interaction with trial number
mTrialInteraction_4 <- glmer(isRight~conditionC*totalTrialNumC+ (1+totalTrialNumC|subject),data=subset(d,experiment=="4"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(mTrialInteraction_4)
confint(mTrialInteraction_4,method="Wald")[4:7,]

#effect size
effsize::cohen.d(subj_4$accuracy,subj_4$condition)

#Low-performing participants in experiment 4
#summarizing subject performance by block
subj_4_block=ddply(filter(d,experiment=="4"),.(subject,experiment,condition,blockNum),summarize,
                   accuracy=mean(isRight))
#participants with lower than 75% accuracy on final block
non_learners <- unique(subj_4_block$subject[subj_4_block$accuracy<0.75&subj_4_block$blockNum==3])
subj_4$non_learner <- ifelse(subj_4$subject %in% non_learners,1,0)
#what conditions do low learners belong to?
table(subj_4$condition[subj_4$non_learner==1])
#15 high, 31 low nameability condition
binom.test(31,46,p=0.5)
#model fit with these participants removed
mLearn_4 <- glmer(isRight~conditionC+ (1|subject),data=subset(filter(d,experiment=="4"),!(subject %in% non_learners)),family=binomial)
summary(mLearn_4)
confint(mLearn_4,method="Wald")[2:3,]

####Experiment 4 - Graph####

subjOverall_byBlock_exp4 <- summarySEwithin(filter(subj_blocks, experiment=="4"),"accuracy",betweenvars=c("condition"),withinvars=c("blockNum"),idvar="subject")

p4 <- ggplot(subjOverall_byBlock_exp4, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block")+
  ylab("Training Accuracy")+
  ylim(c(0.48,1))+
  geom_hline(yintercept=0.5,linetype="dotted",size=1.1)+
  theme(legend.position=c(.25, .8), legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))
p4

####Experiment 4 - Tangram Properties####

#summarize nameability (Simpson's diversity index)
tangram_properties %>%
  summarySE("simpson_diversity",groupvars=c("nameability")) %>%
  mutate(low_ci=simpson_diversity-ci,high_ci=simpson_diversity+ci)
t.test(filter(tangram_properties,nameability=="high")$simpson_diversity,
       filter(tangram_properties,nameability=="low")$simpson_diversity,var.equal=T)

#summarize nameability (Simpson's diversity index) for images used in experiment 4
tangram_properties %>%
  filter(exp4=="yes") %>%
  summarySE("simpson_diversity",groupvars=c("nameability")) %>%
  mutate(low_ci=simpson_diversity-ci,high_ci=simpson_diversity+ci)
#t-test
t.test(filter(tangram_properties,exp4=="yes"&nameability=="high")$simpson_diversity,
       filter(tangram_properties,exp4=="yes"&nameability=="low")$simpson_diversity,var.equal=T)

#####Supplementary Materials: S2 - Self-reported Verbal Strategy####

#add no strategy columns
verbal_strategy$noStrategy <- ifelse(verbal_strategy$singleFeature==0&verbal_strategy$multipleFeature==0&verbal_strategy$holisticStrategy==0,1,0)
verbal_strategy$noStrategy_DimType <- ifelse(verbal_strategy$singleDimType==0&verbal_strategy$combinedDimType==0,1,0)

#convert to long form for tallying strategies used
verbal_strategy_long <- gather(verbal_strategy, codingType, is_used,singleFeature:noStrategy_DimType)

#combine with data (subset to experiments with verbal strategy info)
verbal_strategy_d <- filter(d, experiment %in% c("1A","1B","3A","3B","s4_shape+color"))
verbal_strategy_d <- merge(verbal_strategy_d,verbal_strategy)

####S2: Experiments 1A-1B

#overview over strategies reported
verbalSum1 <- subset(verbal_strategy_long, (experiment=="1A"|experiment=="1B")&codingType %in% c("singleFeature","multipleFeature","holisticStrategy","noStrategy")) %>%
  group_by(condition, codingType) %>%
  summarize(
    avg=mean(as.numeric(is_used)),
    N=sum(!is.na(as.numeric(is_used))),
    n=sum(as.numeric(is_used)),
    se.lower=binom.test(sum(as.numeric(is_used)),sum(!is.na(as.numeric(is_used))))$conf.int[1],
    se.upper=binom.test(sum(as.numeric(is_used)),sum(!is.na(as.numeric(is_used))))$conf.int[2])
verbalSum1

#test for differences in strategy use between conditions
#single feature
m <- glm(is_used~condition,data=subset(verbal_strategy_long, (experiment=="1A"|experiment=="1B")&codingType=="singleFeature"),family=binomial)
summary(m)
Anova(m,type="III")
#multiple feature
m <- glm(is_used~condition,data=subset(verbal_strategy_long, (experiment=="1A"|experiment=="1B")&codingType=="multipleFeature"),family=binomial)
summary(m)
Anova(m,type="III")
#holistic strategy
m <- glm(is_used~condition,data=subset(verbal_strategy_long, (experiment=="1A"|experiment=="1B")&codingType=="holisticStrategy"),family=binomial)
summary(m)
Anova(m,type="III")

#Strategy use and categorization accuracy
m <- glmer(isRight~conditionC+singleFeature+multipleFeature+holisticStrategy+ (1|subject),data=subset(verbal_strategy_d,(experiment=="1A"|experiment=="1B")),family=binomial)
summary(m)

####S2: Experiments 1A-1B - Graphing
p1 <- ggplot(verbalSum1,aes(codingType,avg,group=condition,color=condition,fill=condition))+
  geom_bar(stat="identity",position=position_dodge(.75),size=1,width=0.7)+
  geom_errorbar(aes(ymin=se.lower,ymax=se.upper),position=position_dodge(.75),width=0,color="black")+
  scale_x_discrete(name = "Strategy Type",limits=c("singleFeature","multipleFeature","holisticStrategy","noStrategy"),labels=c("single\nfeature","multiple\nfeature","holistic","none"))+
  ylab("Percent Responses")+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Strategy Type")+
  theme_classic(base_size=18)+
  theme(legend.position=c(0.5,0.8),legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))#+

#summarize subjects
subj_accuracy <- d %>%
  group_by(experiment,subject,condition) %>%
  summarize(accuracy=mean(isRight))
#merge with verbal strategy
subj_accuracy <- merge(subj_accuracy,verbal_strategy)

#summarize accuracy for each strategy type
sumSubjStrategySingleFeature_1 <- filter(subj_accuracy,experiment=="1A"|experiment=="1B") %>%
  group_by(condition,singleFeature) %>%
  summarize(strategy="single feature",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="singleFeature")
sumSubjStrategyMultipleFeature_1 <- filter(subj_accuracy,experiment=="1A"|experiment=="1B") %>%
  group_by(condition,multipleFeature) %>%
  summarize(strategy="multiple feature",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="multipleFeature")
sumSubjStrategyHolistic_1 <- filter(subj_accuracy,experiment=="1A"|experiment=="1B") %>%
  group_by(condition,holisticStrategy) %>%
  summarize(strategy="holistic",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="holisticStrategy")
sumSubjStrategyNoStrategy_1<- filter(subj_accuracy,experiment=="1A"|experiment=="1B") %>%
  group_by(condition,noStrategy) %>%
  summarize(strategy="none",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="noStrategy")
strategySummarized_1 <- bind_rows(sumSubjStrategySingleFeature_1,sumSubjStrategyMultipleFeature_1,sumSubjStrategyHolistic_1,sumSubjStrategyNoStrategy_1)

strategySummarized_1$strategyF <- factor(strategySummarized_1$strategy,levels=c("single feature","multiple feature","holistic","none"))
p1_Acc <- ggplot(strategySummarized_1,aes(condition,acc,group=as.factor(strategyCode),fill=as.factor(strategyCode)))+
  geom_bar(stat="identity",position=position_dodge(.9))+
  geom_errorbar(aes(ymin=se.lower,ymax=se.upper),color="black",position=position_dodge(.9),width=0.1)+
  facet_wrap(~strategyF)+
  scale_fill_brewer(palette="Accent",name="Strategy Used?",direction=-1,breaks=c(1,0),labels=c("yes","no"))+
  theme_classic(base_size=18)+
  ylab("Accuracy")

#plot all together
strat1 <- plot_grid(p1,p1_Acc,labels=c("A","B"))
strat1

####S2: Experiments 3A-3B

#overview over strategies reported
verbalSum3 <- subset(verbal_strategy_long, (experiment=="3A"|experiment=="3B")&codingType %in% c("singleFeature","multipleFeature","holisticStrategy","noStrategy")) %>%
  group_by(condition, codingType) %>%
  summarize(
    avg=mean(as.numeric(is_used)),
    N=sum(!is.na(as.numeric(is_used))),
    n=sum(as.numeric(is_used)),
    se.lower=binom.test(sum(as.numeric(is_used)),sum(!is.na(as.numeric(is_used))))$conf.int[1],
    se.upper=binom.test(sum(as.numeric(is_used)),sum(!is.na(as.numeric(is_used))))$conf.int[2])
verbalSum3

#test for differences in strategy use between conditions
#single feature
m <- glm(is_used~condition,data=subset(verbal_strategy_long, (experiment=="3A"|experiment=="3B")&codingType=="singleFeature"),family=binomial)
summary(m)
Anova(m,type="III")
#multiple feature
m <- glm(is_used~condition,data=subset(verbal_strategy_long, (experiment=="3A"|experiment=="3B")&codingType=="multipleFeature"),family=binomial)
summary(m)
Anova(m,type="III")
#holistic strategy
m <- glm(is_used~condition,data=subset(verbal_strategy_long, (experiment=="3A"|experiment=="3B")&codingType=="holisticStrategy"),family=binomial)
summary(m)
Anova(m,type="III")

#Strategy use and categorization accuracy
m <- glmer(isRight~conditionC+singleFeature+multipleFeature+holisticStrategy+ (1|subject),data=subset(verbal_strategy_d,(experiment=="3A"|experiment=="3B")),family=binomial)
summary(m)

####S2: Experiments 3A-3B - Graphing
p3 <- ggplot(verbalSum3,aes(codingType,avg,group=condition,color=condition,fill=condition))+
  geom_bar(stat="identity",position=position_dodge(.75),size=1,width=0.7)+
  geom_errorbar(aes(ymin=se.lower,ymax=se.upper),position=position_dodge(.75),width=0,color="black")+
  scale_x_discrete(name = "Strategy Type", limits=c("singleFeature","multipleFeature","holisticStrategy","noStrategy"),labels=c("single\nfeature","multiple\nfeature","holistic","none"))+
  ylab("Percent Responses")+
  xlab("Strategy Type")+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  theme_classic(base_size=18)+
  theme(legend.position=c(.8, .8),legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))+
  ylim(0,1)

#summarize accuracy for each strategy type
sumSubjStrategySingleFeature_3 <- filter(subj_accuracy,experiment=="3A"|experiment=="3B") %>%
  group_by(condition,singleFeature) %>%
  summarize(strategy="single feature",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="singleFeature")
sumSubjStrategyMultipleFeature_3 <- filter(subj_accuracy,experiment=="3A"|experiment=="3B") %>%
  group_by(condition,multipleFeature) %>%
  summarize(strategy="multiple feature",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="multipleFeature")
sumSubjStrategyHolistic_3 <- filter(subj_accuracy,experiment=="3A"|experiment=="3B") %>%
  group_by(condition,holisticStrategy) %>%
  summarize(strategy="holistic",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="holisticStrategy")
sumSubjStrategyNoStrategy_3<- filter(subj_accuracy,experiment=="3A"|experiment=="3B") %>%
  group_by(condition,noStrategy) %>%
  summarize(strategy="none",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="noStrategy")
strategySummarized_3 <- bind_rows(sumSubjStrategySingleFeature_3,sumSubjStrategyMultipleFeature_3,sumSubjStrategyHolistic_3,sumSubjStrategyNoStrategy_3)

strategySummarized_3$strategyF <- factor(strategySummarized_3$strategy,levels=c("single feature","multiple feature","holistic","none"))
p3_Acc <- ggplot(strategySummarized_3,aes(condition,acc,group=as.factor(strategyCode),fill=as.factor(strategyCode)))+
  geom_bar(stat="identity",position=position_dodge(.9))+
  geom_errorbar(aes(ymin=se.lower,ymax=se.upper),color="black",position=position_dodge(.9),width=0.1)+
  facet_wrap(~strategyF)+
  scale_fill_brewer(palette="Accent",name="Strategy Used?",direction=-1,breaks=c(1,0),labels=c("yes","no"))+
  theme_classic(base_size=18)+
  ylab("Accuracy")+
  coord_cartesian(ylim=c(0,1.1))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5,0.75,1))

#plot all together
strat3 <- plot_grid(p3,p3_Acc,labels=c("A","B"))
strat3

####Supplementary Materials: S3 - Verbal Interference####

#summarizing subject performance
subj_s3 <- filter(d,experiment=="s3_verbal_interference") %>%
  group_by(subject,condition,Age,Gender,L1) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000&RT>=200]),
    num_excludedTrials_rt=sum(RT>5000|RT<200),
    total_trials = sum(!is.na(RT)))

subj_s3_overall <- filter(d,experiment=="s3_verbal_interference") %>%
  group_by(subject,Age,Gender,L1) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000&RT>=200]),
    num_excludedTrials_rt=sum(RT>5000|RT<200),
    total_trials = sum(!is.na(RT)))

#overall descriptives
exps3_descriptives <- subj_s3_overall %>%
  group_by() %>%
  summarize(N=sum(!is.na(unique(subject))), gender_f=sum(Gender=="Female"), avg_age = round(mean(Age,na.rm=T),2), min_age=round(min(Age,na.rm=T),2),max_age=round(max(Age,na.rm=T),2), native_lang=sum(grepl("English",L1)))
exps3_descriptives

#overall accuracy
#high condition
mean(subj_s3$accuracy[subj_s3$condition=="high"])
confint(lm(accuracy~1,data=subset(subj_s3,condition=="high")))
#low condition
mean(subj_s3$accuracy[subj_s3$condition=="low"])
confint(lm(accuracy~1,data=subset(subj_s3,condition=="low")))

#interaction with experiment 2A
#add column for interference parameter
d$interference <- ifelse(d$experiment=="s3_verbal_interference","verbal interference",
                         ifelse(d$experiment=="2A","no verbal interference - Exp 2A",NA))
d$interferenceC <- ifelse(d$interference=="verbal interference",0.5,
                          ifelse(d$interference=="no verbal interference - Exp 2A",-0.5,NA))
mInterferenceInteraction_s3 <- glmer(isRight~conditionC*interferenceC+(1+conditionC|subject),data=filter(d, experiment %in% c("2A", "s3_verbal_interference")),family=binomial)
summary(mInterferenceInteraction_s3)
confint(mInterferenceInteraction_s3,method="Wald")[4:7,]
#model fit controlling for item type and random item effects
mInterferenceInteraction_item_s3 <- glmer(isRight~conditionC*interferenceC+ imageTypeC+(1+conditionC+imageTypeC|subject)+ (1|imageName),data=filter(d, experiment %in% c("2A", "s3_verbal_interference")),family=binomial)
summary(mInterferenceInteraction_item_s3)
confint(mInterferenceInteraction_item_s3,method="Wald")

#models fit while removing trials with very short (<200 ms) and very long (>5000ms) RTs
mInterferenceInteraction_rtExclusions_s3 <- glmer(isRight~conditionC*interferenceC+(1+conditionC|subject),data=filter(d, experiment %in% c("2A", "s3_verbal_interference")&(RT>=200&RT<=5000)),family=binomial)
summary(mInterferenceInteraction_s3)
confint(mInterferenceInteraction_s3,method="Wald")[4:7,]

#model fit controlling for item type and random item effects
mInterferenceInteraction_rtExclusions_item_s3 <- glmer(isRight~conditionC*interferenceC+ imageTypeC+(1+conditionC+imageTypeC|subject)+ (1|imageName),data=filter(d, experiment %in% c("2A", "s3_verbal_interference")&(RT>=200&RT<=5000)),family=binomial)
summary(mInterferenceInteraction_rtExclusions_item_s3)
confint(mInterferenceInteraction_rtExclusions_item_s3,method="Wald")

####Supplementary Materials: S3 - Graph ####
subj_blocks_s3_2A <- d %>% 
  filter(!is.na(interference)) %>%
  group_by(experiment,interference,subject,condition,blockNum) %>%
  summarize(accuracy=mean(isRight))

subjOverall_byBlock_s3_2A <- summarySEwithin(subj_blocks_s3_2A,"accuracy",betweenvars=c("interference"),withinvars=c("condition","blockNum"),idvar="subject")

#plot by block
ps3 <- ggplot(subjOverall_byBlock_s3_2A, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05),linetype="solid")+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block")+
  ylab("Training Accuracy")+
  scale_y_continuous(breaks=seq(0.4,1,0.1), limits=c(0.38,1))+
  theme(legend.position=c(.8, .3),legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))+
  facet_wrap(~interference)
ps3

####Supplementary Materials: S4 - Combined color and shape feature study####

#summarizing subject performance
subj_s4 <- filter(d,experiment=="s4_shape+color") %>%
  group_by(subject,condition,Age,Gender,L1,totalTimeMinutes) %>%
  summarize(
    accuracy=mean(isRight),
    rt=mean(RT[RT<=5000]),
    num_excludedTrials_rt=sum(RT>5000),
    total_trials = sum(!is.na(RT)))

#overall descriptives
exps4_descriptives <- subj_s4 %>%
  group_by() %>%
  summarize(
    N=sum(!is.na(unique(subject))), 
    gender_f=sum(Gender=="Female"), 
    avg_age = round(mean(Age,na.rm=T),2), 
    min_age=round(min(Age,na.rm=T),2),
    max_age=round(max(Age,na.rm=T),2), 
    native_lang=sum(grepl("English",L1)),
    num_high_condition=sum(condition=="high"),
    mean_time=mean(totalTimeMinutes),
    sd_time=sd(totalTimeMinutes))
exps4_descriptives

#overall accuracy
#high condition
mean(subj_s4$accuracy[subj_s4$condition=="high"])
confint(lm(accuracy~1,data=subset(subj_s4,condition=="high")))
#low condition
mean(subj_s4$accuracy[subj_s4$condition=="low"])
confint(lm(accuracy~1,data=subset(subj_s4,condition=="low")))

#overall RTs
#high condition
mean(subj_s4$rt[subj_s4$condition=="high"])
confint(lm(rt~1,data=subset(subj_s4,condition=="high")))
#low condition
mean(subj_s4$rt[subj_s4$condition=="low"])
confint(lm(rt~1,data=subset(subj_s4,condition=="low")))

#overall model fit
mLearn_s4 <- glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment=="s4_shape+color"),family=binomial)
summary(mLearn_s4)
confint(mLearn_s4,method="Wald")[2:3,]

#overall model fit w/ random intercept for stimulus
mLearn_s4_item <- glmer(isRight~conditionC+ (1|subject)+(1|imageName),data=subset(d,experiment=="s4_shape+color"),family=binomial)
summary(mLearn_s4_item)
confint(mLearn_s4_item,method="Wald")[3:4,]

#interaction with trial number
mTrialInteraction_s4 <- glmer(isRight~conditionC*totalTrialNumC+ (1+totalTrialNumC|subject),data=subset(d,experiment=="s4_shape+color"),family=binomial,glmerControl(optimizer="bobyqa"))
summary(mTrialInteraction_s4)
confint(mTrialInteraction_s4,method="Wald")[4:7,]

#look at non-learners/ low-performers
#create 8
d$blockNum8 <- ifelse(d$totalTrialNum<8,1,
                      ifelse(d$totalTrialNum<17,2,
                             ifelse(d$totalTrialNum<25,3,
                                    ifelse(d$totalTrialNum<33,4,
                                           ifelse(d$totalTrialNum<41,5,6)))))
#summarizing subject performance
subj_s4_block <- subset(d,experiment=="s4_shape+color") %>%
  group_by(subject,experiment,condition,blockNum8) %>%
  summarize(accuracy=mean(isRight))

#participants with lower than 75% accuracy on 8 final trials
non_learners <- unique(subj_s4_block$subject[subj_s4_block$accuracy<0.75&subj_s4_block$blockNum8==6])
subj_s4$non_learner <- ifelse(subj_s4$subject %in% non_learners,1,0)

#what conditions
table(subj_s4$condition[subj_s4$non_learner==1])
#33 high, 32 low nameability condition
binom.test(33,65,p=0.5)

#model fit with high-performing participants only
m <- glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment %in% c("s4_shape+color")&!(subject %in% non_learners)),family=binomial)
summary(m)
confint(m,method="Wald")[2:3,]

m <- glmer(isRight~conditionC+ (1|subject),data=subset(d,experiment %in% c("s4_shape+color")&(subject %in% non_learners)),family=binomial)
summary(m)
confint(m,method="Wald")[2:3,]

#visualization of the distribution of accuracy in high vs. low conditions for high- and low-performing participants
ggplot(subj_s4,aes(accuracy))+geom_histogram()+facet_wrap(~condition+non_learner)

#Verbal Strategy

#overview over strategies reported
verbalSums4 <- subset(verbal_strategy_long, (experiment=="s4_shape+color")&codingType %in% c("singleFeature","multipleFeature","holisticStrategy","noStrategy","singleDimType","combinedDimType","noStrategy_DimType")) %>%
  group_by(condition, codingType) %>%
  summarize(
    avg=mean(as.numeric(is_used)),
    N=sum(!is.na(as.numeric(is_used))),
    n=sum(as.numeric(is_used)),
    se.lower=binom.test(sum(as.numeric(is_used)),sum(!is.na(as.numeric(is_used))))$conf.int[1],
    se.upper=binom.test(sum(as.numeric(is_used)),sum(!is.na(as.numeric(is_used))))$conf.int[2])
verbalSums4

#test for differences in strategy use between conditions
#single feature
m <- glm(is_used~condition,data=subset(verbal_strategy_long, (experiment=="s4_shape+color")&codingType=="singleFeature"),family=binomial)
summary(m)
Anova(m,type="III")
#multiple feature
m <- glm(is_used~condition,data=subset(verbal_strategy_long, (experiment=="s4_shape+color")&codingType=="multipleFeature"),family=binomial)
summary(m)
Anova(m,type="III")
#holistic strategy
m <- glm(is_used~condition,data=subset(verbal_strategy_long, (experiment=="s4_shape+color")&codingType=="holisticStrategy"),family=binomial)
summary(m)
Anova(m,type="III")
#single dimension
m <- glm(is_used~condition,data=subset(verbal_strategy_long, (experiment=="s4_shape+color")&codingType=="singleDimType"),family=binomial)
summary(m)
Anova(m,type="III")
#multiple dimension
m <- glm(is_used~condition,data=subset(verbal_strategy_long, (experiment=="s4_shape+color")&codingType=="combinedDimType"),family=binomial)
summary(m)
Anova(m,type="III")

#Strategy use and categorization accuracy
#features used
m <- glmer(isRight~conditionC+singleFeature+multipleFeature+holisticStrategy+ (1|subject),data=subset(verbal_strategy_d,(experiment=="s4_shape+color")),family=binomial)
summary(m)
#single vs. combined dimension
m <- glmer(isRight~conditionC+singleDimType+combinedDimType+ (1|subject),data=subset(verbal_strategy_d,(experiment=="s4_shape+color")),family=binomial)
summary(m)

##### Supplementary Materials: S4 - Graph ####

subjOverall_byBlock_exps4 <- summarySEwithin(subset(subj_blocks,experiment=="s4_shape+color"),"accuracy",betweenvars=c("condition"),withinvars=c("blockNum"),idvar="subject")

ps4 <- ggplot(subjOverall_byBlock_exps4, aes(blockNum,accuracy,color=condition,group=condition))+
  geom_line(size=2,position=position_dodge(.05))+
  geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),width=0,size=0.75,position=position_dodge(.05))+
  geom_point(aes(fill=condition),size=4,position=position_dodge(.05))+
  theme_classic(base_size=24)+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  xlab("Block")+
  ylab("Training Accuracy")+
  scale_y_continuous(breaks=seq(0.5,1,0.1),limits=c(0.45,1))+
  geom_hline(yintercept=0.5,linetype="dotted",size=1.1)+
  theme(legend.position=c(.2, .75), legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))
ps4

####S4: Verbal Strategy - Graphing
ps4_verbal <- ggplot(verbalSums4,aes(codingType,avg,group=condition,color=condition,fill=condition))+
  geom_bar(stat="identity",position=position_dodge(.75),size=1,width=0.7)+
  geom_errorbar(aes(ymin=se.lower,ymax=se.upper),position=position_dodge(.75),width=0,color="black")+
  scale_x_discrete(name = "Strategy Type", limits=c("singleFeature","multipleFeature","holisticStrategy","noStrategy"),labels=c("single\nfeature","multiple\nfeature","holistic","none"))+
  ylab("Percent Responses")+
  xlab("Strategy Type")+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  theme_classic(base_size=18)+
  theme(legend.position=c(.5, .8),legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))

#summarize accuracy for each strategy type
sumSubjStrategySingleFeature_s4 <- filter(subj_accuracy,experiment=="s4_shape+color") %>%
  group_by(condition,singleFeature) %>%
  summarize(strategy="single feature",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="singleFeature")
sumSubjStrategyMultipleFeature_s4 <- filter(subj_accuracy,experiment=="s4_shape+color") %>%
  group_by(condition,multipleFeature) %>%
  summarize(strategy="multiple feature",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="multipleFeature")
sumSubjStrategyHolistic_s4 <- filter(subj_accuracy,experiment=="s4_shape+color") %>%
  group_by(condition,holisticStrategy) %>%
  summarize(strategy="holistic",
            acc=mean(accuracy),
            se.lower=ifelse(sum(!is.na(accuracy))>1,t.test(accuracy)$conf.int[1],NA),
            se.upper=ifelse(sum(!is.na(accuracy))>1,t.test(accuracy)$conf.int[2],NA)) %>%
  rename(strategyCode="holisticStrategy")
sumSubjStrategyNoStrategy_s4<- filter(subj_accuracy,experiment=="s4_shape+color") %>%
  group_by(condition,noStrategy) %>%
  summarize(strategy="none",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="noStrategy")
strategySummarized_s4 <- bind_rows(sumSubjStrategySingleFeature_s4,sumSubjStrategyMultipleFeature_s4,sumSubjStrategyHolistic_s4,sumSubjStrategyNoStrategy_s4)

strategySummarized_s4$strategyF <- factor(strategySummarized_s4$strategy,levels=c("single feature","multiple feature","holistic","none"))
ps4_Acc <- ggplot(strategySummarized_s4,aes(condition,acc,group=as.factor(strategyCode),fill=as.factor(strategyCode)))+
  geom_bar(stat="identity",position=position_dodge(.9))+
  geom_errorbar(aes(ymin=se.lower,ymax=se.upper),color="black",position=position_dodge(.9),width=0.1)+
  facet_wrap(~strategyF)+
  scale_fill_brewer(palette="Accent",name="Strategy Used?",direction=-1,breaks=c(1,0),labels=c("yes","no"))+
  theme_classic(base_size=18)+
  ylab("Accuracy")+
  coord_cartesian(ylim=c(0,1.1))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5,0.75,1))

#plot all together
strat_s4 <- plot_grid(ps4_verbal,ps4_Acc,labels=c("A","B"))
strat_s4

#plot dimensions used
ps4_dim=ggplot(subset(verbalSums4,codingType %in% c("singleDimType","combinedDimType","noStrategy_DimType")),aes(codingType,avg,group=condition,color=condition,fill=condition))+
  geom_bar(stat="identity",position=position_dodge(.75),size=1,width=0.7)+
  geom_errorbar(aes(ymin=se.lower,ymax=se.upper),position=position_dodge(.75),width=0,color="black")+
  scale_x_discrete(name = "Dimensional Combination", limits=c("singleDimType","combinedDimType","noStrategy_DimType"),labels=c("single\ndimension","combining\ndimensions","none"))+
  ylab("Percent Responses")+
  xlab("Dimensional Combination")+
  scale_color_brewer(palette="Set1",name="Nameability")+
  scale_fill_brewer(palette="Set1",name="Nameability")+
  theme_classic(base_size=18)+
  theme(legend.position=c(.8, .8),legend.text=element_text(size=16),legend.title=element_text(size=16,face="bold"))
ps4_dim


sumSubjStrategySingleDim_s4 <- filter(subj_accuracy,experiment=="s4_shape+color") %>%
  group_by(condition,singleDimType) %>%
  summarize(strategy="single dim.",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="singleDimType")
sumSubjStrategyCombinedDim_s4<- filter(subj_accuracy,experiment=="s4_shape+color") %>%
  group_by(condition,combinedDimType) %>%
  summarize(strategy="combining dim.",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="combinedDimType")
sumSubjStrategynoStrategyDim_s4 <- filter(subj_accuracy,experiment=="s4_shape+color") %>%
  group_by(condition,noStrategy_DimType) %>%
  summarize(strategy="none",
            acc=mean(accuracy),
            se.lower=t.test(accuracy)$conf.int[1],
            se.upper=t.test(accuracy)$conf.int[2]) %>%
  rename(strategyCode="noStrategy_DimType")
strategySummarized_s4 <- bind_rows(sumSubjStrategySingleDim_s4,sumSubjStrategyCombinedDim_s4,sumSubjStrategynoStrategyDim_s4)
#experiment S4
strategySummarized_s4$strategyF=factor(strategySummarized_s4$strategy,levels=c("single dim.","combining dim.","none"))
ps4_DimAcc=ggplot(strategySummarized_s4,aes(condition,acc,group=as.factor(strategyCode),fill=as.factor(strategyCode)))+
  geom_bar(stat="identity",position=position_dodge(.9))+
  geom_errorbar(aes(ymin=se.lower,ymax=se.upper),color="black",position=position_dodge(.9),width=0.1)+
  facet_wrap(~strategyF,ncol=2)+
  scale_fill_brewer(palette="Accent",name="Strategy Used?",direction=-1,breaks=c(1,0),labels=c("yes","no"))+
  theme_classic(base_size=18)+
  ylab("Accuracy")+
  coord_cartesian(ylim=c(0,1.1))+
  scale_y_continuous(breaks=c(0, 0.25, 0.5,0.75,1))
ps4_DimAcc

strat_s4_dim=plot_grid(ps4_dim,ps4_DimAcc,labels=c("A","B"))
strat_s4_dim
