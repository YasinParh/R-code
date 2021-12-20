### mortality
data=read.csv("E:/Study_BSU/FALL_2019/THESIS/data/dataset_thesis.csv", header=TRUE,sep=",")
summary(data)
dim(data)
str(data)
library(tidyverse)
attach(data)
table(ethnicity)
data_cvd=mutate(data,status_cvd=ifelse(mortality_grouping=="CVD", 1,0),status_cancer=ifelse(mortality_grouping=="Cancer", 1,0),
                status_other=ifelse(mortality_grouping=="Other", 1,0)) ## separating different types of mortality status

###removing of variables
## ethnicity most of cases belong to same category
dat=select(data_cvd,-c(ï..ID, test_number,record_date,ethnicity,hdl,ldl,cancer_comment, cancer_yn, heart_comments, 
                        lung_comments,artery_comments,other_disease_comments,
                        peak_rpe,peak_dbp,peak_sbp,death_cause_generic,death_date,death_year, 
                        mortality_grouping, protocol, percent_body_fat, waist))
dim(dat)
#dim(data)
summary(dat)
## distribution of fitness rank
par(mfrow=c(1,1))
hist(fitness_rank, col = 5, main = "Distribution of Fitness Rank")
mean(fitness_rank)
## categories CRF
dat=mutate(dat,CRF=ifelse(fitness_rank<=0.33, c("low"), c("high"))) 
dat=mutate(dat,crf_d=ifelse(fitness_rank<=.33, 0,1))

table(dat$CRF)
attach(dat)
##Data Frame for our analysis
dat=data.frame(sex,age,weight,height,bmi,fitness_rank,glucose,trig,total_cholesterol,resting_sbp,resting_dbp,
               VO2_rel, VO2_abs,max_rer,resting_hr, max_hr,sex_code, mortality_status,inactivity,
               smoker,obesity,dyslipidemia,hypertension,diabetes,follow_up_yrs,med_lipid,med_diabetes,
               med_hypertensives, CRF,crf_d,status_cvd, status_cancer, status_other)
sum(is.na(dat))
dat=na.omit(dat)
dim(dat)
attach(dat)

train=dat[1:2000,]
test=dat[2001:3203,]
dim(test)
tail(test)



library(survival)
library(ggplot2)
library(ranger)
library(ggfortify)
library(dplyr)
### Bi-variate analysis for mortality under three scenarios
## kaplan mier estimate
### All cause mortality

par(mfrow=c(2,4))
km.fit.crf=survfit(Surv(follow_up_yrs,mortality_status)~CRF,data=dat)
plot(km.fit.crf , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time (year)", ylab="Probability ofSurvival", 
     main="CRF", sub = "Log-rank test p value 0.10")
# add a legend
legend(0, .3, c("CRF High", "CRF Low"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.crf,xlab="time", ylab="Survival", main="KM estimate for CRF")

km.fit.somking=survfit(Surv(follow_up_yrs,mortality_status)~smoker,data=dat)
plot(km.fit.somking , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time (year)", ylab="Probability ofSurvival", 
     main="Smoking Status", sub = "Log-rank test p value<.0001 " )
# add a legend
legend(0, .3, c("Non-Smoker", "Smoker"), 
       lty = c(1:2), col=c(5,4))


#autoplot(km.fit.somking)

km.fit.gender=survfit(Surv(follow_up_yrs,mortality_status)~sex,data=dat)
plot(km.fit.gender , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Sex", sub = "Log-rank p value <.0001")
# add a legend
legend(0, .3, c("Female", "Male"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.gender)


km.fit.hp=survfit(Surv(follow_up_yrs,mortality_status)~hypertension,data=dat)
plot(km.fit.hp , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Hypertension", sub = "Log-rank test p value <.0001")

legend(0, .3, c("HP absence", "HP presence"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.hp,xlab="time", ylab="Survival", main="KM estimate for hpertention")

#par(mfrow=c(1,4))
km.fit.obesity=survfit(Surv(follow_up_yrs,mortality_status)~obesity,data=dat)
plot(km.fit.obesity , col=c(5,4), lty=c(1:2),lwd = 1,
     xlab="Time(year)", ylab="Probability of Survival", main="Obesity",
     sub = "Log-rank test p value =.01")
# add a legend
legend(0, .3, c("Not obese", "Obese"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.obesity, xlab="time", ylab="Survival", main="KM estimate for obesity")

km.fit.dys=survfit(Surv(follow_up_yrs,mortality_status)~dyslipidemia,data=dat)

plot(km.fit.dys , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Dyslipidemia", sub = "Log-rank test p value <.0001")
# add a legend
legend(0, .3, c("Dysp absence", "Dysp presence"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.dys, xlab="time", ylab="Survival", main="KM estimate for dyslipidemia")

km.fit.inac=survfit(Surv(follow_up_yrs,mortality_status)~inactivity,data=dat)
plot(km.fit.inac , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", main="
     physical activity level", sub = "Log-rank test p value = .04")
# add a legend
legend(0, .3, c("Active", "Inactive"), 
       lty = c(1:2), col=c(5,4))

#autoplot(km.fit.inac, xlab="time", ylab="Survival", main="KM estimate for inactivity")

km.fit.diab=survfit(Surv(follow_up_yrs,mortality_status)~diabetes,data=dat)
plot(km.fit.diab , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Diabetes", sub = "Log-rank test p value <.0001")
# add a legend
legend(0, .3, c("Diab abcence", "Diab precence"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.diab, xlab="time", ylab="Survival", main="KM estimate for diabetes")


## Log rank test

l.crf=survdiff(Surv(follow_up_yrs,mortality_status)~CRF, rho=0, data=dat)
l.crf

l.smoker=survdiff(Surv(follow_up_yrs,mortality_status)~smoker, rho=0, data=dat)
l.smoker
l.sex=survdiff(Surv(follow_up_yrs,mortality_status)~sex, rho=0, data=dat)
l.sex
l.obes=survdiff(Surv(follow_up_yrs,mortality_status)~obesity, rho=0, data=dat)
l.obes
l.inac=survdiff(Surv(follow_up_yrs,mortality_status)~inactivity, rho=0, data=dat)
l.inac
l.diab=survdiff(Surv(follow_up_yrs,mortality_status)~diabetes, rho=0, data=dat)
l.diab

l.hyp=survdiff(Surv(follow_up_yrs,mortality_status)~hypertension, rho=0, data=dat)
l.hyp
l.dys=survdiff(Surv(follow_up_yrs,mortality_status)~dyslipidemia, rho=0, data=dat)
l.dys



### Mortality due to CVD
par(mfrow=c(2,4))
km.fit.crf=survfit(Surv(follow_up_yrs,status_cvd)~CRF,data=dat)
plot(km.fit.crf , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability ofSurvival", 
     main="CRF ", sub = "Log-rank test p value = 0.70")
# add a legend
legend(0, .3, c("CRF High", "CRF Low"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.crf,xlab="time", ylab="Survival", main="KM estimate for CRF")

km.fit.somking=survfit(Surv(follow_up_yrs,status_cvd)~smoker,data=dat)
plot(km.fit.somking , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability ofSurvival", 
     main="Smoking Status ", sub = "Log-rank test p value = 0.002")
# add a legend
legend(0, .3, c("Non-Smoker", "Smoker"), 
       lty = c(1:2), col=c(5,4))

#autoplot(km.fit.somking)

km.fit.gender=survfit(Surv(follow_up_yrs,status_cvd)~sex,data=dat)
plot(km.fit.gender , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Sex", sub = "Log-rank test p value <.0001")
# add a legend
legend(0, .3, c("Female", "Male"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.gender)


km.fit.hp=survfit(Surv(follow_up_yrs,status_cvd)~hypertension,data=dat)
plot(km.fit.hp , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Hypertension ", sub = "Log-rank test p value < .0001")

legend(0, .3, c("HP absence", "HP presence"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.hp,xlab="time", ylab="Survival", main="KM estimate for hpertention")

#par(mfrow=c(1,4))
km.fit.obesity=survfit(Surv(follow_up_yrs,status_cvd)~obesity,data=dat)
plot(km.fit.obesity , col=c(5,4), lty=c(1:2),lwd = 1,
     xlab="Time(year)", ylab="Probability of Survival", main="Obesity",
     sub = "Log-rank test p value 0.80")
# add a legend
legend(0, .3, c("Not obese", "Obese"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.obesity, xlab="time", ylab="Survival", main="KM estimate for obesity")

km.fit.dys=survfit(Surv(follow_up_yrs,status_cvd)~dyslipidemia,data=dat)

plot(km.fit.dys , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Dyslipidemia ", sub = "Log-rank test p value = 0.02")
# add a legend
legend(0, .3, c("Dysp absence", "Dysp presence"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.dys, xlab="time", ylab="Survival", main="KM estimate for dyslipidemia")

km.fit.inac=survfit(Surv(follow_up_yrs,status_cvd)~inactivity,data=dat)
plot(km.fit.inac , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Inactivity", sub = "Log-rank test p value = 0.80")
# add a legend
legend(0, .3, c("Active", "Inactive"), 
       lty = c(1:2), col=c(5,4))

#autoplot(km.fit.inac, xlab="time", ylab="Survival", main="KM estimate for inactivity")

km.fit.diab=survfit(Surv(follow_up_yrs,status_cvd)~diabetes,data=dat)
plot(km.fit.diab , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Diabetes", sub = "Log-rank test p value <.0001")
# add a legend
legend(0, .3, c("Diab abcence", "Diab precence"), 
       lty = c(1:2), col=c(5,4))
## Log rank test


l.crf=survdiff(Surv(follow_up_yrs,status_cvd)~CRF, rho=0, data=dat)
l.crf

l.smoker=survdiff(Surv(follow_up_yrs,status_cvd)~smoker, rho=0, data=dat)
l.smoker
l.sex=survdiff(Surv(follow_up_yrs,status_cvd)~sex, rho=0, data=dat)
l.sex
l.obes=survdiff(Surv(follow_up_yrs,status_cvd)~obesity, rho=0, data=dat)
l.obes
l.inac=survdiff(Surv(follow_up_yrs,status_cvd)~inactivity, rho=0, data=dat)
l.inac
l.diab=survdiff(Surv(follow_up_yrs,status_cvd)~diabetes, rho=0, data=dat)
l.diab

l.hyp=survdiff(Surv(follow_up_yrs,status_cvd)~hypertension, rho=0, data=dat)
l.hyp
l.dys=survdiff(Surv(follow_up_yrs,status_cvd)~dyslipidemia, rho=0, data=dat)
l.dys


### Mortality due to Cancer

par(mfrow=c(2,4))
km.fit.crf=survfit(Surv(follow_up_yrs,status_cancer)~CRF,data=dat)
plot(km.fit.crf , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability ofSurvival", main="CRF ",
     sub = "Log-rank test p value = 0.30")
# add a legend
legend(0, .3, c("CRF High", "CRF Low"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.crf,xlab="time", ylab="Survival", main="KM estimate for CRF")

km.fit.somking=survfit(Surv(follow_up_yrs,status_cancer)~smoker,data=dat)
plot(km.fit.somking , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability ofSurvival", 
     main="Smoking Status ", sub = "Log-rank test p value = 0.01")
# add a legend
legend(0, .3, c("Non-Smoker", "Smoker"), 
       lty = c(1:2), col=c(5,4))


#autoplot(km.fit.somking)

km.fit.gender=survfit(Surv(follow_up_yrs,status_cancer)~sex,data=dat)
plot(km.fit.gender , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Sex", sub = "Log-rank test p value = 0.001")
# add a legend
legend(0, .3, c("Female", "Male"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.gender)


km.fit.hp=survfit(Surv(follow_up_yrs,status_cancer)~hypertension,data=dat)
plot(km.fit.hp , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival",
     main="Hypertension ", sub = "Log-rank test p value = 0.03")

legend(0, .3, c("HP absence", "HP presence"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.hp,xlab="time", ylab="Survival", main="KM estimate for hpertention")

#par(mfrow=c(1,4))
km.fit.obesity=survfit(Surv(follow_up_yrs,status_cancer)~obesity,data=dat)
plot(km.fit.obesity , col=c(5,4), lty=c(1:2),lwd = 1,
     xlab="Time(year)", ylab="Probability of Survival", main="Obesity",
     sub = "Log-rank test p value = 0.06")
# add a legend
legend(0, .3, c("Not obese", "Obese"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.obesity, xlab="time", ylab="Survival", main="KM estimate for obesity")

km.fit.dys=survfit(Surv(follow_up_yrs,status_cancer)~dyslipidemia,data=dat)

plot(km.fit.dys , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Dyslipidemia ", sub = "Log-rank test p value = 0.02")
# add a legend
legend(0, .3, c("Dysp absence", "Dysp presence"), 
       lty = c(1:2), col=c(5,4))
#autoplot(km.fit.dys, xlab="time", ylab="Survival", main="KM estimate for dyslipidemia")

km.fit.inac=survfit(Surv(follow_up_yrs,status_cancer)~inactivity,data=dat)
plot(km.fit.inac , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Inactivity", sub = "Log-rank test p value = 0.90")
# add a legend
legend(0, .3, c("Active", "Inactive"), 
       lty = c(1:2), col=c(5,4))

#autoplot(km.fit.inac, xlab="time", ylab="Survival", main="KM estimate for inactivity")

km.fit.diab=survfit(Surv(follow_up_yrs,status_cancer)~diabetes,data=dat)
plot(km.fit.diab , col=c(5,4), lty=c(1:2),lwd = 1,
     xmax = 2000,xlab="Time(year)", ylab="Probability of Survival", 
     main="Diabetes", sub = "Log-rank test p value = 0.70")
# add a legend
legend(0, .3, c("Diab abcence", "Diab precence"), 
       lty = c(1:2), col=c(5,4))
## Log rank test


l.crf=survdiff(Surv(follow_up_yrs,status_cancer)~CRF, rho=0, data=dat)
l.crf

l.smoker=survdiff(Surv(follow_up_yrs,status_cancer)~smoker, rho=0, data=dat)
l.smoker
l.sex=survdiff(Surv(follow_up_yrs,status_cancer)~sex, rho=0, data=dat)
l.sex
l.obes=survdiff(Surv(follow_up_yrs,status_cancer)~obesity, rho=0, data=dat)
l.obes
l.inac=survdiff(Surv(follow_up_yrs,status_cancer)~inactivity, rho=0, data=dat)
l.inac
l.diab=survdiff(Surv(follow_up_yrs,status_cancer)~diabetes, rho=0, data=dat)
l.diab

l.hyp=survdiff(Surv(follow_up_yrs,status_cancer)~hypertension, rho=0, data=dat)
l.hyp
l.dys=survdiff(Surv(follow_up_yrs,status_cancer)~dyslipidemia, rho=0, data=dat)
l.dys



## Multivariate semi-parametric regression analysis
## All cause mortality
coxmodel.crf=coxph(Surv(follow_up_yrs,mortality_status)~CRF+age+sex+bmi+trig+obesity+glucose+
                     total_cholesterol+dyslipidemia+hypertension+inactivity+smoker+diabetes)
summary(coxmodel.crf)

## ROC
par(mfrow=c(1,3))
library(risksetROC)
eta=coxmodel.crf$linear.predictors
AUC=NULL
out=CoxWeights(marker = eta,Stime = follow_up_yrs,status = mortality_status,predict.time = 20)
AUC=out$AUC
AUC
ROC30=risksetROC(Stime=follow_up_yrs,status=mortality_status,marker=eta,predict.time = 20,
   
                             method = "Cox", main="ROC for all-cause mortality", sub="AUC =0.78")



## mortality due to cvd
coxmodel.crf.cvd=coxph(Surv(follow_up_yrs,status_cvd)~CRF+age+sex+bmi+trig+obesity+glucose+
                         total_cholesterol+dyslipidemia+hypertension+inactivity+smoker+diabetes)
summary(coxmodel.crf.cvd)

## ROC

eta=coxmodel.crf.cvd$linear.predictors
AUC=NULL
out=CoxWeights(marker = eta,Stime = follow_up_yrs,status =status_cvd,predict.time = 20)
AUC=out$AUC
AUC
ROC20=risksetROC(Stime=follow_up_yrs,status=status_cvd,marker=eta,predict.time = 20,
                 method = "Cox", main="ROC for mortality due to CVD",sub="AUC = 0.81")

## mortality due to cancer

coxmodel.crf.cancer=coxph(Surv(follow_up_yrs,status_cancer)~CRF+age+sex+bmi+trig+obesity+glucose+
                            total_cholesterol+dyslipidemia+hypertension+inactivity+smoker+diabetes)
summary(coxmodel.crf.cancer)

## ROC
eta=coxmodel.crf.cancer$linear.predictors
AUC=NULL
out=CoxWeights(marker = eta,Stime = follow_up_yrs,status =status_cancer,predict.time = 15)
AUC=out$AUC
AUC
ROC20=risksetROC(Stime=follow_up_yrs,status=status_cancer,marker=eta,predict.time = 15,
                 method = "Cox", main="ROC for motality due to cancer", sub="AUC=0.76")

### assessing the suitability of parametric regression model
## All Cause mortality
### Plot of log-cumulative hazard vs log t

cumhazard<- survfit(Surv(follow_up_yrs, mortality_status) ~ 1, data=dat)
par(mfrow=c(1,3))
plot(log(cumhazard$time),log(-log(cumhazard$surv)), lty=2:3, col="blue",
     xlab="log (time)", ylab="log(Cumulative Hazard)", main = " All-cause mortality", sub = "a")
##Weibull

library(SurvRegCensCov)

p.r.s.w.crf=WeibullReg(Surv(follow_up_yrs,mortality_status)~CRF+age+sex+bmi+trig+obesity+glucose+
                         total_cholesterol+dyslipidemia+hypertension+inactivity+smoker+diabetes,data=dat)
p.r.s.w.crf
summary(p.r.s.w.crf)
x=-(coefficients(p.r.s.w.crf))/0.42
expx=exp(-(coefficients(p.r.s.w.crf))/0.42)

## Mortality due to CVD

cumhazard<- survfit(Surv(follow_up_yrs, status_cvd) ~ 1, data=dat)
plot(log(cumhazard$time),log(-log(cumhazard$surv)), lty=2:3, col="green",
     xlab="log (time)", ylab="log(Cumulative Hazard)", main = " Mortality due to CVD",sub = "b")

##Weibull
p.r.s.w.fitness.cvd=WeibullReg(Surv(follow_up_yrs,status_cvd)~fitness_rank+age+sex+bmi+trig+obesity+glucose+
                                 total_cholesterol+dyslipidemia+hypertension+inactivity+smoker+diabetes)

summary(p.r.s.w.fitness.cvd)

p.r.s.w.crf.cvd=WeibullReg(Surv(follow_up_yrs,status_cvd)~CRF+age+sex+bmi+trig+obesity+glucose+
                             total_cholesterol+dyslipidemia+hypertension+inactivity+smoker+diabetes)
p.r.s.w.crf.cvd
summary(p.r.s.w.crf.cvd)
x=-(coefficients(p.r.s.w.crf.cvd))/0.387
expx=exp(-(coefficients(p.r.s.w.crf.cvd))/0.387)
AIC(p.r.s.w.crf.cvd)

### mortality due to Cancer
cumhazard<- survfit(Surv(follow_up_yrs, status_cancer) ~ 1, data=dat)
plot(log(cumhazard$time),log(-log(cumhazard$surv)), lty=2:3, col=2,
     xlab="log (time)", ylab="log(Cumulative Hazard)", main = " Mortality due to cancer", sub = "c")

p.r.s.w.crf.cancer=WeibullReg(Surv(follow_up_yrs,status_cancer)~CRF+age+sex+bmi+trig+obesity+glucose+
                                total_cholesterol+dyslipidemia+hypertension+inactivity+smoker+diabetes)
p.r.s.w.crf.cancer
round(p.r.s.w.crf.cancer$coef,3)
round(p.r.s.w.crf.cancer$HR,2)
summary(p.r.s.w.crf.cancer)

x=-(coefficients(p.r.s.w.crf.cancer))/0.522
expx=exp(-(coefficients(p.r.s.w.crf.cancer))/0.522)
AIC(p.r.s.w.crf.cancer)



##############################################################################
### Subgroup for all cause mortality
library(partykit)
colnames(dat)[colnames(dat)=="follow_up_yrs"] <- "time"
colnames(dat)[colnames(dat)=="mortality_status"] <- "cens"

dat=as.data.frame(dat)
attach(dat)

train=dat[1:2000,]
test=dat[2001:3203,]
dim(test)
tail(test)

#### Training dataset
wbreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  survreg(y ~ 0 + x, weights = weights, dist = "weibull", ...)
}

logLik.survreg <- function(object, ...)
  structure(object$loglik[2], df = sum(object$df), class = "logLik")

mob.survival.reg.train=mob(Surv(time,cens) ~ CRF | age+sex+bmi+trig+glucose+
                       total_cholesterol+hypertension+dyslipidemia+smoker+obesity+inactivity, data=train, fit = wbreg, 
                     control = mob_control(alpha=.1,minsize = 700))


print(mob.survival.reg.train,node=1)

plot(mob.survival.reg.train)

### test

mob.survival.reg.test=mob(Surv(time,cens) ~ CRF | age+sex+bmi+trig+glucose+
                             total_cholesterol+hypertension+dyslipidemia+smoker+obesity+inactivity, data=test, fit = wbreg, 
                           control = mob_control(alpha=.1,minsize = 500))


print(mob.survival.reg.test,node=1)

plot(mob.survival.reg.test)


####Sruvival regression for all cause mortality
wbreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  survreg(y ~ 0 + x, weights = weights, dist = "weibull", ...)
}

logLik.survreg <- function(object, ...)
  structure(object$loglik[2], df = sum(object$df), class = "logLik")

mob.survival.reg=mob(Surv(time,cens) ~ CRF | age+sex+bmi+trig+glucose+
                       total_cholesterol+hypertension+dyslipidemia+smoker+obesity+inactivity, data=dat, fit = wbreg, 
                     control = mob_control(alpha=.1,minsize = 1100))


print(mob.survival.reg,node=1)

plot(mob.survival.reg)

#### Subgroup 1 for all cause death
all.all.cph=coxph(Surv(follow_up_yrs,mortality_status)~ CRF, data=dat)
summary(all.all.cph)


#### Subgroup 1
subgroup.1.all=subset(dat,age<=49)
dim(subgroup.1.all)

cox.subgroup.1.all=coxph(Surv(time,cens)~ CRF, data=subgroup.1.all)
summary(cox.subgroup.1.all)

km.all.1=survfit(Surv(time,cens)~CRF,data=subgroup.1.all)
autoplot(km.all.1)

## Weibull reg
table(subgroup.1.all$cens)

WeibullReg(Surv(time,cens)~CRF,subgroup.1.all)

#### Subgroup 2 for all cause death

subgroup.2.all=subset(dat,age>49)
dim(subgroup.2.all)


cox.subgroup.2.all=coxph(Surv(time,cens)~ CRF, data=subgroup.2.all)
summary(cox.subgroup.2.all)

km2=survfit(Surv(time,cens)~CRF,data=subgroup.2.all)
autoplot(km2)

## Weibull reg
table(subgroup.2.all$cens)

WeibullReg(Surv(time,cens)~CRF,subgroup.2.all)

###### Kaplan Meier curve
par(mfrow=c(1,2))
#l.all.all=survdiff(Surv(time,cens)~CRF, rho=0, data=dat)
#l.all.all
#km.all.all=survfit(Surv(time,cens)~CRF,data=dat)
#plot(km.all.all, col=c(3,6), lty = 1:2,  xlab="Time", ylab="Probability of survival", 
#    main = "All individuals",sub=" Log-rank test p-value =0.1")
#legend(0,.2,c("CRF high","CRF low"),col=c(3,6),lty=1,2)

#l.all.1=survdiff(Surv(time,cens)~CRF, rho=0, data=subgroup.1.all)
#l.all.1

km.all.1=survfit(Surv(time,cens)~CRF,data=subgroup.1.all)
plot(km.all.1, col=c(3,6), lty = 1:2,  xlab="Time", ylab="Probability of survival", 
     main = "Subgroup I age <= 49 years",sub=" Log-rank test p-value <.0001")
legend(0,.3,c("CRF high","CRF low"),col=c(3,6),lty=1,2)

l.all.2=survdiff(Surv(time,cens)~CRF, rho=0, data=subgroup.2.all)
l.all.2
km.all.2=survfit(Surv(time,cens)~CRF,data=subgroup.2.all)
plot(km.all.2, col=c(3,6), lty = 1:2,  xlab="Time", ylab="Probability of survival",
     main = "Subgroup II age >49 years",sub=" Log-rank test p-value =0.05")
legend(0,.3,c("CRF high","CRF low"),col=c(3,6),lty=1,2)


### CVD

####Sruvival regression for cvd
wbreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  survreg(y ~ 0 + x, weights = weights, dist = "weibull", ...)
}

logLik.survreg <- function(object, ...)
  structure(object$loglik[2], df = sum(object$df), class = "logLik")

### TRain
mob.survival.cvd.train=mob(Surv(time,status_cvd) ~ CRF | age+sex+weight+bmi+trig+glucose+
                       total_cholesterol+hypertension+dyslipidemia+smoker+obesity+inactivity, data=train, fit = wbreg, 
                     control = mob_control(alpha=.1,minsize = 300))

plot(mob.survival.cvd.train)

### test
mob.survival.cvd.test=mob(Surv(time,status_cvd) ~ CRF | age+sex+weight+bmi+trig+glucose+
                             total_cholesterol+hypertension+dyslipidemia+smoker+obesity+inactivity, data=test, fit = wbreg, 
                           control = mob_control(alpha=.1,minsize = 400))

plot(mob.survival.cvd.test)

### Whole data set
mob.survival.cvd=mob(Surv(time,status_cvd) ~ CRF | age+sex+weight+bmi+trig+glucose+
                       total_cholesterol+hypertension+dyslipidemia+smoker+obesity+inactivity, data=dat, fit = wbreg, 
                     control = mob_control(alpha=.1,minsize = 700))


plot(mob.survival.cvd)

#### for all CVD death
cvd.all.cph=coxph(Surv(time,status_cvd)~ CRF, data=dat)
summary(cvd.all.cph)


#### Subgroup 1
subgroup.1.cvd=subset(dat,age<=43)
dim(subgroup.1.cvd)

cox.subgroup.1.cvd=coxph(Surv(time,status_cvd)~ CRF, data=subgroup.1.cvd)
summary(cox.subgroup.1.cvd)

km.cvd.1=survfit(Surv(time,status_cvd)~CRF,data=subgroup.1.cvd)
autoplot(km.cvd.1)

### Weibull reg
table(subgroup.1.cvd$status_cvd)

WeibullReg(Surv(time,status_cvd)~CRF,subgroup.1.cvd)

#### Subgroup 2

subgroup.2.cvd=subset(dat,age<=53 & age>43)
dim(subgroup.2.cvd)

cox.subgroup.2.cvd=coxph(Surv(time,status_cvd)~ CRF, data=subgroup.2.cvd)
summary(cox.subgroup.2.cvd)

km.cvd.2=survfit(Surv(time,status_cvd)~CRF,data=subgroup.2.cvd)
autoplot(km.cvd.2)

### Weibull reg
table(subgroup.2.cvd$status_cvd)

WeibullReg(Surv(time,status_cvd)~CRF,subgroup.2.cvd)

#### Subgroup 3

subgroup.3.cvd=subset(dat,age>53)
dim(subgroup.3.cvd)

cox.subgroup.3.cvd=coxph(Surv(time,status_cvd)~ CRF, data=subgroup.3.cvd)
summary(cox.subgroup.3.cvd)

km.cvd.3=survfit(Surv(time,status_cvd)~CRF,data=subgroup.3.cvd)
autoplot(km.cvd.3)
## Weibull
table(subgroup.3.cvd$status_cvd)

WeibullReg(Surv(time,status_cvd)~CRF,subgroup.3.cvd)

### km curve

par(mfrow=c(1,3))
#l.cvd.all=survdiff(Surv(time,status_cvd)~CRF, rho=0, data=dat)
#l.cvd.all
#km.cvd.all=survfit(Surv(time,status_cvd)~CRF,data=dat)
#plot(km.cvd.all, col=c(3,6), lty = 1:2,  xlab="Time", ylab="Probability of survival", 
#    main = "All individuals",sub=" Log-rank test p-value =0.7")
#legend(0,.4,c("CRF high","CRF low"),col=c(3,6),lty=1,2)


#l.cvd.1=survdiff(Surv(time,status_cvd)~CRF, rho=0, data=subgroup.1.cvd)
#l.cvd.1

km.cvd.1=survfit(Surv(time,status_cvd)~CRF,data=subgroup.1.cvd)
plot(km.cvd.1, col=c(3,6), lty = 1:2,  xlab="Time (year)", ylab="Probability of survival", 
     main = "Subgroup I (age <= 43)",sub=" Log-rank test p-value 0.2")
legend(0,.4,c("CRF high","CRF low"),col=c(3,6),lty=1,2)

l.cvd.2=survdiff(Surv(time,status_cvd)~CRF, rho=0, data=subgroup.2.cvd)
l.cvd.2

km.cvd.2=survfit(Surv(time,status_cvd)~CRF,data=subgroup.2.cvd)
plot(km.cvd.2, col=c(3,6), lty = 1:2,  xlab="Time (year)", ylab="Probability of survival", 
     main = "Subgroup II (43< age <=53)",sub=" Log-rank test p-value 0.02")
legend(0,.4,c("CRF high","CRF low"),col=c(3,6),lty=1,2)

l.cvd.3=survdiff(Surv(time,status_cvd)~CRF, rho=0, data=subgroup.3.cvd)
l.cvd.3

km.cvd.3=survfit(Surv(time,status_cvd)~CRF,data=subgroup.3.cvd)
plot(km.cvd.3, col=c(3,6), lty = 1:2,  xlab="Time (year)", ylab="Probability of survival", 
     main = "Subgroup III (age > 53)",sub=" Log-rank test p-value 1.00")
legend(0,.4,c("CRF high","CRF low"),col=c(3,6),lty=1,2)


### Cancer

wbreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
  survreg(y ~ 0 + x, weights = weights, dist = "weibull", ...)
}

logLik.survreg <- function(object, ...)
  structure(object$loglik[2], df = sum(object$df), class = "logLik")

### TRain

mob.survival.cancer.train=mob(Surv(time,status_cancer) ~ CRF | age+sex+weight+bmi+trig+glucose+
                          total_cholesterol+hypertension+dyslipidemia+smoker+obesity+inactivity, data=train, fit = wbreg, 
                        control = mob_control(alpha=.1,minsize = 700))


plot(mob.survival.cancer.train)

### test
mob.survival.cancer.test=mob(Surv(time,status_cancer) ~ CRF | age+sex+weight+bmi+trig+glucose+
                                total_cholesterol+hypertension+dyslipidemia+smoker+obesity+inactivity, data=test, fit = wbreg, 
                              control = mob_control(alpha=.1,minsize = 300))


plot(mob.survival.cancer.test)

mob.survival.cancer=mob(Surv(time,status_cancer) ~ CRF | age+sex+weight+bmi+trig+glucose+
                          total_cholesterol+hypertension+dyslipidemia+smoker+obesity+inactivity, data=dat, fit = wbreg, 
                        control = mob_control(alpha=.1,minsize = 800))


plot(mob.survival.cancer)

print(mob.survival.cancer, node=1)


#### Subgroup 1 for all cause death
cancer.all.cph=coxph(Surv(time,status_cancer)~ CRF, data=dat)
summary(cancer.all.cph)


#### Subgroup 1
subgroup.1.cancer=subset(dat,age<=37)
dim(subgroup.1.cancer)

cox.subgroup.1.cancer=coxph(Surv(time,status_cancer)~ CRF, data=subgroup.1.cancer)
summary(cox.subgroup.1.cancer)

km.cancer.1=survfit(Surv(time,status_cancer)~CRF,data=subgroup.1.cancer)
autoplot(km.cancer.1)

## Weibull reg
table(subgroup.1.cancer$status_cancer)

WeibullReg(Surv(time,status_cancer)~CRF,subgroup.1.cancer)

#### Subgroup 2
subgroup.2.cancer=subset(dat, age<=49 & age>37)
dim(subgroup.2.cancer)

cox.subgroup.2.cancer=coxph(Surv(time,status_cancer)~ CRF, data=subgroup.2.cancer)
summary(cox.subgroup.2.cancer)

km.cancer.2=survfit(Surv(time,status_cancer)~CRF,data=subgroup.2.cancer)
autoplot(km.cancer.2)

## Weibull reg
table(subgroup.2.cancer$status_cancer)

WeibullReg(Surv(time,status_cancer)~CRF,subgroup.2.cancer)


#### Subgroup 3
subgroup.3.cancer=subset(dat, age>49)
dim(subgroup.3.cancer)

cox.subgroup.3.cancer=coxph(Surv(time,status_cancer)~ CRF, data=subgroup.3.cancer)
summary(cox.subgroup.3.cancer)

km.cancer.3=survfit(Surv(time,status_cancer)~CRF,data=subgroup.3.cancer)
autoplot(km.cancer.1)

## Weibull reg
table(subgroup.3.cancer$status_cancer)

WeibullReg(Surv(time,status_cancer)~CRF,subgroup.3.cancer)


### km curve
## All
par(mfrow=c(1,3))
#l.cancer.all=survdiff(Surv(time,status_cancer)~CRF, rho=0, data=dat)
#l.cancer.all
#km.cancer.all=survfit(Surv(time,status_cancer)~CRF,data=dat)
#plot(km.cancer.all, col=c(3,6), lty = 1:2,  xlab="Time", ylab="Probability of survival", 
#    main = "All individuals",sub=" Log-rank test p-value =0.3")
#legend(0,.4,c("CRF high","CRF low"),col=c(3,6),lty=1,2)


#l.cvd.1=survdiff(Surv(time,status_cancer)~CRF, rho=0, data=subgroup.1.cancer)
#l.cvd.1

km.cancer.1=survfit(Surv(time,status_cancer)~CRF,data=subgroup.1.cancer)
plot(km.cancer.1, col=c(3,6), lty = 1:2,  xlab="Time (year)", ylab="Probability of survival", 
     main = "Subgroup I (age <= 37)",sub=" Log-rank test p-value 0.2")
legend(0,.4,c("CRF high","CRF low"),col=c(3,6),lty=1,2)


l.cvd.2=survdiff(Surv(time,status_cancer)~CRF, rho=0, data=subgroup.2.cancer)
l.cvd.2

km.cancer.2=survfit(Surv(time,status_cancer)~CRF,data=subgroup.2.cancer)
plot(km.cancer.2, col=c(3,6), lty = 1:2,  xlab="Time (year)", ylab="Probability of survival", 
     main = "Subgroup II (37<age <= 49)",sub=" Log-rank test p-value 0.2")
legend(0,.4,c("CRF high","CRF low"),col=c(3,6),lty=1,2)

l.cvd.3=survdiff(Surv(time,status_cancer)~CRF, rho=0, data=subgroup.3.cancer)
l.cvd.3

km.cancer.3=survfit(Surv(time,status_cancer)~CRF,data=subgroup.3.cancer)
plot(km.cancer.3, col=c(3,6), lty = 1:2,  xlab="Time (year)", ylab="Probability of survival", 
     main = "Subgroup III (age >49)",sub=" Log-rank test p-value 0.03")
legend(0,.4,c("CRF high","CRF low"),col=c(3,6),lty=1,2)
