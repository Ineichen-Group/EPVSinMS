#load libraries
library(ggplot2)
library(dplyr)
library(bestNormalize)
library(car)
library(faraway)
library(MASS)
library(AUC)
library(DescTools)
library(lme4)
library(Matrix)
library(irr)
library(readxl)
library(writexl)
library(pdftools)
library(stringr)
library(tidyverse)
library(forcats)
library(rlang)

#Data adjustment
#load file
PVS$MS_type <- as.factor(PVS$MS_type)
PVS$OCB <- as.factor(PVS$OCB)
PVS$BrainSegVolNotVent <- as.numeric(PVS$BrainSegVolNotVent)
PVS$CSOT <- as.numeric(PVS$CSOT)
PVS$BGCT <- as.numeric(PVS$BGCT)
PVS$BGAT <- as.numeric(PVS$BGAT)
PVS$BS <- as.numeric(PVS$BS)
PVS$Potter_CSO <- as.numeric(PVS$Potter_CSO)
PVS$Potter_BGC <- as.numeric(PVS$Potter_BGC)
PVS$Potter_BGA <- as.numeric(PVS$Potter_BGA)
PVS$Potter_BS <- as.numeric(PVS$Potter_BS)
PVS$T2_Count_FS <- as.numeric(PVS$T2_Count_FS)
PVS$Relapses <- as.numeric(PVS$Relapses)

Long$Relapses <- as.numeric(Long$Relapses)
Long$T2_Count_FS <- as.numeric(Long$T2_Count_FS)
Long$MS_type <- as.factor(Long$MS_type)
Long$Deltavol_CSO_bin <- as.factor(Long$Deltavol_CSO_bin)
Long$Deltavol_CSO_bin2 <- as.factor(Long$Deltavol_CSO_bin2)
Long$Deltavol_BG_bin <- as.factor(Long$Deltavol_BG_bin)
Long$Deltavol_BG_bin2 <- as.factor(Long$Deltavol_BG_bin2)
Long$Deltavol_tot_bin <- as.factor(Long$Deltavol_tot_bin)
Long$Deltavol_tot2 <- as.factor(Long$Deltavol_tot2)
#no log transformations in int type possible, only in numeric

Check:
summary(PVS$MS_type)
summary(PVS$OCB)
summary(PVS$T2_Volume_FS)
summary(PVS$T2_Count_FS)

#Descriptive statistics
str(PVS)
summary(PVS)

summary(PVS$Sex)
summary(PVS$MS_type)
summary(PVS$Age_Scan)
sd(PVS$Age_Scan)
summary(PVS$EDSS)
sd(PVS$EDSS, na.rm=TRUE)


summary(PVS$CSOT)
ggplot(PVS, aes(x=CSOT)) + 
  geom_histogram(binwidth = 0.5) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

summary(PVS$BGAT)
ggplot(PVS, aes(x=BGAT)) + 
  geom_histogram(binwidth = 0.5) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

summary(PVS$BGCT)
ggplot(PVS, aes(x=BGCT)) + 
  geom_histogram(binwidth = 0.5) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

summary(PVS$BS)
ggplot(PVS, aes(x=BS)) + 
  geom_histogram(binwidth = 0.5) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold"))

#Are BGCT and BGAT correlated with each other?
cor.test(PVS$BGCT, PVS$BGAT)
#yes, they are

#check for normal distribution
qqnorm(PVS$CSOT) #data points should be on a straight line under Gauss dis
shapiro.test(PVS$CSOT)
shapiro.test(PVS$BGAT) #p<0.05 indicates non-normal distributed dataset
install.packages("e1071")
library("e1071")
kurtosis(PVS$CSOT) #should be close to 0 for normal distribution
skewness(PVS$CSOT) #should be close to 0 for normal distribution
mean(PVS$CSOT)
median(PVS$CSOT) 
#mean and median should be close together for normal distribution
#Conclusion: CSOT is not normally distributed
logCSOT <- log(PVS$CSOT) #log-transform data
is.na(logCSOT) <- sapply(logCSOT, is.infinite) #remove inf values
logCSOT[!is.finite(logCSOT)] <- 0 #I do not want to remove the 0, thus I replace NA by 0
#similar for BG and BS
logBGCT <- log(PVS$BGCT)
logBGAT <- log(PVS$BGAT)
#no need to remove NA/Inf from BG because no patient w/o PVS in BG. Do not use log transformation for BS
#BoxCox transformation does not work since data contains 0s
#Yeo-Johnson Power transformation, CRAN package (=BoxCox transformation for data sets with 0 or negative values) https://rdrr.io/cran/bestNormalize/man/yeojohnson.html
yeojohnson(PVS$CSOT, eps = 0.001, standardize = TRUE)
YeoCSOT <- (PVS$CSOT)^0.17

qqnorm(Long$CSO_D) #data points should be on a straight line under Gauss dis
shapiro.test(Long$CSO_D) #p<0.05 indicates non-normal distributed dataset


#explorative statistics
#1.) age
ggplot(PVS, aes(x=Age_Scan, y=CSOT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#CSOT
cor.test(PVS$CSOT, PVS$Age_Scan)
lm(CSOT ~ Age_Scan, data = PVS) #http://r-statistics.co/Linear-Regression.html
lm(PVS$CSOT ~ PVS$Age_Scan) #same as above
summary(lm(CSOT ~ Age_Scan, data = PVS)) #dist=Intercept+(β∗Age_Scan)
summary(lm(logCSOT ~ Age_Scan, data = PVS))
summary(lm(logCSOT ~ Sex, data = PVS))
summary(lm(logCSOT ~ MS_type, data = PVS))
summary(lm(logCSOT ~ MS_type+Age_Scan, data = PVS))
summary(lm(logCSOT ~ Age_Scan+Sex, data=PVS))
summary(lm(logCSOT ~ MS_type+Sex, data=PVS))
summary(lm(logCSOT ~ Age_Scan+Sex+MS_type, data = PVS))

summary(lm(CSOT ~ Age_Scan+Sex, data=PVS))

summary(PVS$MS_type)

summary(lm(SqrtCSOT ~ Age_Scan, data = PVS))
summary(lm(SqrtCSOT ~ Sex, data = PVS))
summary(lm(SqrtCSOT ~ MS_type, data = PVS))
summary(lm(SqrtCSOT ~ MS_type+Age_Scan, data = PVS))
summary(lm(SqrtCSOT ~ Age_Scan+Sex, data=PVS))
summary(lm(SqrtCSOT ~ MS_type+Sex, data=PVS))
summary(lm(SqrtCSOT ~ Age_Scan+Sex+MS_type, data = PVS))

summary(lm(YeoCSOT ~ Age_Scan, data = PVS))
summary(lm(YeoCSOT ~ Sex, data = PVS))
summary(lm(YeoCSOT ~ MS_type, data = PVS))
summary(lm(YeoCSOT ~ MS_type+Age_Scan, data = PVS))
summary(lm(YeoCSOT ~ Age_Scan+Sex, data=PVS))
summary(lm(YeoCSOT ~ MS_type+Sex, data=PVS))
summary(lm(YeoCSOT ~ Age_Scan+Sex+MS_type, data = PVS))
# Y = intercept + B1*X1 + B2*X2
# CSOT = intercept + Estimate_Age_Scan*Age_Scan + Estimate_SexM*SexM
# The magnitude of the t statistics provides a means to judge relative importance of the independent variables (the higher, the more important)

#BGAT
summary(lm(logBGAT ~ Age_Scan, data = PVS))
summary(lm(logBGAT ~ Sex, data = PVS))
summary(lm(logBGAT ~ MS_type, data = PVS))
summary(lm(logBGAT ~ MS_type+Age_Scan, data = PVS))
summary(lm(logBGAT ~ Age_Scan+Sex, data=PVS))
summary(lm(logBGAT ~ MS_type+Sex, data=PVS))
summary(lm(logBGAT ~ Age_Scan+Sex+MS_type, data = PVS))

#BGCT
summary(lm(logBGCT ~ Sex, data = PVS))
summary(lm(logBGCT ~ Age_Scan, data = PVS))
summary(lm(logBGCT ~ MS_type, data = PVS))
summary(lm(logBGCT ~ Age_Scan+Sex, data=PVS))
summary(lm(logBGCT ~ MS_type+Sex, data=PVS))
summary(lm(logBGCT ~ MS_type+Age_Scan, data = PVS))
summary(lm(logBGCT ~ Age_Scan+Sex+MS_type, data = PVS))

#BS
summary(lm(BS ~ Sex, data = PVS))
summary(lm(BS ~ Age_Scan, data = PVS))
summary(lm(BS ~ MS_type, data = PVS))
summary(lm(BS ~ Age_Scan+Sex, data=PVS))
summary(lm(BS ~ MS_type+Sex, data=PVS))
summary(lm(BS ~ MS_type+Age_Scan, data = PVS))
summary(lm(BS ~ Age_Scan+Sex+MS_type, data = PVS))

summary(fit <- lm(log(CSOT) ~ Age_Scan+Sex, data=PVS))
fitted <- fitted.values(fit)

#fitted graph
ggplot(PVS, aes(x=fitted, y=CSOT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#How good does the model fit
#R-squared should be >0.7, F-statistics: the higher the better

#2.) sex
ggplot(PVS, aes(x=PVS$Sex, y=PVS$CSOT)) +
  geom_point(size=3, shape=20)

wilcox.test(PVS$CSOT~PVS$Sex) #non-parametric t-test (Mann Whitney U test)


#3.) MS phenotype
ggplot(PVS, aes(x=PVS$MS_type, y=PVS$CSOT)) +
  geom_point(size=3, shape=20)

kruskal.test(PVS$CSOT ~ PVS$MS_type) #non-paramatric ANOVA
pairwise.wilcox.test(PVS$CSOT, PVS$MS_type, p.adjust.method = "BH")
#non-parametric test with pairwise comparisons

#4.) SDMT
ggplot(PVS, aes(x=SDMT_Z, y=logCSOT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(PVS, aes(x=SDMT_Z, y=logBGCT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(PVS, aes(x=SDMT_Z, y=Duration)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

plot(PVS$SDMT_raw, PVS$SDMT_Z)
logDuration <- log(PVS$Duration)

#CSOT
summary(lm(SDMT_Z ~ logCSOT+Age_Scan+CortexVol, data=PVS))
summary(lm(SDMT_Z ~ logCSOT+Age_Scan+TotalGrayVol, data=PVS))

#BGCT
summary(lm(SDMT_Z ~ logBGCT+Duration+BrainSegVolNotVent, data=PVS))
summary(lm(SDMT_Z ~ logBGCT+Duration+CortexVol, data=PVS))

#BGAT
summary(lm(SDMT_Z ~ logBGAT+Duration+BrainSegVolNotVent, data=PVS))
summary(lm(SDMT_Z ~ logBGAT+Duration+CortexVol, data=PVS))

#BS
summary(lm(SDMT_Z ~ BS+Duration+BrainSegVolNotVent, data=PVS))
summary(lm(SDMT_Z ~ BS+Duration+CortexVol, data=PVS))





#5.) EDSS
ggplot(PVS, aes(x=PVS$EDSS, y=PVS$CSOT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(EDSS ~ logCSOT+Age_Scan+Sex+MS_type, data=PVS))

summary(lm(EDSS ~ logBGCT+Age_Scan+Sex+MS_type, data=PVS))

summary(lm(EDSS ~ logBGAT+Age_Scan+Sex+MS_type, data=PVS))

summary(lm(EDSS ~ BS+Age_Scan+Sex+MS_type, data=PVS))



#6.) Brain atrophy
#CSOT
summary(lm(logCSOT ~ SupraTentorialVolNotVent+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(CSOT ~ SupraTentorialVolNotVent+Sex, data=PVS))
summary(lm(CSOT ~ SupraTentorialVolNotVent+Age_Scan, data=PVS))
summary(lm(CSOT ~ SupraTentorialVolNotVent+Sex, data=PVS))
summary(lm(logCSOT ~ TotalGrayVol+Age_Scan+Sex, data=PVS))
summary(lm(CSOT ~ TotalGrayVol+Age_Scan+Sex, data=PVS))
summary(lm(logCSOT ~ CortexVol+Age_Scan+Sex, data=PVS))
summary(lm(CSOT ~ CortexVol+Age_Scan+Sex, data=PVS))
summary(lm(logCSOT ~ BrainSegVolNotVent+Age_Scan+Sex, data=PVS))
summary(lm(CSOT ~ BrainSegVolNotVent+Age_Scan+Sex, data=PVS))

summary(lm(CSOT ~ CortexVol+Age_Scan+Sex, data=PVS))
summary(lm(logCSOT ~ BrainSegVolNotVent+Age_Scan+Sex, data=PVS))
summary(lm(CSOT ~ BrainSegVolNotVent+Age_Scan+Sex, data=PVS))

summary(lm(SupraTentorialVolNotVent ~ logCSOT+Age_Scan+Sex+MS_type+Duration, data=PVS))
summary(lm(SupraTentorialVolNotVent ~ logCSOT+Age_Scan+Sex+Duration, data=PVS))
summary(lm(SupraTentorialVolNotVent ~ logCSOT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(SupraTentorialVolNotVent ~ Age_Scan+Duration, data=PVS))
summary(lm(SupraTentorialVolNotVent ~ Age_Scan+Sex+MS_type, data=PVS))

summary(lm(BrainSegVolNotVent ~ logCSOT+Age_Scan+Sex+MS_type+Duration, data=PVS))
summary(lm(BrainSegVolNotVent ~ logCSOT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(BrainSegVolNotVent ~ Age_Scan+Sex+MS_type, data=PVS))

summary(lm(CortexVol ~ logCSOT+Age_Scan+Sex+MS_type+Duration, data=PVS))
summary(lm(CortexVol ~ logCSOT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(CortexVol ~ Age_Scan+Sex+MS_type, data=PVS))

summary(lm(TotalGrayVol ~ logCSOT+Age_Scan+Sex+MS_type+Duration, data=PVS))
summary(lm(TotalGrayVol ~ logCSOT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(TotalGrayVol ~ Age_Scan+Sex+MS_type, data=PVS))

ggplot(PVS, aes(x=SupraTentorialVolNotVent, y=logCSOT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#New analysis august with longitudinal freesurfer data
ggplot(Long, aes(x=CSOT, y=WMF1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(WMF1 ~ CSOT+Age_Scan+Sex+MS_type+Duration, data=Long))

ggplot(Long, aes(x=CSOT, y=GMF1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(GMF1 ~ CSOT+Age_Scan+Sex+MS_type+Duration, data=Long))

ggplot(Long, aes(x=CSOT, y=BPF1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(BPF1 ~ CSOT+Age_Scan+Sex+MS_type+Duration, data=Long))

#BGCT
summary(lm(logBGCT ~ SupraTentorialVolNotVent+Age_Scan+Sex, data=PVS))
summary(lm(BGCT ~ SupraTentorialVolNotVent+Sex, data=PVS))
summary(lm(BGCT ~ SupraTentorialVolNotVent+Age_Scan, data=PVS))#
summary(lm(BGCT ~ SupraTentorialVolNotVent+Sex, data=PVS))
summary(lm(logBGCT ~ TotalGrayVol+Age_Scan+Sex+MS_type, data=PVS))##
summary(lm(BGCT ~ TotalGrayVol+Age_Scan+Sex, data=PVS))###
summary(lm(logBGCT ~ CortexVol+Age_Scan+Sex, data=PVS))##
summary(lm(BGCT ~ CortexVol+Age_Scan+Sex, data=PVS))###
summary(lm(logBGCT ~ BrainSegVolNotVent+Age_Scan+Sex, data=PVS))
summary(lm(BGCT ~ BrainSegVolNotVent+Age_Scan+Sex, data=PVS))

summary(lm(SupraTentorialVolNotVent ~ logBGCT+Age_Scan+Sex+MS_type+Duration, data=PVS))
summary(lm(SupraTentorialVolNotVent ~ logBGCT+Age_Scan+Sex+Duration, data=PVS))
summary(lm(SupraTentorialVolNotVent ~ logBGCT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(SupraTentorialVolNotVent ~ Age_Scan+Duration, data=PVS))
summary(lm(SupraTentorialVolNotVent ~ Age_Scan+Sex+MS_type, data=PVS))

summary(lm(BrainSegVolNotVent ~ logBGCT+Age_Scan+Sex+MS_type+Duration, data=PVS))
summary(lm(BrainSegVolNotVent ~ logBGCT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(BrainSegVolNotVent ~ Age_Scan+Sex+MS_type, data=PVS))

summary(lm(CortexVol ~ logBGCT+Age_Scan+Sex+MS_type+Duration, data=PVS))
summary(lm(CortexVol ~ logBGCT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(CortexVol ~ Age_Scan+Sex+MS_type, data=PVS))

summary(lm(TotalGrayVol ~ logBGCT+Age_Scan+Sex+MS_type+Duration, data=PVS))
summary(lm(TotalGrayVol ~ logBGCT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(TotalGrayVol ~ Age_Scan+Sex+MS_type, data=PVS))

ggplot(PVS, aes(x=TotalGrayVol, y=logBGCT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#BGAT
summary(lm(logBGAT ~ SupraTentorialVolNotVent+Age_Scan+Sex, data=PVS))##
summary(lm(BGAT ~ SupraTentorialVolNotVent+Sex, data=PVS))##
summary(lm(BGAT ~ SupraTentorialVolNotVent+Age_Scan, data=PVS))
summary(lm(logBGAT ~ TotalGrayVol+Age_Scan+Sex, data=PVS))##
summary(lm(BGAT ~ TotalGrayVol+Age_Scan+Sex, data=PVS))##
summary(lm(logBGAT ~ CortexVol+Age_Scan+Sex, data=PVS))
summary(lm(BGAT ~ CortexVol+Age_Scan+Sex, data=PVS))
summary(lm(logBGAT ~ BrainSegVolNotVent+Age_Scan+Sex+MS_type, data=PVS))##
summary(PVS$MS_type)

summary(lm(BGAT ~ MS_type, data=PVS))
str(PVS$MS_type)

summary(lm(SupraTentorialVolNotVent ~ logBGAT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(BrainSegVolNotVent ~ logBGAT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(CortexVol ~ logBGAT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(TotalGrayVol ~ logBGAT+Age_Scan+Sex+MS_type, data=PVS))

ggplot(PVS, aes(x=BrainSegVolNotVent, y=logBGAT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#BS
summary(lm(BS ~ SupraTentorialVolNotVent+Sex, data=PVS))
summary(lm(BS ~ SupraTentorialVolNotVent+Age_Scan, data=PVS))
summary(lm(BS ~ TotalGrayVol+Age_Scan+Sex, data=PVS))
summary(lm(BS ~ CortexVol+Age_Scan+Sex, data=PVS))
summary(lm(BS ~ BrainSegVolNotVent+Age_Scan+Sex, data=PVS))
#-> no sig results

summary(lm(SupraTentorialVolNotVent ~ BS+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(BrainSegVolNotVent ~ BS+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(CortexVol ~ BS+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(TotalGrayVol ~ BS+Age_Scan+Sex+MS_type, data=PVS))

#7.) T1/T2 lesions
#CSOT
summary(lm(T2_Count_FS ~ logCSOT+Age_Scan+Sex, data=PVS))
summary(lm(T2_Volume_FS ~ logCSOT+Age_Scan+Sex, data=PVS))

ggplot(Long, aes(x=CSOT, y=WMH1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(WMH1 ~ CSOT+Age_Scan+Sex+MS_type, data=Long))

#BGCT
summary(lm(T2_Count_FS ~ logBGCT+Age_Scan+Sex, data=PVS))
summary(lm(T2_Volume_FS ~ logBGCT+Age_Scan+Sex, data=PVS))

#BGAT
summary(lm(T2_Count_FS ~ logBGAT+Age_Scan+Sex, data=PVS))
summary(lm(T2_Volume_FS ~ logBGAT+Age_Scan+Sex, data=PVS))

#BS
summary(lm(T2_Count_FS ~ BS+Age_Scan+Sex, data=PVS))
summary(lm(T2_Volume_FS ~ BS+Age_Scan+Sex, data=PVS))

#7.) Relapses
summary(lm(Relapses ~ logCSOT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(Relapses ~ logBGCT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(Relapses ~ logBGAT+Age_Scan+Sex+MS_type, data=PVS))

ggplot(PVS, aes(x=Relapses, y=logCSOT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#8.) OCB
ggplot(PVS, aes(x=PVS$OCB, y=logCSOT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

t.test(logCSOT~PVS$OCB)
t.test(PVS$CSOT~PVS$OCB)
#no sig difference in OCB prevalence upon PVS increase

#9.) PVS Diameter (surrogate marker for PVS volume)
summary(lm(EDSS ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
summary(lm(SDMT_Z ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))

ggplot(PVS, aes(x=logCSO_D, y=SDMT_Z)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#same test as above with removal of score 3 and 4
Ribi <- log(Diametertest$CSO_D)
summary(lm(SDMT_Z ~ Ribi+Age_Scan+Sex+MS_type, data=Long))
#not any more significant after removing 3 and 4, what about only removing 4?
Ribi2 <- log(Diametertest2$CSO_D)
summary(lm(SDMT_Z ~ Ribi2+Age_Scan+Sex+MS_type, data=Long))

ggplot(PVS, aes(x=Ribi2, y=SDMT_Z)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(SDMT_Z ~ CSO_D_2+Sex+MS_type, data=Long)) #age can be removed from the model since SDMT_Z is already normalized to age
#is SDMT_Z normally distributed?
qqnorm(Long$SDMT_Z)
shapiro.test(Long$SDMT_Z) #p<0.05 indicates non-normal distributed dataset
install.packages("e1071")
library("e1071")
kurtosis(Long$SDMT_Z) #should be close to 0 for normal distribution
skewness(Long$SDMT_Z, na.rm=TRUE) #should be close to 0 for normal distribution
mean(Long$SDMT_Z, na.rm=TRUE)
median(Long$SDMT_Z, na.rm=TRUE) 
#tests suggest normal distribution of SDMT_Z

t.test(SDMT_Z~CSO_D_2, data=Long)
#There seems to be an outlier in the dilated PVS group, rerun analysis with excluding it
t.test(SDMT_Z_E~CSO_D_2, data=Long)
#What happens if you also exclude the top "outlier" in group 1?
t.test(SDMT_Z_E2~CSO_D_2, data=Long)
summary(lm(SDMT_Z_E2 ~ CSO_D_2+Sex+MS_type, data=Long))
summary(lm(SDMT_Z_E2 ~ CSO_D_2+Sex+MS_type+CortexVol, data=Long))
summary(lm(SDMT_Z_E2 ~ CSO_D_2+Sex+MS_type+TotalGrayVol, data=Long))
summary(lm(SDMT_Z_E2 ~ CSO_D_2+Sex+MS_type+SupraTentorialVolNotVent, data=Long))

wilcox.test(Long$SDMT_Z~Long$CSO_D_2) 

ggplot(Long, aes(x=CSO_D_2, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)


#what if you choose 3 mm as cut.off? only 4 patients in the upper group remain
ggplot(Long, aes(x=CSO_D_3, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

t.test(SDMT_Z_E2~CSO_D_3, data=Long)
wilcox.test(Long$SDMT_Z~Long$CSO_D_3) 


#dose-dependency?
summary(lm(SDMT_Z_E2 ~ CSO_D+Sex+MS_type, data=Long))
ggplot(Long, aes(x=CSO_D, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)
#problem: too few observations!

#longitudinal analysis
plot(Long$AUCnorm) #spot outliers
summary(lm(AUCnorm ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
#after excluding again 2 outliers (same as above)
summary(lm(AUCnorm_E ~ CSO_D_2+Age_Scan+Sex+MS_type, data=Long))
ggplot(Long, aes(x=CSO_D_2, y=AUCnorm_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(AUCnorm_E2 ~ CSO_D_2+Age_Scan+Sex+MS_type, data=Long))
ggplot(Long, aes(x=CSO_D_2, y=AUCnorm_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

t.test(AUCnorm_E~CSO_D_2, data=Long)
t.test(AUCnorm_E2~CSO_D_2, data=Long)

summary(lm(T2_Volume_FS ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
summary(lm(T2_Count_FS ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
summary(lm(T2_volume_D ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
summary(lm(T2_count_D ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))

summary(lm(BrainSegVolNotVent ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
summary(lm(CortexVol ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
summary(lm(TotalGrayVol ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
summary(lm(SupraTentorialVolNotVent ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Supratent_D ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=logCSO_D, y=Supratent_D)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Subcort_D ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=logCSO_D, y=Subcort_D)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Gray_D ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Cortex_D ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Brain_D ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))

summary(lm(Relapses ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))
summary(lm(OCB ~ logCSO_D+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=logCSO_D, y=OCB)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#10.) Potter versus counting





#Graphs
#1.) Age vs. logCSOT
ggplot(PVS, aes(x=Age_Scan, y=logCSOT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#2.) Sex vs. logCSOT
ggplot(PVS, aes(x=Sex, y=logCSOT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#3.) Sex vs. logBGCT
ggplot(PVS, aes(x=Sex, y=logBGCT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

Duration <- c(9.9, 21.3, 13.2, 11.2, 13, 15.2, 12.6, 15, 11.8, 14.3, 9.9, 7.6, 7.9, 11.4, 10.6, 11.2, 12.7, 9.9, 8.6, 9.2, 13.4, 6.8, 7.5, 14.6, 7.3, 10.6, 21.6, 7.7, 11.4, 10.1, 9.6, 22.8, 43, 13.7, 6.9, 7.7, 5.9, 17.8, 6.8, 9.4, 23.6, 7.2, 6.7, 15.9, 8.9, 15.5, 1.4, 6.6, 6.3, 12.6, 16.6, 11.5, 8.5, 6.2, 11.2, 9.1, 6.7, 7.1, 6.4, 16.9, 6.4, 4.9, 4.6, 28.7, 3.1, 6.6, 5.9, 4.5, 4.3, 4.2 , 2.6, 0.3, 3.6, 3, 3.3, 2.1, 2.4, 3.1, 0.6, 2.2, 4.3, 3.4, 5.1, 7.5, 7.3, 1.6, 4.6, 2.2, 2.3, 3.6, 10.8, 2.1, 2.6, 2.5, 19.3, 2.8, 1.6, 2.2, 2.9, 0.8, 0.1, 4, 3.2, 1.4, 0, 0.1, 1.9, 7.6, 0.4, 0.6, 13.4, 0, 2.1, 1.2, 0.6, 0.4, 2.2, 9.4, 1.3, 2.4, 0.5, 1.7, 0.7, 0.8, 1.3, 34.3, 0.1, 0.2, 1.6, 0.3, 1.5, 0.6, 14.7, 0.9, 0.1, 0.6, 0.1, 0, 2.7, 0, 0.3, 7.9)
PVS <- cbind(PVS, Duration)


#longitudinal analysis
#1.) Test for calculating AUC
AUC(x=c(1,3), y=c(1,1)) #DescTool package
plot(x=c(1,2,2.5), y=c(1,2,4), type="l", col="blue", ylim=c(0,4))
SDMTtest <- c(0, 0, 0, 1, 0, 0)
Timetest <- c(0, 5, 7, 12, 28, 35)
SDMT <- data.frame(SDMTtest, Timetest)

list(SDMT)
ggplot(SDMT, aes(x=SDMT$Timetest, y=SDMT$SDMTtest)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

AUC(x=SDMT$Timetest, y=SDMT$SDMTtest)
#conclusions: Did AUC calculation in excel

#2.) BL PVS -> SDMT AUC
summary(lm(AUCnorm ~ logCSOT+Age_SDMT+SDMT_number, data=Long))
summary(lm(AUCnorm ~ logBGCT+Age_SDMT+SDMT_number, data=Long))
summary(lm(AUCnorm ~ logBGAT+Age_SDMT+SDMT_number, data=Long))
summary(lm(AUCnorm ~ BS+Age_SDMT+SDMT_number, data=Long))

#age
ggplot(Long, aes(x=Age_SDMT, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#number of SDMT tests
ggplot(Long, aes(x=SDMT_number, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#logCSOT
ggplot(Long, aes(x=logCSOT, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#3.) BL PVS -> Delta atrophy
summary(lm(Supratent_D ~ logCSOT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Subcort_D ~ logCSOT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Gray_D ~ logCSOT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Cortex_D ~ logCSOT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Brain_D ~ logCSOT+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=logCSOT, y=Supratent_D)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)


ggplot(Long, aes(x=CSO_D, y=Supratent_D)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Supratent_D ~ logBGCT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Subcort_D ~ logBGCT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Gray_D ~ logBGCT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Cortex_D ~ logBGCT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Brain_D ~ logBGCT+Age_Scan+Sex+MS_type, data=Long))

summary(lm(Supratent_D ~ logBGAT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Subcort_D ~ logBGAT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Gray_D ~ logBGAT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Cortex_D ~ logBGAT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Brain_D ~ logBGAT+Age_Scan+Sex+MS_type, data=Long))

summary(lm(Supratent_D ~ BS+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Subcort_D ~ BS+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Gray_D ~ BS+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Cortex_D ~ BS+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Brain_D ~ BS+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=BS, y=Cortex_D)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#4.) BL PVS -> Delta lesions
summary(lm(T2_volume_D ~ logCSOT+Age_Scan+Scan_interval, data=Long))
summary(lm(T2_count_D ~ logCSOT+Age_Scan+Scan_interval, data=Long))

summary(lm(T2_volume_D ~ logBGCT+Age_Scan+Scan_interval, data=Long))
summary(lm(T2_count_D ~ logBGCT+Age_Scan+Scan_interval, data=Long))

summary(lm(T2_volume_D ~ logBGAT+Age_Scan+Scan_interval, data=Long))
summary(lm(T2_count_D ~ logBGAT+Age_Scan+Scan_interval, data=Long))

summary(lm(T2_volume_D ~ BS+Age_Scan+Scan_interval, data=Long))
summary(lm(T2_count_D ~ BS+Age_Scan+Scan_interval, data=Long))

ggplot(Long, aes(x=logCSOT, y=T2_volume_D)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=logCSOT, y=T2_count_D)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)


#Secondary analysis with re-measuring PVS dilation
str(Long$CSO_D_Redone)
str(Long$CSO_D_Number)

  #Descriptive part (CSO)
ggplot(data=subset(Long, !is.na(CSO_D_Redone)),
       aes(x=Sex, y=CSO_D_Redone)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(CSO_D_Number)),
       aes(x=Sex, y=CSO_D_Number)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

t.test(CSO_D_Redone~Sex, data=Long)
t.test(CSO_D_Number~Sex, data=Long)

ggplot(data=subset(Long, !is.na(CSO_D_Redone)),
       aes(x=Age_Scan, y=CSO_D_Redone)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(CSO_D_Number)),
       aes(x=Age_Scan, y=CSO_D_Number)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(CSO_D_Redone ~ Sex+Age_Scan, data=Long))
summary(lm(CSO_D_Number ~ Sex+Age_Scan, data=Long))


#1.) SDMT
#a.) Cross-sectional
ggplot(Long, aes(x=CSO_D_Redone, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(SDMT_Z_E2 ~ CSO_D_Redone+Sex+MS_type, data=Long))

Long3 <- Long %>% filter(!is.na(SDMT_Z_E2)) %>%
  filter(!is.na(CSO_D_Redone)) %>%
  filter(!is.na(Sex)) %>%
  filter(!is.na(MS_type))

SDMT_diameter <- lm(SDMT_Z_E2 ~ CSO_D_Redone+Sex+MS_type, data=Long)

SDMT_diameter_fitted <- fitted.values(SDMT_diameter)

ggplot(Long3, aes(x=CSO_D_Redone, y=SDMT_diameter_fitted)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 14),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 14)) +
        xlab("EPVS diameter") +
        ylab("SDMT z-score fitted") +
        theme_minimal(base_size = 15)

summary(lm(SDMT_Z_E2 ~ CSO_D_Redone+Sex+MS_type+SupraTentorialVolNotVent, data=Long))

summary(lm(SDMT_Z_E2 ~ CSO_D_Redone+Sex+MS_type+T2_Count_FS_E+SupraTentorialVolNotVent, data=Long))
summary(lm(SDMT_Z_E2 ~ CSO_D_Redone+Sex+MS_type+T2_Volume_FS_E+SupraTentorialVolNotVent, data=Long))
summary(lm(SDMT_Z_E2 ~ CSO_D_Redone+Sex+MS_type+CortexVol, data=Long))
summary(lm(SDMT_Z_E2 ~ CSO_D_Redone+Sex+MS_type+TotalGrayVol, data=Long))
summary(lm(SDMT_Z_E2 ~ CSO_D_Redone+Sex+MS_type+SupraTentorialVolNotVent, data=Long))

#with 1.5 mm as cutoff
ggplot(Long, aes(x=CSO_D_Bin, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

t.test(SDMT_Z_E2~CSO_D_Bin, data=Long)
wilcox.test(Long$SDMT_Z~Long$CSO_D_Bin) 

#with 2 mm as cutoff
ggplot(Long, aes(x=CSO_D_Bin2, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

t.test(SDMT_Z_E2~CSO_D_Bin2, data=Long)
wilcox.test(Long$SDMT_Z~Long$CSO_D_Bin) 

#number of dilated PVS
ggplot(Long, aes(x=CSO_D_Number, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(SDMT_Z_E2 ~ CSO_D_Number+Sex+MS_type, data=Long))
summary(lm(SDMT_Z_E2 ~ CSO_D_Number+Sex+MS_type+CortexVol, data=Long))
summary(lm(SDMT_Z_E2 ~ CSO_D_Number+Sex+MS_type+TotalGrayVol, data=Long))
summary(lm(SDMT_Z_E2 ~ CSO_D_Number+Sex+MS_type+SupraTentorialVolNotVent, data=Long))

summary(lm(SDMT_Z_E2 ~ CSO_D_Bin3+Sex+MS_type, data=Long))

Long4 <- Long %>% filter(!is.na(SDMT_Z_E2)) %>%
  filter(!is.na(CSO_D_Bin3)) %>%
  filter(!is.na(Sex)) %>%
  filter(!is.na(MS_type))

SDMT_binned <- lm(SDMT_Z_E2 ~ CSO_D_Bin3+Sex+MS_type, data=Long)

SDMT_binned_fitted <- fitted.values(SDMT_binned)

ggplot(Long4, aes(x=CSO_D_Bin3, y=SDMT_binned_fitted)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 14),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 14)) +
  xlab("EPVS diameter") +
  ylab("SDMT z-score fitted") +
  theme_minimal(base_size = 15)

#b.) SDMT longitudinal
summary(lm(AUCnorm ~ CSO_D_Number+Age_SDMT+SDMT_number, data=Long))
summary(lm(AUCnorm ~ CSO_D_Redone+Age_SDMT+SDMT_number, data=Long))

ggplot(Long, aes(x=CSO_D_Redone, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#2.) Atrophy
#a.) Cross-sectional
summary(lm(BrainSegVolNotVent ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))
summary(lm(BrainSegVolNotVent ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))

summary(lm(CortexVol ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))
summary(lm(CortexVol ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))

summary(lm(TotalGrayVol ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))
summary(lm(TotalGrayVol ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))

summary(lm(SupraTentorialVolNotVent ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))
summary(lm(SupraTentorialVolNotVent ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Redone, y=BrainSegVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

 #with new freesurfer data
ggplot(Long, aes(x=CSO_D_Redone, y=BPF1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=CSO_D_Number, y=BPF1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(BPF1 ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))
summary(lm(BPF1 ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Number, y=WMF1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=CSO_D_Number, y=GMF1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#b.) longitudinal
ggplot(Long, aes(x=CSOT, y=Delta_WMF)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=CSOT, y=Delta_WMF_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_WMF ~ CSOT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Delta_WMF_E1 ~ CSOT+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSOT, y=Delta_GMF)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=CSOT, y=Delta_GMF_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_GMF_E1 ~ CSOT+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSOT, y=Delta_BPF)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_GMF_E1 ~ CSOT+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSOT, y=Delta_BPF_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_BPF_E1 ~ CSOT+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Number, y=Delta_BPF_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_BPF_E1 ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Number, y=Delta_WMF_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=CSO_D_Number, y=Delta_GMF_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_GMF_E1 ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Redone, y=Delta_BPF_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_BPF_E1 ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))

#3.) Lesions
#a.) Cross-sectional

Long5 <- merge(Long, CVRF, by="HIVE")

summary(lm(T2_Volume_FS_E ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))
summary(lm(T2_Count_FS_E ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))
summary(lm(T2_Volume_FS_E ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))
summary(lm(T2_Count_FS_E ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))

summary(lm(T2_Volume_FS_E ~ CSO_D_Number+Age_Scan+Sex+MS_type+Sum_CVRF, data=Long5))#also adjusted for CVRF


ggplot(Long, aes(x=CSO_D_Redone, y=T2_Volume_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=CSO_D_Redone, y=T2_Count_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=CSO_D_Number, y=T2_Volume_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=CSO_D_Number, y=T2_Count_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

CVRF <- read_excel("G:/MS/Experiments/Karolinska/perivascular spaces in MS/Statistics/PVS/CVRF.xlsx", sheet = 1)#CVRF data set load
names(CVRF) <- c("HIVE", "Individual_Timepoint", "Diabetes", "Hypertension", "Dyslipidemia", "Antiplatelet", "Sum_CVRF")

Long2 <- Long %>% filter(!is.na(T2_Volume_FS_E)) %>%
  filter(!is.na(T2_Count_FS_E)) %>%
  filter(!is.na(CSO_D_Redone)) %>%
  filter(!is.na(CSO_D_Number)) %>%
  filter(!is.na(Age_Scan)) %>%
  filter(!is.na(Sex)) %>%
  filter(!is.na(MS_type))

Long2 <- merge(Long2, CVRF, by="HIVE")#merge Long2 with CVRF

Long4 <- Long %>% filter(!is.na(WMH1_E1)) %>%
  filter(!is.na(T2_Count_FS_E)) %>%
  filter(!is.na(CSO_D_Redone)) %>%
  filter(!is.na(CSO_D_Number)) %>%
  filter(!is.na(Age_Scan)) %>%
  filter(!is.na(Sex)) %>%
  filter(!is.na(MS_type))

Long4 <- merge(Long4, CVRF, by="HIVE")#merge Long4 with CVRF

Vol_diameter <- lm(T2_Volume_FS_E ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long2)
Cou_diameter <- lm(T2_Count_FS_E ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long)
Vol_number <- lm(T2_Volume_FS_E ~ CSO_D_Number+Age_Scan+Sex+MS_type+Sum_CVRF, data=Long2)
Cou_number <- lm(T2_Count_FS_E ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long)

Vol_diameter_fitted <- fitted.values(Vol_diameter)
Cou_diameter_fitted <- fitted.values(Cou_diameter)
Vol_number_fitted <- fitted.values(Vol_number)
Cou_number_fitted <- fitted.values(Cou_number)

WMH_number <- lm(WMH1_E1 ~ CSO_D_Number+Age_Scan+Sex+MS_type+Sum_CVRF, data=Long4)
WMH_fitted <- fitted.values(WMH_number)
WMH_fitted1000 <- WMH_fitted/1000

ggplot(Long2, aes(x=CSO_D_Redone, y=Vol_diameter_fitted)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long2, aes(x=CSO_D_Redone, y=Cou_diameter_fitted)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long2, aes(x=CSO_D_Number, y=Vol_number_fitted)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long2, aes(x=CSO_D_Number, y=Cou_number_fitted)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#  -->unadjusted analysis
summary(lm(T2_Volume_FS_E ~ CSO_D_Redone, data=Long))
summary(lm(T2_Count_FS_E ~ CSO_D_Redone, data=Long))
summary(lm(T2_Volume_FS_E ~ CSO_D_Number, data=Long))
summary(lm(T2_Count_FS_E ~ CSO_D_Number, data=Long))

#With new freesurfer data
ggplot(Long, aes(x=CSO_D_Number, y=WMH1_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(WMH1_E1 ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))
summary(lm(WMH1_E1 ~ CSO_D_Number, data=Long))

ggplot(Long, aes(x=CSO_D_Redone, y=WMH1_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(WMH1_E1 ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))
summary(lm(WMH1_E1 ~ CSO_D_Redone, data=Long))

ggplot(Long, aes(x=T2_Volume_FS_E, y=WMH1_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

cor.test(Long$T2_Volume_FS_E, Long$WMH1_E1)

#b.) Longitudinal
ggplot(Long, aes(x=CSOT, y=Delta_T2_Volume)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_T2_Volume ~ CSOT+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSOT, y=Delta_WMH)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_WMH ~ CSOT+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSOT, y=Delta_WMH)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_WMH ~ CSOT+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Number, y=Delta_WMH)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_WMH ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Number, y=Delta_T2_Volume)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_WMH ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Redone, y=Delta_T2_Volume)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#4.) EDSS
ggplot(Long, aes(x=CSO_D_Redone, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(EDSS ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Number, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(EDSS ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))

#5.) Relapses
ggplot(Long, aes(x=CSO_D_Redone, y=Relapses)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Relapses ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Number, y=Relapses)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Relapses ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))

#6.) MS type
ggplot(Long, aes(x=MS_type, y=CSO_D_Redone)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(MS_type ~ CSO_D_Redone+Age_Scan+Sex, data=Long))

ggplot(Long, aes(x=MS_type, y=CSO_D_Number)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Relapses ~ CSO_D_Number+Age_Scan+Sex, data=Long))


#Longitudinal changes of EPVS dilation
#1.) descriptive
summary(Long$CSO_D_Number)

#2.) SDMT
ggplot(Long, aes(x=Bin_Diam, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Bin_Diam2, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

t.test(SDMT_Z_E2~Bin_Diam, data=Long)

ggplot(Long, aes(x=Bin_Diam, y=AUCnorm_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#3.) Brain volume
ggplot(Long, aes(x=Bin_Diam, y=SupraTentorialVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Bin_Diam, y=CortexVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Bin_Diam, y=TotalGrayVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#4.) T2 lesion
ggplot(Long, aes(x=Bin_Diam, y=T2_Volume_FS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Bin_Diam, y=T2_Count_FS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#5.) EDSS
ggplot(Long, aes(x=Bin_Diam, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

t.test(EDSS~Bin_Diam, data=Long)

#6.) MS type
ggplot(Long, aes(x=MS_type, y=Bin_Diam)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#7.) Relapses
ggplot(Long, aes(x=Bin_Diam, y=Relapses)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

t.test(Relapses~Bin_Diam, data=Long)

#Exploration Tysabri
#Hypothesis: MS patients who had Tysabri have less dilated PVS

ggplot(Long, aes(x=Tysabri, y=CSO_D_Redone)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Tysabri, y=CSO_D_Number)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)


#manual segmentation for PVS volume assessment (CSO and BG)
Str(Long$Vol_PVS_CSO)
Str(Long$Vol_PVS_BG)
summary(Long$Vol_PVS_CSO)
summary(Long$Vol_PVS_BG)

  #Descriptive part
Vol_men <- MS_sub <- subset(Long, Sex=="M")
Vol_women <- MS_sub <- subset(Long, Sex=="F")

    #PVS CSO Vol
ggplot(data=subset(Long, !is.na(Vol_PVS_CSO)),
       aes(x=Sex, y=Vol_PVS_CSO)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

t.test(Vol_PVS_CSO~Sex, data=Long)

ggplot(data=subset(Long, !is.na(Vol_PVS_CSO)),
       aes(x=Age_Scan, y=Vol_PVS_CSO)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Vol_PVS_CSO ~ Sex+Age_Scan, data=Long))

    #PVS BG Vol
ggplot(data=subset(Long, !is.na(Vol_PVS_BG)),
       aes(x=Sex, y=Vol_PVS_BG)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

t.test(Vol_PVS_BG~Sex, data=Long)

ggplot(data=subset(Long, !is.na(Vol_PVS_BG)),
       aes(x=Age_Scan, y=Vol_PVS_BG)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Vol_PVS_BG ~ Sex+Age_Scan, data=Long))

#1.) SDMT
str(Long$Deltavol_CSO)
str(Long$Deltavol_BG_bin)
str(Long$Vol_PVS_CSO)

#Cross-sectional PVS volume and cross-sectional SDMT
  #CSO
ggplot(Long, aes(x=Vol_PVS_CSO, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(SDMT_Z_E2 ~ Vol_PVS_CSO+Sex+MS_type+SupraTentorialVolNotVent, data=Long))
summary(lm(SDMT_Z_E2 ~ Vol_PVS_CSO+Sex+SupraTentorialVolNotVent, data=Long))

  #BG
ggplot(Long, aes(x=Vol_PVS_BG, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(SDMT_Z_E2 ~ Vol_PVS_BG+Sex+MS_type+SupraTentorialVolNotVent, data=Long))

  #Total
ggplot(Long, aes(x=Vol_PVS_tot, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(SDMT_Z_E2 ~ Vol_PVS_tot+Sex+MS_type+SupraTentorialVolNotVent, data=Long))

#Cross-sectional PVS volume and longitudinal SDMT (AUC)
ggplot(Long, aes(x=Vol_PVS_CSO, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(AUCnorm ~ Vol_PVS_CSO+Age_SDMT+SDMT_number, data=Long))

ggplot(Long, aes(x=Vol_PVS_BG, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(AUCnorm ~ Vol_PVS_BG+Age_SDMT+SDMT_number, data=Long))

ggplot(Long, aes(x=Vol_PVS_tot, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(AUCnorm ~ Vol_PVS_tot+Age_SDMT+SDMT_number, data=Long))

#Longitudinal PVS volume change and cross-sectional SDMT
  #with absolute volume changes
ggplot(Long, aes(x=Deltavol_CSO, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Deltavol_BG, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Deltavol_tot, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

  #with binned results
ggplot(data=subset(Long, !is.na(Deltavol_CSO_bin)),
       aes(x=Deltavol_CSO_bin, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_BG_bin)),
       aes(x=Deltavol_BG_bin, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_tot_bin)),
       aes(x=Deltavol_tot_bin, y=SDMT_Z_E2)) +
    geom_point(size=3, shape=20) +
    geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(aov(SDMT_Z_E2 ~ Deltavol_BG_bin, data = Long))

#Longitudinal PVS volume change and longitudinal SDMT change (AUC)
  #with absolute volume changes
ggplot(Long, aes(x=Deltavol_CSO, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Deltavol_BG, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Deltavol_tot, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

  #binned
ggplot(data=subset(Long, !is.na(Deltavol_CSO_bin)),
       aes(x=Deltavol_CSO_bin, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_BG_bin)),
       aes(x=Deltavol_BG_bin, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_tot_bin)),
       aes(x=Deltavol_tot_bin, y=AUCnorm)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(aov(AUCnorm ~ Deltavol_CSO_bin, data = Long))
summary(aov(AUCnorm ~ Deltavol_BG_bin, data = Long))

#2.) brain volume
#cross-sectional PVS volume with brain volume
  #CSO
ggplot(Long, aes(x=Vol_PVS_CSO, y=BrainSegVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_CSO, y=CortexVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_CSO, y=TotalGrayVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_CSO, y=SupraTentorialVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(BrainSegVolNotVent ~ Vol_PVS_CSO+Age_Scan+Sex+MS_type, data=Long))
summary(lm(CortexVol ~ Vol_PVS_CSO+Age_Scan+Sex+MS_type, data=Long))
summary(lm(TotalGrayVol ~ Vol_PVS_CSO+Age_Scan+Sex+MS_type, data=Long))
summary(lm(SupraTentorialVolNotVent ~ Vol_PVS_CSO+Age_Scan+Sex+MS_type, data=Long))

  #after excluding subject with high PVS volume and brain volume
ggplot(Long, aes(x=Vol_PVS_CSO_E1, y=BrainSegVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(BrainSegVolNotVent ~ Vol_PVS_CSO_E1+Age_Scan+Sex+MS_type, data=Long))
summary(lm(CortexVol ~ Vol_PVS_CSO_E1+Age_Scan+Sex+MS_type, data=Long))
summary(lm(TotalGrayVol ~ Vol_PVS_CSO_E1+Age_Scan+Sex+MS_type, data=Long))
summary(lm(SupraTentorialVolNotVent ~ Vol_PVS_CSO_E1+Age_Scan+Sex+MS_type, data=Long))

  #BG
ggplot(Long, aes(x=Vol_PVS_BG, y=BrainSegVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_BG, y=CortexVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_BG, y=TotalGrayVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_BG, y=SupraTentorialVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(BrainSegVolNotVent ~ Vol_PVS_BG+Age_Scan+Sex+MS_type, data=Long))
summary(lm(CortexVol ~ Vol_PVS_BG+Age_Scan+Sex+MS_type, data=Long))
summary(lm(TotalGrayVol ~ Vol_PVS_BG+Age_Scan+Sex+MS_type, data=Long))
summary(lm(SupraTentorialVolNotVent ~ Vol_PVS_BG+Age_Scan+Sex+MS_type, data=Long))

  #Total
ggplot(Long, aes(x=Vol_PVS_tot, y=BrainSegVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_tot, y=CortexVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_tot, y=TotalGrayVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_tot, y=SupraTentorialVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#Longitudinal changes in PVS volume with brain volume (binned)
  #CSO
ggplot(data=subset(Long, !is.na(Deltavol_CSO_bin)),
       aes(x=Deltavol_CSO_bin, y=BrainSegVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_CSO_bin)),
       aes(x=Deltavol_CSO_bin, y=CortexVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_CSO_bin)),
       aes(x=Deltavol_CSO_bin, y=TotalGrayVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_CSO_bin)),
       aes(x=Deltavol_CSO_bin, y=SupraTentorialVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

  #BG
ggplot(data=subset(Long, !is.na(Deltavol_BG_bin)),
       aes(x=Deltavol_BG_bin, y=BrainSegVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_BG_bin)),
       aes(x=Deltavol_BG_bin, y=CortexVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_BG_bin)),
       aes(x=Deltavol_BG_bin, y=TotalGrayVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_BG_bin)),
       aes(x=Deltavol_BG_bin, y=SupraTentorialVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(aov(BrainSegVolNotVent ~ Deltavol_BG_bin, data = Long))
summary(aov(CortexVol ~ Deltavol_BG_bin, data = Long))
summary(aov(TotalGrayVol ~ Deltavol_BG_bin, data = Long))
summary(aov(SupraTentorialVolNotVent ~ Deltavol_BG_bin, data = Long))

  #total
ggplot(data=subset(Long, !is.na(Deltavol_tot_bin)),
       aes(x=Deltavol_tot_bin, y=BrainSegVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_tot_bin)),
       aes(x=Deltavol_tot_bin, y=CortexVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_tot_bin)),
       aes(x=Deltavol_tot_bin, y=TotalGrayVol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_tot_bin)),
       aes(x=Deltavol_tot_bin, y=SupraTentorialVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(aov(BrainSegVolNotVent ~ Deltavol_tot_bin, data = Long))
summary(aov(CortexVol ~ Deltavol_tot_bin, data = Long))
summary(aov(TotalGrayVol ~ Deltavol_tot_bin, data = Long))
summary(aov(SupraTentorialVolNotVent ~ Deltavol_tot_bin, data = Long))

  #absolute changes (not binned)
ggplot(Long, aes(x=Deltavol_CSO, y=BrainSegVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Deltavol_BG, y=BrainSegVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Deltavol_tot, y=BrainSegVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#with new longitudinal freesurfer data
ggplot(Long, aes(x=Vol_PVS_CSO_E1, y=BPF1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(BPF1 ~ Vol_PVS_CSO_E1+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=Vol_PVS_CSO_E1, y=WMF1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_CSO_E1, y=GMF1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(GMF1 ~ Vol_PVS_CSO_E1+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=Vol_PVS_CSO_E1, y=Delta_BPF_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_BPF_E1 ~ Vol_PVS_CSO_E1+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=Vol_PVS_CSO_E1, y=Delta_WMF_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_CSO_E1, y=Delta_GMF_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#3.) Lesions
#Cross-sectional PVS volume with lesion measures
  #CSO
ggplot(Long, aes(x=Vol_PVS_CSO, y=T2_Volume_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_CSO, y=T2_Count_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(T2_Volume_FS_E ~ Vol_PVS_CSO+Age_Scan+Sex+MS_type, data=Long))
summary(lm(T2_Count_FS_E ~ Vol_PVS_CSO+Age_Scan+Sex+MS_type, data=Long))

  #BG
ggplot(Long, aes(x=Vol_PVS_BG, y=T2_Volume_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_BG, y=T2_Count_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(T2_Volume_FS_E ~ Vol_PVS_BG+Age_Scan+Sex+MS_type, data=Long))
summary(lm(T2_Count_FS_E ~ Vol_PVS_BG+Age_Scan+Sex+MS_type, data=Long))

  #total
ggplot(Long, aes(x=Vol_PVS_tot, y=T2_Volume_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_tot, y=T2_Count_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(T2_Volume_FS_E ~ Vol_PVS_tot+Age_Scan+Sex+MS_type, data=Long))
summary(lm(T2_Count_FS_E ~ Vol_PVS_tot+Age_Scan+Sex+MS_type, data=Long))

#Longitudinal PVS volume changes
  #CSO
ggplot(data=subset(Long, !is.na(Deltavol_CSO_bin)),
       aes(x=Deltavol_CSO_bin, y=T2_Volume_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_CSO_bin)),
       aes(x=Deltavol_CSO_bin, y=T2_Count_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

  #BG
ggplot(data=subset(Long, !is.na(Deltavol_BG_bin)),
       aes(x=Deltavol_BG_bin, y=T2_Volume_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_BG_bin)),
       aes(x=Deltavol_BG_bin, y=T2_Count_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

  #total
ggplot(data=subset(Long, !is.na(Deltavol_tot_bin)),
       aes(x=Deltavol_tot_bin, y=T2_Volume_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_tot_bin)),
       aes(x=Deltavol_tot_bin, y=T2_Count_FS_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#with new freesurfer data from longitudinal assessment
ggplot(Long, aes(x=Vol_PVS_CSO_E1, y=WMH1_E1)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(WMH1_E1 ~ Vol_PVS_CSO_E1+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=Vol_PVS_CSO_E1, y=Delta_T2_Volume)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Delta_T2_Volume ~ Vol_PVS_CSO_E1+Age_Scan+Sex+MS_type, data=Long))

#4.) EDSS
#cross-sectional PVS
ggplot(Long, aes(x=Vol_PVS_CSO, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(EDSS ~ Vol_PVS_CSO+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=Vol_PVS_BG, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(EDSS ~ Vol_PVS_BG+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=Vol_PVS_tot, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(EDSS ~ Vol_PVS_tot+Age_Scan+Sex+MS_type, data=Long))

  #longitudinal PVS changes
ggplot(data=subset(Long, !is.na(Deltavol_CSO_bin)),
       aes(x=Deltavol_CSO_bin, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_BG_bin)),
       aes(x=Deltavol_BG_bin, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_tot_bin)),
       aes(x=Deltavol_tot_bin, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#5.) Relapses
  #cross-sectional PVS
ggplot(Long, aes(x=Vol_PVS_CSO, y=Relapses)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_BG, y=Relapses)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=Vol_PVS_tot, y=Relapses)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

  #Longitudinal PVS changes
ggplot(data=subset(Long, !is.na(Deltavol_CSO_bin)),
       aes(x=Deltavol_CSO_bin, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_BG_bin)),
       aes(x=Deltavol_BG_bin, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_tot_bin)),
       aes(x=Deltavol_tot_bin, y=EDSS)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#6.) MS type
  #Cross-sectional
ggplot(Long, aes(x=MS_type, y=Vol_PVS_CSO)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=MS_type, y=Vol_PVS_BG)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=MS_type, y=Vol_PVS_tot)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(aov(Vol_PVS_CSO ~ MS_type, data = Long))
summary(aov(Vol_PVS_BG ~ MS_type, data = Long))
summary(aov(Vol_PVS_tot ~ MS_type, data = Long))

  #Longitudinal PVS changes
ggplot(data=subset(Long, !is.na(Deltavol_CSO_bin)),
       aes(x=Deltavol_CSO_bin, y=MS_type)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_BG_bin)),
       aes(x=Deltavol_BG_bin, y=MS_type)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(data=subset(Long, !is.na(Deltavol_tot_bin)),
       aes(x=Deltavol_tot_bin, y=MS_type)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

MS_sub <- subset(Long, MS_type==1|MS_type==2)

#Inter-readout correlation
  #Count vs. Volume
ggplot(Long, aes(x=CSOT, y=Vol_PVS_CSO)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=BGCT, y=Vol_PVS_BG)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=BGAT, y=Vol_PVS_BG)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=BGACT, y=Vol_PVS_BG)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=BGCT+BGAT, y=Vol_PVS_BG)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

cor.test(Long$CSOT, Long$Vol_PVS_CSO)
cor.test(Long$BGCT, Long$Vol_PVS_BG)
cor.test(Long$BGAT, Long$Vol_PVS_BG)
cor.test(Long$BGACT, Long$Vol_PVS_BG)

  #Dilation versus volume
ggplot(Long, aes(x=CSO_D_Number, y=Vol_PVS_CSO)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

cor.test(Long$CSO_D_Number, Long$Vol_PVS_CSO)

ggplot(Long, aes(x=BGA_D, y=Vol_PVS_BG)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

cor.test(Long$BGA_D, Long$Vol_PVS_BG)

ggplot(Long, aes(x=BGC_D, y=Vol_PVS_BG)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

cor.test(Long$BGC_D, Long$Vol_PVS_BG)


summary(Long$Deltavol_CSO_bin)
summary(Long$Deltavol_BG_bin)

summary(Long$Vol_PVS_CSO)
summary(Long$Vol_PVS_CSO_2)
summary(Long$Vol_PVS_CSO_3)

summary(Long$Vol_PVS_BG)
summary(Long$Vol_PVS_BG_2)
summary(Long$Vol_PVS_BG_3)

#Association of PVS to Gd/Gad lesions
#1.) PVS numbers to Gd lesions
ggplot(Long, aes(x=CSOT, y=Gd1_vol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Gd1_vol ~ CSOT, data=Long))
summary(lm(Gd1_vol ~ CSOT+Age_Scan+Sex+MS_type, data=Long))

#2.) PVS volume to Gd lesions
ggplot(Long, aes(x=Vol_PVS_CSO, y=Gd1_vol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Gd1_vol ~ Vol_PVS_CSO+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=Vol_PVS_CSO, y=Gd1_bin)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

#3.) PVS diameter to Gd lesions
ggplot(Long, aes(x=CSO_D_Redone, y=Gd1_vol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Gd1_vol ~ CSO_D_Redone+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Number, y=Gd1_vol)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Gd1_vol ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))

ggplot(Long, aes(x=CSO_D_Number, y=Gd1_bin)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

ggplot(Long, aes(x=CSO_D_Number, y=Gd1_number)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)


ggplot(Long, aes(x=T2_Volume_FS_E, y=Gd1_vol_E)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(Gd1_vol_E ~ T2_Volume_FS_E+Age_Scan+Sex+MS_type, data=Long))

#additional figures
  #EPVS number vs T2
Number_T2 <- ggplot(Long2, aes(x=CSOT, y=Vol_number_fitted)) +
  geom_point(size=2, shape=16, color="blue") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 14),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 14)) +
  ylab("T2 lesion volume (ml)") +
  xlab("EPVS Number") +
  scale_y_continuous(name="T2 lesion volume (ml)", limits=c(0, 15))

  #EPVS volume vs T2
Volume_T2 <- ggplot(Long2, aes(x=Vol_PVS_CSO, y=Vol_number_fitted)) +
  geom_point(size=2, shape=16, color="blue") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 14),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 14)) +
  ylab("T2 lesion volume (ml)") +
  xlab("EPVS Volume (ul)") +
  scale_y_continuous(name="T2 lesion volume (ml)", limits=c(0, 15))

  #EPVS Number vs brain Volume
Number_brainvolume <- ggplot(Long2, aes(x=CSOT, y=BPF1)) +
  geom_point(size=2, shape=16, color="blue") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 14),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 14), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 14)) +
  xlab("EPVS Number") +
  scale_y_continuous(name="Brain parenchymal fraction")

  #EPVS volume vs brain volume
Volume_brainvolume <- ggplot(Long2, aes(x=Vol_PVS_CSO, y=BPF1)) +
  geom_point(size=2, shape=16, color="blue") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 14),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 14), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 14)) +
  xlab("EPVS Volume (ul)") +
  scale_y_continuous(name="Brain parenchymal fraction")

  #EPVS diameter vs brain volume
Diameter_brainvolume <- ggplot(Long2, aes(x=CSO_D_Number, y=BPF1)) +
  geom_point(size=2, shape=16, color="blue") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 14),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 14), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 14)) +
  xlab("#of dilated EPVS") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_continuous(name="Brain parenchymal fraction")

##########Validation cohort#############3

###Freesurfer output
#load and merge files
Seg_MultipleMS <- read.table("aseg_stats_multiplems.txt", sep="\t", header=TRUE)
Seg_MyelinMS <- read.table("aseg_stats_REMY.txt", sep="\t", header=TRUE)
Seg_Validation <- rbind(Seg_MultipleMS, Seg_MyelinMS)
#prepare variables
Seg_Validation <- subset(Seg_Validation, select=c("Measure.volume", "Left.Cerebellum.White.Matter", "Right.Cerebellum.White.Matter", 
                                                  "CerebralWhiteMatterVol", "Brain.Stem",
                                                  "CC_Posterior", "CC_Mid_Posterior", "CC_Central", "CC_Mid_Anterior", "CC_Anterior", 
                                                  "BrainSegVolNotVent", "TotalGrayVol", "EstimatedTotalIntraCranialVol", 
                                                  "WM.hypointensities"))
Seg_Validation$WM <- Seg_Validation$Left.Cerebellum.White.Matter+
                      Seg_Validation$Right.Cerebellum.White.Matter+
                      Seg_Validation$CerebralWhiteMatterVol+
                      Seg_Validation$Brain.Stem+
                      Seg_Validation$CC_Posterior+Seg_Validation$CC_Mid_Posterior+
                      Seg_Validation$CC_Central+Seg_Validation$CC_Mid_Anterior+Seg_Validation$CC_Anterior
Seg_Validation$BPF <- Seg_Validation$BrainSegVolNotVent/Seg_Validation$EstimatedTotalIntraCranialVol
Seg_Validation$GMF <- Seg_Validation$TotalGrayVol/Seg_Validation$EstimatedTotalIntraCranialVol
Seg_Validation$WMF <- Seg_Validation$WM/Seg_Validation$EstimatedTotalIntraCranialVol
Seg_Validation <- subset(Seg_Validation, select=c("Measure.volume", "BPF", "GMF", "WMF", "WM.hypointensities", "BrainSegVolNotVent"))
names(Seg_Validation) <- c("ID", "BPF", "GMF", "WMF", "T1_lesions", "Brain_volume_ul")
Seg_Validation

#load validation cohort excel file
EPVS_validation <- read.table("EPVS_quantification_validationcohort_sorted.txt", sep="\t", header=TRUE, stringsAsFactors = F)
EPVS_validation <- merge(EPVS_validation, Seg_Validation, by="ID")

#merge table with T2 lesion output from LST
LST_T2 <- read_excel("LST_output.xlsx")
EPVS_validation <- merge(EPVS_validation, LST_T2, by="ID", all.x = T)

####T2_PDF report results extraction Lesion segmentation toolbox SPM
LST_output <- pdf_text("HTML_report_LST_LPA.pdf")
  #LST_output <- str_split(LST_output, pattern = "\n|\r")
  #LST_output <- paste(unlist(LST_output), sep = " ")
  #LST_output <- paste(LST_output, collapse = " ")

  #str_count(" ples_lpa_mSW1-1875_M12_FLAIR.nii")#=33 characters
  #str_count(" mREMY040_M0_FLAIR.nii")#=22 characters
  #str_count("Lesion map ples_lpa_mSW1-1875_M12_FLAIR.nii Algorithm used for segmentation LPA Lesion volume")#93

  #str_view_all(LST_output, regex("Lesion map.{1,}\\b", ignore_case = T))
  #str_view_all(LST_output, regex("Lesion map.{1,200} ml", ignore_case = T, dotall = T))
LST_output <- str_extract_all(LST_output, regex("Lesion map.{1,200} ml", ignore_case = T, dotall = T))
LST_output <- Filter(length, LST_output)
LST_output <- paste(unlist(LST_output), collapse = ";")

  #str_view_all(LST_output, regex("SW1-\\d{4}_M\\d{1,2}|REMY\\d{3}_M\\d{1,2}", ignore_case = T))
  #str_view_all(LST_output, regex("\\d{1,2}\\.\\d{1,}|\\b0\\b", ignore_case = T))
Subject.temp <- str_extract_all(LST_output, regex("SW1-\\d{4}_M\\d{1,2}|REMY\\d{3}_M\\d{1,2}", ignore_case = T))
Volume.temp <- str_extract_all(LST_output, regex("\\d{1,2}\\.\\d{1,}|\\b0\\b", ignore_case = T))
Subject.temp.df <- data.frame(matrix(unlist(Subject.temp), nrow=65, byrow=T),stringsAsFactors=FALSE)
Volume.temp.df <- data.frame(matrix(unlist(Volume.temp), nrow=65, byrow=T),stringsAsFactors=FALSE)
LST_output <- cbind(Subject.temp.df, Volume.temp.df)
names(LST_output) <- c("ID", "T2_volume_ml")
write_xlsx(LST_output, "LST_output.xlsx")

###proximity analysis data extraction (validation cohort)
Prox_Remy <- read.delim(file = "EPVSvolumesREMY.txt")
Prox_Remy1 <- as.data.frame(Prox_Remy[seq(1, nrow(Prox_Remy), 3), ])
Prox_Remy2 <- as.data.frame(Prox_Remy[seq(2, nrow(Prox_Remy), 3), ])
Prox_Remy3 <- as.data.frame(Prox_Remy[seq(3, nrow(Prox_Remy), 3), ])
Prox_Remy <- cbind(Prox_Remy1, Prox_Remy2, Prox_Remy3)
Prox_Remy$`Prox_Remy[seq(1, nrow(Prox_Remy), 3), ]` <- gsub("\\ \\d{1,}\\.000000", "", Prox_Remy$`Prox_Remy[seq(1, nrow(Prox_Remy), 3), ]`)
BGminusT1maskbin <- Prox_Remy[1:38,1:3]
names(BGminusT1maskbin) <- c("BGminusT1maskbin", "ID", "n")
BGminusT2maskbin <- Prox_Remy[39:76,1:3]
names(BGminusT2maskbin) <- c("BGminusT2maskbin", "ID", "n")
CSOminusT1maskbin <- Prox_Remy[77:114,1:3]
names(CSOminusT1maskbin) <- c("CSOminusT1maskbin", "ID", "n")
CSOminusT2maskbin <- Prox_Remy[115:152,1:3]
names(CSOminusT2maskbin) <- c("CSOminusT2maskbin", "ID", "n")
CSO1minusT1maskbin <- Prox_Remy[153:190,1:3]
names(CSO1minusT1maskbin) <- c("CSO1minusT1maskbin", "ID", "n")
CSO1minusT2maskbin <- Prox_Remy[191:228,1:3]
names(CSO1minusT2maskbin) <- c("CSO1minusT2maskbin", "ID", "n")
CSO2minusT1maskbin <- Prox_Remy[229:266,1:3]
names(CSO2minusT1maskbin) <- c("CSO2minusT1maskbin", "ID", "n")
CSO2minusT2maskbin <- Prox_Remy[267:304,1:3]
names(CSO2minusT2maskbin) <- c("CSO2minusT2maskbin", "ID", "n")
BGmaskbin <- Prox_Remy[305:342,1:3]
names(BGmaskbin) <- c("BGmaskbin", "ID", "n")
BGmaskbinexp <- Prox_Remy[343:380,1:3]
names(BGmaskbinexp) <- c("BGmaskbinexp", "ID", "n")
CSO1maskbin <- Prox_Remy[381:418,1:3]
names(CSO1maskbin) <- c("CSO1maskbin", "ID", "n")
CSO1maskbinexp <- Prox_Remy[419:456,1:3]
names(CSO1maskbinexp) <- c("CSO1maskbinexp", "ID", "n")
CSO2maskbin <- Prox_Remy[457:494,1:3]
names(CSO2maskbin) <- c("CSO2maskbin", "ID", "n")
CSO2maskbinexp <- Prox_Remy[495:532,1:3]
names(CSO2maskbinexp) <- c("CSO2maskbinexp", "ID", "n")
CSOmaskbin <- Prox_Remy[533:570,1:3]
names(CSOmaskbin) <- c("CSOmaskbin", "ID", "n")
CSOmaskbinexp <- Prox_Remy[571:608,1:3]
names(CSOmaskbinexp) <- c("CSOmaskbinexp", "ID", "n")
T1maskbin <- Prox_Remy[609:646,1:3]
names(T1maskbin) <- c("T1maskbin", "ID", "n")
T2maskbin <- Prox_Remy[647:684,1:3]
names(T2maskbin) <- c("T2maskbin", "ID", "n")
Prox_Remy <- Reduce(function(...) merge(..., by = "ID", all=TRUE), list(BGminusT1maskbin[,1:2], BGminusT2maskbin[,1:2],
                                                                   CSOminusT1maskbin[,1:2], CSOminusT2maskbin[,1:2],
                                                                   CSO1minusT1maskbin[,1:2], CSO1minusT2maskbin[,1:2],
                                                                   CSO2minusT1maskbin[,1:2], CSO2minusT2maskbin[,1:2],
                                                                   BGmaskbin[,1:2], BGmaskbinexp[,1:2],
                                                                   CSO1maskbin[,1:2], CSO1maskbinexp[,1:2],
                                                                   CSO2maskbin[,1:2], CSO2maskbinexp[,1:2],
                                                                   CSOmaskbin[,1:2], CSOmaskbinexp[,1:2],
                                                                   T1maskbin[,1:2], T2maskbin[,1:2]))

Prox_Mul  <- read.delim(file = "EPVSvolumesMULT.txt")
Prox_Mul1 <- as.data.frame(Prox_Mul[seq(1, nrow(Prox_Mul), 3), ])
Prox_Mul2 <- as.data.frame(Prox_Mul[seq(2, nrow(Prox_Mul), 3), ])
Prox_Mul3 <- as.data.frame(Prox_Mul[seq(3, nrow(Prox_Mul), 3), ])
Prox_Mul <- cbind(Prox_Mul1, Prox_Mul2, Prox_Mul3)
Prox_Mul$`Prox_Mul[seq(1, nrow(Prox_Mul), 3), ]` <- gsub("\\ \\d{1,}\\.000000", "", Prox_Mul$`Prox_Mul[seq(1, nrow(Prox_Mul), 3), ]`)
BGminusT1maskbin <- Prox_Mul[1:25,1:3]
names(BGminusT1maskbin) <- c("BGminusT1maskbin", "ID", "n")
BGminusT2maskbin <- Prox_Mul[26:50,1:3]
names(BGminusT2maskbin) <- c("BGminusT2maskbin", "ID", "n")
CSOminusT1maskbin <- Prox_Mul[51:75,1:3]
names(CSOminusT1maskbin) <- c("CSOminusT1maskbin", "ID", "n")
CSOminusT2maskbin <- Prox_Mul[76:100,1:3]
names(CSOminusT2maskbin) <- c("CSOminusT2maskbin", "ID", "n")
CSO1minusT1maskbin <- Prox_Mul[101:125,1:3]
names(CSO1minusT1maskbin) <- c("CSO1minusT1maskbin", "ID", "n")
CSO1minusT2maskbin <- Prox_Mul[126:150,1:3]
names(CSO1minusT2maskbin) <- c("CSO1minusT2maskbin", "ID", "n")
CSO2minusT1maskbin <- Prox_Mul[151:175,1:3]
names(CSO2minusT1maskbin) <- c("CSO2minusT1maskbin", "ID", "n")
CSO2minusT2maskbin <- Prox_Mul[176:200,1:3]
names(CSO2minusT2maskbin) <- c("CSO2minusT2maskbin", "ID", "n")
BGmaskbin <- Prox_Mul[201:225,1:3]
names(BGmaskbin) <- c("BGmaskbin", "ID", "n")
BGmaskbinexp <- Prox_Mul[226:250,1:3]
names(BGmaskbinexp) <- c("BGmaskbinexp", "ID", "n")
CSO1maskbin <- Prox_Mul[251:275,1:3]
names(CSO1maskbin) <- c("CSO1maskbin", "ID", "n")
CSO1maskbinexp <- Prox_Mul[276:300,1:3]
names(CSO1maskbinexp) <- c("CSO1maskbinexp", "ID", "n")
CSO2maskbin <- Prox_Mul[301:325,1:3]
names(CSO2maskbin) <- c("CSO2maskbin", "ID", "n")
CSO2maskbinexp <- Prox_Mul[326:350,1:3]
names(CSO2maskbinexp) <- c("CSO2maskbinexp", "ID", "n")
CSOmaskbin <- Prox_Mul[351:375,1:3]
names(CSOmaskbin) <- c("CSOmaskbin", "ID", "n")
CSOmaskbinexp <- Prox_Mul[376:400,1:3]
names(CSOmaskbinexp) <- c("CSOmaskbinexp", "ID", "n")
T1maskbin <- Prox_Mul[401:425,1:3]
names(T1maskbin) <- c("T1maskbin", "ID", "n")
T2maskbin <- Prox_Mul[426:450,1:3]
names(T2maskbin) <- c("T2maskbin", "ID", "n")
Prox_Mul <- Reduce(function(...) merge(..., by = "ID", all=TRUE), list(BGminusT1maskbin[,1:2], BGminusT2maskbin[,1:2],
                                                                        CSOminusT1maskbin[,1:2], CSOminusT2maskbin[,1:2],
                                                                        CSO1minusT1maskbin[,1:2], CSO1minusT2maskbin[,1:2],
                                                                        CSO2minusT1maskbin[,1:2], CSO2minusT2maskbin[,1:2],
                                                                        BGmaskbin[,1:2], BGmaskbinexp[,1:2],
                                                                        CSO1maskbin[,1:2], CSO1maskbinexp[,1:2],
                                                                        CSO2maskbin[,1:2], CSO2maskbinexp[,1:2],
                                                                        CSOmaskbin[,1:2], CSOmaskbinexp[,1:2],
                                                                        T1maskbin[,1:2], T2maskbin[,1:2]))

Prox_Remy$cohort <- "REMYDI"
Prox_Mul$cohort <- "MULTIPLEMS"
Prox = rbind(Prox_Remy, Prox_Mul)
Prox[,2:19] <- sapply(Prox[2:19],as.numeric)#convert numbers to numeric (columns 2:19)
Prox$ID <- as.character(Prox$ID)
Prox$cohort <- as.factor(Prox$cohort)
Prox$ID <- str_replace_all(Prox$ID, "_M0", "")

#secondary analysis with expansion only by one voxel
Prox_Mul_Vox1  <- read.delim(file = "EPVSvolumes_MMS_Vox1.txt")
Prox_Mul1_vox1 <- as.data.frame(Prox_Mul_Vox1[seq(1, nrow(Prox_Mul_Vox1), 3), ])
Prox_Mul2_vox1 <- as.data.frame(Prox_Mul_Vox1[seq(2, nrow(Prox_Mul_Vox1), 3), ])
Prox_Mul3_vox1 <- as.data.frame(Prox_Mul_Vox1[seq(3, nrow(Prox_Mul_Vox1), 3), ])
Prox_Mul_Vox1 <- cbind(Prox_Mul1_vox1, Prox_Mul2_vox1, Prox_Mul3_vox1)
Prox_Mul_Vox1$`Prox_Mul_Vox1[seq(1, nrow(Prox_Mul_Vox1), 3), ]` <- gsub("\\ \\d{1,}\\.000000", "", Prox_Mul_Vox1$`Prox_Mul_Vox1[seq(1, nrow(Prox_Mul_Vox1), 3), ]`)
BGmaskbinexpVox1 <- Prox_Mul_Vox1[1:25,1:3]
names(BGmaskbinexpVox1) <- c("BGmaskbinexpVox1", "ID", "n")
CSO1maskbinexpVox1 <- Prox_Mul_Vox1[26:50,1:3]
names(CSO1maskbinexpVox1) <- c("CSO1maskbinexpVox1", "ID", "n")
CSO2maskbinexpVox1 <- Prox_Mul_Vox1[51:75,1:3]
names(CSO2maskbinexpVox1) <- c("CSO2maskbinexpVox1", "ID", "n")
CSOmaskbinexpVox1 <- Prox_Mul_Vox1[76:100,1:3]
names(CSOmaskbinexpVox1) <- c("CSOmaskbinexpVox1", "ID", "n")
CSOminusT1maskbinVox1 <- Prox_Mul_Vox1[101:125,1:3]
names(CSOminusT1maskbinVox1) <- c("CSOminusT1maskbinVox1", "ID", "n")
CSO1minusT1maskbinVox1 <- Prox_Mul_Vox1[126:150,1:3]
names(CSO1minusT1maskbinVox1) <- c("CSO1minusT1maskbinVox1", "ID", "n")
CSO2minusT1maskbinVox1 <- Prox_Mul_Vox1[151:175,1:3]
names(CSO2minusT1maskbinVox1) <- c("CSO2minusT1maskbinVox1", "ID", "n")
BGminusT1maskbinVox1 <- Prox_Mul_Vox1[176:200,1:3]
names(BGminusT1maskbinVox1) <- c("BGminusT1maskbinVox1", "ID", "n")
CSOminusT2maskbinVox1 <- Prox_Mul_Vox1[201:225,1:3]
names(CSOminusT2maskbinVox1) <- c("CSOminusT2maskbinVox1", "ID", "n")
CSO1minusT2maskbinVox1 <- Prox_Mul_Vox1[226:250,1:3]
names(CSO1minusT2maskbinVox1) <- c("CSO1minusT2maskbinVox1", "ID", "n")
CSO2minusT2maskbinVox1 <- Prox_Mul_Vox1[251:275,1:3]
names(CSO2minusT2maskbinVox1) <- c("CSO2minusT2maskbinVox1", "ID", "n")
BGminusT2maskbinVox1 <- Prox_Mul_Vox1[276:300,1:3]
names(BGminusT2maskbinVox1) <- c("BGminusT2maskbinVox1", "ID", "n")
Prox_Mul_Vox1 <- Reduce(function(...) merge(..., by = "ID", all=TRUE), list(BGmaskbinexpVox1[,1:2], CSO1maskbinexpVox1[,1:2],
                                                                            CSO2maskbinexpVox1[,1:2], CSOmaskbinexpVox1[,1:2],
                                                                            CSOminusT1maskbinVox1[,1:2], CSO1minusT1maskbinVox1[,1:2],
                                                                            CSO2minusT1maskbinVox1[,1:2], BGminusT1maskbinVox1[,1:2],
                                                                            CSOminusT2maskbinVox1[,1:2], CSO1minusT2maskbinVox1[,1:2],
                                                                            CSO2minusT2maskbinVox1[,1:2], BGminusT2maskbinVox1[,1:2]))

Prox_Rem_Vox1  <- read.delim(file = "EPVSvolumes_REMY_Vox1.txt")
Prox_Rem1_vox1 <- as.data.frame(Prox_Rem_Vox1[seq(1, nrow(Prox_Rem_Vox1), 3), ])
Prox_Rem2_vox1 <- as.data.frame(Prox_Rem_Vox1[seq(2, nrow(Prox_Rem_Vox1), 3), ])
Prox_Rem3_vox1 <- as.data.frame(Prox_Rem_Vox1[seq(3, nrow(Prox_Rem_Vox1), 3), ])
Prox_Rem_Vox1 <- cbind(Prox_Rem1_vox1, Prox_Rem2_vox1, Prox_Rem3_vox1)
Prox_Rem_Vox1$`Prox_Rem_Vox1[seq(1, nrow(Prox_Rem_Vox1), 3), ]` <- gsub("\\ \\d{1,}\\.000000", "", Prox_Rem_Vox1$`Prox_Rem_Vox1[seq(1, nrow(Prox_Rem_Vox1), 3), ]`)
BGmaskbinexpVox1 <- Prox_Rem_Vox1[1:38,1:3]
names(BGmaskbinexpVox1) <- c("BGmaskbinexpVox1", "ID", "n")
CSO1maskbinexpVox1 <- Prox_Rem_Vox1[39:76,1:3]
names(CSO1maskbinexpVox1) <- c("CSO1maskbinexpVox1", "ID", "n")
CSO2maskbinexpVox1 <- Prox_Rem_Vox1[77:114,1:3]
names(CSO2maskbinexpVox1) <- c("CSO2maskbinexpVox1", "ID", "n")
CSOmaskbinexpVox1 <- Prox_Rem_Vox1[115:152,1:3]
names(CSOmaskbinexpVox1) <- c("CSOmaskbinexpVox1", "ID", "n")
CSOminusT1maskbinVox1 <- Prox_Rem_Vox1[153:190,1:3]
names(CSOminusT1maskbinVox1) <- c("CSOminusT1maskbinVox1", "ID", "n")
CSO1minusT1maskbinVox1 <- Prox_Rem_Vox1[191:228,1:3]
names(CSO1minusT1maskbinVox1) <- c("CSO1minusT1maskbinVox1", "ID", "n")
CSO2minusT1maskbinVox1 <- Prox_Rem_Vox1[229:266,1:3]
names(CSO2minusT1maskbinVox1) <- c("CSO2minusT1maskbinVox1", "ID", "n")
BGminusT1maskbinVox1 <- Prox_Rem_Vox1[267:304,1:3]
names(BGminusT1maskbinVox1) <- c("BGminusT1maskbinVox1", "ID", "n")
CSOminusT2maskbinVox1 <- Prox_Rem_Vox1[305:342,1:3]
names(CSOminusT2maskbinVox1) <- c("CSOminusT2maskbinVox1", "ID", "n")
CSO1minusT2maskbinVox1 <- Prox_Rem_Vox1[343:380,1:3]
names(CSO1minusT2maskbinVox1) <- c("CSO1minusT2maskbinVox1", "ID", "n")
CSO2minusT2maskbinVox1 <- Prox_Rem_Vox1[381:418,1:3]
names(CSO2minusT2maskbinVox1) <- c("CSO2minusT2maskbinVox1", "ID", "n")
BGminusT2maskbinVox1 <- Prox_Rem_Vox1[419:456,1:3]
names(BGminusT2maskbinVox1) <- c("BGminusT2maskbinVox1", "ID", "n")
Prox_Rem_Vox1 <- Reduce(function(...) merge(..., by = "ID", all=TRUE), list(BGmaskbinexpVox1[,1:2], CSO1maskbinexpVox1[,1:2],
                                                                            CSO2maskbinexpVox1[,1:2], CSOmaskbinexpVox1[,1:2],
                                                                            CSOminusT1maskbinVox1[,1:2], CSO1minusT1maskbinVox1[,1:2],
                                                                            CSO2minusT1maskbinVox1[,1:2], BGminusT1maskbinVox1[,1:2],
                                                                            CSOminusT2maskbinVox1[,1:2], CSO1minusT2maskbinVox1[,1:2],
                                                                            CSO2minusT2maskbinVox1[,1:2], BGminusT2maskbinVox1[,1:2]))

Prox_Vox1 = rbind(Prox_Rem_Vox1, Prox_Mul_Vox1)
Prox_Vox1[,2:13] <- sapply(Prox_Vox1[2:13],as.numeric)#convert numbers to numeric (columns 2:19)
Prox_Vox1$ID <- as.character(Prox_Vox1$ID)
Prox_Vox1$ID <- str_replace_all(Prox_Vox1$ID, "_M0", "")

Prox <- merge(Prox, Prox_Vox1, by = "ID", all = T)

Sex_Val <- EPVS_validation
Sex_Val$ID <- str_replace_all(Sex_Val$ID, "_M\\d{1,2}", "")
Sex_Val <- subset(Sex_Val, select = c("ID", "Sex", "age", "Cohort", "T2_volume_ml", "Brain_volume_ul"))
Prox <- merge(Prox, Sex_Val, by = "ID", all.x = T)
###proximity analysis data extraction (control cohort)




###proximity analysis data extraction (primary cohort STOPMS)
STOPMS <- read.delim(file = "EPVSvolumes_STOPMS.txt")
STOPMS1 <- as.data.frame(STOPMS[seq(1, nrow(STOPMS), 3), ])
STOPMS2 <- as.data.frame(STOPMS[seq(2, nrow(STOPMS), 3), ])
STOPMS3 <- as.data.frame(STOPMS[seq(3, nrow(STOPMS), 3), ])
STOPMS <- cbind(STOPMS1, STOPMS2, STOPMS3)
STOPMS$`STOPMS[seq(1, nrow(STOPMS), 3), ]` <- gsub("\\ \\d{1,}\\.000000", "", STOPMS$`STOPMS[seq(1, nrow(STOPMS), 3), ]`)
BGminusT1maskbin <- STOPMS[1:134,1:3]
names(BGminusT1maskbin) <- c("BGminusT1maskbin", "ID", "n")
CSOminusT1maskbin <- STOPMS[135:268,1:3]
names(CSOminusT1maskbin) <- c("CSOminusT1maskbin", "ID", "n")
BGmaskbin <- STOPMS[269:402,1:3]
names(BGmaskbin) <- c("BGmaskbin", "ID", "n")
BGmaskbinexp <- STOPMS[403:536,1:3]
names(BGmaskbinexp) <- c("BGmaskbinexp", "ID", "n")
CSOmaskbin <- STOPMS[537:670,1:3]
names(CSOmaskbin) <- c("CSOmaskbin", "ID", "n")
CSOmaskbinexp <- STOPMS[671:804,1:3]
names(CSOmaskbinexp) <- c("CSOmaskbinexp", "ID", "n")
T1maskbin <- STOPMS[805:938,1:3]
names(T1maskbin) <- c("T1maskbin", "ID", "n")
STOPMS <- Reduce(function(...) merge(..., by = "ID", all=TRUE), list(BGminusT1maskbin[,1:2],CSOminusT1maskbin[,1:2],
                                                                          BGmaskbin[,1:2], BGmaskbinexp[,1:2],
                                                                        CSOmaskbin[,1:2], CSOmaskbinexp[,1:2],
                                                                        T1maskbin[,1:2]))

STOPMS[,2:8] <- sapply(STOPMS[2:8],as.numeric)#convert numbers to numeric (columns 2:7)

#### Proximity analysis #####
#thoughts:  #BGmaskbinexp = total volume of EPVS
            #BGminusT1 = volume of EPVS NOT adjacent to lesion
            #BGmaskbinexp - BGminusT1 = EPVS volume adjacent to lesion
Prox$prox_BG_T1 <- (Prox$BGmaskbinexp - Prox$BGminusT1maskbin)/(Prox$BGmaskbinexp+0.0001)
Prox$prox_CSO1_T1 <- (Prox$CSO1maskbinexp - Prox$CSO1minusT1maskbin)/(Prox$CSO1maskbinexp+0.0001)
Prox$prox_CSO2_T1 <- (Prox$CSO2maskbinexp - Prox$CSO2minusT1maskbin)/(Prox$CSO2maskbinexp+0.0001)
Prox$prox_CSO_T1 <- (Prox$CSOmaskbinexp - Prox$CSOminusT1maskbin)/(Prox$CSOmaskbinexp+0.0001)

Prox$prox_BG_T2 <- (Prox$BGmaskbinexp - Prox$BGminusT2maskbin)/(Prox$BGmaskbinexp+0.0001)
Prox$prox_CSO1_T2 <- (Prox$CSO1maskbinexp - Prox$CSO1minusT2maskbin)/(Prox$CSO1maskbinexp+0.0001)
Prox$prox_CSO2_T2 <- (Prox$CSO2maskbinexp - Prox$CSO2minusT2maskbin)/(Prox$CSO2maskbinexp+0.0001)
Prox$prox_CSO_T2 <- (Prox$CSOmaskbinexp - Prox$CSOminusT2maskbin)/(Prox$CSOmaskbinexp+0.0001)

Prox$CSO1toCSO1exp <- Prox$CSO1maskbinexp/Prox$CSO1maskbin
Prox$CSO2toCSO2exp <- Prox$CSO2maskbinexp/Prox$CSO2maskbin
Prox$CSOtoCSOexp <- Prox$CSOmaskbinexp/Prox$CSOmaskbin
Prox$BGtoBGexp <- Prox$BGmaskbinexp/Prox$BGmaskbin
plot(Prox$CSOtoCSOexp)

Prox$CSO1_T1_corr <- Prox$prox_CSO1_T1*Prox$CSO1toCSO1exp
Prox$CSO2_T1_corr <- Prox$prox_CSO2_T1*Prox$CSO2toCSO2exp
Prox$CSO_T1_corr <- Prox$prox_CSO_T1*Prox$CSOtoCSOexp
Prox$BG_T1_corr <- Prox$prox_BG_T1*Prox$BGtoBGexp

Prox$CSO1_T2_corr <- Prox$prox_CSO1_T2*Prox$CSO1toCSO1exp
Prox$CSO2_T2_corr <- Prox$prox_CSO2_T2*Prox$CSO2toCSO2exp
Prox$CSO_T2_corr <- Prox$prox_CSO_T2*Prox$CSOtoCSOexp
Prox$BG_T2_corr <- Prox$prox_BG_T2*Prox$BGtoBGexp

hist(Prox$prox_CSO1_T1, breaks = 50)
hist(Prox$prox_CSO1_T2, breaks = 50)
hist(Prox$prox_CSO2_T1, breaks = 50)
hist(Prox$prox_CSO2_T2, breaks = 50)

summary(Prox$prox_CSO_T1*100)
sd(Prox$prox_CSO_T1*100)
summary(Prox$prox_CSO_T2*100)
sd(Prox$prox_CSO_T2*100)

summary(Prox$prox_CSO1_T1*100)
sd(Prox$prox_CSO1_T1*100)
summary(Prox$prox_CSO1_T2*100)
sd(Prox$prox_CSO2_T1*100)

summary(Prox$prox_CSO2_T1*100)
sd(Prox$prox_CSO1_T2*100)
summary(Prox$prox_CSO2_T2*100)
sd(Prox$prox_CSO2_T2*100)

summary(Prox$prox_BG_T1*100)
sd(Prox$prox_BG_T1*100, na.rm = T)
summary(Prox$prox_BG_T2*100)
sd(Prox$prox_BG_T2*100, na.rm = T)

summary(Prox$prox_CSO_T1[Prox$prox_CSO_T1 > 0])
plot(Prox$prox_CSO_T1[Prox$prox_CSO_T1 > 0])

summary(Prox$prox_CSO_T2[Prox$prox_CSO_T1 > 0])
plot(Prox$prox_CSO_T1[Prox$prox_CSO_T1 > 0])

summary(Prox$prox_CSO_T1[Prox$prox_CSO_T1 > 0])
plot(Prox_CSOT1$prox_CSO_T1[Prox$prox_CSO_T1 > 0])

summary(Prox$prox_CSO_T1[Prox$prox_CSO_T1 > 0])
plot(Prox_CSOT1$prox_CSO_T1[Prox$prox_CSO_T1 > 0])


    #how much volume do EPVS constitute of whole brain volume?
Segmentation <- rbind(Seg_MultipleMS, Seg_MyelinMS)

summary(Segmentation$CerebralWhiteMatterVol)#median cerebral white matter volume = 469369 mm3 = 469369 ul
summary(Segmentation$BrainSegVolNotVent)#median total cerebral volume = 1163953 mm3 = 1163953 ul

summary(Prox$CSOmaskbin)#median CSO EPVS volume = 326 ul
(326/1163953)*100#CSO EPVS make around 0.028% of  cerebral volume
summary(Prox$BGmaskbin)#median BG EPVS volume = 147 ul
(147/1163953)*100#BG EPVS make around 0.013% of  cerebral volume

#summary(Segmentation$WM.hypointensities))#median total cerebral volume = 1308 mm3 = 1308 ul
(1308/1163953)*100#WM hypointensities make around 0.112% of  cerebral volume (median T1 lesion volume 10 times higher compared
                    #to BG EPVS volume and 3 times higher compared to CSO EPVS volume)

summary(Prox$T2_volume_ml)#median T2 lesion volume = 297 mm3 = 297 ul
(297/1163953)*100#T2 lesions make around 0.026% of cerebral volume (median T2 lesion volume approx= CSO EPVS volume snd twice BG EPVS volume

summary(Long$BrainSegVolNotVent)#median total cerebral volume = 1107658 mm3 = 1107658 ul
summary(Long$Vol_PVS_CSO)#median CSO volume 168 ul
(168/1107658)*100#CSO EPVS make around 0.015% of cerebral volume

summary(Long$Vol_PVS_BG)#median BG volume 136 ul
(136/1107658)*100#BG EPVS make around 0.012% of  cerebral volume

summary(Long$T2_Volume_FS)#median T2 lesion volume 2100 ul
(2100/1107658)*100#T2 makes around 0.189% of cerebral volume

summary(Long$WMH1)#median T1 lesion volume 2084 ul
(2084/1107658)*100#T1 makes around 0.188% of cerebral volume



summary(EPVS_controls$Brain_volume_ul)#median total brain volume 1153608 ul
summary(EPVS_controls$CSO_Total_Vol)#median CSO volume 145 ul
(145/1153608)*100#CSO EPVS make around 0.013% of cerebral volume

summary(EPVS_controls$BG_Vol)#median BG volume 159 ul
(159/1153608)*100#BG EPVS make around 0.014% of  cerebral volume



summary(lm(prox_CSO_T1 ~ age+Sex+Cohort, data=Prox))
summary(lm(prox_CSO_T2 ~ age+Sex+Cohort, data=Prox))

summary(lm(prox_CSO1_T1-prox_CSO2_T1 ~ age+Sex+Cohort, data=Prox))
summary(lm(prox_CSO1_T2-prox_CSO2_T2 ~ age+Sex+Cohort, data=Prox))

summary(lm(prox_CSO1_T1 ~ age+Sex+Cohort, data=Prox))
summary(lm(prox_CSO2_T1 ~ age+Sex+Cohort, data=Prox))

summary(lm(prox_CSO1_T2 ~ age+Sex+Cohort, data=Prox))
summary(lm(prox_CSO2_T2 ~ age+Sex+Cohort, data=Prox))

#parametric tests
t.test(Prox$prox_CSO1_T1, Prox$prox_CSO2_T1, alternative = c("two.sided"), paired = F)#non-sig
t.test(Prox$prox_CSO1_T2, Prox$prox_CSO2_T2, alternative = c("two.sided"), paired = F)#non-sig

t.test(Prox$prox_CSO1_T1[Prox$cohort == "REMYDI"], Prox$prox_CSO2_T1[Prox$cohort == "REMYDI"], alternative = c("two.sided"), paired = F)#non-sig
t.test(Prox$prox_CSO1_T2[Prox$cohort == "REMYDI"], Prox$prox_CSO2_T2[Prox$cohort == "REMYDI"], alternative = c("two.sided"), paired = F)#non-sig

t.test(Prox$prox_CSO1_T1[Prox$cohort == "MULTIPLEMS"], Prox$prox_CSO2_T1[Prox$cohort == "MULTIPLEMS"], alternative = c("two.sided"), paired = F)#non-sig
t.test(Prox$prox_CSO1_T2[Prox$cohort == "MULTIPLEMS"], Prox$prox_CSO2_T2[Prox$cohort == "MULTIPLEMS"], alternative = c("two.sided"), paired = F)#non-sig

#non-parametric tests (data w/o Gaussian distribution)
wilcox.test(Prox$prox_CSO1_T1,Prox$prox_CSO2_T1)
wilcox.test(Prox$prox_CSO1_T2,Prox$prox_CSO2_T2)

wilcox.test(Prox$prox_CSO1_T1[Prox$cohort == "REMYDI"],Prox$prox_CSO2_T1[Prox$cohort == "REMYDI"])
wilcox.test(Prox$prox_CSO1_T2[Prox$cohort == "REMYDI"],Prox$prox_CSO2_T2[Prox$cohort == "REMYDI"])

wilcox.test(Prox$prox_CSO1_T1[Prox$cohort == "MULTIPLEMS"],Prox$prox_CSO2_T1[Prox$cohort == "MULTIPLEMS"])
wilcox.test(Prox$prox_CSO1_T2[Prox$cohort == "MULTIPLEMS"],Prox$prox_CSO2_T2[Prox$cohort == "MULTIPLEMS"])

prox_CSO1_T1 <- lm(prox_CSO1_T1 ~ age+Sex+Cohort, data=Prox)
prox_CSO1_T1 <- fitted.values(prox_CSO1_T1)
prox_CSO2_T1 <- lm(prox_CSO2_T1 ~ age+Sex+Cohort, data=Prox)
prox_CSO2_T1 <- fitted.values(prox_CSO2_T1)
prox_CSO1_T2 <- lm(prox_CSO1_T2 ~ age+Sex+Cohort, data=Prox)
prox_CSO1_T2 <- fitted.values(prox_CSO1_T2)
prox_CSO2_T2 <- lm(prox_CSO2_T2 ~ age+Sex+Cohort, data=Prox)
prox_CSO2_T2 <- fitted.values(prox_CSO2_T2)

wilcox.test(prox_CSO1_T1,prox_CSO2_T1)
wilcox.test(prox_CSO1_T2,prox_CSO2_T2)

summary(prox_CSO1_T1*100)
summary(prox_CSO2_T1*100)
summary(prox_CSO1_T2*100)
summary(prox_CSO2_T2*100)


#subgroup analysis (only patients which have some CSO2)
Prox_sub <- Prox %>% filter(CSO2maskbin > 0)

summary(Prox_sub$prox_CSO1_T1)
sd(Prox_sub$prox_CSO1_T1)
summary(Prox_sub$prox_CSO2_T1)
sd(Prox_sub$prox_CSO2_T1)

summary(Prox_sub$prox_CSO1_T2)
sd(Prox_sub$prox_CSO1_T2)
summary(Prox_sub$prox_CSO2_T2)
sd(Prox_sub$prox_CSO2_T2)

wilcox.test(Prox_sub$prox_CSO1_T1,Prox_sub$prox_CSO2_T1)
wilcox.test(Prox_sub$prox_CSO1_T2,Prox_sub$prox_CSO2_T2)

prox_CSO1_T1 <- lm(prox_CSO1_T1 ~ age+Sex+Cohort, data=Prox_sub)
prox_CSO1_T1 <- fitted.values(prox_CSO1_T1)
prox_CSO2_T1 <- lm(prox_CSO2_T1 ~ age+Sex+Cohort, data=Prox_sub)
prox_CSO2_T1 <- fitted.values(prox_CSO2_T1)
prox_CSO1_T2 <- lm(prox_CSO1_T2 ~ age+Sex+Cohort, data=Prox_sub)
prox_CSO1_T2 <- fitted.values(prox_CSO1_T2)
prox_CSO2_T2 <- lm(prox_CSO2_T2 ~ age+Sex+Cohort, data=Prox_sub)
prox_CSO2_T2 <- fitted.values(prox_CSO2_T2)

wilcox.test(prox_CSO1_T1,prox_CSO2_T1)
wilcox.test(prox_CSO1_T2,prox_CSO2_T2)

summary(prox_CSO1_T1*100)
summary(prox_CSO2_T1*100)
summary(prox_CSO1_T2*100)
summary(prox_CSO2_T2*100)

#secondary analysis with masks only expanded by 1 voxel
Prox$prox_BG_T1_Vox1 <- (Prox$BGmaskbinexpVox1 - Prox$BGminusT1maskbinVox1)/(Prox$BGmaskbinexpVox1+0.0001)
Prox$prox_CSO1_T1_Vox1 <- (Prox$CSO1maskbinexpVox1 - Prox$CSO1minusT1maskbinVox1)/(Prox$CSO1maskbinexpVox1+0.0001)
Prox$prox_CSO2_T1_Vox1 <- (Prox$CSO2maskbinexpVox1 - Prox$CSO2minusT1maskbinVox1)/(Prox$CSO2maskbinexpVox1+0.0001)
Prox$prox_CSO_T1_Vox1 <- (Prox$CSOmaskbinexpVox1 - Prox$CSOminusT1maskbinVox1)/(Prox$CSOmaskbinexpVox1+0.0001)

Prox$prox_BG_T2_Vox1 <- (Prox$BGmaskbinexpVox1 - Prox$BGminusT2maskbinVox1)/(Prox$BGmaskbinexpVox1+0.0001)
Prox$prox_CSO1_T2_Vox1 <- (Prox$CSO1maskbinexpVox1 - Prox$CSO1minusT2maskbinVox1)/(Prox$CSO1maskbinexpVox1+0.0001)
Prox$prox_CSO2_T2_Vox1 <- (Prox$CSO2maskbinexpVox1 - Prox$CSO2minusT2maskbinVox1)/(Prox$CSO2maskbinexpVox1+0.0001)
Prox$prox_CSO_T2_Vox1 <- (Prox$CSOmaskbinexpVox1 - Prox$CSOminusT2maskbinVox1)/(Prox$CSOmaskbinexpVox1+0.0001)

summary(Prox$prox_CSO1_T1_Vox1)
sd(Prox$prox_CSO1_T1_Vox1)
summary(Prox$prox_CSO2_T1_Vox1)
sd(Prox$prox_CSO2_T1_Vox1)

summary(Prox$prox_CSO1_T2_Vox1)
sd(Prox$prox_CSO1_T2_Vox1)
summary(Prox$prox_CSO2_T2_Vox1)
sd(Prox$prox_CSO2_T2_Vox1)

summary(Prox$prox_CSO_T1_Vox1)
sd(Prox$prox_CSO_T1_Vox1)
summary(Prox$prox_CSO_T2_Vox1)
sd(Prox$prox_CSO_T2_Vox1)

summary(Prox$prox_BG_T1_Vox1)
sd(Prox$prox_BG_T1_Vox1)
summary(Prox$prox_BG_T2_Vox1)
sd(Prox$prox_BG_T2_Vox1)

wilcox.test(Prox$prox_CSO1_T1_Vox1,Prox$prox_CSO2_T1_Vox1)
wilcox.test(Prox$prox_CSO1_T2_Vox1,Prox$prox_CSO2_T2_Vox1)

summary(Prox$prox_CSO1_T1_Vox1[Prox$prox_CSO1_T1_Vox1 > 0])
sd(Prox$prox_CSO1_T1_Vox1[Prox$prox_CSO1_T1_Vox1 > 0])
summary(Prox$prox_CSO2_T1_Vox1[Prox$prox_CSO2_T1_Vox1 > 0])
sd(Prox$prox_CSO2_T1_Vox1[Prox$prox_CSO2_T1_Vox1 > 0])

summary(Prox$prox_CSO1_T2_Vox1[Prox$prox_CSO1_T2_Vox1 > 0])
sd(Prox$prox_CSO1_T2_Vox1[Prox$prox_CSO1_T2_Vox1 > 0])
summary(Prox$prox_CSO2_T2_Vox1[Prox$prox_CSO2_T2_Vox1 > 0])
sd(Prox$prox_CSO2_T2_Vox1[Prox$prox_CSO2_T2_Vox1 > 0])

summary(Prox$prox_CSO_T1_Vox1[Prox$prox_CSO_T1_Vox1 > 0])
sd(Prox$prox_CSO_T1_Vox1[Prox$prox_CSO_T1_Vox1 > 0])
summary(Prox$prox_CSO_T2_Vox1[Prox$prox_CSO_T2_Vox1 > 0])
sd(Prox$prox_CSO_T2_Vox1[Prox$prox_CSO_T2_Vox1 > 0])

summary(Prox$prox_BG_T1_Vox1[Prox$prox_BG_T1_Vox1 > 0])
sd(Prox$prox_BG_T1_Vox1[Prox$prox_BG_T1_Vox1 > 0])
summary(Prox$prox_BG_T2_Vox1[Prox$prox_BG_T2_Vox1 > 0])
sd(Prox$prox_BG_T2_Vox1[Prox$prox_BG_T2_Vox1 > 0])

summary(Prox_sub$prox_CSO1_T1_Vox1)
sd(Prox_sub$prox_CSO1_T1_Vox1)
summary(Prox_sub$prox_CSO2_T1_Vox1)
sd(Prox_sub$prox_CSO2_T1_Vox1)

summary(Prox_sub$prox_CSO1_T2_Vox1)
sd(Prox_sub$prox_CSO1_T2_Vox1)
summary(Prox_sub$prox_CSO2_T2_Vox1)
sd(Prox_sub$prox_CSO2_T2_Vox1)

t.test(Prox_sub$prox_CSO1_T1, Prox_sub$prox_CSO2_T1, alternative = c("two.sided"), paired = F)
t.test(Prox_sub$prox_CSO1_T2, Prox_sub$prox_CSO2_T2, alternative = c("two.sided"), paired = F)



#proximity analysis STOPMS
STOPMS$prox_BG_T1 <- (STOPMS$BGmaskbinexp - STOPMS$BGminusT1maskbin)/(STOPMS$BGmaskbinexp+0.0001)
STOPMS$prox_CSO_T1 <- (STOPMS$CSOmaskbinexp - STOPMS$CSOminusT1maskbin)/(STOPMS$CSOmaskbinexp+0.0001)

hist(STOPMS$prox_BG_T1, breaks = 100)
hist(STOPMS$prox_CSO_T1, breaks = 100)

summary(STOPMS$prox_BG_T1)
sd(STOPMS$prox_BG_T1, na.rm = T)
summary(STOPMS$prox_CSO_T1)
sd(STOPMS$prox_CSO_T1, na.rm = T)








#### Validation cohort analysis ####

#data preparation
str(EPVS_validation$T1_lesions)
EPVS_validation$CSOT <- as.numeric(EPVS_validation$CSOT)
EPVS_validation$BGCT <- as.numeric(EPVS_validation$BGCT)
EPVS_validation$BGAT <- as.numeric(EPVS_validation$BGAT)
EPVS_validation$Max_D_CSOT <- as.numeric(EPVS_validation$Max_D_CSOT)
EPVS_validation$CSOT2plus <- as.numeric(EPVS_validation$CSOT2plus)
EPVS_validation$BST <- as.numeric(EPVS_validation$BST)
EPVS_validation$BG_Vol <- as.numeric(EPVS_validation$BG_Vol)
EPVS_validation$CSO1_Vol <- as.numeric(EPVS_validation$CSO1_Vol)
EPVS_validation$CSO2_Vol <- as.numeric(EPVS_validation$CSO2_Vol)
EPVS_validation$CSO_Total_Vol <- as.numeric(EPVS_validation$CSO_Total_Vol)
EPVS_validation$BPF <- as.numeric(EPVS_validation$BPF)
EPVS_validation$WMF <- as.numeric(EPVS_validation$WMF)
EPVS_validation$GMF <- as.numeric(EPVS_validation$GMF)
EPVS_validation$T1_lesions <- as.numeric(EPVS_validation$T1_lesions)
EPVS_validation$T2_volume_ml <- as.numeric(EPVS_validation$T2_volume_ml)
EPVS_validation$Sex <- as.factor(EPVS_validation$Sex)
EPVS_validation$Cohort <- as.factor(EPVS_validation$Cohort)
EPVS_validation$age <- as.numeric(EPVS_validation$age)
EPVS_validation$T1_lesions <- EPVS_validation$T1_lesions/1000

ggplot(EPVS_validation, aes(x=T1_lesions, y=T2_volume_ml)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)
cor.test(EPVS_validation$T1_lesions, EPVS_validation$T2_volume_ml)#T1 and T2 lesions are highly correlated

plot(EPVS_validation$CSOT)

#create new variable for total EPVS
EPVS_validation$EPVS_vol <- EPVS_validation$CSO_Total_Vol+EPVS_validation$BG_Vol
EPVS_validation$BGtoBrainvolume <- EPVS_validation$BG_Vol/EPVS_validation$Brain_volume_ul
EPVS_validation$CSOtoBrainvolume <- EPVS_validation$CSO_Total_Vol/EPVS_validation$Brain_volume_ul

plot(EPVS_validation$BGtoBrainvolume)
hist(EPVS_validation$BGtoBrainvolume, breaks = 10)

plot(EPVS_validation$CSOtoBrainvolume)
hist(EPVS_validation$CSOtoBrainvolume, breaks = 10)

summary(lm(CSOtoBrainvolume ~ age+Sex+Cohort, data=EPVS_validation))#

#descriptive statistics
table(EPVS_validation$BST)

summary(EPVS_validation$CSOT)
summary(EPVS_validation$BGCT)
summary(EPVS_validation$BGAT)
summary(EPVS_validation$BST)

summary(EPVS_validation$CSO_Total_Vol)
summary(EPVS_validation$BG_Vol)

summary(Long$Vol_PVS_CSO)

table(EPVS_validation$Potter_CSO)
table(EPVS_validation$Potter_BGC)
table(EPVS_validation$Potter_BGA)
table(EPVS_validation$Potter_BS)

table(EPVS_validation$CSOT2plus)
summary(EPVS_validation$CSOT2plus)
sd(EPVS_validation$CSOT2plus)

#Sex -> no significant hits
summary(lm(CSOT ~ age+Sex+Cohort, data=EPVS_validation))#CSO EPVS numbers borderline sig higher when using log(CSOT)
summary(lm(BGAT ~ age+Sex+Cohort, data=EPVS_validation))
summary(lm(BGCT ~ age+Sex+Cohort, data=EPVS_validation))
summary(lm(BST ~ age+Sex+Cohort, data=EPVS_validation))

summary(lm(EPVS_vol ~ age+Sex, data=EPVS_validation))#CSO EPVS volume borderline sig higher when using log(CSOT)
summary(lm(CSO_Total_Vol ~ age+Sex, data=EPVS_validation))#CSO EPVS volume borderline sig higher when using log(CSOT)

summary(lm(Max_D_CSOT ~ age+Sex, data=EPVS_validation))

#EPVS numbers
  #to brain volumetry
  summary(lm(BPF ~ CSOT+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(WMF ~ CSOT+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(GMF ~ CSOT+age+Sex+Cohort, data=EPVS_validation))#ns
  
  summary(lm(BPF ~ BGAT+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(WMF ~ BGAT+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(GMF ~ BGAT+age+Sex+Cohort, data=EPVS_validation))#ns
  
  summary(lm(BPF ~ BGCT+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(WMF ~ BGCT+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(GMF ~ BGCT+age+Sex+Cohort, data=EPVS_validation))#ns
  
  summary(lm(BPF ~ BST+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(WMF ~ BST+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(GMF ~ BST+age+Sex+Cohort, data=EPVS_validation))#ns

  #to lesion measures
  summary(lm(T1_lesions ~ CSOT+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(T2_volume_ml ~ CSOT+age+Sex+Cohort, data=EPVS_validation))#ns
  
  summary(lm(T1_lesions ~ BGAT+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(T2_volume_ml ~ BGAT+age+Sex+Cohort, data=EPVS_validation))#ns

  summary(lm(T1_lesions ~ BGCT+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(T2_volume_ml ~ BGCT+age+Sex+Cohort, data=EPVS_validation))#ns

  summary(lm(T1_lesions ~ BST+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(T2_volume_ml ~ BST+age+Sex+Cohort, data=EPVS_validation))#ns

#EPVS volume
  #to brain volumetry
  summary(lm(BPF ~ CSO1_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(WMF ~ CSO1_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(GMF ~ CSO1_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  
  summary(lm(BPF ~ CSO2_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(WMF ~ CSO2_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(GMF ~ CSO2_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  
  summary(lm(BPF ~ CSO_Total_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(WMF ~ CSO_Total_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(GMF ~ CSO_Total_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  
  summary(lm(BPF ~ EPVS_vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(WMF ~ EPVS_vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(GMF ~ EPVS_vol+age+Sex+Cohort, data=EPVS_validation))#ns
  
  #to lesion measures
  summary(lm(T1_lesions ~ CSO1_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(T2_volume_ml ~ CSO1_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  
  summary(lm(T1_lesions ~ CSO2_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(T2_volume_ml ~ CSO2_Vol+age+Sex+Cohort, data=EPVS_validation))#ns

  summary(lm(T1_lesions ~ CSO_Total_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(T2_volume_ml ~ CSO_Total_Vol+age+Sex+Cohort, data=EPVS_validation))#ns
  
  summary(lm(T1_lesions ~ EPVS_vol+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(T2_volume_ml ~ EPVS_vol+age+Sex+Cohort, data=EPVS_validation))#ns
  
#EPVS diameter
  #to brain volumetry
  summary(lm(BPF ~ CSOT2plus+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(WMF ~ CSOT2plus+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(GMF ~ CSOT2plus+age+Sex+Cohort, data=EPVS_validation))#ns
  
  summary(lm(BPF ~ Max_D_CSOT+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(WMF ~ Max_D_CSOT+age+Sex+Cohort, data=EPVS_validation))#ns
  summary(lm(GMF ~ Max_D_CSOT+age+Sex+Cohort, data=EPVS_validation))#ns
  
  #to lesion measures
  summary(lm(T1_lesions ~ CSOT2plus+age+Sex+Cohort, data=EPVS_validation))#sig*
  summary(lm(T2_volume_ml ~ CSOT2plus+age+Sex+Cohort, data=EPVS_validation))#sig**

  summary(lm(T1_lesions ~ Max_D_CSOT+age+Sex, data=EPVS_validation))#ns
  summary(lm(T2_volume_ml ~ Max_D_CSOT+age+Sex, data=EPVS_validation))#ns

  #subgroup analysis
  summary(lm(T1_lesions ~ CSOT2plus+age+Sex, data=EPVS_validation[EPVS_validation$Cohort == "REMYDI",]))
  summary(lm(T1_lesions ~ CSOT2plus+age+Sex, data=EPVS_validation[EPVS_validation$Cohort == "MULTIPLEMS",]))
  
  summary(lm(T2_volume_ml ~ CSOT2plus+age+Sex, data=EPVS_validation[EPVS_validation$Cohort == "REMYDI",]))
  summary(lm(T2_volume_ml ~ CSOT2plus+age+Sex, data=EPVS_validation[EPVS_validation$Cohort == "MULTIPLEMS",]))
  
  ggplot(EPVS_validation, aes(x=CSOT2plus, y=T2_volume_ml)) +
    geom_point(size=3, shape=20) +
    geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)
  
  ggplot(EPVS_validation[EPVS_validation$Cohort == "REMYDI",], aes(x=CSOT2plus, y=T2_volume_ml)) +
    geom_point(size=3, shape=20) +
    geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)
  
  ggplot(EPVS_validation[EPVS_validation$Cohort == "MULTIPLEMS",], aes(x=CSOT2plus, y=T2_volume_ml)) +
    geom_point(size=3, shape=20) +
    geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)
  
  ggplot(EPVS_validation[EPVS_validation$Cohort == "REMYDI",], aes(x=CSOT2plus, y=T1_lesions)) +
    geom_point(size=3, shape=20) +
    geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)
  
  ggplot(EPVS_validation[EPVS_validation$Cohort == "MULTIPLEMS",], aes(x=CSOT2plus, y=T1_lesions)) +
    geom_point(size=3, shape=20) +
    geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)
  
  EPVS_validation_REM_E <- subset(EPVS_validation_REM, ID != "REMY005_M0")
  summary(lm(T1_lesions ~ CSOT2plus+age+Sex, data=EPVS_validation_REM_E))#when removing REMY005 --> higher sig
  summary(lm(T2_volume_ml ~ CSOT2plus+age+Sex, data=EPVS_validation_REM_E))#when removing REMY005 --> higher sig
  
  EPVS_validation_REM_E2 <- EPVS_validation_REM %>% filter(CSOT2plus < 4)
  summary(lm(T1_lesions ~ CSOT2plus+age+Sex, data=EPVS_validation_REM_E2))#when removing #EPVS >4 --> higher sig
  summary(lm(T2_volume_ml ~ CSOT2plus+age+Sex, data=EPVS_validation_REM_E2))#when removing #EPVS >4 --> higher sig
  
  #all
  T1_diameter <- lm(T1_lesions ~ CSOT2plus+age+Sex, data=EPVS_validation)#without excluding patients
  T1_diameter_fitted <- fitted.values(T1_diameter)
  T2_diameter <- lm(T2_volume_ml ~ CSOT2plus+age+Sex, data=EPVS_validation)#without excluding patients
  T2_diameter_fitted <- fitted.values(T2_diameter)
  
  #REM
  T1_diameter_REM <- lm(T1_lesions ~ CSOT2plus+age+Sex, data=EPVS_validation_REM)#without excluding patients
  T1_diameter_REM_fitted <- fitted.values(T1_diameter_REM)
  T2_diameter_REM <- lm(T2_volume_ml ~ CSOT2plus+age+Sex, data=EPVS_validation_REM)#without excluding patients
  T2_diameter_REM_fitted <- fitted.values(T2_diameter_REM)
  
  #MUL
  T1_diameter_MUL <- lm(T1_lesions ~ CSOT2plus+age+Sex, data=EPVS_validation_MUL)#without excluding patients
  T1_diameter_MUL_fitted <- fitted.values(T1_diameter_MUL)
  T2_diameter_MUL <- lm(T2_volume_ml ~ CSOT2plus+age+Sex, data=EPVS_validation_MUL)#without excluding patients
  T2_diameter_MUL_fitted <- fitted.values(T2_diameter_MUL)
  
#### control cohort analysis ####
EPVS_controls <- read.delim(file = "EPVS_controls.txt", stringsAsFactors = F)

#merge with segmentation
#load and merge files
Seg_controls <- read.table("aseg_stats_STOPMScontrols.txt", sep="\t", header=TRUE)

#prepare variables
Seg_controls <- subset(Seg_controls, select=c("Measure.volume", "Left.Cerebellum.White.Matter", "Right.Cerebellum.White.Matter", 
                                                  "CerebralWhiteMatterVol", "Brain.Stem",
                                                  "CC_Posterior", "CC_Mid_Posterior", "CC_Central", "CC_Mid_Anterior", "CC_Anterior", 
                                                  "BrainSegVolNotVent", "TotalGrayVol", "EstimatedTotalIntraCranialVol", 
                                                  "WM.hypointensities"))
Seg_controls$WM <- Seg_controls$Left.Cerebellum.White.Matter+
  Seg_controls$Right.Cerebellum.White.Matter+
  Seg_controls$CerebralWhiteMatterVol+
  Seg_controls$Brain.Stem+
  Seg_controls$CC_Posterior+Seg_controls$CC_Mid_Posterior+
  Seg_controls$CC_Central+Seg_controls$CC_Mid_Anterior+Seg_controls$CC_Anterior
Seg_controls$BPF <- Seg_controls$BrainSegVolNotVent/Seg_controls$EstimatedTotalIntraCranialVol
Seg_controls$GMF <- Seg_controls$TotalGrayVol/Seg_controls$EstimatedTotalIntraCranialVol
Seg_controls$WMF <- Seg_controls$WM/Seg_controls$EstimatedTotalIntraCranialVol
Seg_controls <- subset(Seg_controls, select=c("Measure.volume", "BPF", "GMF", "WMF", "WM.hypointensities", "BrainSegVolNotVent"))
names(Seg_controls) <- c("ID", "BPF", "GMF", "WMF", "T1_lesions", "Brain_volume_ul")
Seg_controls
  
#load validation cohort excel file
EPVS_controls <- merge(EPVS_controls, Seg_controls, by="ID")
  
#adjust variables
EPVS_controls$ID <- as.factor(EPVS_controls$ID)
EPVS_controls$Sex <- as.factor(EPVS_controls$Sex)
EPVS_controls$Cohort <- as.factor(EPVS_controls$Cohort)
EPVS_controls$CSOR <- as.numeric(EPVS_controls$CSOR)
EPVS_controls$CSOL <- as.numeric(EPVS_controls$CSOL)
EPVS_controls$CSOT <- as.numeric(EPVS_controls$CSOT)
EPVS_controls$BGCR <- as.numeric(EPVS_controls$BGCR)
EPVS_controls$BGCL <- as.numeric(EPVS_controls$BGCL)
EPVS_controls$BGCT <- as.numeric(EPVS_controls$BGCT)
EPVS_controls$BGAR <- as.numeric(EPVS_controls$BGAR)
EPVS_controls$BGAL <- as.numeric(EPVS_controls$BGAL)
EPVS_controls$BGAT <- as.numeric(EPVS_controls$BGAT)
EPVS_controls$BST <- as.numeric(EPVS_controls$BST)
EPVS_controls$CSOT2plus <- as.numeric(EPVS_controls$CSOT2plus)
EPVS_controls$Potter_CSO <- as.factor(EPVS_controls$Potter_CSO)
EPVS_controls$Potter_BGA <- as.factor(EPVS_controls$Potter_BGA)
EPVS_controls$Potter_BGC <- as.factor(EPVS_controls$Potter_BGC)
EPVS_controls$Potter_BST <- as.factor(EPVS_controls$Potter_BST)

#descriptive statistics
table(EPVS_controls$Sex)
summary(EPVS_controls$age)
sd(EPVS_controls$age)#mean age 45 --> very high, have to filter controls for age
EPVS_controls <- EPVS_controls %>% filter(age < 57)
table(EPVS_controls$Sex)
summary(EPVS_controls$age)
sd(EPVS_controls$age)
t.test(EPVS_controls$age, EPVS_validation$age)#not significant different from validation cohort
t.test(EPVS_controls$age, Long$Age_Scan)#not significant different from validation cohort

summary(EPVS_controls$BPF)
IQR(EPVS_controls$BPF)

summary(EPVS_controls$CSOT)
summary(EPVS_controls$BGCT)
summary(EPVS_controls$BGAT)
summary(EPVS_controls$BST)

summary(EPVS_controls$BG_Vol)
summary(EPVS_controls$CSO_Total_Vol)

table(EPVS_controls$Potter_CSO)
table(EPVS_controls$Potter_BGC)
table(EPVS_controls$Potter_BGA)
table(EPVS_controls$Potter_BST)

table(EPVS_controls$CSOT2plus)
summary(EPVS_controls$CSOT2plus)
sd(EPVS_controls$CSOT2plus)

summary(Long$CSO_D_Number)
sd(Long$CSO_D_Number, na.rm = T)

#statistics
Primary <- subset(Long, select = c("HIVE", "CSOT", "BGAT", "BGCT", "BS", "Age_Scan", "Sex", "Vol_PVS_CSO", "Vol_PVS_BG", "BPF1", "CSO_D_Number"))
Primary <- Primary %>% mutate(cohort = "Primary")
names(Primary) <- c("ID", "CSOT", "BGAT", "BGCT", "BST", "Age", "Sex", "CSO_Vol", "BG_Vol", "BPF", "Num_diameter", "Cohort")
Primary <- Primary %>% filter(!is.na(CSOT))
Primary$Sex <- gsub("F", "Female", Primary$Sex)
Primary$Sex <- gsub("M", "Male", Primary$Sex)

Validation <- subset(EPVS_validation, select = c("ID", "CSOT", "BGAT", "BGCT", "BST", "age", "Sex", "CSO_Total_Vol", "BG_Vol", "BPF", "CSOT2plus"))
Validation <- Validation %>% mutate(cohort = "Validation")
names(Validation) <- c("ID", "CSOT", "BGAT", "BGCT", "BST", "Age", "Sex", "CSO_Vol", "BG_Vol", "BPF", "Num_diameter", "Cohort")
Validation <- Validation %>% filter(!is.na(ID))

Control <- subset(EPVS_controls, select = c("ID", "CSOT", "BGAT", "BGCT", "BST", "age", "Sex", "CSO_Total_Vol", "BG_Vol", "BPF", "CSOT2plus"))
Control <- Control %>% mutate(cohort = "Control")
names(Control) <- c("ID", "CSOT", "BGAT", "BGCT", "BST", "Age", "Sex", "CSO_Vol", "BG_Vol", "BPF", "Num_diameter", "Cohort")
Control <- Control %>% filter(!is.na(ID))

Prim_Ctrl<- rbind(Primary, Control)
Val_Ctrl <- rbind(Validation, Control)
Prim_Val <- rbind(Primary, Validation)

summary(lm(CSO_Vol ~ Cohort+Age+Sex, data=Prim_Ctrl))#ns
summary(lm(CSO_Vol ~ Cohort+Age+Sex, data=Val_Ctrl))#sig***

summary(lm(BG_Vol ~ Cohort+Age+Sex, data=Prim_Ctrl))#ns
summary(lm(BG_Vol ~ Cohort+Age+Sex, data=Val_Ctrl))#ns

summary(lm(CSOT ~ Cohort+Age+Sex, data=Prim_Ctrl))#sig**
summary(lm(CSOT ~ Cohort+Age+Sex, data=Val_Ctrl))#sig***

summary(lm(BGAT ~ Cohort+Age+Sex, data=Prim_Ctrl))#sig***
summary(lm(BGAT ~ Cohort+Age+Sex, data=Val_Ctrl))#sig***

summary(lm(BGCT ~ Cohort+Age+Sex, data=Prim_Ctrl))#sig***
summary(lm(BGCT ~ Cohort+Age+Sex, data=Val_Ctrl))#ns

summary(lm(BST ~ Cohort+Age+Sex, data=Prim_Ctrl))#ns
summary(lm(BST ~ Cohort+Age+Sex, data=Val_Ctrl))#sig* but minus

summary(lm(Num_diameter ~ Cohort+Age+Sex, data=Prim_Ctrl))#sig*
summary(lm(Num_diameter ~ Cohort+Age+Sex, data=Val_Ctrl))#sig**

summary(lm(CSOT ~ Age+Sex, data=Prim_Val))#sig*
summary(lm(CSOT ~ Age+Sex, data=Primary))#
summary(lm(CSOT ~ Age+Sex, data=Validation))#
summary(lm(CSOT ~ Age+Sex, data=Control))#

summary(lm(BGAT ~ Age+Sex, data=Prim_Val))#ns
summary(lm(BGAT ~ Age+Sex, data=Primary))#ns
summary(lm(BGAT ~ Age+Sex, data=Validation))#ns
summary(lm(BGAT ~ Age+Sex, data=Control))#ns

summary(lm(BGCT ~ Age+Sex, data=Prim_Val))#
summary(lm(BGCT ~ Age+Sex, data=Primary))#
summary(lm(BGCT ~ Age+Sex, data=Validation))#
summary(lm(BGCT ~ Age+Sex, data=Control))#

summary(lm(CSO_Vol ~ Age+Sex, data=Prim_Val))#Trend
summary(lm(CSO_Vol ~ Age+Sex, data=Primary))#ns
summary(lm(CSO_Vol ~ Age+Sex, data=Validation))#sig**
summary(lm(CSO_Vol ~ Age+Sex, data=Control))#ns


cor.test(Long$CSOT, Long$Vol_PVS_CSO)
cor.test(Long$BGCT, Long$Vol_PVS_BG)
cor.test(Long$BGAT, Long$Vol_PVS_BG)
cor.test(Long$BGACT, Long$Vol_PVS_BG)

cor.test(Prim_Val$BGCT, Prim_Val$BG_Vol)
cor.test(Prim_Val$BGAT, Prim_Val$BG_Vol)

cor.test(Prim_Val$CSOT, Prim_Val$CSO_Vol)
cor.test(Prim_Val$Num_diameter, Prim_Val$CSO_Vol)
cor.test(Prim_Val$CSOT, Prim_Val$Num_diameter)


  ###Graphs
  #all
  Plot_T1 <- ggplot(EPVS_validation, aes(x=CSOT2plus, y=T1_diameter_fitted)) +
    geom_point(size=2, shape=16, color="blue") +
    geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
                linetype="dashed", color="darkred") +
    theme_minimal(base_size = 18) +
    theme(axis.text.x = element_text(color = "black", 
                                     size = 24),
          axis.text.y = element_text(color = "black", 
                                     size = 24), 
          axis.title = element_text(face = "bold", color = "black", 
                                    size = 24)) +
    ylab("T1 lesion volume (ml)") +
    xlab("#of dilated EPVS") +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
    scale_y_continuous(name="T1 lesion volume (ml)")
  
  Plot_T2 <- ggplot(EPVS_validation, aes(x=CSOT2plus, y=T2_diameter_fitted)) +
    geom_point(size=2, shape=16, color="blue") +
    geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
                linetype="dashed", color="darkred") +
    theme_minimal(base_size = 18) +
    theme(axis.text.x = element_text(color = "black", 
                                     size = 24),
          axis.text.y = element_text(color = "black", 
                                     size = 24), 
          axis.title = element_text(face = "bold", color = "black", 
                                    size = 24)) +
    ylab("T2 lesion volume (ml)") +
    xlab("#of dilated EPVS") +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
    scale_y_continuous(name="T2 lesion volume (ml)")
  
  
  #REMYDI
  Plot_T1_REM <- ggplot(EPVS_validation_REM, aes(x=CSOT2plus, y=T1_diameter_REM_fitted)) +
    geom_point(size=2, shape=16, color="blue") +
    geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
                linetype="dashed", color="darkred") +
    theme_minimal(base_size = 18) +
    theme(axis.text.x = element_text(face = "bold", color = "black", 
                                     size = 14),
          axis.text.y = element_text(face = "bold", color = "black", 
                                     size = 14)) +
    ylab("T1 lesion volume (ml)") +
    xlab("#of dilated EPVS") +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
    scale_y_continuous(name="T1 lesion volume (ml)")

  Plot_T2_REM <- ggplot(EPVS_validation_REM, aes(x=CSOT2plus, y=T2_diameter_REM_fitted)) +
    geom_point(size=2, shape=16, color="blue") +
    geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
                linetype="dashed", color="darkred") +
    theme_minimal(base_size = 18) +
    theme(axis.text.x = element_text(face = "bold", color = "black", 
                                     size = 14),
          axis.text.y = element_text(face = "bold", color = "black", 
                                     size = 14)) +
    ylab("T2 lesion volume (ml)") +
    xlab("#of dilated EPVS") +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
    scale_y_continuous(name="T2 lesion volume (ml)")
  
  #MULTIPLEMS
  Plot_T1_MUL <- ggplot(EPVS_validation_MUL, aes(x=CSOT2plus, y=T1_diameter_MUL_fitted)) +
    geom_point(size=2, shape=16, color="blue") +
    theme_minimal(base_size = 18) +
    theme(axis.text.x = element_text(face = "bold", color = "black", 
                                     size = 14),
          axis.text.y = element_text(face = "bold", color = "black", 
                                     size = 14)) +
    ylab("T1 lesion volume (ml)") +
    xlab("#of dilated EPVS") +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
    scale_y_continuous(name="T1 lesion volume (ml)")
  
  Plot_T2_MUL <- ggplot(EPVS_validation_MUL, aes(x=CSOT2plus, y=T2_diameter_MUL_fitted)) +
    geom_point(size=2, shape=16, color="blue") +
    theme_minimal(base_size = 18) +
    theme(axis.text.x = element_text(face = "bold", color = "black", 
                                     size = 14),
          axis.text.y = element_text(face = "bold", color = "black", 
                                     size = 14)) +
    ylab("T2 lesion volume (ml)") +
    xlab("#of dilated EPVS") +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
    scale_y_continuous(name="T2 lesion volume (ml)")
  
  
  
  
  
  
#######################
########Paper##########
#######################

#Methods
citation()
R.Version()

#Results
summary(Long$CSO_D)
summary(Long$CSOT)
summary(Long$BGCT)
summary(Long$BGAT)
summary(Long$BS)

summary(as.factor(Long$Potter_CSO))
summary(as.factor(Long$Potter_BGC))
summary(as.factor(Long$Potter_BGA))

summary(lm(logCSOT ~ Age_Scan+Sex, data=PVS))
summary(lm(logBGCT ~ Age_Scan+Sex, data=PVS))

ggplot(Long, aes(x=Sex, y=BGCT)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(lm(SupraTentorialVolNotVent ~ logCSOT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(TotalGrayVol ~ logCSOT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(CortexVol ~ logCSOT+Age_Scan+Sex+MS_type, data=PVS))

summary(lm(Supratent_D ~ logCSOT+Age_Scan+Sex+MS_type, data=Long))
summary(lm(Gray_D ~ logCSOT+Age_Scan+Sex+MS_type, data=Long))

summary(lm(T2_Count_FS ~ logCSOT+Age_Scan+Sex+MS_type, data=PVS))
summary(lm(T2_Count_FS ~ logCSOT+Age_Scan+Sex+MS_type, data=PVS))

ggplot(Long, aes(x=logCSOT, y=SupraTentorialVolNotVent)) +
  geom_point(size=3, shape=20) +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95)

summary(Long$TotalGrayVol)
SD(Long$TotalGrayVol, na.rm=TRUE)
summary(Long$SupraTentorialVolNotVent)
SD(Long$SupraTentorialVolNotVent, na.rm=TRUE)
summary(Long$Supratent_D)
summary(Long$Gray_D)

t.test(SDMT_Z_E2~CSO_D_2, data=Long)
summary(lm(SDMT_Z_E2 ~ CSO_D_2+Sex+MS_type+Duration, data=Long))
summary(lm(SDMT_Z_E2 ~ CSOT+Sex+MS_type+Duration, data=Long))
wilcox.test(Long$SDMT_Z~Long$CSO_D_2) 

ggplot(Long, aes(x=CSO_D_2, y=SDMT_Z_E2)) +
  geom_point(size=3, shape=20) +
  stat_summary(geom = "point", fun.y = "mean", col = "black",
    size = 3, shape = 24, fill = "red")

summary(Long$SDMT_Z_E2)
SD(Long$SDMT_Z_E2, na.rm=TRUE)

summary(Long4$WMH1_E1)
IQR(Long4$WMH1_E1)
summary(Long4$T2_Volume_FS_E)
IQR(Long4$T2_Volume_FS_E)
summary(Long4$BPF1)
IQR(Long4$BPF1)

summary(Long4$EDSS)

#Validation cohort tables
summary(EPVS_validation$Sex)
summary(EPVS_validation$Cohort)
summary(EPVS_validation$age)
sd(EPVS_validation$age)
summary(EPVS_validation$T1_lesions)
IQR(EPVS_validation$T1_lesions)
summary(EPVS_validation$T2_volume_ml)
IQR(EPVS_validation$T2_volume_ml)
summary(EPVS_validation$BPF)
IQR(EPVS_validation$BPF)



EPVS_validation_REM <- subset(EPVS_validation, Cohort == "REMYDI")
EPVS_validation_MUL <- subset(EPVS_validation, Cohort == "MULTIPLEMS")

summary(EPVS_validation_REM$Sex)
summary(EPVS_validation_REM$age)
sd(EPVS_validation_REM$age)
summary(EPVS_validation_REM$T1_lesions)
IQR(EPVS_validation_REM$T1_lesions)
summary(EPVS_validation_REM$T2_volume_ml)
IQR(EPVS_validation_REM$T2_volume_ml)
summary(EPVS_validation_REM$BPF)
IQR(EPVS_validation_REM$BPF)

summary(EPVS_validation_MUL$Sex)
summary(EPVS_validation_MUL$age)
sd(EPVS_validation_MUL$age)
summary(EPVS_validation_MUL$T1_lesions)
IQR(EPVS_validation_MUL$T1_lesions)
summary(EPVS_validation_MUL$T2_volume_ml)
IQR(EPVS_validation_MUL$T2_volume_ml)
summary(EPVS_validation_MUL$BPF)
IQR(EPVS_validation_MUL$BPF)

summary(EPVS_validation$CSOT)
IQR(EPVS_validation$CSOT)
summary(EPVS_validation$BGCT)
IQR(EPVS_validation$BGCT)
summary(EPVS_validation$BGAT)
IQR(EPVS_validation$BGAT)
summary(EPVS_validation$BST)
IQR(EPVS_validation$BST)

summary(EPVS_validation$CSO_Total_Vol)
IQR(EPVS_validation$CSO_Total_Vol)
summary(EPVS_validation$BG_Vol)
IQR(EPVS_validation$BG_Vol)

#Graphs
  #SDMT graph
Plot_SDMT <- ggplot(Long4, aes(x=CSO_D_Bin3, y=SDMT_binned_fitted)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 14),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 14)) +
  xlab("EPVS diameter") +
  ylab("SDMT z-score fitted")

  #Lesion graphs
    #T2 volume
Plot_volume <- ggplot(Long2, aes(x=CSO_D_Number, y=Vol_number_fitted)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 24),
        axis.text.y = element_text(color = "black", 
                                   size = 24), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 24)) +
  ylab("T2 lesion volume (ml)") +
  xlab("#of dilated EPVS") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_continuous(name="T2 lesion volume (ml)", limits=c(0, 15))


    #T2 count
Plot_count <- ggplot(Long2, aes(x=CSO_D_Number, y=Cou_number_fitted)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 14),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 14)) +
  ylab("T2 lesion count fitted") +
  xlab("#of dilated EPVS") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_continuous(name="T2 lesion count fitted", limits=c(10, 30))

Plot_WMH <- ggplot(Long4, aes(x=CSO_D_Number, y=WMH_fitted1000)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 24),
        axis.text.y = element_text(color = "black", 
                                   size = 24), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 24)) +
  ylab("T1 lesion volume (ml)") +
  xlab("#of dilated EPVS") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_continuous(name="T1 lesion volume (ml)")


### EPVS prevalence ###
Primary <- subset(Long, select = c("HIVE", "CSOT", "BGAT", "BGCT", "BS", "Age_Scan", "Sex", "Vol_PVS_CSO", "Vol_PVS_BG", "BPF1", "CSO_D_Number"))
Primary <- Primary %>% mutate(cohort = "Primary")
names(Primary) <- c("ID", "CSOT", "BGAT", "BGCT", "BST", "Age", "Sex", "CSO_Vol", "BG_Vol", "BPF", "Num_diameter", "Cohort")
Primary <- Primary %>% filter(!is.na(ID))

Validation <- subset(EPVS_validation, select = c("ID", "CSOT", "BGAT", "BGCT", "BST", "age", "Sex", "CSO_Total_Vol", "BG_Vol", "BPF", "CSOT2plus"))
Validation <- Validation %>% mutate(cohort = "Validation")
names(Validation) <- c("ID", "CSOT", "BGAT", "BGCT", "BST", "Age", "Sex", "CSO_Vol", "BG_Vol", "BPF", "Num_diameter", "Cohort")
Validation <- Validation %>% filter(!is.na(ID))

Control <- subset(EPVS_controls, select = c("ID", "CSOT", "BGAT", "BGCT", "BST", "age", "Sex", "CSO_Total_Vol", "BG_Vol", "BPF", "CSOT2plus"))
Control <- Control %>% mutate(cohort = "Control")
names(Control) <- c("ID", "CSOT", "BGAT", "BGCT", "BST", "Age", "Sex", "CSO_Vol", "BG_Vol", "BPF", "Num_diameter", "Cohort")
Control <- Control %>% filter(!is.na(ID))

boxplot_df <- rbind(Primary, Control, Validation)

boxplot_df1 <- boxplot_df %>% mutate(Class = "SOC")#Semioval center
boxplot_df1 <- subset(boxplot_df1, select = c("CSOT", "Cohort", "Class"))
names(boxplot_df1) <- c("EPVScount", "Cohort", "Class")

boxplot_df2 <- boxplot_df %>% mutate(Class = "BGA")#Anterior commissure
boxplot_df2 <- subset(boxplot_df2, select = c("BGAT", "Cohort", "Class"))
names(boxplot_df2) <- c("EPVScount", "Cohort", "Class")

boxplot_df3 <- boxplot_df %>% mutate(Class = "BGC")#Crus anterius
boxplot_df3 <- subset(boxplot_df3, select = c("BGCT", "Cohort", "Class"))
names(boxplot_df3) <- c("EPVScount", "Cohort", "Class")

boxplot_df4 <- boxplot_df %>% mutate(Class = "BST")
boxplot_df4 <- subset(boxplot_df4, select = c("BST", "Cohort", "Class"))
names(boxplot_df4) <- c("EPVScount", "Cohort", "Class")

boxplot_df <- rbind(boxplot_df1, boxplot_df2, boxplot_df3)

Plot_EPVScount <- boxplot_df %>% 
  ggplot(aes(x=Class, y=EPVScount, fill=Cohort)) +
  geom_boxplot(outlier.alpha = F, size = 0.65) + 
  stat_boxplot(geom = "errorbar", colour = "black") +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), alpha=0.2, size = 1.25) +
  labs(x="", y="EPVS count per location") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 22),
        axis.text.y = element_text(color = "black", 
                                   size = 22), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 22)) +
  scale_fill_manual(values=rep(c("#3399FF", "#FF9900", "#FF3300"), 3))

Primary2 <- Primary
Primary2 <- Primary2 %>% 
  mutate(Cohort = str_replace(Cohort, "Primary", "Exploratory"))
Primary2$CSO_Vol <- Primary2$CSO_Vol+80
boxplot_df <- rbind(Primary2, Control, Validation)

boxplot2_df1 <- boxplot_df %>% mutate(Class = "CSO")#semioval center
boxplot2_df1 <- subset(boxplot2_df1, select = c("CSO_Vol", "Cohort", "Class"))
names(boxplot2_df1) <- c("EPVSvolume", "Cohort", "Class")

boxplot2_df2 <- boxplot_df %>% mutate(Class = "BG")#Basal ganglia
boxplot2_df2 <- subset(boxplot2_df2, select = c("BG_Vol", "Cohort", "Class"))
names(boxplot2_df2) <- c("EPVSvolume", "Cohort", "Class")

boxplot2_df <- rbind(boxplot2_df1, boxplot2_df2)

Plot_EPVSvolume <- boxplot2_df %>% 
  ggplot(aes(x=Class, y=EPVSvolume, fill=Cohort)) +
  geom_boxplot(outlier.alpha = F, size = 0.65) + 
  stat_boxplot(geom = "errorbar", colour = "black") +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), alpha=0.2, size = 1.25) +
  labs(x="", y="EPVS volume per location (ul)") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 22),
        axis.text.y = element_text(color = "black", 
                                   size = 22), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 22)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_fill_manual(values=rep(c("#3399FF", "#FF9900", "#FF3300"), 2))


boxplot_df <- rbind(Primary, Control, Validation)
boxplot_df <- boxplot_df %>% filter(!is.na(Num_diameter))

Plot_EPVSdiameter <- ggplot(data=boxplot_df, aes(x=Cohort, y=Num_diameter, fill=Cohort)) +
  geom_bar(stat="summary", fun.y = "mean", position="dodge", colour="black", size=0.6) +
  geom_point(position=position_jitterdodge(jitter.width = 0.4), alpha=0.2, size = 1.25) +
  scale_fill_brewer(palette="Set1")+
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 22),
        axis.text.y = element_text(color = "black", 
                                   size = 22), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 22)) +
  ylab("Number of dilated EPVS") +
  xlab("") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_fill_manual(values=c("#3399FF", "#FF9900", "#FF3300")) +
  theme(axis.text.x = element_blank())
  

ggsave(plot = Plot_EPVScount, width = 9, height = 5.5, filename = "Plot_EPVScount.jpeg")
ggsave(plot = Plot_EPVSvolume, width = 7, height = 5.5, filename = "Plot_EPVSvolume.jpeg")
ggsave(plot = Plot_EPVSdiameter, width = 4, height = 5.2, filename = "Plot_EPVSdiameter.jpeg")

median(boxplot_df$CSOT[boxplot_df$Cohort == "Primary"])
range(boxplot_df$CSOT[boxplot_df$Cohort == "Primary"])
median(boxplot_df$CSOT[boxplot_df$Cohort == "Validation"])
range(boxplot_df$CSOT[boxplot_df$Cohort == "Validation"])
median(boxplot_df$CSOT[boxplot_df$Cohort == "Control"])
range(boxplot_df$CSOT[boxplot_df$Cohort == "Control"])
median(boxplot_df$BGAT[boxplot_df$Cohort == "Primary"])
range(boxplot_df$BGAT[boxplot_df$Cohort == "Primary"])
median(boxplot_df$BGAT[boxplot_df$Cohort == "Validation"])
range(boxplot_df$BGAT[boxplot_df$Cohort == "Validation"])
median(boxplot_df$BGAT[boxplot_df$Cohort == "Control"])
range(boxplot_df$BGAT[boxplot_df$Cohort == "Control"])

STOPMS_prev <- c(96, 100, 100, 31, 100, 90, 100, 24, 100, 100, 100, 47)
STOPMS_EPVS <- rep(c("CSO", "BGC", "BGA", "BS"), 3)
Cohort <- c(rep("Primary", 4), rep("Validation", 4), rep("Control", 4))
STOPMS_prevalence <- data.frame(STOPMS_prev, STOPMS_EPVS, Cohort)

ggplot(data=STOPMS_prevalence, aes(x=STOPMS_EPVS, y=STOPMS_prev, fill=Cohort)) +
  geom_boxplot(stat="identity", position=position_dodge(), colour="black", size=1) +
  geom_errorbar(aes(ymin=Prox_percentage, ymax=Prox_sd), width=.2,
                position=position_dodge(.9), size = 1) +
  scale_fill_brewer(palette="Set2")+
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 18),
        axis.text.y = element_text(color = "black", 
                                   size = 18), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 18)) +
  ylab("%of EPVS volume in vicinity to lesion") +
  xlab("")

#### bar plots for proximity ####
Prox_percentage <- c(0.3, 0.15, 0.26, 0.16, 0.59, 0.08, 0.05, 0.05)
Prox_sd <- c(0.67, 0.39, 0.75, 0.38, 1.1, 0.42, 0.28, 0.30)
Prox_EPVS <- c(rep("CSO", 2), rep("CSO1", 2), rep("CSO2", 2), rep("BG", 2))
Lesions <- rep(c("T1 lesions", "T2 lesions"), 4)
Prox_bar <- data.frame(Prox_percentage, Prox_sd, Prox_EPVS, Lesions)

Validation <- ggplot(data=Prox_bar, aes(x=Prox_EPVS, y=Prox_percentage, fill=Lesions)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", size=1) +
  geom_errorbar(aes(ymin=Prox_percentage, ymax=Prox_sd), width=.2,
                position=position_dodge(.9), size = 1) +
  scale_fill_brewer(palette="Paired")+
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 18),
        axis.text.y = element_text(color = "black", 
                                   size = 18), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 18)) +
  ylab("%of EPVS volume in vicinity to lesion") +
  xlab("")

ggsave(plot = Validation, width = 10, height = 5.5, filename = "Validation.jpeg")

STOPMS_percentage <- c(0.5, 0.28, 0.67, 0.63)
STOPMS_sd <- c(1.68, 0.54, 3.29, 2.28)
STOPMS_EPVS <- c(rep("CSO", 2), rep("BG", 2))
Lesions <- rep(c("T1 lesions", "T2 lesions"), 2)
STOPMS_bar <- data.frame(STOPMS_percentage, STOPMS_sd, STOPMS_EPVS, Lesions)

STOPMS <- ggplot(data=STOPMS_bar, aes(x=STOPMS_EPVS, y=STOPMS_percentage, fill=Lesions)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", size=1) +
  geom_errorbar(aes(ymin=STOPMS_percentage, ymax=STOPMS_sd), width=.2,
                position=position_dodge(.9), size = 1) +
  scale_fill_brewer(palette="Paired")+
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 18),
        axis.text.y = element_text(color = "black", 
                                   size = 18), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 18)) +
  ylab("%of EPVS volume in vicinity to lesion") +
  xlab("")

ggsave(plot = STOPMS, width = 5, height = 5.5, filename = "STOPMS.jpeg")

### boxplot proximity
Prox_T1_SOC <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                   "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO_T1_Vox1"))
names(Prox_T1_SOC) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                    "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_Vox1")
Prox_T1_SOC$type <- c("SOC")
Prox_T1_SOC1 <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO1_T1_Vox1"))
names(Prox_T1_SOC1) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                    "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_Vox1")
Prox_T1_SOC1$type <- c("SOC1")
Prox_T1_SOC2 <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO2_T1_Vox1"))
names(Prox_T1_SOC2) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                    "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_Vox1")
Prox_T1_SOC2$type <- c("SOC2")
Prox_T1_BG <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_BG_T1_Vox1"))
names(Prox_T1_BG) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                    "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_Vox1")
Prox_T1_BG$type <- c("BG")
Prox_T1 <- rbind(Prox_T1_SOC, Prox_T1_SOC1, Prox_T1_SOC2, Prox_T1_BG)
Prox_T1$lesiontype <- c("T1")

Prox_T2_SOC <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO_T2_Vox1"))
names(Prox_T2_SOC) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_Vox1")
Prox_T2_SOC$type <- c("SOC")
Prox_T2_SOC1 <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                        "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO1_T2_Vox1"))
names(Prox_T2_SOC1) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_Vox1")
Prox_T2_SOC1$type <- c("SOC1")
Prox_T2_SOC2 <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                        "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO2_T2_Vox1"))
names(Prox_T2_SOC2) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_Vox1")
Prox_T2_SOC2$type <- c("SOC2")
Prox_T2_BG <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                      "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_BG_T2_Vox1"))
names(Prox_T2_BG) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_Vox1")
Prox_T2_BG$type <- c("BG")
Prox_T2 <- rbind(Prox_T2_SOC, Prox_T2_SOC1, Prox_T2_SOC2, Prox_T2_BG)
Prox_T2$lesiontype <- c("T2")

Prox_box <- rbind(Prox_T1, Prox_T2)
Prox_box$prox_Vox1 <- Prox_box$prox_Vox1*100

Prox_box %>% 
  ggplot(aes(x=type, y=prox_Vox1, fill=lesiontype)) +
  geom_boxplot(outlier.alpha = F, size = 0.65) + 
  stat_boxplot(geom = "errorbar", colour = "black") +
  scale_fill_brewer(palette="Paired") +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), alpha=0.2, size = 1.25, color= "black") +
  labs(x="", y="%of EPVS in vicinity to lesion") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 22),
        axis.text.y = element_text(color = "black", 
                                   size = 22), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 22))


Prox_T1_SOC <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO_T1"))
names(Prox_T1_SOC) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                        "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox")
Prox_T1_SOC$type <- c("SOC")
Prox_T1_SOC1 <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                        "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO1_T1"))
names(Prox_T1_SOC1) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                         "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox")
Prox_T1_SOC1$type <- c("SOC1")
Prox_T1_SOC2 <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                        "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO2_T1"))
names(Prox_T1_SOC2) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                         "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox")
Prox_T1_SOC2$type <- c("SOC2")
Prox_T1_BG <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                      "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_BG_T1"))
names(Prox_T1_BG) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox")
Prox_T1_BG$type <- c("BG")
Prox_T1 <- rbind(Prox_T1_SOC, Prox_T1_SOC1, Prox_T1_SOC2, Prox_T1_BG)
Prox_T1$lesiontype <- c("T1")

Prox_T2_SOC <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO_T2"))
names(Prox_T2_SOC) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                        "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox")
Prox_T2_SOC$type <- c("SOC")
Prox_T2_SOC1 <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                        "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO1_T2"))
names(Prox_T2_SOC1) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                         "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox")
Prox_T2_SOC1$type <- c("SOC1")
Prox_T2_SOC2 <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                        "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_CSO2_T2"))
names(Prox_T2_SOC2) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                         "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox")
Prox_T2_SOC2$type <- c("SOC2")
Prox_T2_BG <- subset(Prox, select = c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                                      "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox_BG_T2"))
names(Prox_T2_BG) <- c("ID", "cohort", "CSOmaskbin", "CSOmaskbinexp", "CSO1maskbin", "CSO1maskbinexp", "CSO2maskbin", "CSO2maskbinexp",
                       "BGmaskbin", "T1maskbin", "T2maskbin", "BGmaskbinexp", "prox")
Prox_T2_BG$type <- c("BG")
Prox_T2 <- rbind(Prox_T2_SOC, Prox_T2_SOC1, Prox_T2_SOC2, Prox_T2_BG)
Prox_T2$lesiontype <- c("T2")

Prox_box <- rbind(Prox_T1, Prox_T2)
Prox_box$prox <- Prox_box$prox*100

Prox_box %>% 
  ggplot(aes(x=type, y=prox, fill=lesiontype)) +
  geom_boxplot(outlier.alpha = F, size = 0.65) + 
  stat_boxplot(geom = "errorbar", colour = "black") +
  scale_fill_brewer(palette="Paired") +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), alpha=0.2, size = 1.25, color= "black") +
  labs(x="", y="%of EPVS volume in vicinity to lesion") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 20),
        axis.text.y = element_text(color = "black", 
                                   size = 20), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 20)) +
        scale_y_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6),
        label = c("0", "1", "2", "3", "4", "5", "6"))

Prox_box_sub <- Prox_box %>% filter(prox < 4)

Prox_box <- Prox_box_sub %>% 
  ggplot(aes(x=type, y=prox, fill=lesiontype)) +
  geom_boxplot(outlier.alpha = F, size = 0.65) + 
  stat_boxplot(geom = "errorbar", colour = "black") +
  scale_fill_brewer(palette="Paired") +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), alpha=0.2, size = 1.25, color= "black") +
  labs(x="", y="% of VRS volume in vicinity to lesion") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 19),
        axis.text.y = element_text(color = "black", 
                                   size = 19), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 19)) +
  scale_y_continuous(breaks = c(0, 1, 2, 3, 4),
                     label = c("0", "1", "2", "3", "4"))


ggsave(plot = Prox_box, width = 7, height = 5.5, filename = "Prox_box_VRS.jpeg")

#For STOPMS
Prox_T1_SOC_STOPMS <- subset(STOPMS, select = c("ID", "CSOmaskbin", "CSOmaskbinexp",
                                       "BGmaskbin", "T1maskbin", "BGmaskbinexp", "prox_CSO_T1"))
names(Prox_T1_SOC_STOPMS) <- c("ID", "CSOmaskbin", "CSOmaskbinexp",
                                "BGmaskbin", "T1maskbin", "BGmaskbinexp", "prox")
Prox_T1_SOC_STOPMS$type <- c("SOC")

Prox_T1_BG_STOPMS <- subset(STOPMS, select = c("ID", "CSOmaskbin", "CSOmaskbinexp",
                                                "BGmaskbin", "T1maskbin", "BGmaskbinexp", "prox_BG_T1"))
names(Prox_T1_BG_STOPMS) <- c("ID", "CSOmaskbin", "CSOmaskbinexp",
                               "BGmaskbin", "T1maskbin", "BGmaskbinexp", "prox")
Prox_T1_BG_STOPMS$type <- c("BG")

Prox_T1_STOPMS <- rbind(Prox_T1_BG_STOPMS, Prox_T1_SOC_STOPMS)
Prox_T1_STOPMS$prox <- Prox_T1_STOPMS$prox*100

Prox_T1_STOPMS %>% 
  ggplot(aes(x=type, y=prox)) +
  geom_boxplot(outlier.alpha = F, size = 0.65) + 
  stat_boxplot(geom = "errorbar", colour = "black") +
  scale_fill_brewer(palette="Paired") +
  geom_point(position=position_identity(), alpha=0.2, size = 1.25, color= "black") +
  labs(x="", y="% of VRS volume in vicinity to lesion") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 19),
        axis.text.y = element_text(color = "black", 
                                   size = 19), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 19))


Prox_T1_STOPMS_sub <- Prox_T1_STOPMS %>% filter(prox < 10)

Prox_box_STOPMS <- Prox_T1_STOPMS_sub %>% 
  ggplot(aes(x=type, y=prox, fill="blue")) +
  geom_boxplot(outlier.alpha = F, size = 0.65) + 
  stat_boxplot(geom = "errorbar", colour = "black") +
  scale_fill_brewer(palette="Paired") +
  geom_point(position=position_jitter(width = 0.15), alpha=0.2, size = 1.25, color= "black") +
  labs(x="", y="% of VRS volume in vicinity to lesion") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 19),
        axis.text.y = element_text(color = "black", 
                                   size = 19), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 19)) +
  theme(legend.position="none")

ggsave(plot = Prox_box_STOPMS, width = 3.5, height = 5.5, filename = "Prox_box_VRS_STOPMS.jpeg")




install.packages("gridExtra")
library("gridExtra")
grid.arrange(Plot_volume, Plot_count)




multiplot <- function(..., plotlist = NULL, file, cols = 2, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


multiplot(Plot_WMH, Plot_volume)
multiplot(Plot_T1, Plot_T2)
multiplot(Plot_T1_REM, Plot_T2_REM)
multiplot(Plot_T1_MUL, Plot_T2_MUL)
multiplot(Number_brainvolume, Volume_brainvolume, Diameter_brainvolume)
multiplot(Number_T2, Volume_T2)



#Inter-/intra-rater agreement
#Randomization
sample(Long$HIVE, 25, replace=F)

#merge files
Carmen <- read.csv2("F:\\MS\\Experiments\\Karolinska\\perivascular spaces in MS\\Carmen.csv", stringsAsFactors=FALSE)
Ben <- read.csv2("F:\\MS\\Experiments\\Karolinska\\perivascular spaces in MS\\Ben.csv", stringsAsFactors=FALSE)
Rater_merge <- merge(Carmen, Ben, by="HIVE", all=T)
Rater_merge <- Rater_merge %>% filter(!is.na(HIVE.ID))
Rater_merge <- Rater_merge[c(1,3,11,4,14,5,15,6,12,7,13,7,8,9,2)]
Rater_CSOT <- subset(Rater_merge, select=c(CSOT_C, CSOT_B))
Rater_CAT <- subset(Rater_merge, select=c(CAT_C, CAT_B))
Rater_BS <- subset(Rater_merge, select=c(BS_C, BS_B))
Rater_Diameter <- subset(Rater_merge, select=c(Diameter_max_C, Diameter_max_B))
Rater_Number <- subset(Rater_merge, select=c(Number_D_C, Number_D_B))

write_xlsx(Rater_merge, "F:\\MS\\Experiments\\Karolinska\\perivascular spaces in MS\\Rater_merge_kappa.xlsx")

Kappa_CSOT <- read.csv2("G:\\MS\\Experiments\\Karolinska\\perivascular spaces in MS\\rater agreement\\Rater_merge_kappa_CSOT.csv", stringsAsFactors=FALSE)
Kappa_CAT <- read.csv2("G:\\MS\\Experiments\\Karolinska\\perivascular spaces in MS\\rater agreement\\Rater_merge_kappa_CAT.csv", stringsAsFactors=FALSE)
Kappa_BS <- read.csv2("G:\\MS\\Experiments\\Karolinska\\perivascular spaces in MS\\rater agreement\\Rater_merge_kappa_BS.csv", stringsAsFactors=FALSE)
Kappa_Diameter <- read.csv2("G:\\MS\\Experiments\\Karolinska\\perivascular spaces in MS\\rater agreement\\Rater_merge_kappa_Diameter.csv", stringsAsFactors=FALSE)
Kappa_Number <- read.csv2("G:\\MS\\Experiments\\Karolinska\\perivascular spaces in MS\\rater agreement\\Rater_merge_kappa_Number.csv", stringsAsFactors=FALSE)

#inter-rater agreement
#the number of times raters agree divided by the number of things being rated
agree(Kappa_CSOT)
agree(Kappa_CAT)
agree(Kappa_BS)
agree(Kappa_Diameter)
agree(Kappa_Number)

#Cohen's kappa
#Cohens kappa also takes into account agreement by chance
kappa2(Kappa_CSOT)
kappa2(Kappa_CAT)
kappa2(Kappa_BS)
kappa2(Kappa_Diameter)
kappa2(Kappa_Number)


#multiple comparison correction (Benjamini-Hochburg)
p.adjust(c(0.02, 0.31, 0.005, 0.02, 0.0001, 0.01, 0.0001, 0.0001, 0.0001, 0.0001,
           0.01, 0.01, 0.0003, 0.001, 0.0001, 0.03, 0.6), "BH")

#p values for paper
summary(lm(T2_Volume_FS_E ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))#reported in paper
summary(lm(T2_Volume_FS_E ~ CSO_D_Number, data=Long))

summary(lm(WMH1_E1 ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=Long))#reported in paper
summary(lm(WMH1_E1 ~ CSO_D_Number, data=Long))

cor.test(Long$T2_Volume_FS_E, Long$WMH1_E1)





##########################################
########## Postmortem analysis ###########
##########################################


Pre_post_EPVS <- read_excel("G:/MS/Experiments/Karolinska/perivascular spaces in MS/Statistics/PVS/Pre_post_MRI.xlsx", sheet = 3)
Pre_post_EPVS$MRI <- factor(Pre_post_EPVS$MRI,levels = c("Premortem", "Postmortem"))#order bars

ggplot(data=Pre_post_EPVS, aes(x=Mode, y=Volume, fill=MRI)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_minimal(base_size = 18) +
  geom_errorbar(aes(ymin=Volume-SD, ymax=Volume+SD), width=.2, position=position_dodge(.9)) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 19),
        axis.text.y = element_text(color = "black", 
                                   size = 19), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 19))

#t-test
Pre_post_EPVS_t <- read_excel("G:/MS/Experiments/Karolinska/perivascular spaces in MS/Statistics/PVS/Pre_post_MRI.xlsx", sheet = 2)

t.test(Pre_post_EPVS_t$Volume[Pre_post_EPVS_t$MRI=="Pre"&Pre_post_EPVS_t$Mode=="EPVS"],
       Pre_post_EPVS_t$Volume[Pre_post_EPVS_t$MRI=="Post"&Pre_post_EPVS_t$Mode=="EPVS"],
       alternative = "two.sided", var.equal = TRUE)

t.test(Pre_post_EPVS_t$Volume[Pre_post_EPVS_t$MRI=="Pre"&Pre_post_EPVS_t$Mode=="MS"],
       Pre_post_EPVS_t$Volume[Pre_post_EPVS_t$MRI=="Post"&Pre_post_EPVS_t$Mode=="MS"],
       alternative = "two.sided", var.equal = TRUE)




###Graphs
#Exploratory
Plot_T1_Exploratory <- ggplot(Long4, aes(x=CSO_D_Number, y=WMH_fitted1000)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 32),
        axis.text.y = element_text(color = "black", 
                                   size = 32), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 32)) +
  ylab("T1 lesion volume (ml)") +
  xlab("#of dilated VRS") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_continuous(limits = c(0, 15),
                     breaks = c(2, 4, 6, 8, 10, 12, 14),
                     name="T1 lesion volume (ml)")

Plot_T2_Exploratory <- ggplot(Long2, aes(x=CSO_D_Number, y=Vol_number_fitted)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 32),
        axis.text.y = element_text(color = "black", 
                                   size = 32), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 32)) +
  ylab("T2 lesion volume (ml)") +
  xlab("#of dilated VRS") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_continuous(limits = c(0, 15),
                     breaks = c(2, 4, 6, 8, 10, 12, 14),
                     name="T2 lesion volume (ml)")


#Validation
Plot_T1_Validation <- ggplot(EPVS_validation, aes(x=CSOT2plus, y=T1_diameter_fitted)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 32),
        axis.text.y = element_text(color = "black", 
                                   size = 32), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 32)) +
  ylab("T1 lesion volume (ml)") +
  xlab("#of dilated VRS") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_continuous(limits = c(0, 5),
                     breaks = c(1, 2, 3, 4, 5),
                     name="T1 lesion volume (ml)")

Plot_T2_Validation <- ggplot(EPVS_validation, aes(x=CSOT2plus, y=T2_diameter_fitted)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 32),
        axis.text.y = element_text(color = "black", 
                                   size = 32), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 32)) +
  ylab("T2 lesion volume (ml)") +
  xlab("#of dilated VRS") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_continuous(limits = c(0, 5),
                     breaks = c(1, 2, 3, 4, 5),
                     name="T2 lesion volume (ml)")


ggsave(plot = Plot_T1_Exploratory, width = 7, height = 5.5, filename = "Plot_T1_Exploratory.jpeg")
ggsave(plot = Plot_T2_Exploratory, width = 7, height = 5.5, filename = "Plot_T2_Exploratory.jpeg")
ggsave(plot = Plot_T1_Validation, width = 7, height = 5.5, filename = "Plot_T1_Validation.jpeg")
ggsave(plot = Plot_T2_Validation, width = 7, height = 5.5, filename = "Plot_T2_Validation.jpeg")


####Revisions EBioMedicine
#pool exploratory and validation cohorts
#T1 lesions exploratory cohort
df1 <- as.data.frame(Long4[, c("HIVE", "CSO_D_Number", "Age_Scan", "Sex", "MS_type")])
df1$T1_lesions <- WMH_fitted1000
names(df1) <- c("ID", "CSO_D_Number", "Age_Scan", "Sex", "MS_type", "T1_lesions")

#T2 lesions exploratory cohort
df2 <- as.data.frame(Long2[, c("HIVE", "CSO_D_Number", "Age_Scan", "Sex", "MS_type")])
df2$T2_lesions <- Vol_number_fitted
names(df2) <- c("ID", "CSO_D_Number", "Age_Scan", "Sex", "MS_type", "T2_lesions")

#Validation cohort
EPVS_validation$Age_scan <- EPVS_validation$Scan_year - EPVS_validation$YOB
EPVS_validation$MS_type <- 1
df3 <- as.data.frame(EPVS_validation[, c("ID", "CSOT2plus", "Age_scan", "Sex", "MS_type")])
df3$Sex <- str_replace_all(df3$Sex, "Male", "M")
df3$Sex <- str_replace_all(df3$Sex, "Female", "F")
df3$T1_lesions <- T1_diameter_fitted
df3$T2_lesions <- T2_diameter_fitted
names(df3) <- c("ID", "CSO_D_Number", "Age_Scan", "Sex", "MS_type", "T1_lesions", "T2_lesions")

#pool T1 dataset
df_T1_pooled <- bind_rows(df1, df3)

#pool T2 dataset
df_T2_pooled <- bind_rows(df2, df3)

#pooled plot T1
Plot_T1_pooled <- ggplot(df_T1_pooled, aes(x=CSO_D_Number, y=T1_lesions)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 32),
        axis.text.y = element_text(color = "black", 
                                   size = 32), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 32)) +
  ylab("T1 lesion volume (ml)") +
  xlab("#of dilated VRS") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_continuous(limits = c(0, 15),
                     breaks = c(2, 4, 6, 8, 10, 12, 14),
                     name="T1 lesion volume (ml)")

#pooled plot T2
Plot_T2_pooled <- ggplot(df_T2_pooled, aes(x=CSO_D_Number, y=T2_lesions)) +
  geom_point(size=2, shape=16, color="blue") +
  geom_smooth(method=lm, se=TRUE, fullrange=FALSE, level=0.95,
              linetype="dashed", color="darkred") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(color = "black", 
                                   size = 32),
        axis.text.y = element_text(color = "black", 
                                   size = 32), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 32)) +
  ylab("T2 lesion volume (ml)") +
  xlab("#of dilated VRS") +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                     label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")) +
  scale_y_continuous(limits = c(0, 15),
                     breaks = c(2, 4, 6, 8, 10, 12, 14),
                     name="T2 lesion volume (ml)")

ggsave(plot = Plot_T1_pooled, width = 7, height = 5.5, filename = "Plot_T1_pooled_R1.jpeg")
ggsave(plot = Plot_T2_pooled, width = 7, height = 5.5, filename = "Plot_T2_pooled_R1.jpeg")


#linear models
summary(lm(T1_lesions ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=df_T1_pooled))

summary(lm(T2_lesions ~ CSO_D_Number+Age_Scan+Sex+MS_type, data=df_T2_pooled))

#pooled plot EPVS measures
#boxplot_df_pooled <- boxplot_df
boxplot_df_pooled$Cohort <- str_replace_all(boxplot_df_pooled$Cohort, "Primary", "Pooled")
boxplot_df_pooled$Cohort <- str_replace_all(boxplot_df_pooled$Cohort, "Validation", "Pooled")
boxplot_df_pooled$Class <- str_replace_all(boxplot_df_pooled$Class, "BGA", "BG")
boxplot_df_pooled$Class <- str_replace_all(boxplot_df_pooled$Class, "SOC", "CSO")
boxplot_df_pooled <- subset(boxplot_df_pooled, Class !=  "BGC")

Plot_EPVScount_pooled <- boxplot_df_pooled %>% 
  ggplot(aes(x=Class, y=EPVScount, fill=Cohort)) +
  geom_boxplot(outlier.alpha = F, size = 0.65) + 
  stat_boxplot(geom = "errorbar", colour = "black") +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), alpha=0.2, size = 1.25) +
  labs(x="", y="EPVS count per location") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 22),
        axis.text.y = element_text(color = "black", 
                                   size = 22), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 22)) +
  scale_fill_manual(values=rep(c("#3399FF", "#33CC00"), 3))

#boxplot2_df_pooled <- boxplot2_df
boxplot2_df_pooled$Cohort <- str_replace_all(boxplot2_df_pooled$Cohort, "Exploratory", "Pooled")
boxplot2_df_pooled$Cohort <- str_replace_all(boxplot2_df_pooled$Cohort, "Validation", "Pooled")
boxplot2_df_pooled$Class <- str_replace_all(boxplot2_df_pooled$Class, "BGA", "BG")
boxplot2_df_pooled <- subset(boxplot2_df_pooled, Class !=  "BGC")

Plot_EPVSvolume_pooled <- boxplot2_df_pooled %>% 
  ggplot(aes(x=Class, y=EPVSvolume, fill=Cohort)) +
  geom_boxplot(outlier.alpha = F, size = 0.65) + 
  stat_boxplot(geom = "errorbar", colour = "black") +
  geom_point(position=position_jitterdodge(jitter.width = 0.1), alpha=0.2, size = 1.25) +
  labs(x="", y="EPVS volume per location (ul)") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 22),
        axis.text.y = element_text(color = "black", 
                                   size = 22), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 22)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_fill_manual(values=rep(c("#3399FF", "#33CC00"), 2))

#boxplot_df_pooled2 <- boxplot_df
boxplot_df_pooled2$Cohort <- str_replace_all(boxplot_df_pooled2$Cohort, "Primary", "Pooled")
boxplot_df_pooled2$Cohort <- str_replace_all(boxplot_df_pooled2$Cohort, "Validation", "Pooled")
boxplot_df_pooled2$Class <- str_replace_all(boxplot_df_pooled2$Class, "BGA", "BG")


Plot_EPVSdiameter_pooled <- ggplot(data=boxplot_df_pooled2, aes(x=Cohort, y=Num_diameter, fill=Cohort)) +
  geom_bar(stat="summary", fun.y = "mean", position="dodge", colour="black", size=0.6) +
  geom_point(position=position_jitterdodge(jitter.width = 0.4), alpha=0.2, size = 1.25) +
  scale_fill_brewer(palette="Set1")+
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 22),
        axis.text.y = element_text(color = "black", 
                                   size = 22), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 22)) +
  ylab("Number of dilated EPVS") +
  xlab("") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_fill_manual(values=c("#3399FF", "#33CC00")) +
  theme(axis.text.x = element_blank())

ggsave(plot = Plot_EPVScount_pooled, width = 4.5, height = 5.5, filename = "Plot_EPVScount_pooled.jpeg")
ggsave(plot = Plot_EPVSvolume_pooled, width = 4.5, height = 5.5, filename = "Plot_EPVSvolume_pooled.jpeg")
ggsave(plot = Plot_EPVSdiameter_pooled, width = 4.5, height = 5.5, filename = "Plot_EPVSdiameter_pooled.jpeg")

median(boxplot_df_pooled2$CSOT[boxplot_df_pooled2$Cohort == "Pooled"])
range(boxplot_df_pooled2$CSOT[boxplot_df_pooled2$Cohort == "Pooled"])
median(boxplot_df_pooled2$BGAT[boxplot_df_pooled2$Cohort == "Pooled"])
range(boxplot_df_pooled2$BGAT[boxplot_df_pooled2$Cohort == "Pooled"])
median(boxplot_df_pooled2$CSO_Vol[boxplot_df_pooled2$Cohort == "Pooled"], na.rm = T)
range(boxplot_df_pooled2$CSO_Vol[boxplot_df_pooled2$Cohort == "Pooled"], na.rm = T)


#Proximity graphs
df1 <- as.data.frame(Prox_T1_STOPMS_sub[, c("ID", "prox", "type")])
df2 <- as.data.frame(Prox_box_sub[, c("ID", "prox", "type")])
df_prox_pooled <- bind_rows(df1, df2)
df_prox_pooled$type <- str_replace_all(df_prox_pooled$type, "SOC1", "CSO")
df_prox_pooled$type <- str_replace_all(df_prox_pooled$type, "SOC2", "CSO")
df_prox_pooled$type <- str_replace_all(df_prox_pooled$type, "SOC", "CSO")

Prox_box_pooled <- df_prox_pooled %>% 
  ggplot(aes(x=type, y=prox, fill="blue")) +
  geom_boxplot(outlier.alpha = F, size = 0.65) + 
  stat_boxplot(geom = "errorbar", colour = "black") +
  scale_fill_brewer(palette="Paired") +
  geom_point(position=position_jitter(width = 0.15), alpha=0.2, size = 1.25, color= "black") +
  labs(x="", y="% of VRS volume in vicinity to lesion") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 19),
        axis.text.y = element_text(color = "black", 
                                   size = 19), 
        axis.title = element_text(face = "bold", color = "black", 
                                  size = 19)) +
  theme(legend.position="none")

ggsave(plot = Prox_box_pooled, width = 3.5, height = 5.5, filename = "Prox_box_pooled.jpeg")

#Longitudinal analysis
df_longitudinal <- subset(Long5, CSO_D_Redone_2 !="F")
median(df_longitudinal$CSOT, na.rm = T)
range(df_longitudinal$CSOT, na.rm = T)
median(df_longitudinal$BGCT, na.rm = T)
range(df_longitudinal$BGCT, na.rm = T)
median(df_longitudinal$BS, na.rm = T)
range(df_longitudinal$BS, na.rm = T)
median(df_longitudinal$Vol_PVS_CSO, na.rm = T)
range(df_longitudinal$Vol_PVS_CSO, na.rm = T)
median(df_longitudinal$Vol_PVS_BG, na.rm = T)
range(df_longitudinal$Vol_PVS_BG, na.rm = T)

median(df_longitudinal$EDSS, na.rm = T)
range(df_longitudinal$EDSS, na.rm = T)
median(df_longitudinal$SDMT_Z, na.rm = T)
range(df_longitudinal$SDMT_Z, na.rm = T)
median(df_longitudinal$SDMT_Z, na.rm = T)
range(df_longitudinal$SDMT_Z, na.rm = T)
median(df_longitudinal$T2_Volume_FS, na.rm = T)
range(df_longitudinal$T2_Volume_FS, na.rm = T)
median(df_longitudinal$WMH1, na.rm = T)
range(df_longitudinal$WMH1, na.rm = T)

cor.test(df_longitudinal$Deltavol_CSO, df_longitudinal$T2_volume_D)

median(Long5$Duration)
range(Long5$Duration)
median(Long5$Age_Scan)
range(Long5$Age_Scan)

##Correlation pre-postmortem
VRS_premortem <- c(9.375, 10.38, 2, 21.62, 10.5, 7.625, 5.375, 10.25,  12.55611111)
VRS_postmortem <- c(1.807, 1.908, 0.502, 3.966, 1.104, 0.9036, 3.464, 0.502, 1.729177778)
Correlation_prepost <- data.frame(VRS_premortem, VRS_postmortem)
cor.test(Correlation_prepost$VRS_premortem, Correlation_prepost$VRS_postmortem)
plot(Correlation_prepost$VRS_premortem, Correlation_prepost$VRS_postmortem)

range(Long5$T2_Volume_FS)
range(Long5$WMH1)
range(Long5$BPF1)

range(EPVS_validation$T2_volume_ml)
range(EPVS_validation$T1_lesions)
range(EPVS_validation$BPF)

range(Control$BPF)
