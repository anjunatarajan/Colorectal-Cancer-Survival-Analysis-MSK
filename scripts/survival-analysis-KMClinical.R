library(dplyr)
library(tidyverse)
library(ggplot2)
library(survival)
library(survminer)
library(readxl)

#Kaplan Meier Curves for Clinical data

alldata <- read_excel("X1CRC_FINAL_GENEDATA.xls")

# Primary Tumor Location
KM_Tumor_Location <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_Tumor_Location

KM_TL_FIT <- survfit(KM_Tumor_Location ~ alldata$Primary_Tumor_Location, data = alldata)
summary(KM_TL_FIT)
ggsurvplot(KM_TL_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# Sex
KM_Sex <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_Sex

KM_Sex_FIT <- survfit(KM_Sex~ alldata$Sex, data = alldata)
summary(KM_Sex_FIT)
ggsurvplot(KM_Sex_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# Age at Diagnosis above or below Median Age of Cohort
alldata$Above_or_Below_Median <- ifelse(alldata$`Age at Diagnosis` <= 54,0,1)
ab_median <- alldata$Above_or_Below_Median
View(ab_median)
km_ab_median <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
km_ab_median

km_ab_median_fit <- survfit(km_ab_median~ alldata$Above_or_Below_Median, data = alldata)
summary(km_ab_median_fit)
ggsurvplot(km_ab_median_fit, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE, legend.labs = c("Below Median Age", "Above Median Age"))

# Metastasectomy
KM_Metsy <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_Metsy

KM_Metsy_FIT <- survfit(KM_Metsy~ alldata$Metastasectomy, data = alldata)
summary(KM_Metsy_FIT)
ggsurvplot(KM_Metsy_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)

# Stage at Diagnosis

sad_res <- pairwise_survdiff(Surv(OS_Months, OS_Status ~ Stage_At_Diagnosis), data = sad_data)
sad_res

KM_Stage <- Surv(time = alldata$`Overall Survival (Months)`, event = alldata$`Overall Survival Status (Num)`)
KM_Stage

KM_Stage_FIT <- survfit(KM_Stage~ alldata$Stage_At_Diagnosis, data = alldata)
summary(KM_Stage_FIT)
ggsurvplot(KM_Stage_FIT, data = alldata, pval = TRUE, risk.table = TRUE, conf.int = TRUE)
