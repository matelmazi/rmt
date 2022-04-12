#' ---
#' title: "Personalized Dynamic predictions for Longitudinal Outcomes"
#' subtitle: "An application to stroke patients"
#' author: "Matea Elmazi, Research Master at NIHES"
#' date: "`r Sys.setenv(LANG = 'en_US.UTF-8'); format(Sys.Date(), '%d %B %Y')`"
#' output: 
#'   html_document:
#'     toc: true
#'     toc_float:
#'       collapsed: false
#'     df_print: paged
#' ---
#' 

#' ## Load packages 
library(lattice)
library(ggplot2)
library(gtrendsR)
library(nlme)
library(splines)
library(memisc)
library(latticeExtra)
library(moments)
library(shiny)
library(gridExtra)
library(grid)
library(JMbayes)


#' Load the data.
load('BIactNew2.RData')

#' Save dataset in data frame type.
BIactNew2 <- data.frame(BIactNew2)

#' Summary of the missing values
table(is.na(BIactNew2))

#' Transform numeric variables to factor variables
BIactNew2 <- BIactNew2[BIactNew2$WRITINGHAND %in% c(1,2),]
BIactNew2$WRITINGHAND <- factor(x=BIactNew2$WRITINGHAND, levels = c(1,2), labels = c("left", "right")) 
BIactNew2$Geslacht <- factor(x=BIactNew2$Geslacht, levels = c(1,2), labels = c("male", "female"))
BIactNew2$PARTNER <- factor(x=BIactNew2$PARTNER, levels = c(0,1), labels = c("no", "yes"))
BIactNew2$BAMFORD <- factor(x=BIactNew2$BAMFORD, levels = c(1,2,3), labels = c("1", "2", "3"))

#' Order the dataset
BIactNew2 <- BIactNew2[order(BIactNew2$Patient_nummer, BIactNew2$Days, BIactNew2$FMBETOT),]

#' ##  Number of the repeated measurements
#' Number of repeated measurements per patient:
length.no.NA <- function (x) sum(!is.na(x))
ns <- with(BIactNew2, tapply(BITOT, Patient_nummer, length.no.NA))

#' We produce numerical summaries and a visualization using a boxplot:
c(summary(ns), sd = sd(ns))
boxplot(ns, col = "lightgrey", ylab = "Number Repeated Measurements per Patient")

#' ## Descriptive Statistics
#' Create baseline dataset
base_BIactNew2 <- BIactNew2[!duplicated(BIactNew2[, "Patient_nummer"]), ]
#'Summary of the baseline characteristics
summary(base_BIactNew2)

#' ## Baseline measurements
#' Here we repeat the baseline values as many times as the repeated measurements
BIactNew2$BI_base <- rep(base_BIactNew2$BITOT, table(BIactNew2$Patient_nummer))
BIactNew2$FM_base <- rep(base_BIactNew2$FMBETOT, table(BIactNew2$Patient_nummer))
BIactNew2$MI_UE_base <- rep(base_BIactNew2$MITOTBE, table(BIactNew2$Patient_nummer))
BIactNew2$MI_LE_base <- rep(base_BIactNew2$MITOTOE, table(BIactNew2$Patient_nummer))

#' ## Last measurements
#' Create a data set containing only the last measurement per patient
tail_BIactNew2 <- BIactNew2[tapply(row.names(BIactNew2), BIactNew2$Patient_nummer, tail, 1), ]

#' Here we repeat the last value as many times as the repeated measurements
BIactNew2$BI_tail <- rep(tail_BIactNew2$BITOT, table(BIactNew2$Patient_nummer))
BIactNew2$FM_tail <- rep(tail_BIactNew2$FMBETOT, table(BIactNew2$Patient_nummer))
BIactNew2$MI_UE_tail <- rep(tail_BIactNew2$MITOTBE, table(BIactNew2$Patient_nummer))
BIactNew2$MI_LE_tail <- rep(tail_BIactNew2$MITOTOE, table(BIactNew2$Patient_nummer))

#' Exclude patients with 1 measurement
table(table(BIactNew2$Patient_nummer))

dat1mes <- lapply(split(BIactNew2, BIactNew2$Patient_nummer), function(x) as.data.frame(x))
iDs <- lapply(dat1mes, function(x) if (dim(x)[1] == 1) print(x$Patient_nummer))
vecIDS <- do.call(c, iDs)

#' Create a dataset that contains patients with more than one measurement
data_morethanonerow <-  BIactNew2[!BIactNew2$Patient_nummer %in% vecIDS, ]


#' ## Model building
#' Definition of different models structures and evaluation with AIC to determine first the most suitable number of splines and after the most suitable model structure
#' No interaction assumed - different splines (3,4,5)
m1.3 <-  lme(BITOT ~ ns(Days, 3) + Geslacht + StrokeLeeftijd+ MI_UE_base + MI_LE_base + FM_base, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 3))),
             data = BIactNew2, method = "ML")

m1.4 <-  lme(BITOT ~ ns(Days, 4) + Geslacht + StrokeLeeftijd+ MI_UE_base + MI_LE_base + FM_base, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 4))),
             data = BIactNew2, method = "ML")

m1.5 <-  lme(BITOT ~ ns(Days, 5) + Geslacht + StrokeLeeftijd+ MI_UE_base + MI_LE_base + FM_base, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
             data = BIactNew2, method = "ML")
AIC(m1.3, m1.4, m1.5)

#' Interaction between time, age and gender - different splines (3,4,5)
m2.3 <-  lme(BITOT ~ ns(Days, 3) *(Geslacht + StrokeLeeftijd) + MI_UE_base + MI_LE_base + FM_base, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 3))),
             data = BIactNew2, method = "ML")

m2.4 <-  lme(BITOT ~ ns(Days, 4) *(Geslacht + StrokeLeeftijd) + MI_UE_base + MI_LE_base + FM_base, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 4))),
             data = BIactNew2, method = "ML")

m2.5 <-  lme(BITOT ~ ns(Days, 5) *(Geslacht + StrokeLeeftijd) + MI_UE_base + MI_LE_base + FM_base, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
             data = BIactNew2, method = "ML")
AIC(m2.3, m2.4, m2.5)

#' Interaction between time, motricity index(lower & upper extremity) and fugl meyer
m3.3 <-  lme(BITOT ~ ns(Days, 3) *(MI_UE_base + MI_LE_base + FM_base)+ Geslacht + StrokeLeeftijd ,
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 3))),
             data = BIactNew2, method = "ML")

m3.4 <-  lme(BITOT ~ ns(Days, 4) *(MI_UE_base + MI_LE_base + FM_base)+ Geslacht + StrokeLeeftijd ,
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 4))),
             data = BIactNew2, method = "ML")

m3.5 <-  lme(BITOT ~ ns(Days, 5) *(MI_UE_base + MI_LE_base + FM_base)+ Geslacht + StrokeLeeftijd ,
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
             data = BIactNew2, method = "ML")
AIC(m3.3, m3.4, m3.5)

#' Interaction between time, motricity index(lower & upper extremity) and fugl meyer &  between time Gender and age
m4.3 <-  lme(BITOT ~ ns(Days, 3)*(Geslacht + StrokeLeeftijd) + ns(Days, 3) *(MI_UE_base + MI_LE_base + FM_base),
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 3))),
             data = BIactNew2, method = "ML")

m4.4 <-  lme(BITOT ~ ns(Days, 4)*(Geslacht + StrokeLeeftijd) + ns(Days, 4) *(MI_UE_base + MI_LE_base + FM_base),
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 4))),
             data = BIactNew2, method = "ML")

m4.5 <-  lme(BITOT ~ ns(Days, 5)*(Geslacht + StrokeLeeftijd) + ns(Days, 5) *(MI_UE_base + MI_LE_base + FM_base),
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
             data = BIactNew2, method = "ML")
AIC(m4.3, m4.4, m4.5)

#' Interaction between time, motricity index(lower & upper extremity) and fugl meyer &  between time Gender and age
#' Main effects of NIHSSTOT and Bamford
m5.3 <-  lme(BITOT ~ ns(Days, 3)*(Geslacht + StrokeLeeftijd) + 
                     ns(Days, 3)*(MI_UE_base + MI_LE_base + FM_base) +
                     NIHSSTOT + BAMFORD, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 3))),
             data = BIactNew2, method = "ML")

m5.4 <-  lme(BITOT ~ ns(Days, 4)*(Geslacht + StrokeLeeftijd) + 
                     ns(Days, 4)*(MI_UE_base + MI_LE_base + FM_base) +
                     NIHSSTOT + BAMFORD, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 4))),
             data = BIactNew2, method = "ML")

m5.5 <-  lme(BITOT ~ ns(Days, 5)*(Geslacht + StrokeLeeftijd) + 
                     ns(Days, 5)*(MI_UE_base + MI_LE_base + FM_base) +
                     NIHSSTOT + BAMFORD, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
             data = BIactNew2, method = "ML")
AIC(m5.3, m5.4, m5.5)

#' Interaction between time, motricity index(lower & upper extremity) and fugl meyer &  between time Gender and age
#' Interaction between  time NIHSSTOT and Bamford
m6.3 <- lme(BITOT ~ ns(Days, 3)*(Geslacht + StrokeLeeftijd) + 
                     ns(Days, 3)*(MI_UE_base + MI_LE_base + FM_base) +
                     ns(Days, 3)*(NIHSSTOT + BAMFORD), 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 3))),
             data = BIactNew2, method = "ML")

m6.4 <- lme(BITOT ~ ns(Days, 4)*(Geslacht + StrokeLeeftijd) + 
                     ns(Days, 4)*(MI_UE_base + MI_LE_base + FM_base) +
                     ns(Days, 4)*(NIHSSTOT + BAMFORD), 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 4))),
             data = BIactNew2, method = "ML")

m6.5 <- lme(BITOT ~ ns(Days, 5)*(Geslacht + StrokeLeeftijd) + 
                     ns(Days, 5)*(MI_UE_base + MI_LE_base + FM_base) +
                     ns(Days, 5)*(NIHSSTOT + BAMFORD), 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
             data = BIactNew2, method = "ML")

m6.5_no_FM <- lme(BITOT ~ ns(Days, 5)*(Geslacht + StrokeLeeftijd) + 
                          ns(Days, 5)*(MI_UE_base + MI_LE_base) +
                          ns(Days, 5)*(NIHSSTOT + BAMFORD), 
                  random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
                  data = BIactNew2, method = "ML")

m6.5_no_MI <-  lme(BITOT ~ ns(Days, 5)*(Geslacht + StrokeLeeftijd) + 
                     ns(Days, 5)*(FM_base) +
                     ns(Days, 5)*(NIHSSTOT + BAMFORD), 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
             data = BIactNew2, method = "ML")

m6.5_no_NIHSS <-  lme(BITOT ~ ns(Days, 5)*(Geslacht + StrokeLeeftijd) + 
                     ns(Days, 5)*(MI_UE_base + MI_LE_base + FM_base) +
                     ns(Days, 5)*(BAMFORD), 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
             data = BIactNew2, method = "ML")
AIC(m6.5, m6.5_no_FM, m6.5_no_MI, m6.5_no_NIHSS)


#' Interaction between time, motricity index(lower & upper extremity) and fugl meyer &  between time and Gender
#' Interaction between  time NIHSSTOT and Bamford
##### NO AGE included
m7.3 <-  lme(BITOT ~ ns(Days, 3)*(Geslacht) + ns(Days, 3)*(MI_UE_base + MI_LE_base + FM_base) +
                     ns(Days, 3)*(NIHSSTOT + BAMFORD), 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 3))),
             data = BIactNew2, method = "ML")

m7.4 <-  lme(BITOT ~ ns(Days, 4)*(Geslacht) + ns(Days, 4)*(MI_UE_base + MI_LE_base + FM_base) +
                     ns(Days, 4)*(NIHSSTOT + BAMFORD), 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 4))),
             data = BIactNew2, method = "ML")

m7.5 <-  lme(BITOT ~ ns(Days, 5)*(Geslacht) + ns(Days, 5)*(MI_UE_base + MI_LE_base + FM_base) +
                     ns(Days, 5)*(NIHSSTOT + BAMFORD), 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
             data = BIactNew2, method = "ML")
AIC(m7.3, m7.4, m7.5)


#' Main effects of all covariates
#' ##### NO AGE included
m8.3 <-  lme(BITOT ~ ns(Days, 3) + Geslacht + MI_UE_base + MI_LE_base + FM_base + NIHSSTOT + BAMFORD, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 3))),
             data = BIactNew2, method = "ML")

m8.4 <-  lme(BITOT ~ ns(Days, 4) + Geslacht + MI_UE_base + MI_LE_base + FM_base + NIHSSTOT + BAMFORD, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 4))),
             data = BIactNew2, method = "ML")

m8.5 <-  lme(BITOT ~ ns(Days, 5) + Geslacht + MI_UE_base + MI_LE_base + FM_base + NIHSSTOT + BAMFORD, 
             random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
             data = BIactNew2, method = "ML")

AIC(m8.3, m8.4, m8.5)

AIC(m1.5,m2.5,m3.5,m4.5, m5.5, m6.5, m7.5, m8.5)



