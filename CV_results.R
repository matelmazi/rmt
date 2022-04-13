#' Load packages
library(nlme)
library(splines)
library(ggstatsplot)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridisLite)
library(viridis)
library(grid)
library(ggpubr)

###################################################################
#' Load the dataset
load('BIactNew2.RData')

#' Transform variables from numeric to factor 
BIactNew2 <- BIactNew2[BIactNew2$WRITINGHAND %in% c(1,2),]
BIactNew2$WRITINGHAND <- factor(x=BIactNew2$WRITINGHAND, levels = c(1,2), labels = c("left", "right"))
BIactNew2$Geslacht <- factor(x=BIactNew2$Geslacht, levels = c(1,2), labels = c("male", "female"))
BIactNew2$BAMFORD <- factor(x=BIactNew2$BAMFORD, levels = c(1,2,3), labels = c("1", "2", "3"))
BIactNew2$PARTNER <- factor(x=BIactNew2$PARTNER, levels = c(0,1), labels = c("no", "yes"))

#' Order the dataset
BIactNew2 <- BIactNew2[order(BIactNew2$Patient_nummer, BIactNew2$Days, BIactNew2$FMBETOT),]
#' Exclude patients with 1 measurement
table(table(BIactNew2$Patient_nummer))
dat1mes <- lapply(split(BIactNew2, BIactNew2$Patient_nummer), function(x) as.data.frame(x))
iDs <- lapply(dat1mes, function(x) if (dim(x)[1] == 1) print(x$Patient_nummer))
vecIDS <- do.call(c, iDs)
#' Create a dataset that contains patients with more than one measurement
data_morethanonerow <-  BIactNew2[!BIactNew2$Patient_nummer %in% vecIDS, ]

#' ## Baseline measurements
base_BIactNew2 <- BIactNew2[!duplicated(BIactNew2[, "Patient_nummer"]), ]
#' Here we repeat the baseline value as many times as the repeated measurements
#' We start with the Barthel Index total score
BIactNew2$BI_base <- rep(base_BIactNew2$BITOT, table(BIactNew2$Patient_nummer))
BIactNew2$FM_base <- rep(base_BIactNew2$FMBETOT, table(BIactNew2$Patient_nummer))
BIactNew2$MI_UE_base <- rep(base_BIactNew2$MITOTBE, table(BIactNew2$Patient_nummer))
BIactNew2$MI_LE_base <- rep(base_BIactNew2$MITOTOE, table(BIactNew2$Patient_nummer))


###################################################################
#' Load Cross validation results
load("results_m5.5.RData")
load("results_m6.5.RData")
load("results_m7.5.RData")
load("res_m6_without_FM.RData")
load("res_m6_without_MI.RData")


#' Check if the columns id and time are the same in the files to ensure all patients were validated
cbind(final_5.5$Patient_nummer,final_6.5$Patient_nummer, final_7.5$Patient_nummer,
      final_m6_without_FM$Patient_nummer, final_m6_without_MI$Patient_nummer)

##################################################################
#'Remove lines with NA - keep only the predictions and then transform the number of visits from numeric to factor 
final_5.5 <- final_5.5[which(!is.na(final_5.5$low)),]
final_5.5$no_visits <- factor(x=final_5.5$no_visits , levels = c(1,2,3,4,5,6,7), labels = c("1", "2", "3","4", "5", "6", "7"))

final_6.5 <- final_6.5[which(!is.na(final_6.5$low)),]
final_6.5$no_visits <- factor(x=final_6.5$no_visits , levels = c(1,2,3,4,5,6,7), labels = c("1", "2", "3","4", "5", "6", "7"))

final_7.5 <- final_7.5[which(!is.na(final_7.5$low)),]

final_m6_without_MI <- final_m6_without_MI[which(!is.na(final_m6_without_MI$low)),]
final_m6_without_MI$no_visits <- factor(x=final_m6_without_MI$no_visits , levels = c(1,2,3,4,5,6,7), labels = c("1", "2", "3","4", "5", "6", "7"))

final_m6_without_FM <- final_m6_without_FM[which(!is.na(final_m6_without_FM$low)),]
final_m6_without_FM$no_visits <- factor(x=final_m6_without_FM$no_visits , levels = c(1,2,3,4,5,6,7), labels = c("1", "2", "3","4", "5", "6", "7"))

##################################################################
#' Create a data frame that contains the errors of all models and the days 
All_errors <- data.frame(Gender=final_5.5$Geslacht, Error_m5=final_5.5$Error, Error_m7=final_7.5$Error,
                        Error_m6=final_6.5$Error, Error_m6_noFM=final_m6_without_FM$Error,
                        Error_m6_noMI=final_m6_without_MI$Error)
##################################################################
#' Number of serial measurements used and absolute error per model 
par(mfrow = c(1,3))

boxplot(final_6.5$Error~final_6.5$no_visits, horizontal=FALSE,
        cex.axis = 0.9, ylab = "Absolute Error",
        cex.lab = 1.1, xlab = "Number of serial measurements used", main = "Mixed-Effects Model m6.5", 
        cex.main = 1.3, ylim = c(0, 20),
        outline = FALSE)

boxplot(final_m6_without_FM $Error~final_m6_without_FM $no_visits, horizontal=FALSE,
        cex.axis = 0.9, ylab = "Absolute Error",
        cex.lab = 1.1, xlab = "Number of serial measurements used", main = "Mixed-Effects Model m6.5 without FM", 
        cex.main = 1.3, ylim = c(0, 20),
        outline = FALSE)

boxplot(final_m6_without_MI$Error~final_m6_without_MI$no_visits, horizontal=FALSE,
        cex.axis = 0.9, ylab = "Absolute Error",
        cex.lab = 1.1, xlab = "Number of serial measurements used", main = "Mixed-Effects Model m6.5 without MI", 
        cex.main = 1.3, ylim = c(0, 20),
        outline = FALSE)

boxplot(final_5.5$Error~final_5.5$no_visits, horizontal=FALSE,
        cex.axis = 0.9, ylab = "Absolute Error",
        cex.lab = 1.1, xlab = "Number of serial measurements used", main = "Mixed-Effects Model m5.5", 
        cex.main = 1.3, ylim = c(0, 20),
        outline = FALSE)

boxplot(final_7.5$Error~final_7.5$no_visits, horizontal=FALSE,
        cex.axis = 0.9, ylab = "Absolute Error",
        cex.lab = 1.1, xlab = "Number of serial measurements used", main = "Mixed-Effects Model m7.5", 
        cex.main = 1.3, ylim = c(0, 20),
        outline = FALSE)

##################################################################
#' Create violin error plots to compare the predicting performance of the models
p1 <- ggplot(final_6.5, aes(x=no_visits, y=Error, fill=no_visits)) + 
  geom_violin(trim=FALSE) +
  stat_summary(fun=median) +
  scale_fill_brewer(palette="RdBu") +
  labs(y = "Absolute Error", x = "Number of serial measurements used")


p2 <- ggplot(final_m6_without_FM, aes(x=no_visits, y=Error, fill=no_visits)) + 
  geom_violin(trim=FALSE) +
  scale_fill_brewer(palette="RdBu") +
  stat_summary(fun=median) +
  labs(y = "Absolute Error", x = "Number of serial measurements used")


p3 <- ggplot(final_m6_without_MI, aes(x=no_visits, y=Error, fill=no_visits)) + 
  geom_violin(trim=FALSE) +
  stat_summary(fun=median) +
  scale_fill_brewer(palette="RdBu") +
  labs(y = "Absolute Error", x = "Number of serial measurements used")


p4 <- ggplot(final_5.5, aes(x=no_visits, y=Error, fill=no_visits)) + 
  geom_violin(trim=FALSE)+
  stat_summary(fun=median) +
  scale_fill_brewer(palette="RdBu")+
  labs(y = "Absolute Error", x = "Number of serial measurements used")

#' Combine violin plots into a single figure
ggarrange(
  p1, p2, p3, p4, labels = c("Model m6.5", "Model m6.5 without FM", "Model m6.5 without MI", "Model m7.5"),
  common.legend = FALSE, legend = "none", ncol = 2, nrow = 2
)

ggarrange(
  p1, p2, labels = c("Model m6.5", "Model m6.5 without FM"),
  common.legend = FALSE, legend = "none", ncol = 2, nrow = 1
)