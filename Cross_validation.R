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
library(nlme)
library(splines)
library(memisc)

#' Load the data.
load('BIactNew2.RData')

#' Save dataset in data frame type.
BIactNew2 <- data.frame(BIactNew2)

#' Transform variables from numeric to factor 
BIactNew2 <- BIactNew2[BIactNew2$WRITINGHAND %in% c(1,2),]
BIactNew2$WRITINGHAND <- factor(x=BIactNew2$WRITINGHAND, levels = c(1,2), labels = c("left", "right"))
BIactNew2$Geslacht <- factor(x=BIactNew2$Geslacht, levels = c(1,2), labels = c("male", "female"))
BIactNew2$PARTNER <- factor(x=BIactNew2$PARTNER, levels = c(0,1), labels = c("no", "yes"))
BIactNew2$BAMFORD <- factor(x=BIactNew2$BAMFORD, levels = c(1,2,3), labels = c("1", "2", "3"))

#' Order the dataset based on the patient number, the days since stroke and the FM score
BIactNew2 <- BIactNew2[order(BIactNew2$Patient_nummer, BIactNew2$Days, BIactNew2$FMBETOT),]

#' Create baseline data set
base_BIactNew2 <- BIactNew2[!duplicated(BIactNew2[, "Patient_nummer"]), ]

#' Create baseline covariates
#' Here we repeat the baseline value as many times as the repeated measurements
BIactNew2$BI_base <- rep(base_BIactNew2$BITOT, table(BIactNew2$Patient_nummer))
BIactNew2$FM_base <- rep(base_BIactNew2$FMBETOT, table(BIactNew2$Patient_nummer))
BIactNew2$MI_UE_base <- rep(base_BIactNew2$MITOTBE, table(BIactNew2$Patient_nummer))
BIactNew2$MI_LE_base <- rep(base_BIactNew2$MITOTOE, table(BIactNew2$Patient_nummer))

#'load functions for prediction
source("Functions.R")

#' Exclude patients with 1 measurement
table(table(BIactNew2$Patient_nummer))
dat1mes <- lapply(split(BIactNew2, BIactNew2$Patient_nummer), function(x) as.data.frame(x))
iDs <- lapply(dat1mes, function(x) if (dim(x)[1] == 1) print(x$Patient_nummer))
vecIDS <- do.call(c, iDs)
#' Create a dataset that contains patients with more than one measurement
data_morethanonerow <-  BIactNew2[!BIactNew2$Patient_nummer %in% vecIDS, ]


#'## Cross - Validation for Model: m6.5 without Fugl Meyer
Ind_pred_noFM <- list()
#' create a vector that includes unique patients
vector_id <- unique(data_morethanonerow$Patient_nummer)

for (i in vector_id) {
  trainingData <- data_morethanonerow[!data_morethanonerow$Patient_nummer %in% i, ]
  testingData <- data_morethanonerow[data_morethanonerow$Patient_nummer %in% i, ]
  
  #' Fit the model that you want to validate under the training data set
  m6_without_FM <-  lme(BITOT ~ ns(Days, 5)*(Geslacht+StrokeLeeftijd) + ns(Days, 5)*(MI_UE_base+MI_LE_base) + ns(Days, 5)*(NIHSSTOT+BAMFORD), 
                        random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
                        data = trainingData)
  
  #' Create a list to save the predictions
  IndvPrediction_noFM<- list()
  #' Loop to make predictions for every extra measurement per patient
  for (j in 1:(dim(testingData)[1]-1)) {
    
    testingData_without <- testingData[1:j, ]
    
    # Select last prediction time from the testingData to predict on the last time point for each patient
    pred_time <- tail(testingData$Days, n=1)
    # Obtain predictions
    IndvPrediction_noFM[[j]] <- IndvPred_lme_2(lmeObject = m6_without_FM, newdata = testingData_without, "Days", times = pred_time, M = 500, 
                                              interval = "prediction", return_data = TRUE)
    # Obtain number of visit for every prediction
    IndvPrediction_noFM[[j]]$no_visits<- dim(testingData_without)[1]
  }
  Ind_pred_noFM[[i]]<- do.call(rbind,IndvPrediction_noFM)
  print(i) #check patient     
}
final_m6_without_FM <- do.call(rbind,Ind_pred_noFM)
# calculate the error
final_m6_without_FM$Error <- abs(final_m6_without_FM$BITOT - final_m6_without_FM$pred)
# save the dataset with the predictions
save(final_m6_without_FM, file = "res_m6_without_FM.RData")

#' The same procedure can be performed for each model separately, for simplicity there is only code for only one model