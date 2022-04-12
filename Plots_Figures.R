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

#' Transform numeric variables to factor variables
BIactNew2 <- BIactNew2[BIactNew2$WRITINGHAND %in% c(1,2),]
BIactNew2$WRITINGHAND <- factor(x=BIactNew2$WRITINGHAND, levels = c(1,2), labels = c("left", "right")) 
BIactNew2$Geslacht <- factor(x=BIactNew2$Geslacht, levels = c(1,2), labels = c("male", "female"))
BIactNew2$PARTNER <- factor(x=BIactNew2$PARTNER, levels = c(0,1), labels = c("no", "yes"))
BIactNew2$BAMFORD <- factor(x=BIactNew2$BAMFORD, levels = c(1,2,3), labels = c("1", "2", "3"))

#' Order the dataset
BIactNew2 <- BIactNew2[order(BIactNew2$Patient_nummer, BIactNew2$Days, BIactNew2$FMBETOT),]

#' Create baseline dataset
base_BIactNew2 <- BIactNew2[!duplicated(BIactNew2[, "Patient_nummer"]), ]

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

#' ## Appendix plots 
#' Plots of baseline characteristics (BI, Age, Gender, Writing hand)

plot(x = base_BIactNew2$StrokeLeeftijd, y = base_BIactNew2$BITOT, ylab = "Barthel Index", xlab = "Age", 
     cex.axis = 1, cex.lab = 1, col = base_BIactNew2$Geslacht)
legend(30, 25, legend = c("male", "female"), col = c(1,2)) 

xyplot(x = BITOT ~ StrokeLeeftijd, group = Geslacht, data = base_BIactNew2, type = "smooth", 
       lwd = 2, col = c("red", "blue"))

xyplot(x = BITOT ~ StrokeLeeftijd | Geslacht, data = base_BIactNew2, type = c("p", "smooth"), 
       lwd = 2, col = c("red"))   

xyplot(x = BITOT ~ StrokeLeeftijd | WRITINGHAND, data = base_BIactNew2, type = c("p", "smooth"), 
       lwd = 2, col = c("red")) 

#' ## Individual patient plots
#' Subject-specific trajectories for all patients
xyplot(x = BITOT ~ Days, group = Patient_nummer, data = BIactNew2, type = "l", col = "black")

#' Subject-specific trajectories per gender
xyplot(BITOT ~ Days | Geslacht, group = Patient_nummer, data = BIactNew2, 
       panel = function (x, y, ...) {
         panel.xyplot(x, y, type = "l", col = 1, ...)
         panel.loess(x, y, col = 2, lwd = 2)
       }, xlab = "Days since stroke", ylab = "Barthel Index")

#' Subject-specific trajectories per writing hand
xyplot(BITOT ~ Days | WRITINGHAND, group = Patient_nummer, data = BIactNew2, 
       panel = function (x, y, ...) {
         panel.xyplot(x, y, type = "l", col = 1, ...)
         panel.loess(x, y, col = 2, lwd = 2)
       }, xlab = "Days since stroke", ylab = "Barthel Index")

#' ## Scales of baseline scores
#' Creation of a new categorical variable with the different scales of baseline Barthel Index
BIactNew2$BIcategory[BIactNew2$BI_base < 4] <- "BI: 0-3"
BIactNew2$BIcategory[BIactNew2$BI_base < 13 & BIactNew2$BI_base >=4] <- "BI: 4-12"
BIactNew2$BIcategory[BIactNew2$BI_base >=13] <- "BI: 13-20"


#' Creation of a new categorical variable with the different scales of Fugl Mayer
BIactNew2$FMcat[BIactNew2$FM_base < 4] <- "FM: 0-3"
BIactNew2$FMcat[BIactNew2$FM_base < 15 & BIactNew2$FM_base >=4] <- "FM: 4-15"
BIactNew2$FMcat[BIactNew2$FM_base >=15] <- "FM: 16-66"

#' Creation of a new categorical variable with the different scales of Motricity index - Lower extremity
BIactNew2$MI_LEcat[BIactNew2$MI_LE_base < 23] <- "MI_LE: 0-22"
BIactNew2$MI_LEcat[BIactNew2$MI_LE_base < 75 & BIactNew2$MI_LE_base >=23] <- "MI_LE: 23-74"
BIactNew2$MI_LEcat[BIactNew2$MI_LE_base >=75] <- "MI_LE: 75-100"

#' Creation of a new categorical variable with the different scales of Motricity index - Upper extremity
BIactNew2$MI_UEcat[BIactNew2$MI_UE_base < 31 ] <- "MI_UE: 0-30"
BIactNew2$MI_UEcat[BIactNew2$MI_UE_base < 70  & BIactNew2$MI_UE_base >= 31] <- "MI_UE: 31-69"
BIactNew2$MI_UEcat[BIactNew2$MI_UE_base >=70] <- "MI_UE: 70-100"

#' ## Plots exploration of splines with different df
# To do this step model m2.3 is used but it should be loaded first
# Any other model under study of course can also be used 
set.seed(123)
ind <- sample(unique(data_morethanonerow$Patient_nummer), 49, replace = F)

data_morethanonerow$fitted_marg <- fitted(m2.3, level = 0) 
data_morethanonerow$fitted_subj <- fitted(m2.3, level = 1) # subject-specific
p1 <- xyplot(BITOT + fitted_marg + fitted_subj ~ Days | Patient_nummer, data = data_morethanonerow,
             panel = function (x, y, ...) {
                     x.mat <- matrix(x, ncol = 3)
                     y.mat <- matrix(y, ncol = 3)
                     panel.xyplot(x.mat[, 1], y.mat[, 1], type = "p", col = "black")
                     panel.xyplot(x.mat[, 3], y.mat[, 3], type = "l", lwd = 2, col = "blue")
             }, subset = Patient_nummer %in% ind, layout = c(7, 7), as.table = TRUE, 
             xlab = "Time (Days)", ylab = "Barthel index")

data_morethanonerow$fitted_marg <- fitted(m2.4, level = 0) 
data_morethanonerow$fitted_subj <- fitted(m2.4, level = 1) # subject-specific
p2 <- xyplot(BITOT + fitted_marg + fitted_subj ~ Days | Patient_nummer, data = data_morethanonerow,
             panel = function (x, y, ...) {
                     x.mat <- matrix(x, ncol = 3)
                     y.mat <- matrix(y, ncol = 3)
                     panel.xyplot(x.mat[, 1], y.mat[, 1], type = "p", col = "black")
                     panel.xyplot(x.mat[, 3], y.mat[, 3], type = "l", lwd = 2, col = "red")
             }, subset = Patient_nummer %in% ind, layout = c(7, 7), as.table = TRUE, 
             xlab = "Time (Days)", ylab = "Barthel index")

data_morethanonerow$fitted_marg <- fitted(m2.5, level = 0) 
data_morethanonerow$fitted_subj <- fitted(m2.5, level = 1) # subject-specific
p3 <- xyplot(BITOT + fitted_marg + fitted_subj ~ Days | Patient_nummer, data = data_morethanonerow,
             panel = function (x, y, ...) {
                     x.mat <- matrix(x, ncol = 3)
                     y.mat <- matrix(y, ncol = 3)
                     panel.xyplot(x.mat[, 1], y.mat[, 1], type = "p", col = "black")
                     panel.xyplot(x.mat[, 3], y.mat[, 3], type = "l", lwd = 2, col = "green")
             }, subset = Patient_nummer %in% ind, layout = c(7, 7), as.table = TRUE, 
             xlab = "Time (Days)", ylab = "Barthel index")

p1 + as.layer(p2) + as.layer(p3)


############
# Figure 1 #
############
#500 x 500
par(mfrow = c(1, 3))
###  BI: 0-3  ###
BI03 <- BIactNew2[BIactNew2$BIcategory == "BI: 0-3",]

BI_1 <- xyplot(BITOT ~ Days, group = Patient_nummer, data = BI03, type = c("l"), cex = 0.6,col = "grey", 
               strip = FALSE, pch = 16, grid = FALSE, lwd = 2, ylim = c(0,21), xlim = c(0,230), cex.main =0.5,
               xlab = list("Days since stroke"), ylab = list("Barthel Index"), main = "Baseline BI: 0-3", 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata1 <- BI03[BI03$Patient_nummer == 7,]
BI_2 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "blue4", lwd = 4, 
               data = subdata1, ylim = c(0,21), xlim = c(0,230), 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata2 <- BI03[BI03$Patient_nummer == 150,]
BI_3 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "blue", lwd = 4, 
               data = subdata2, ylim = c(0,21), xlim = c(0,230), 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata3 <- BI03[BI03$Patient_nummer == 215,]
BI_4 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "darkred", lwd = 4,
               data = subdata3, ylim = c(0,21), xlim = c(0,230), 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata4 <- BI03[BI03$Patient_nummer == 350,]
BI_5 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "darkcyan", lwd = 4, 
               data = subdata4, ylim = c(0,21), xlim = c(0,230), 
               par.settings = list(fontsize = list(text = 13, points = 18)))


plot1<- BI_1 + as.layer(BI_2) + as.layer(BI_3) + as.layer(BI_4) + as.layer(BI_5) 

###  BI: 4-12  ###
BI0412 <- BIactNew2[BIactNew2$BIcategory == "BI: 4-12",]

BI_6 <- xyplot(BITOT ~ Days, group = Patient_nummer, data = BI0412, type = c("l"), cex = 0.6,col = "grey", 
               strip = FALSE, pch = 16, grid = FALSE, lwd = 2, ylim = c(0,21), xlim = c(0,230), cex.main =0.5,
               xlab = list("Days since stroke"), ylab = list("Barthel Index"), main = "Baseline BI: 4-12", 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata1 <- BI0412[BI0412$Patient_nummer ==124,]
BI_7 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "pink", lwd = 4, 
               data = subdata1, ylim = c(0,21), xlim = c(0,230), 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata2 <- BI0412[BI0412$Patient_nummer == 145,]
BI_8 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "darkcyan", lwd = 4, 
               data = subdata2, ylim = c(0,21), xlim = c(0,230), 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata3 <- BI0412[BI0412$Patient_nummer == 201,]
BI_9 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "darkred", lwd = 4,
               data = subdata3, ylim = c(0,21), xlim = c(0,230), 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata4 <- BI0412[BI0412$Patient_nummer == 322,]
BI_10 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "purple", lwd = 4, 
               data = subdata4, ylim = c(0,21), xlim = c(0,230), 
               par.settings = list(fontsize = list(text = 13, points = 18)))


plot2<- BI_6 + as.layer(BI_7) + as.layer(BI_8) + as.layer(BI_9) + as.layer(BI_10) 

###  BI: 13-20  ###
BI1320 <- BIactNew2[BIactNew2$BIcategory == "BI: 13-20",]

BI_11 <- xyplot(BITOT ~ Days, group = Patient_nummer, data = BI1320, type = c("l"), cex = 0.6,col = "grey", 
               strip = FALSE, pch = 16, grid = FALSE, lwd = 2, ylim = c(0,21), xlim = c(0,230), cex.main =0.5,
               xlab = list("Days since stroke"), ylab = list("Barthel Index"), main = "Baseline BI: 13-20", 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata1 <- BI1320[BI1320$Patient_nummer ==124,]
BI_12 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "blue4", lwd = 4, 
               data = subdata1, ylim = c(0,21), xlim = c(0,230), 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata2 <- BI1320[BI1320$Patient_nummer == 145,]
BI_13 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "darkcyan", lwd = 4, 
               data = subdata2, ylim = c(0,21), xlim = c(0,230), 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata3 <- BI1320[BI1320$Patient_nummer == 201,]
BI_14 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "pink", lwd = 4,
               data = subdata3, ylim = c(0,21), xlim = c(0,230), 
               par.settings = list(fontsize = list(text = 13, points = 18)))

subdata4 <- BI1320[BI1320$Patient_nummer == 322,]
BI_15 <- xyplot(BITOT ~ Days, groups = Patient_nummer, type = c("p", "l"), pch = 16, cex = 0.6, col = "purple", lwd = 4, 
                data = subdata4, ylim = c(0,21), xlim = c(0,230), 
                par.settings = list(fontsize = list(text = 13, points = 18)))


plot3<- BI_11 + as.layer(BI_12) + as.layer(BI_13) + as.layer(BI_14) + as.layer(BI_15) 

########## SAVE THE FIGURES########## 
tiff("Figure_1.tiff", width = 4, height = 4, units = 'in', res = 300)
plot1
dev.off()

tiff("Figure_2.tiff", width = 4, height = 4, units = 'in', res = 300)
plot2
dev.off()

tiff("Figure_3.tiff", width = 4, height = 4, units = 'in', res = 300)
plot3
dev.off()


#' ## Dynamic Prediction
#'load functions for prediction
source("Functions.R")

m6.5_no_FM <- lme(BITOT ~ ns(Days, 5)*(Geslacht + StrokeLeeftijd) + 
                    ns(Days, 5)*(MI_UE_base + MI_LE_base) +
                    ns(Days, 5)*(NIHSSTOT + BAMFORD), 
                  random = list(Patient_nummer = pdDiag(form = ~ ns(Days, df = 5))),
                  data = BIactNew2, method = "ML")



DynPlots <- function(model.output = model.output, newdata, timeVar, 
                     main_title = "Dynamic predictions"){
  
# Load individual prediction ------------------------------------
data <- model.output$data
formYx <- formula(model.output)
yOutcome <- formYx[[2]]

IndvPrediction95 <- IndvPred_lme_2(lmeObject = model.output, newdata, timeVar, times = NULL, M = 500, 
                                 interval = "prediction", return_data = TRUE)

IndvPrediction68 <- IndvPred_lme_2(lmeObject = model.output, newdata, timeVar, times = NULL, M = 500, 
                                 interval = "prediction", return_data = TRUE, level = 0.68)

pred95 <- IndvPrediction95[which(!is.na(IndvPrediction95$low)),]
pred68 <- IndvPrediction68[which(!is.na(IndvPrediction68$low)),]

nopred <- IndvPrediction95[which(is.na(IndvPrediction95$low)),]

timeVariable <- pred95[[timeVar]]

xyplot(pred ~ timeVariable , main = main_title, data = pred95,
       type = "l", col = rgb(0.6769,0.4447,0.7114, alpha = 1), lty = c(1, 2, 2), lwd = 3,
       ylim = c(0,30), xlim = c(0,230), ylab = list("BI", cex = 1.5), xlab = list("Days since stroke", cex = 1.5),
       scales = list(x = list(cex = 1.3) , y = list(cex = 1.3)),
       panel = function(x, y,  ...) {
         panel.xyplot(x, y, ...)
         panel.polygon(c(pred95[,"Days"], rev(pred95[,"Days"])), 
                       c(pred95[,"upp"], rev(pred95[,"low"])),
                       col = "bisque", border=NA)
         panel.polygon(c(pred68[,"Days"], rev(pred68[,"Days"])), 
                       c(pred68[,"upp"], rev(pred68[,"low"])),
                       border = NA,
                       col =rgb(0.6769,0.4447,0.7114, alpha = 0.4))
         panel.points(x = nopred[[timeVar]], y = nopred[[yOutcome]], cex = 1.2, pch = 16, col = "black");
         panel.lines(x = rep(tail(nopred[[timeVar]], n = 1), 200), y = seq(-100, 100, length = 200), col = "grey", lty = 3, lwd = 2)
         panel.lines(x = pred95$Days, y = pred95$pred, col = "black", lty = 1, lwd = 2)
       })
}

#' Illustration for Patient 
newPatient <- BIactNew2[BIactNew2$Patient_nummer == 45, ]
newPatient0 <- newPatient[1:1, ]
newPatient1 <- newPatient[1:2, ]
newPatient2 <- newPatient[1:3, ]
newPatient3 <- newPatient[1:4, ]
newPatient4 <- newPatient[1:5, ]
newPatient5 <- newPatient[1:6, ]


p40 <- DynPlots(model.output = m6.5_no_FM, newdata = newPatient0, 
                timeVar = "Days",
                main_title = "")
p41 <- DynPlots(model.output = m6.5_no_FM, newdata = newPatient1, 
                timeVar = "Days",
                main_title = "")
p42 <- DynPlots(model.output = m6.5_no_FM, newdata = newPatient2, 
                timeVar = "Days",
                main_title = "")
p43 <- DynPlots(model.output = m6.5_no_FM, newdata = newPatient3, 
                timeVar = "Days",
                main_title = "")
p44 <- DynPlots(model.output = m6.5_no_FM, newdata = newPatient4, 
                timeVar = "Days",
                main_title = "")
p45 <- DynPlots(model.output = m6.5_no_FM, newdata = newPatient4, 
                timeVar = "Days",
                main_title = "")

grid.arrange(p40, p41, p42, p43, p44, p45, nrow = 3, widths=c(3,3))

grid.arrange(p41, p40, p43, p42, p45, p44, nrow = 3, widths=c(3,3))
