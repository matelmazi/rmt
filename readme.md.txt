**Personalized Dynamic Prediction for Barthel Index**

This repository includes the code to obtain dynamic predictions and to perform cross-validation.

The data set consists of 13 variables as shown below.

Patient_nummer: discrete numerical variable 
StrokeLeeftijd: discrete numerical variable 
Geslacht: numerical as male(1), female (2)
Bamford: Laci(1), Paci(2), Taci(3)
MITOTBE: discrete numerical variable from 0 to 100
MITOTOE: discrete numerical variable from 0 to 100
FMBETOT: discrete numerical variable from 0 to 66
NIHSSTOT: discrete numerical variable from 0 to 24
WRITINGHAND: right(1), left(2)
Days: discrete numerical variable from 0 to 296
Partner: yes(1), no(0)
BITOT: discrete numerical variable from 0 to 20
BITOTact: discrete numerical variable from 0 to 16

The different scripts can be runned with the specific order:
1. Model_building: to explore the data(descriptives, existence of missing values,etc), define and analyze different model structures.
2. Cross_validation: code to perform cross validation
3. CV results: mostly error plots to check lowest model error from the cross-validation
4. Plots_Figures: R code that was used to obtain the figures
Functions should be loaded with the command source() at all previous scripts to use the IndvPred_lme_2 function were needed. 
