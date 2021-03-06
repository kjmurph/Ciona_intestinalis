#install.packages("nlme")
library(nlme)
#install.packages("boot")
library(boot)
#install.packages("hydroGOF")
library(hydroGOF)
#install.packages("JBTools")
library(JBTools)
#install.packages("Metrics")
#library(Metrics)
#install.packages("Hmisc")
library(gridExtra)
#install.packages("data.table")
#library(data.table)

#######################

# Reset working directory

setwd("")

# Combined years dataset
DataComb <- read.csv("Data_Comb_Slim.csv")

# Create 2014 standalone dataset
Data14 <- DataComb[!DataComb$Year == "15",]

Data14_Edit <- Data14[!Data14$SiteCode =="LR" & !Data14$SiteCode =="SA",]

# Create 2015 standalone dataset
Data15 <- DataComb[!DataComb$Year == "14",]

# Edit the 2015 dataset so that the missing values for VarAccel are removed, as well as removing all the sites that are not included
Data_Edit <- Data15[!Data15$SiteCode =="SP" & !Data15$SiteCode =="IN" & !Data15$SiteCode =="LB" & !Data15$SiteCode =="SA" & !Data15$SiteCode =="LR",]

# Edit the combined dataset so that all missing values are removed. This is used for predictions later
Data_Comb_Edit <- DataComb[!is.na(DataComb$MaxSal) & !is.na(DataComb$MinpH) & !DataComb$SiteCode =="SP" & !DataComb$SiteCode =="IN" & !DataComb$SiteCode =="LB" & !DataComb$SiteCode =="SA" & !DataComb$SiteCode =="LR",]

# Edit 2014 dataset so that all missing values are removed
Data_14_Best <- Data14[!is.na(Data14$AvgTemp) & !is.na(Data14$AvgSal) & !is.na(Data14$MaxVarAccel) & !Data14$SiteCode =="SA" & !Data14$SiteCode =="LR",]


#############################
#         Functions         #
#############################

# SiteCodeRaneff1 <- function(Ranef, Data){
#   names(Ranef) <- "Var1"
#   Ranef$SiteCode <- row.names(Ranef)
#   Data <- merge(Data, Ranef, by="SiteCode", all.y=TRUE)
#   Data
# }


#############################
#     Top-down Model        #
#############################

Best_2014.nlme <- nlme(PercentCiona  ~ I(100 *((inv.logit(exp(MaxPop) *Daysafterdeploy^Lag/(exp(Days50)^Lag + Daysafterdeploy^Lag))-0.5)*2)),
                     data=Data_14_Best,
                     fixed = list(  MaxPop ~ I(AvgSal-31), 
                                    Days50 ~ I(AvgSal-31), 
                                    Lag ~I(AvgSal-31)),
                     random = MaxPop ~1|SiteCode,
                     start=c(-0.555676 ,  0,
                             4.058143, 0,  
                             12.101252,  0), control=nlmeControl(maxIter=50), verbose=T)


#############################
#      Model Diagnostics    #
#############################

# Basic model diagnostics for the Top-down Model

plot(augPred(Best_2014.nlme, primary = "Daysafterdeploy",minimum = 0, level = 0:1), main = "Best_2014.nlme")

# Create plot of the fitted values to investigate proportion of random and fixed effects

plot(Best_2014.nlme)
temp <- Data_14_Best
temp$resid <- resid(Best_2014.nlme)
temp$fitted0 <- fitted(Best_2014.nlme, level=0)
temp$fitted1 <- fitted(Best_2014.nlme, level=1)
plot(temp$fitted0, temp$resid, xlab="Population Fitted", ylab="Residual")
abline(h=0)
plot(temp$PercentCiona, temp$fitted0, xlab="Observed", ylab="Population Fitted")
abline(0,1)
lapply(1:length(unique(temp$SiteCode)), function(i, temp){
  temp1 <- temp[temp$SiteCode==unique(temp$SiteCode)[i],]
  plot(temp1$fitted0, temp1$resid, main=unique(temp$SiteCode)[i], xlab="Population Fitted", ylab="Residual")
  abline(h=0)
  plot(temp1$fitted1, temp1$resid, main=unique(temp$SiteCode)[i], xlab="Per Site Fitted", ylab="Residual")
  abline(h=0)
  plot(temp1$PercentCiona, temp1$fitted0, main=unique(temp$SiteCode)[i], xlab="Observed", ylab="Population Fitted")
  abline(0,1)
  plot(temp1$PercentCiona, temp1$fitted1, main=unique(temp$SiteCode)[i], xlab="Observed",ylab="Per Site Fitted")
  abline(0,1)
}, temp=temp)



##############################
#       Model Validation     #
##############################

# Assign name to fitted values from the model
prediction <- fitted(Best_2014.nlme)

# Assign name to fixed-effects fitted values only from the model
prediction_Fixed <- fitted(Best_2014.nlme, level = 0)

# Assign name to the observed values from the data
observation <- Data_14_Best$PercentCiona


# Calculate Mean Absolute Error (MAE) for the full model 
MAE_Full <- mae(prediction, observation)

# Calculate Mean Absolute Error (MAE) for the fixed-effects of model 
MAE_Fixed <- mae(prediction_Fixed, observation)


##############

# Calculate Root Mean Square Error (RMSE) for the random & fixed (full) effects
RMSE_Full <- rmse(prediction, observation)

# Calculate Root Mean Square Error (RMSE) for the fixed-effects only
RMSE_Fixed <- rmse(prediction_Fixed, observation)


##############

# Calculate Model Efficiency for the full model 
MEF_Full <- MEF(prediction, observation)

# Calculate Model Efficiency  for the fixed-effects of model 
MEF_Fixed <- MEF(prediction_Fixed, observation)


##############

# Calculate Spearman rho for the full model 
S_rho_Full <- cor(prediction, observation, method = "spearman")

# Calculate Spearman rho for the fixed-effects of model 
S_rho_Fixed <- cor(prediction_Fixed, observation, method = "spearman")


Mod_Val_Results <- list( MAE_Full, MAE_Fixed, RMSE_Full, RMSE_Fixed, MEF_Full, MEF_Fixed, S_rho_Full, S_rho_Fixed)


  
#####################################################################################################
#  Combine observed values and fitted (with fixed-effect separate as well) values into single file  #
#####################################################################################################

# Create sub-dataset that matches that used in the given model 

datasub <- Data_14_Best

# requires datasub to be the appropriate subset used in the model

# Create copy of datasub to avoid altering original

datasubcopy <- datasub

# Add new column "OF" (Obsered/Fitted) and fill in with "observed" values corresponding to the observed "PercentCiona" values 
datasubcopy$OF <- "observed"

# Assign the fixed-effect fitted values a name 

fixfittedvals <- fitted(Best_2014.nlme, level = 0)

# Assign a new datasubcopy that will be for the fitted values from Mod6.nlme

datasubwFR <- datasubcopy


# Replace the observed "PercentCiona" values with Global_Model1 fixed effect values

datasubwFR$PercentCiona <- fixfittedvals

# Add new column "FR" (Fixed/Random) and fill in with "Fixed" values corresponding to the Fixed-effect fitted "PercentCiona" values 

datasubwFR$FR <- "Fixed"

# Assign a new datasubcopy that will be for the Random effect values from Global_Model1

datasubwRF <- datasubcopy

# Assign the fixed-effect fitted values a name 

Ranfittedvals <- fitted(Best_2014.nlme, level = 1)


# Replace the observed "PercentCiona" values with Global_Model1 fixed effect values

datasubwRF$PercentCiona <- Ranfittedvals


# Add new column "FR" (Fixed/Random) and fill in with "Random" values corresponding to the Random-effect fitted "PercentCiona" values 

datasubwRF$FR <- "Random"

# Merge the two versions of datasubcopy that have the fitted and observed values.

# Add new column "OF" (Obsered/Fitted) and fill in with "observed" values corresponding to the observed "PercentCiona" values 
datasubcopy$FR <- "observed"


datasubcomb <- rbindlist(list(datasubwFR,datasubwRF, datasubcopy))


datasubcomb1 <- rbindlist(list(datasubwFR,datasubwRF))


############################################################
#     Plot the Fixed and Random values on the same plot    #
############################################################

Model_Fit_Plot_2014 <- xyplot(PercentCiona~Daysafterdeploy |factor(SiteCode, levels(reorder(SiteCode, MeanPercentCiona2014))),
                              groups = FR,
                              layout = c(1,11),
                              xlab = "",
                              ylab = "",
                              ylab.right = list(
                                label=substitute(paste(italic('C. intestinalis'), " relative abundance (%)")), cex = 1.25),
                              data = datasubcomb1,
                              par.settings = list(
                                layout.widths = list(
                                  right.padding = 1, left.padding = -2)),
                              type = c("a"), col = c("black"), lty = c("dotted", "solid"), lwd = c(2, 2),
                              scales = list(
                                y=list(
                                  at=seq(0,100,50),labels = c("0", "50", "100"), cex=0.75, font=2, tck =c(0,1), alternating = c(2, 0)),
                                x=list(
                                  at=seq(30,150,30), cex=.85, tck =c(1,0), alternating = c(1, 1), font=1)),
                              xlim = c(5, 165),
                              strip = FALSE)


Ciona_Plot_2014 <- xyplot(PercentCiona~Daysafterdeploy |factor(SiteCode, levels(reorder(SiteCode, MeanPercentCiona2014))),
                          layout = c(1,11),
                          xlab = "",
                          ylab= list('Site', cex = 1.25),
                          par.settings = list(
                            layout.widths = list(
                              right.padding = -2, left.padding = 1)),
                          data = Data14_Edit,
                          type = c("a","p"), col = c("black"),
                          scales = list(
                            y=list(
                              at=seq(0,100,50),labels = c("0", "50", "100"), cex=0.6, font=2, tck =c(0,1), alternating = c(0, 0)),
                            x=list(
                              at=seq(30,150,30), cex=.85, tck =c(1,0), alternating = c(1, 1), font=1)),
                          xlim = c(5, 165),
                          strip = FALSE,
                          strip.left = strip.custom(horizontal = FALSE, bg = 'white'),
                          par.strip.text=list(cex=.75, font=2))


grid.arrange(Ciona_Plot_2014, Model_Fit_Plot_2014,
             widths=c(1, 1), 
             nrow=1, ncol=2,
             bottom = "Time (Days after deployment)")



########################################################
#     2014 Best Model  MEF Nested Cross Validation     #
########################################################


# 10-repeated 10-fold cross-validation
runs = 10
nfold = 10
results = c()

for (run in 1:runs)
{
  split = sample(rep(1:nfold, length = 330), 330) # Split n in 10 sets
  resample = lapply(1:nfold, function(x,spl) list(cal=which(spl!=x), val=which(spl==x)), spl=split) # Create values for 10 cal/val sets. Don't worry about the details in this line, just check the result
  
  # Here we loop through all cal/val combinations
  for(fold in 1:nfold)
  {
    # Create the cal/val data sets for the current 'fold'
    cal = Data_14_Best[resample[[fold]]$cal,]
    val  = Data_14_Best[resample[[fold]]$val,]
    
    # Fit a model with the current calibration data set
    model = Best_2014.nlme # Note that we are reusing the formula already defined
    
    # Calculate MEF and put it in results
    observations = val$PercentCiona
    predictions = predict(model, newdata=val)
    
    
    #MAE = mae(observations, predictions)
    #RMSE = rmse(observations, predictions)
    MEF = MEF(predictions, observations)
    #cor = cor(predictions, observations, method = "spearman")
    
    #results = c(results, MAE) # This attaches 'MAE' to any previous content of 'results'
    #results = c(results, RMSE) # This attaches 'RMSE' to any previous content of 'results'
    results = c(results, MEF) # This attaches 'MEF' to any previous content of 'results'
    #results = c(results, cor) # This attaches 'cor' to any previous content of 'results'
  }
}

Cross_Val_Results <- mean(results)


##########################################################
#                     Predictions                        #
##########################################################


Predict_2015 <- predict(Best_2014.nlme, Data_Edit[!Data_Edit$SiteCode =="SH",])
Predict_2015_Fixed <- predict(Best_2014.nlme, level = 0, Data_Edit[!Data_Edit$SiteCode =="SH",])


RMSE_2015 <- rmse(Data_Edit[!Data_Edit$SiteCode =="SH",]$PercentCiona, Predict_2015)
RMSE_2015_Fixed <- rmse(Data_Edit[!Data_Edit$SiteCode =="SH",]$PercentCiona, Predict_2015_Fixed)


MAE_2015 <- mae(Data_Edit[!Data_Edit$SiteCode =="SH",]$PercentCiona, Predict_2015)
MAE_2015_Fixed <- mae(Data_Edit[!Data_Edit$SiteCode =="SH",]$PercentCiona, Predict_2015_Fixed)


MEF_2015 <- MEF(Predict_2015, Data_Edit[!Data_Edit$SiteCode =="SH",]$PercentCiona)
MEF_2015_Fixed <- MEF(Predict_2015_Fixed, Data_Edit[!Data_Edit$SiteCode =="SH",]$PercentCiona)


S_rho_2015 <- cor(Predict_2015, Data_Edit[!Data_Edit$SiteCode =="SH",]$PercentCiona, method = "spearman")
S_rho_2015_Fixed <- cor(Predict_2015_Fixed, Data_Edit[!Data_Edit$SiteCode =="SH",]$PercentCiona, method = "spearman")

#####################################################
  
Predict_Comb <- predict(Best_2014.nlme, Data_Comb_Edit[!Data_Comb_Edit$SiteCode =="SH",])
Predict_Comb_Fixed <- predict(Best_2014.nlme, level = 0, Data_Comb_Edit[!Data_Comb_Edit$SiteCode =="SH",])


RMSE_Comb <- rmse(Data_Comb_Edit[!Data_Comb_Edit$SiteCode =="SH",]$PercentCiona, Predict_Comb)
RMSE_Comb_Fixed <- rmse(Data_Comb_Edit[!Data_Comb_Edit$SiteCode =="SH",]$PercentCiona, Predict_Comb_Fixed)


MAE_Comb <- mae(Data_Comb_Edit[!Data_Comb_Edit$SiteCode =="SH",]$PercentCiona, Predict_Comb)
MAE_Comb_Fixed <- mae(Data_Comb_Edit[!Data_Comb_Edit$SiteCode =="SH",]$PercentCiona, Predict_Comb_Fixed)


MEF_Comb <- MEF(Predict_Comb, Data_Comb_Edit[!Data_Comb_Edit$SiteCode =="SH",]$PercentCiona)
MEF_Comb_Fixed <- MEF(Predict_Comb_Fixed, Data_Comb_Edit[!Data_Comb_Edit$SiteCode =="SH",]$PercentCiona)


S_rho_Comb <- cor(Predict_Comb, Data_Comb_Edit[!Data_Comb_Edit$SiteCode =="SH",]$PercentCiona, method = "spearman")
S_rho_Comb_Fixed <- cor(Predict_Comb_Fixed, Data_Comb_Edit[!Data_Comb_Edit$SiteCode =="SH",]$PercentCiona, method = "spearman")


