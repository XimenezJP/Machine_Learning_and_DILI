library(tidyverse)
library(dbdataset)
library(caret)
library(mlbench)
library(glmnet)
library(Boruta)
library(pROC)
library(Hmisc)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)

## Regularized Random Forest (RRF) -------------------------------------
set.seed(100)
RRF_model <- train(DILLI_Status ~ ., 
                   data=data_full_model, 
                   method="RRF")

RRF_importance <- varImp(RRF_model, scale=FALSE)

print(RRF_importance)

plot(RRF_importance, main='Variable Importance')

## Least Absolute Shrinkage and Selection Operator (LASSO) --------------
set.seed(100)

x <- as.matrix(data_full_model[,-1]) # all X vars
y <- as.double(as.matrix(ifelse(data_full_model[,1] == 'NO_DILI_CONCERN', 0, 1)))

cv.lasso <- cv.glmnet(x, 
                      y, 
                      family='binomial', 
                      alpha=1, 
                      parallel=TRUE, 
                      standardize=TRUE, 
                      type.measure='auc')
plot(cv.lasso)

plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE)

cat('Min Lambda: ', cv.lasso$lambda.min, '\n 1Sd Lambda: ', cv.lasso$lambda.1se)

df_coef <- round(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.min)), 2)
df_coef[df_coef[, 1] != 0, ]


## Step wise Forward and Backward Selection ---------------------------
base.mod <- glm(DILLI_Status ~ 1, data = data_full_model, family = "binomial")  

all.mod <- glm(DILLI_Status ~ .,  data = data_full_model, family = "binomial") 

stepMod <- step(base.mod, 
                scope = list(lower = base.mod, upper = all.mod), 
                direction = "both", 
                trace = 0, 
                steps = 1000) 


shortlistedVars <- names(unlist(stepMod[[1]])) 
shortlistedVars <- shortlistedVars[!shortlistedVars %in% "(Intercept)"] # remove intercept

print(shortlistedVars)

## Genetic Algorithm ---------------------------------------------------
ga_ctrl <- gafsControl(functions = rfGA,
                       method = "cv",
                       repeats = 3)

set.seed(100)
ga_obj <- gafs(x=data_full_model[,-1], 
               y=data_full_model$DILLI_Status, 
               iters = 500, 
               gafsControl = ga_ctrl)

ga_obj
print(ga_obj$optVariables)

## Simulated Annealing -------------------------------------------------
sa_ctrl <- safsControl(functions = rfSA,
                       method = "repeatedcv",
                       repeats = 3,
                       improve = 5)


set.seed(100)
sa_obj <- safs(x=data_full_model[,-1], 
               y=data_full_model$DILLI_Status,
               safsControl = sa_ctrl)

sa_obj
print(sa_obj$optVariables)
