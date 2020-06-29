library(tidyverse)
library(dbdataset)
library(caret)
library(Boruta)
library(pROC)
library(Hmisc)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)

data_base <- read.delim(file = "fulldataset_last version.csv", 
                        header = TRUE, sep = ",", row.names = NULL)

CYP_action   <- as.data.frame(dbdataset::Actions_Enzyme_Drug) %>%
                filter(action == 'inhibitor')
Enzymes_drug <- as.data.frame(dbdataset::Enzymes_Drug) %>% 
                select(parent_key, id, name) %>%
                filter(id == 'BE0002638' | id == 'BE0002362' | id == 'BE0002363') %>%
                rename(drug_id = parent_key, enzyme_id = id, cyp_name = name)

Drugs        <- as.data.frame(dbdataset::Drugs) %>% 
                select(primary_key, name) %>%
                rename(drug_id = primary_key, drug_name = name)

res_ids <- inner_join(Drugs, Enzymes_drug, by=c("drug_id"))

res_ids <- inner_join(res_ids, CYP_action, by=c("enzyme_id"))

res_ids <-res_ids %>% mutate_all(funs(toupper))

data_model <- inner_join(data_base, res_ids, by=c("DrugName" = "drug_name"))


data_model <- data_model %>%
  mutate(CYP2D6_inhibitor = ifelse(cyp_name == "CYTOCHROME P450 2D6", 1, 0),
         CYP3A4_inhibitor = ifelse(cyp_name == "CYTOCHROME P450 3A4", 1, 0),
         CYP3A5_inhibitor = ifelse(cyp_name == "CYTOCHROME P450 3A5", 1, 0))  %>%
  select(-c(drug_id, enzyme_id, cyp_name, action))

data_sum <- data_model %>%
            group_by(DrugName) %>% 
            mutate(Frequency_CYP2D6_inhibitor = sum(CYP2D6_inhibitor),
                   Frequency_CYP3A4_inhibitor = sum(CYP3A4_inhibitor),
                   Frequency_CYP3A5_inhibitor = sum(CYP3A5_inhibitor)) %>%
            mutate(vDILIConcern = recode(vDILIConcern, `VNO-DILI-CONCERN`       = "NO_DILI_CONCERN",
                                                       `AMBIGUOUS DILI-CONCERN` = "NO_DILI_CONCERN",
                                                       `VMOST-DILI-CONCERN`     = "DILI_CONCERN",
                                                       `VLESS-DILI-CONCERN`     = "DILI_CONCERN" )) %>%
            rename(DILLI_Status = vDILIConcern) %>%
            filter(!duplicated(DrugName)) %>%
            select(-c(CYP2D6_inhibitor, CYP3A4_inhibitor, CYP3A5_inhibitor))

data_sum <- data_sum[!(data_sum$Molar_dose_mmol == 0),]

data_sum[,c(3:9,12)] <- scale(log(data_sum[,c(3:9,12)]))
data_sum[,10:11] <- scale(data_sum[,10:11])

write.csv(data_sum, file = 'data_model.csv', row.names = FALSE)

data_full_model <- data_sum[,-1]

boruta_output <- Boruta(DILLI_Status ~ ., 
                        data_full_model, 
                        doTrace = TRUE,
                        mcAdj = TRUE)

plot(boruta_output, cex.axis=.7, las=2, xlab="", main="Variable Importance")

boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
print(boruta_signif)

ind <- sample(2, nrow(data_full_model), replace = TRUE, prob=c(0.7, 0.3))

trainset <-  data_full_model[ind == 1,]
testset  <-  data_full_model[ind == 2,]

## SVM

grid <- expand.grid(sigma = c(0,0.01, 0.02, 0.025, 0.03, 0.04,
                              0.05, 0.06, 0.07,0.08, 0.09, 0.1,
                              0.25, 0.5, 0.75,0.9),
                    C     = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
                              1, 1.5, 2,5))

ctrl <- trainControl(method  = "repeatedcv", # 10fold cross validation
                     number  = 10,
                     repeats = 5, # do 5 repititions of cv
                     summaryFunction = twoClassSummary, # Use AUC to pick the best model
                     classProbs = TRUE)


svm.tune_1 <- train(as.factor(DILLI_Status) ~ Dose_mg + Molar_dose_mmol + Frequency_CYP2D6_inhibitor,
                    data = trainset,
                    method = "svmRadial",
                    tuneLength = 10, 
                    tuneGrid = grid,
                    preProc = c("center","scale"),
                    metric="ROC",
                    trControl=ctrl)

svm.pred <- predict(svm.tune_1,testset)

confusionMatrix(svm.pred,as.factor(testset$DILLI_Status))

svm.probs_1<- predict(svm.tune_1,testset, type="prob")

svm.ROC_1 <- roc(predictor=svm.probs_1$DILI_CONCERN,
                  response=testset$DILLI_Status,
                  levels=rev(levels(as.factor(testset$DILLI_Status))))


plot(svm.ROC_1,main="ROC for SVM", legacy.axes = TRUE)
