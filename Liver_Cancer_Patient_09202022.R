###################################
##title: "Liver Cancer Analysis"##
##Author: yfu2015               ##
##date: "2022-09-20"            ##
##################################



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### Install packages ######

if(!require(tidyverse)) 
  install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(caret)) 
  install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(data.table)) 
  install.packages("data.table", repos = "http://cran.us.r-project.org")
if(!require(rpart)) 
  install.packages("rpart", repos = "http://cran.us.r-project.org")
if(!require(matrixStats)) 
  install.packages("matrixStats", repos = "http://cran.us.r-project.org")
if(!require(dslabs)) 
  install.packages("dslabs", repos = "http://cran.us.r-project.org")
#if(!require(genefilter)) install.packages("genefilter", repos = "http://cran.us.r-project.org")
if(!require(gam)) 
  install.packages("gam", repos = "http://cran.us.r-project.org")
if(!require(gridExtra)) 
  install.packages("gridExtra", repos = "http://cran.us.r-project.org")
if(!require(randomForest)) 
  install.packages("randomForest", repos = "http://cran.us.r-project.org")
if(!require(tinytex)) 
  install.packages("tinytex", repos = "http://cran.us.r-project.org")

##### loading packages #####

library(tidyverse)
library(caret)
library(data.table)
library(rpart)
library(matrixStats)
library(dslabs)
#library(genefilter)
library(gam)
library(gridExtra)
library(randomForest)
library(tinytex)

#tinytex::install_tinytex()
#update.packages(ask = FALSE, checkBuilt = TRUE)
#update.packages("tidyverse")
#update.packages()
#installed.packages()
#remove.packages(pkgs=row.names(x=installed.packages(priority="NA")))
#old.packages()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#liver.df <- read.csv("C:/Users/yfu/Desktop/LiverPatients_download.csv")
liver.df <- read.csv(url(
           "https://raw.githubusercontent.com/yfu2021/liver_cancer/master/LiverPatients_download.csv"
              )
            )
summary(liver.df)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
liver <- liver.df %>% 
         mutate(Diagnosis = factor(ifelse(Dataset==1, 1, 0)),
                Gender = as.numeric(ifelse(Gender=="Female", 0, 1)),
                Albumin_and_Globulin_Ratio = ifelse(is.na(Albumin_and_Globulin_Ratio),
                                        median(Albumin_and_Globulin_Ratio,na.rm=TRUE),
                                        Albumin_and_Globulin_Ratio) # Replace NA with median value
              ) %>% select(-Dataset)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------

str(liver)

sum(liver$Diagnosis==1)
sum(liver$Diagnosis==0)

All_Patients_Age <- liver$Age
hist(All_Patients_Age, main="Histogram of Age")   
Cancer <- liver %>% filter(Diagnosis==1)
Liver_Cancer_Patients_Age <- Cancer$Age
hist(Liver_Cancer_Patients_Age, main="Histogram of Age for Cancer Patient Only")


#Distribution of all predictors vs target value

p1 <- liver %>% ggplot(aes(Diagnosis, Age)) + 
                geom_jitter(width = 0.1, alpha = 0.2)
p2 <- liver %>% ggplot(aes(Diagnosis,Gender)) + 
                geom_jitter(width = 0.1, alpha = 0.2)
p3 <- liver %>% ggplot(aes(Diagnosis,Total_Bilirubin)) + 
                geom_jitter(width = 0.1, alpha = 0.2) 
p4 <- liver %>% ggplot(aes(Diagnosis,Direct_Bilirubin)) + 
                geom_jitter(width = 0.1, alpha = 0.2) 
p5 <- liver %>% ggplot(aes(Diagnosis,Alkaline_Phosphotase)) + 
                geom_jitter(width = 0.1, alpha = 0.2) 
p6 <- liver %>% ggplot(aes(Diagnosis,Alamine_Aminotransferase)) + 
                geom_jitter(width = 0.1, alpha = 0.2) 
p7 <- liver %>% ggplot(aes(Diagnosis,Aspartate_Aminotransferase)) + 
                geom_jitter(width = 0.1, alpha = 0.2)
p8 <- liver %>% ggplot(aes(Diagnosis,Total_Protiens)) + 
                geom_jitter(width = 0.1, alpha = 0.2)
p9 <- liver %>% ggplot(aes(Diagnosis,Albumin)) + 
                geom_jitter(width = 0.1, alpha = 0.2) 
p10 <- liver %>% ggplot(aes(Diagnosis,Albumin_and_Globulin_Ratio)) + 
                geom_jitter(width = 0.1, alpha = 0.2)


grid.arrange(p1, p2,p3,p4,p5,p6,p7,p8,p9,p10, nrow=2, ncol = 5)


#Check if there is any data element with very few non-unique values or close to zero variation

nearZeroVar(liver)



## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
liver.df <- liver.df %>% mutate(Gender=as.numeric( ifelse(Gender=="Female", 0, 1) ))

cor(liver.df, use="pairwise.complete")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# centering and scaling on all predictors

options(digits = 3)

#set.seed(1) # if using R 3.5 or earlier
set.seed(1, sample.kind = "Rounding") # if using R 3.6 or later


# setp 1 - center and scale 10 predictors

Predictors <- liver %>%  select(-Diagnosis)
Diagnosis <- liver$Diagnosis  ###%>%  select(Diagnosis)  

x_centered <- sweep(Predictors, 2, colMeans(Predictors))
x_scaled <- sweep(x_centered, 2, colSds(as.matrix(Predictors)), FUN = "/")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
cor(x_scaled, use="pairwise.complete") ##%>% knitr::kable()
image(as.matrix(cor(x_scaled, use="pairwise.complete")), axes = TRUE, 
      main = "Correlation of All Predictors")  


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# principal component analysis

pca <- prcomp(x_scaled)
summary(pca) 


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Some basic analysis on principal component analysis 

# 1. By looking at the distribution of PC1 and PC2, Liver cancer patients 
#    tend to have larger values of PC2 than non-liver cancer patient

data.frame(pc_1 = pca$x[,1], pc_2 = pca$x[,2], 
           type=Diagnosis ) %>%
  ggplot(aes(pc_1, pc_2, color = type)) +
  geom_point()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#2. Liver cancer patients tend to have larger variance of PC3 
#   than non-liver cancer patient.

data.frame(pc_3 = pca$x[,3], pc_4 = pca$x[,4], 
           type=Diagnosis ) %>%
  ggplot(aes(pc_3, pc_4, color = type)) +
  geom_point()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#3. Distribution of IQRs from PC 1 through PC 10

data.frame(type = Diagnosis, pca$x[,1:10]) %>%        
  gather(key = "PC", value = "value", -type) %>%
  ggplot(aes(PC, value, fill = type)) +
  geom_boxplot()


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Split the Scaled data to train set and test set

set.seed(1, sample.kind="Rounding") # if using R 3.5 or earlier, use `set.seed(1)`
test_index <- createDataPartition(Diagnosis, times = 1, p = 0.1, list = FALSE)
test_x <- x_scaled[test_index,]   ## str(test_x)
test_y <- Diagnosis[test_index]    ## length(test_y)

train_x <- x_scaled[-test_index,]   ## str(train_x)
train_y <- Diagnosis[-test_index]   ## length(train_y)


## ----setup, warning=FALSE--------------------------------------------------------------------------------------------------------------------------------------------
#Logistic regression model   

# set.seed(1) if using R 3.5 or earlier
set.seed(1, sample.kind = "Rounding") # if using R 3.6 or later
train_glm <- train(train_x, train_y, method = "glm")
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
glm_pred <- predict(train_glm, test_x)
logistic_Accuracy2 <-  mean(glm_pred == test_y)

glm_sensitivity <- rbind(sensitivity(glm_pred, test_y), specificity(glm_pred, test_y))
rownames( glm_sensitivity )  <-  c("sensitivity","specificity") 
colnames( glm_sensitivity )  <-  c("Logistic Regression") 

# Top 5 important Predictors

a <- varImp(train_glm)[[1]]
varImp_df_glm <- data.frame(matrix(c(rownames(a), as.numeric(a[,1])), nrow=10, ncol=2, 
                                   dimnames=list(c(seq(1:10)), c("Predictor", "Value"))))
varImp_df_glm <- varImp_df_glm %>% mutate(Value=as.numeric(Value)) %>% arrange(desc(Value)) 
rownames(varImp_df_glm) <- seq(1:10)
Top_5_Glm_Predictors <- varImp_df_glm[1:5,1]   
Top_5_GLM_predictors <- data.frame(Top_5_Glm_Predictors)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#K-nearest neighbors model K-nearest neighbors model    

# set.seed(1)
set.seed(1, sample.kind = "Rounding") # simulate R 3.5
control <- trainControl(method = "cv", number = 10, p = .9)
train_knn <- train(train_x, train_y, method = "knn", 
                   tuneGrid = data.frame(k = seq(9, 71, 2)),
                   trControl = control)
ggplot(train_knn, highlight=TRUE, main="k nearst neighbor Model accuracy vs k")
train_knn$bestTune
knn_pred <- predict( train_knn, test_x)
K_nearst_Accuracy2 <- mean(knn_pred== test_y)

knn_sensitivity <- rbind(sensitivity(knn_pred, test_y), specificity(knn_pred, test_y))
rownames( knn_sensitivity )  <-  c("sensitivity","specificity") 
colnames( knn_sensitivity )  <-  c("K_Nearst_Neighbor") 

# Top 5 important Predictors
a <- varImp(train_knn)[[1]]
varImp_df_knn <- data.frame(matrix(c(rownames(a), as.numeric(a[,1])), nrow=10, ncol=2, 
                                   dimnames=list(c(seq(1:10)), c("Predictor", "Value"))))
varImp_df_knn <- varImp_df_knn %>% mutate(Value=as.numeric(Value)) %>% arrange(desc(Value)) 
rownames(varImp_df_knn) <- seq(1:10)
Top_5_knn_Predictors <- varImp_df_knn[1:5,1]   
Top_5_knn_predictors <- data.frame(Top_5_knn_Predictors)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#random Forest model
                  
# set.seed(1)
set.seed(1, sample.kind = "Rounding") # simulate R 3.5   

tuning <- data.frame(mtry = c(1,2,3,4,5,6,7,8))
train_rf <- train(train_x, train_y,
                  method = "rf",
                  tuneGrid = tuning,
                  importance = TRUE)

ggplot(train_rf,highlight=TRUE, title="random Forest Model Acuracy distribution")

rf_pred <- predict(train_rf, test_x)
rf_Accuracy2 <- mean(rf_pred == test_y)

rf_sensitivity <- rbind(sensitivity(rf_pred, test_y), specificity(rf_pred, test_y))
rownames( rf_sensitivity )  <-  c("sensitivity","specificity") 
colnames( rf_sensitivity )  <-  c("random Forest") 

# Top 5 important Predictors
a <- varImp(train_rf)[[1]]  
varImp_df_rf <- data.frame(matrix(c(rownames(a), as.numeric(a[,1])), nrow=10, ncol=2, 
                                  dimnames=list(c(seq(1:10)), c("Predictor", "Value"))))
varImp_df_rf <- varImp_df_knn %>% mutate(Value=as.numeric(Value)) %>% arrange(desc(Value)) 
rownames(varImp_df_rf) <- seq(1:10)
Top_5_rf_Predictors <- varImp_df_rf[1:5,1]   
Top_5_rf_predictors <- data.frame(Top_5_rf_Predictors)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Local Polynomial Regression Model

# set.seed(1)
set.seed(1, sample.kind = "Rounding") # simulate R 3.5
grid <- expand.grid(span = seq(0.15, 0.65, len = 10), degree = 1)
train_loess <- train(train_x, train_y, method = "gamLoess", tuneGrid=grid)
knitr::opts_chunk$set(warning = FALSE, message = FALSE) 
ggplot(train_loess, highlight=TRUE)
loess_pred <- predict(train_loess, test_x)
Loess_Accuracy2 <- mean(loess_pred == test_y)   

Loess_sensitivity <- rbind(sensitivity(loess_pred, test_y), specificity(loess_pred, test_y))
rownames( Loess_sensitivity )  <-  c("sensitivity","specificity") 
colnames( Loess_sensitivity )  <-  c("Local Polynomial Regression") 

#### Top 5 important Predictors ####
a <- varImp(train_loess)[[1]]
varImp_df_loess <- data.frame(matrix(c(rownames(a), as.numeric(a[,1])), nrow=10, ncol=2, 
                                     dimnames=list(c(seq(1:10)), c("Predictor", "Value"))))
varImp_df_loess <- varImp_df_loess %>% mutate(Value=as.numeric(Value)) %>% arrange(desc(Value)) 
rownames(varImp_df_loess) <- seq(1:10)
Top_5_Loess_Predictors <- varImp_df_loess[1:5,1]   
Top_5_Loess_predictors <- data.frame(Top_5_Loess_Predictors)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Format/combine the results from all ML models

Accuracy_results <- data_frame(Method = "Logistic Regression Model", Accuracy = logistic_Accuracy2)
Accuracy_results <- bind_rows(Accuracy_results,
                              data_frame(Method="K nearst neighbors Model",
                                         Accuracy = K_nearst_Accuracy2 ))
Accuracy_results <- bind_rows(Accuracy_results,
                              data_frame(Method="Random Forest",  
                                         Accuracy =  rf_Accuracy2 ))
Accuracy_results <- bind_rows(Accuracy_results,
                              data_frame(Method="Local Polynomial Regression Model",  
                                         Accuracy =  Loess_Accuracy2 ))
Accuracy_results %>% knitr::kable(align='c')


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Format Sensitivity and Specificity Results of all Models

cbind(glm_sensitivity , knn_sensitivity, rf_sensitivity, Loess_sensitivity) %>% 
      knitr::kable(align='c')


## --------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Format top 5 important predictors/features 

cbind(Top_5_GLM_predictors, Top_5_knn_predictors, Top_5_rf_predictors, Top_5_Loess_predictors)%>% 
      knitr::kable()

