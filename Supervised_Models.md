# ML-KRAS-Data
KRAS healthcare data
library(readxl)
library(caret)
library(e1071)
library(randomForest)
library(mass)
library(doSNOW)
dat <- read_excel("~/Desktop/KRAS_features_normalbyBoth.xlsx")
dat1<- read_excel("~/Desktop/Mayo Clinic data/Key_NoPHI.xlsx")

write.csv(a, file = "modified.csv")

dat2<- dat1[-c(1,101,102),]
a<- as.vector(NULL)
for(i in 1:nrow(dat2)) {
  if(str_detect(dat2[i,3], "Muta") == TRUE) {
    a <- c(a, "Mutant")
  } else if(str_detect(dat2[i,3], "Wild Type") == TRUE) {
    a<- c(a, "Wild")
  } else {
    a<- c(a,"Mutant")
  }
}

clean <- data.frame(dat[,-1], a)
colnames(clean)[43] <- "PType"

## remove extra columns not used for analysis
a <- clean[, -c(1,2,3,4,5)]

## LOOCV with random forest not tuned
library(caret)
cv.train<- trainControl(method = "LOOCV")

first.rf <- train(x = a[,-38], y = a$PType, method = "rf", ntree = 1000, trControl = cv.train)
first.rf

## random forest not tuned
rff<- randomForest(a[,-38], a[,38], importance = TRUE, ntree = 1000)
rff
varImpPlot(rff)

## random forest with reduced variables(acc to variable importance plot) -----------> acc = 63.34%
reducedModel <- a[c(26,31,29,33,35,25, 16,22,15,14,17,20,19, 4,2,3,8,11,7,10, 38)]
rf_reduced1 <- train(x = reducedModel[,-21], y = reducedModel[,21], method = "rf", trControl = cv.train, ntree = 500)
rf_reduced1$finalModel


## tune random forest with reduced variables(acc to variable importance plot) and tuning parameters
reducedModel <- a[c(26,31,29,33,35,25, 16,22,15,14,17,20,19, 4,2,3,8,11,7,10, 38)]
tngrid <- expand.grid(.mtry = c(1:20))
modellist <- list()
for(ntree in c(500,750,1000,1200)){
  rf_reduced <- train(x = reducedModel[,-21], y = reducedModel[,21], method = "rf", ntree = ntree, tuneGrid = tngrid, trControl = cv.train)
  key <- toString(ntree)
  modellist[[key]] <- rf_reduced
}
plot(modellist$`500`)
plot(modellist$`750`)
plot(modellist$`1000`)
plot(modellist$`1200`)  ### accuracy didnt increase, remains same (63.63%)

## Using 20 fold cross validation on tuned values of ntree = 1500 ------> 62.64%
set.seed(100)
cv.train1 <- trainControl(method = "repeatedcv", number = 20, repeats = 2)
rf_cv<- train(x = reducedModel[,-21], y = reducedModel[,21], method = "rf", trControl = cv.train1, ntree = 1500)
rf_cv$finalModel

## Using 30 fold cross validation on tuned values of ntree = 1500 --------> 64.65 %
set.seed(100)
cv.train1 <- trainControl(method = "repeatedcv", number = 30, repeats = 4)
rf_cv<- train(x = reducedModel[,-21], y = reducedModel[,21], method = "rf", trControl = cv.train1, ntree = 1500)
rf_cv$finalModel
rf_cv

#######-----------------------------------------------#################

## SVM MODELS

library(kernlab)
library(e1071)

cv.train2 <- trainControl(method = "LOOCV")
svm1 <- train(x = a[,-38], y = a[,38], method = "svmRadial", sigma = 0.05, cost = 1.5, trControl = cv.train2)
svm1
svm1$finalModel

## svm to tune parameters -------> result gamma = 0.01, cost = 1.2
cv.train2<- tune.control(sampling = "boot")
svm_tune <- tune(svm, PType ~ ., data = a, ranges = list(gamma = c(0.01,0.02, 0.03, 0.04,0.05,0.07,0.1,0.15), cost = c(0.8,0.9,1,1.2,1.5,1.7,2)), tunecontrol = cv.train2)
svm_tune$best.parameters
svm_tune$performances
plot(svm_tune$performances[,1],svm_tune$performances[,3])

## svm using tuned parameters

cv.train2 <- trainControl(method = "LOOCV")
svm3 <- train(x = a[,-38], y = a[,38], method = "svmRadial", sigma = 0.03, cost = 1, trControl = cv.train2)
svm3$finalModel

svm3$results

### variable selection using correlation
library(mlbench)
co <- cor(a[,-38])
hig_corr <- findCorrelation(co, cutoff = 0.5)
print(hig_corr)
c1<- a[,c(2 ,27  ,4  ,8  ,9  ,1 ,29 ,10 ,11 ,12 ,33  ,7 ,35 ,28 ,31 ,25  ,3 ,17 ,24,  5, 22, 14, 19,23,6,30, 38)]
View(c1)

## variable selection using var importance
tr<- trainControl(method = "repeatedcv", number = 20, repeats = 4)
importance_estimation <- train(PType ~., data = c, method = "rf", ntree = 1000, preProcess = "scale", trControl = tr)
importanc <- varImp(importance_estimation, scale = FALSE)
plot(importanc)

## applying SVM with tuned parameters and Variables with importance
cv.train2 <- trainControl(method = "LOOCV")
svm4 <- train(x = c[,-27], y = c[,27], method = "svmRadial", sigma = 0.03, cost = 1, trControl = cv.train2)
svm4$finalModel

## using test and train dataset, applying SVM with tuned parameters and Variables with importance
set.seed(100)
vec <- sample(1:nrow(c), floor(nrow(c)*0.75))
c_train <- c[vec, ]
c_test <- c[-vec, ]

cv.train2 <- trainControl(method = "LOOCV")
svm5 <- train(x = c_train[,-27], y = c_train[,27], method = "svmRadial", sigma = 0.04, cost = 0.5, trControl = cv.train2)
svm5

pred <- predict(svm5, c_test[,-27])
table(pred, c_test[,27])
svm5$bestTune

## using svm e1071 package
svm6 <- svm(PType ~ ., data = c_train, kernel = "radial", cost = 1.2, gamma = 0.1)
svm6
pred1<- predict(svm6, c_test[,-27])
table(true = c_test[,27], predicted = pred1)


### Recursive feature Elimination method to select variables ---------> 4 features are efficient in predicting the results
### RFE == It eliminates features with less weight and constructs model iteratively.

tr.rfe <- rfeControl(functions = rfFuncs, method = 'repeatedcv', number = 20, repeats = 2)
varByRFE <- rfe(x = c1[,-27], y = c1[,27], sizes = c(4:20), rfeControl = tr.rfe)
varByRFE$variables[,c(4,5)]
predictors(varByRFE)
plot(varByRFE)


## ensemble models & Stacking

library(caretEnsemble)
cv.control <- trainControl(method = "repeatedcv", number = 10, repeats = 3, savePredictions=TRUE, classProbs=TRUE  )
algolist <- c('knn', 'svmRadial','rpart','glm', 'nnet')


mod <- caretList(PType ~ ., data = c_train, trControl = cv.control, methodList = algolist)

result <- resamples(mod)
summary(result)
dotplot(result)


## stacking model
stack_control <- trainControl(method = "repeatedcv", number = 20, repeats = 4, savePredictions = TRUE, classProbs = TRUE)

## stack model on randomforest
stack_rf <- caretStack(mod, method = 'rf', metric = 'Accuracy', trControl = stack_control)
stack_rf

pred2<- predict(stack_rf, c_test[,-27], se = FALSE, level = 0.95)
table(pred2, c_test[,27])  ### -----------> accuracy = 50 % in new test dataset


## PCA implementation for the data
## divide data into test and train wrt three Texture features
glcm <- a[ ,c(1:13,38)]
lbp <- a[ ,c(14:23,38)]
gabor <- a[ ,24:38]

set.seed(100)
vec1<- sample(1:nrow(glcm), floor(nrow(glcm)*0.75))
glcm.train <- glcm[vec1, ]
glcm.test <- glcm[-vec1, ]

lbp.train <- lbp[vec1, ]
lbp.test <- lbp[-vec1, ]

gabor.train <- gabor[vec1, ]
gabor.test <- gabor[-vec1, ]


## apply PCA function on training data of glcm texture feature

p1_comp <- prcomp(glcm.train[,-14], scale. = T)
View(p1_comp$rotation) ### ---------> rotated PCA is stored here
dim(p1_comp$x) ### to know how many data is there in PCA (rows * Columns)

st1_dev <- p1_comp$sdev
pr1_var <- st1_dev^2

prop1_var <- pr1_var/sum(pr1_var)

plot(cumsum(prop1_var), xlab = "Principal Component", ylab = "cummulative measure", type = 'b')

train.glcm_pca_comp <- data.frame(Ptype = glcm.train[,14], p1_comp$x[,1:4])

## apply PCA function on training data of lbp texture feature

p2_comp <- prcomp(lbp.train[,-11], scale. = T)

st2_dev <- p2_comp$sdev
pr2_var <- st2_dev^2

prop2_var <- pr2_var/sum(pr2_var)

plot(cumsum(prop2_var), xlab = "Principal Component", ylab = "cummulative measure", type = 'b')

train.lbp_pca_comp <- data.frame(Ptype = lbp.train[,11], p2_comp$x[,1:4])

## apply PCA function on training data of gabor texture feature

p3_comp <- prcomp(gabor.train[,-15], scale. = T)

st3_dev <- p3_comp$sdev
pr3_var <- st3_dev^2

prop3_var <- pr3_var/sum(pr3_var)

plot(cumsum(prop3_var), xlab = "Principal Component", ylab = "cummulative measure", type = 'b')

train.gabor_pca_comp <- data.frame(Ptype = gabor.train[,15], p3_comp$x[,1:7])

## combine all pca components to one dataframe
train.pca <- data.frame(train.glcm_pca_comp, train.lbp_pca_comp[,-1], train.gabor_pca_comp[,-1])

## implement random forest model to the data generated
t.control<- trainControl(method = "repeatedcv", repeats = 2, number = 10)
cart.mod <- train(x = train.pca[,-1], y = train.pca[,1], method = "rpart",cp = 0.1, trControl = t.control)
cart.mod$finalModel

## transform test data to pca

test.glcm_pca_comp <- predict(p1_comp, newdata = glcm.test[,-14])
test.glcm_pca_comp <- as.data.frame(test.glcm_pca_comp)
pred1<- predict(rf.mod, test.glcm_pca_comp)

table(pred1, glcm.test[,14])

trgrid <- expand.grid(cp = c(seq(0.01,0.1,0.01)))
tr <- trainControl(method = "repeatedcv")
mod_cart <- train(x = a[,-38], y = a[,38], method = "rpart", cp = 0.1, trControl = tr, tuneGrid = trgrid)
