###################################################################################
### RUN Leave One Out Cross Validation (RUN_LOOCV.R)
### Author and Maintainer: Elisa Warner (elisawa@umich.edu)
### Last Update : 3/3/2023
###
### Description: This code takes in the residuals saved from RUN_RESD.R and performs
###  leave one out cross-validation in a logit model to predict pseudoprogression from
###  MRI files. Note that geometric PCA is also conducted on the residual densities
###################################################################################

library('stringr')
library(pROC)
library(randomForest)
library("caret")

## SET UP
## SET WORKING DIRECTORY
setwd(".../1Pseudoprogression/UPenn-UMich Code Sharing")

folder.name = getwd()

packages = c("WhiteStripe", "oro.nifti", "doParallel", "stringr", "R.utils", "matlabr",
             "R.matlab", "sp", "RNifti")
package.check = lapply(packages,
                       function(x){
                         if(!x %in% rownames(installed.packages())) {
                           install.packages(x, dependencies = TRUE)
                         }
                         library(x, character.only = TRUE)
                       })
rm(package.check,packages)

pairList = list(c("t1","flair"),
                c("t1post","flair"),
                c("t2","flair"),
                c("t1","t1post"),
                c("t1","t2"),
                c("t1","t1postpre"),
                c("t1post","t2"),
                c("t1post","t1postpre"),
                c("t1","t1post"),
                c("t1postpre","flair"),
                c("t1postpre","t2"),
                c("adc","flair"),
                c("adc","t1"),
                c("adc","t1post"),
                c("adc","t1postpre"),
                c("adc","t2"))


for (i in 1:length(pairList)){
  xname = pairList[[i]][1]
  yname = pairList[[i]][2]
  #print(xname)
  #print(yname)
  
  # LOCATION OF RESIDUAL DATA
  setwd(paste0("/.../1Pseudoprogression/UPenn-UMich Code Sharing/Data_Full42_2023/", xname, "_", yname))
  load(paste0("resd_", xname, "_", yname, ".RData"))

  Combined = resd
  CombinedID = names(resd)

  #### Get labels
  temp = c()
  for (i in 1:length(CombinedID)){
    id = CombinedID[i]
    n = nchar(id)
    if (substr(CombinedID[i],n,n) != "2" && substr(CombinedID[i],n-1,n) != "+1"){
      temp = c(temp, id)
    }
  }

  Combined = Combined[temp]
  CombinedID = names(Combined)

  ### CREATE LABELS
  names = CombinedID
  pseudo_labels = c("3-1","6-1","11-1","17-1","21-1","33-1","38-1","7-1","41-1","47-1","48-1","49-1","50-1")
  for (i in 1:length(pseudo_labels)){
    id = pseudo_labels[i]
    pseudo_labels = c(pseudo_labels, paste0(id,"+1"), paste0(id,"-1"))
    #pseudo_labels = c(pseudo_labels, paste0(id,"+2"), paste0(id,"-2"))
  }

  labels = matrix(0, nrow = 1, ncol = length(names))

  for (i in 1:length(pseudo_labels)) {
    if (pseudo_labels[i] %in% names) {
      idx = which(names == pseudo_labels[i])
      labels[idx] = 1
    }
  }

  ### CREATE PATIENT LIST
  unique = str_count(names,"[-+]")
  patList = c()
  for (i in 1:length(names)){
    if (unique[i] == 1){
      patList = c(patList, names[i])
    }
  }

  #Then loop through with density creation+ classifier for LOOCV
  ##what classifier are we going to use? probit maybe?or a simple RF as well as a first pass?
  ### CODE LOCATION
  source('.../1Pseudoprogression/UPenn-UMich Code Sharing/density_featuresOG.R')

  # prediction = numeric()
  prediction = c()
  gt = c()
  for (i in 1:length(patList)){
    pat = patList[i]
  
    # isolate relevant indices
    idx = c()
    for (n in 1:length(CombinedID)){
      if (substr(CombinedID[n], 1, nchar(pat)) == pat){
        idx = c(idx, n)
      }
    }
  
    x.trn = Combined[-idx]
    names_train = CombinedID[-idx]
    labels_train = labels[-idx]
    
    x.tst = Combined[idx] # pat
    names_test = CombinedID[idx] # pat
    labels_test = labels[idx] #[1] # idx
    ftrs_pdacp = density_features(x.trn, x.tst, perc.var = 99.9)
    
    ftrs_train = as.data.frame(ftrs_pdacp$pScore)
    ftrs_train = ftrs_train #[,c(1:3)]
    ftrs_train$V30=as.factor(as.numeric(labels_train)) # this is where you set the label
    # ftrs_test = ftrs_pdacp$pScore.new
    ftrs_test = as.data.frame(ftrs_pdacp$pScore.new)
    ftrs_test = ftrs_test # [,c(1:3)]
  
    #train the classifier using only the xtrain data

    # Perform training: as.factor(V13) ~ .,
    denyprobit <- glm(V30 ~ ., 
                    family = binomial(link = "logit"), 
                    data = ftrs_train)
    confint(denyprobit)
  
    #Now try to use the trained classifier to classify xtest
    predictions <- predict.glm(denyprobit, 
                             newdata = as.data.frame(ftrs_test),
                             type = "response")
    #prediction[i]=predictions
    prediction = c(prediction, predictions)
    gt = c(gt, labels_test)
  }

  ### NOTE: This commented out code was from the original code. It works
  ## but only in cases where you get a binary output. So I scratch it for now
  ## because roc() function enforces binary output.
  # auc(cases=round(prediction[which(gt==1)]),controls=round(prediction[which(gt==0)]))
    
  ### SAVE LOG FILE
  sink(file = paste0(".../",xname,"_",yname,"_results.txt"))
  myroc <- roc(gt, round(prediction), direction="<")
  print("AUC:")
  x <- auc(myroc)
  print(x)
  print("Sensitivity:")
  print(sensitivity(as.factor(round(prediction)), as.factor(gt)))
  print("Specificity:")
  print(specificity(as.factor(round(prediction)), as.factor(gt)))
  print("Predictions:")
  print(unname(prediction))
  print("Ground Truth")
  print(gt)
  sink()
  rm(labels_train, labels_test, ftrs_pdacp, ftrs_train, ftrs_test)
}
