### Example 5: Cross Validation using Pedigree information and mixed effects model.

In this example we will be using cross validation to illustrate the use of pedigree information in fitting a mixed effects model. 

#### Section A : Import the data and create pedigree objects 
First we read in the data and remove redundant rows so we can convert the object into a pedigree class.

```R
wheat <- read.table("QuantGen/Data/599_yield_raw-1.prn")
wheat_ped <- read.csv("QuantGen/Data/WHEAT599PROGENIE.csv")
wheat_ped <- unique(wheat_ped)
```

Renaming the columns, deleting the original names and creating and rescaling the dependent variable
```R
colnames(wheat) <- c("obs","env","rep","id","gen1","GY")
wheat <- wheat[-1,]
wheat$GY <- as.numeric(as.character(wheat$GY))
wheat$sdGY <- wheat$GY/sd(wheat$GY)
```

Creating the final pedigree objects 

```R
library(pedigreeR)
wheat_ped_edit <- editPed(sire=wheat_ped$gpid1,dam=wheat_ped$gpid2,label=wheat_ped$progenie)
wheat_ped_final <- with(wheat_ped_edit,pedigree(label=label,sire=sire,dam=dam))
```

#### Section B: Split data into training and testing lots

Before splitting the dataset, we truncate it to make it divisible by 5
```R
wheat <- wheat[1:(dim(wheat)[1]-2),]
```

Getting the lot indices for training and testing lots. 
```R
s <- list()
wheatls <- list()

lotsize <- trunc(nrow(wheat)/5)
indices <- 1:nrow(wheat)
for (i in 1:5){
  s[[i]] <- sample(indices,lotsize)
  indices <- indices[!(indices %in% s[[i]])]
  wheatls[[i]] <- wheat[s[[i]],]
}
```

#### Section C: Getting Yhat from testing sets and computing the overall correlation between Yhat and Y. 
Fit data excluding ith lot ( Training sets )  

```R
fm <- list()
wheat_exclude <- list()
for (j in 1:5){
  i <- c(1,2,3,4,5)
  isel <- i[!(i %in% j)]  
  wheat_exclude[[j]] <- do.call("rbind",wheatls[isel])
  fm[[j]] <- pedigreemm(sdGY~env+(1|gen1),data=wheat_exclude[[j]],pedigree=list(gen1=wheat_ped_final))
}
```
Get Yhat on ith lot (using BetaHat computed on fit excluding ith lot) and put them together. ( Testing sets )

```R
predict <- vector()
for (i in 1:5){
  
  predict <- append(predict,predict(fm[[i]],newdata=wheatls[[i]]))
  
}
```
Sort the results and compute the correlation with original dataset

```R
predict_sorted <- predict[order(as.integer(names(predict)))]
result <- cor(predict_sorted,wheat$sdGY) 
```

#### OUTPUT

```R
> result
[1] 0.6999063 
```
[Home](https://github.com/Rpedigree/pedigreeR)
 
