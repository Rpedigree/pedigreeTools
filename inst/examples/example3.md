### Example 3: Subsetting a Pedigree using ```prunePed``` 

#### Section A : Creating the example pedigree and Scrambling it.

Following the previous examples, we start with a scarmbled version of the example data frame pedigree ```pedFrame```

```R
#Creating the example pedigree
pedFrame<-data.frame(sire=as.character(c(NA,NA,NA,NA,NA,1,3,5,6,4,8,1,10,8)),
                  dam= as.character(c(NA,NA,NA,NA,NA,2,2,NA,7,7,NA,9,9,13)),
                  label=as.character(1:14))
#Scrambling the example pedigree:
pedScram<- pedFrame[sample(replace=FALSE, 1:14),] 

```

#### Section B: Using ```prunePed``` to subset the pedigree.

The ```prunePED``` function works by setting a baseline population through a vector of ids, and the number of previous generations ( sires , dams, grandsires and granddams etc. of the baseline population )  the user wishes to select. Here the number of previous generations to select is set through the ```ngen``` argument of prunePED. This selection ends where the parents are unknown regardless of ```ngen```. 

Here we start with an example baseline population vector (```selectVector```) : ```c(12,9,11)``` 
And an ```ngen``` value ```2```

It can be readily seen from the following pictoral representation of the pedigree that the ids that should be returned along with ```c(12,9,11)``` are ```c(1,2,3,5,6,7,8)```

```R
# Applying prunePED to subset the scrambled pedigree 
library(pedigreeR)
pedSelect <- prunePED(pedScram,selectVector=c(12,9,11),ngen=2) 

pedSelect
    sire dam label
 1:   NA  NA     1
 2:    6   7     9
 3:    3   2     7
 4:    1   2     6
 5:    5  NA     8
 6:   NA  NA     3
 7:   NA  NA     2
 8:   NA  NA     5
 9:    1   9    12
10:    8  NA    11

```
#### Section C: Producing a valid pedigree object from the output of ```prunePED```

Following Example 1, we create a sorted and complete pedigree object ```pedFinal``` that can be used for fitting models from the output of prunePED.

```R
pedEdited <- editPed(sire=pedSelect$sire,dam=pedSelect$dam,label=pedSelect$label)
pedFinal <- with(pedEdited, pedigree(label=label,sire=sire,dam=dam))

 pedFinal
   sire  dam
1  <NA> <NA>
3  <NA> <NA>
2  <NA> <NA>
5  <NA> <NA>
7     3    2
6     1    2
8     5 <NA>
9     6    7
11    8 <NA>
12    1    9
```

[Home](https://github.com/Rpedigree/pedigreeR)
 
