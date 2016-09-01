### Example 2: Computing inbreeding, additive genetic relationship matrix A and its inverse.

Following the last example, we create a valid pediree object from a ```data.frame``` representation of the example pedigree:
```R
library(pedigreeR)

pedFrame <-data.frame(sire=as.character(c(NA,NA,NA,NA,NA,1,3,5,6,4,8,1,10,8)),
dam= as.character(c(NA,NA,NA,NA,NA,2,2,NA,7,7,NA,9,9,13)),
label=as.character(1:14))

ped <- with(pedFrame, pedigree(label=label, sire=sire, dam=dam))

```
#### Section A : Obtaining A, Ainv and inbreeding coefficients 

It must be noted here that it is computationally cheaper to set up the inverse of the A matrix than to set up A and then invert it. Hence we compute A or Ainv from the functions ```getA``` or ```getAInv``` as needed. 

Compute the inbreeding coefficients of the pedigree
```R
inb <- inbreeding(ped)
```
Get the additive genetic relationship matrix
```R
 A <- getA(ped)
```
Get the inverse of A
 ```R
 Ainv <- getAInv(ped)
```

Below we show a visualization for the Ainv matrix for the example pedigree:

<img src="https://github.com/siddharth51292/pedigreeR/blob/patch-2/inst/examples/Ainv.png" width="500">

[Home](https://github.com/Rpedigree/pedigreeR)
 
