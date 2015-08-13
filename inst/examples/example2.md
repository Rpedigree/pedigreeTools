#### Example: Computing inbreeding, additive relationships and functions of it.


```R
pede<-data.frame(sire=as.character(c(NA,NA,NA,NA,NA,1,3,5,6,4,8,1,10,8)),
                  dam= as.character(c(NA,NA,NA,NA,NA,2,2,NA,7,7,NA,9,9,13)),
                  label=as.character(1:14))
        #scrambled original pedigree:
        (pede<- pede[sample(replace=FALSE, 1:14),]  )
        (pede<- editPed(sire=pede$sire, dam= pede$dam, label=pede$label)) 
        ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))

```

[Home](https://github.com/Rpedigree/pedigreeR)
 
