#### Example: Pedigree analyses for self-polinated species.

The function ```getASelfing``` computes an additive relationship matrix (A) from a pedigree. The 
relationship matrix is computed using the pedigree functions in the package and 
R-functions that we have created to accommodate selfing. 

```R

library(pedigreeR)

pedigree_info=system.file("data/sample_pedigree_selfing.csv",package="pedigreeR")

pedigree_data=read.csv(file=pedigree_info,header=TRUE)


nCycles=substr(pedigree_data$Generation,start=2,stop=2)
nCycles=as.integer(nCycles)
ID=pedigree_data$id
Par1=pedigree_data$Par1
Par2=pedigree_data$Par2

A=getASelfing(ID=ID,Par1=Par1,Par2=Par2,nCycles=nCycles,nCyclesDefault=6)


```

[Home](https://github.com/Rpedigree/pedigreeR)
 
