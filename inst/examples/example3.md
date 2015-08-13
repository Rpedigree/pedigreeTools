#### Example: Fiting pedigree mixed effects models using pedigreeR and lmer.

```R
library(pedigreeR)

ped_file = system.file("data/mice_pedigree.txt")
ped_info = read.table(file=ped_file,header=TRUE,
                      na.strings="-9",check.names=FALSE)
                        
pat=as.character(ped_info$PAT)
mat=as.character(ped_info$MAT)
id=as.character(ped_info$IID)

#Complete the pedigree
tmp=unique(c(as.character(geno$PAT),
             as.character(geno$MAT)))
            
pat=c(rep(NA,length(tmp)),pat)
mat=c(rep(NA,length(tmp)),mat)
id=c(tmp,id)

tmp=pedigree(pat,mat,id)

#Using the new function, version 0.2.5
A=as.matrix(getA(tmp))
rownames(A)=colnames(A)=id

#index=colnames(A)%in%as.character(rownames(Markers))
#A=A[index,index]


```

[Home](https://github.com/Rpedigree/pedigreeR)
 
