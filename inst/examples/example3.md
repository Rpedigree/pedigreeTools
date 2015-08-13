#### Example: Fiting pedigree mixed effects models using pedigreeR and lmer.

```R
library(pedigreemm)
geno=read.table(file="Markers-chr1-chr19-x.mrk",header=TRUE,
                        na.strings="-9",check.names=FALSE)
                        
                        
pat=as.character(geno$PAT)
mat=as.character(geno$MAT)
id=as.character(geno$IID)

#Complete the pedigree
tmp=unique(c(as.character(geno$PAT),
             as.character(geno$MAT)))
            
pat=c(rep(NA,length(tmp)),pat)
mat=c(rep(NA,length(tmp)),mat)
id=c(tmp,id)

tmp=pedigree(pat,mat,id)

#A=crossprod(as.matrix(relfactor(tmp,labs=id[345:2284])))
#rownames(A)=colnames(A)=id[345:2284]
#A=crossprod(as.matrix(relfactor(tmp)))
#rownames(A)=colnames(A)=id

#Using the new function, version 0.2.5
A=as.matrix(getA(tmp))
rownames(A)=colnames(A)=id

index=colnames(A)%in%as.character(rownames(Markers))
A=A[index,index]


```

[Home](https://github.com/Rpedigree/pedigreeR)
 
