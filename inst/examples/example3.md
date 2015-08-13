#### Example: Fiting pedigree mixed effects models using pedigreeR and lmer.

```R
library(pedigreeR)

#####################################################################
#Reading the information
#####################################################################

ped_file = system.file("data/mice_pedigree.txt",package="pedigreeR")
obesity_file=system.file("data/Obesity.txt", package="pedigreeR")
cage_file=system.file("data/Cage.txt",package="pedigreeR")

obesity=read.table(file=obesity_file,header=TRUE)
cage=read.table(file=cage_file,header=TRUE)

#pedigree info
ped_info = read.table(file=ped_file,header=TRUE,
                      na.strings="-9",check.names=FALSE)

#phenotypes info
obesity=read.table(file=obesity_file, header=TRUE)

#cage info
cage=read.table(file=cage_file,header=TRUE)

#Processing the pedigree
pat=as.character(ped_info$PAT)
mat=as.character(ped_info$MAT)
id=as.character(ped_info$IID)

#Complete the pedigree
tmp=unique(c(as.character(ped_info$PAT),
             as.character(ped_info$MAT)))
            
pat=c(rep(NA,length(tmp)),pat)
mat=c(rep(NA,length(tmp)),mat)
id=c(tmp,id)

tmp=pedigree(pat,mat,id)

#Using the new function, version 0.2.5
A=as.matrix(getA(tmp))
rownames(A)=colnames(A)=id

#####################################################################
#Phenotypes
#####################################################################
out=merge(x=obesity,y=cage,by.x="SUBJECT.NAME",by.y="IID")

#Remove records with missing 
index=rowSums(is.na(out))<1
pheno=out[index,] 

#Common individuals with phenotypes and pedigree info
common=intersect(as.character(pheno$SUBJECT.NAME),rownames(A))

index=as.character(pheno$SUBJECT.NAME)%in%common
pheno=pheno[index,]

index=rownames(A)%in%common

A=A[index,index]

#Sort the pedigree information and pheno information so that they match
index=order(as.character(pheno$SUBJECT.NAME))
pheno=pheno[index,]

index=order(rownames(A))
A=A[index,index]

#Check if every thing matches
if(any(colnames(A)!=as.character(pheno$SUBJECT.NAME))) stop("Ordering problem\n")

#Up to here we have phenotypic and genotypic information

```
[Home](https://github.com/Rpedigree/pedigreeR)
 
