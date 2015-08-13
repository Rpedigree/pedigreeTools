#### Example: Fiting pedigree mixed effects models using pedigreeR and lmer.

```R
library(pedigreeR)

#####################################################################
#Reading the information
#####################################################################

ped_file = system.file("data/mice_pedigree.txt",package="pedigreeR")
obesity_file=system.file("data/Obesity.txt", package="pedigreeR")
cage_file=system.file("data/Cage.txt",package="pedigreeR")

Obesity=read.table(file="Obesity.txt",header=TRUE)
cage=read.table(file="Cage.txt",header=TRUE)

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

Obesity=read.table(file="Obesity.txt",header=TRUE)
cage=read.table(file="Cage.txt",header=TRUE)
out=merge(x=Obesity,y=cage,by.x="SUBJECT.NAME",by.y="IID")

#Remove records with missing values
index=!is.na(out$Obesity.BMI) & !is.na(out$cage) & !is.na(out$Date.Year) & !is.na(out$Date.Season) & !is.na(out$Date.Month) & !is.na(out$EndNormalBW)
pheno=out[index,] 

#Remove individuals without pheno info
index=as.character(pheno$SUBJECT.NAME)%in%as.character(rownames(A))
pheno=pheno[index,]

```
[Home](https://github.com/Rpedigree/pedigreeR)
 
