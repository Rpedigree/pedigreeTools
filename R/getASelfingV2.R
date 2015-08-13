#R Code

#Load routines for pedigree
source("~/Desktop/Bayer/getASelfing.R")



getASelfing=function(ID,Par1,Par2,nCycles,nCyclesDefault, sepChar='-F', verbose=FALSE, fileNewPed=NULL, computeA=TRUE){
	library(pedigreemm)

## Completing and sorting pedigree
TMP<-as.data.frame(editPed(sire=Par1,dam=Par2,label=ID))
 PED<-merge(x=TMP,y=data.frame(ID=ID,nCycles=nCycles,stringsAsFactors=FALSE), all.x=TRUE,by.x=1,by.y=1,sort=FALSE)
 PED=PED[order(PED$gene),]
 
## Adding number of selfing cycles for lines that miss that information 
 index=is.na(PED$nCycles)
 if(length(index)>0){ PED$nCycles[index]=nCyclesDefault}

## Extending the pedigree by adding selfing
  newPED<-PED[integer(),]

  row<-0

  for(i in 1:nrow(PED)){

		id<-PED$label[i]
		
		generation<-PED$gene[i]
		cycles<-PED$nCycles[i]
		row<-row+1	
		newPED[row,1]<-paste(id,0,sep=sepChar)
		newPED[row,2:3]<-PED[i,2:3]
		newPED[row,4]<-paste(PED[i,4],0,sep=sepChar)
	
		if(cycles>1){
			for(j in 1:cycles){
				row<-row+1
				newPED[row,1]<-paste(id,j,sep=sepChar)
				newPED[row,2:3]<-rep(newPED[(row-1),1],2)
				newPED[row,4]<-paste(generation,j,sep=sepChar)
			}
		
			# now we update the id
			newID<-paste(id,cycles,sep=sepChar)
			#tmp<-grep(PED$sire,pattern=id)#*#
			tmp<-which(PED$sire==id)
			if(length(tmp)>0){
				PED$sire[tmp]<-newID
			}
				
			#tmp<-grep(PED$dam,pattern=id)
			tmp<-which(PED$dam==id)
			if(length(tmp)>0){
				PED$dam[tmp]<-newID
			}
		
		}
		if(verbose){print(i)}
	}
		
	if(!is.null(fileNewPed)){
		write.table(newPED,file=fileNewPed,quote=FALSE,row.names=FALSE,col.names=TRUE,sep=",")
	}
	
	if(computeA){
		newPED<-pedigree(sire=newPED$sire,dam=newPED$dam,label=newPED$label)
		A<-getA(newPED)
		tmp<-rownames(A)
		tmp2<-paste(ID,nCycles,sep=sepChar)
		tmpIndex<-which(tmp%in%tmp2)
		A<-A[tmpIndex,tmpIndex]
		TMP<-matrix(unlist(strsplit(rownames(A),split=sepChar)),ncol=2,byrow=TRUE)[,1]
		rownames(A)<-TMP
		colnames(A)<-TMP
		return(A)
	}
}   
   



pedigree_data=read.csv(file="~/Desktop/Bayer/MCP2009_pedigree.csv",header=TRUE)


nCycles=substr(pedigree_data$Generation,start=2,stop=2)
nCycles=as.integer(nCycles)
ID=pedigree_data$id
Par1=pedigree_data$Par1
Par2=pedigree_data$Par2
A=getASelfing(ID=ID,Par1=Par1,Par2=Par2,nCycles=nCycles,nCyclesDefault=6)
 
sum(cov2cor(A)[row(A)>col(A)]>.5)
	
A[c('57','62'),c('57','62')]



plot(diag(A),ylim=c(0,2))
 image(A)
 A<-as.matrix(A)

 SVD<-svd(A,nu=2)
 eigenValues<-SVD$d^2
 plot(eigenValues,ylab='Eeigen-value (A)',col=4,cex=.7)
 
 plot(cumsum(eigenValues)/sum(eigenValues),ylab='Proportion of variance',col=4,cex=.7,type='o');abline(h=1,lty=2)
 
 plot(SVD$u,col=2,xlab='1st PC',ylab='2nd PC')
