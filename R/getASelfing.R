library(pedigreeR)

### Function to get A using pedigreemm

getASelfing<-function(ID,Par1,Par2,nCycles,sepChar='x-x-x',verbose=FALSE,fileNewPed=NULL,computeA=TRUE)
{
	# completes and sorts pedigree
	TMP<-as.data.frame(editPed(sire=Par1,dam=Par2,label=ID))
	PED<-merge(x=TMP,y=data.frame(ID=ID,nCycles=nCycles,stringsAsFactors=FALSE), all.x=TRUE,by.x=1,by.y=1,sort=FALSE)
	
	#I modified this, so that If we do not know 
	#The number of selfing generations for founders we assume we have F6
	
	index=is.na(PED$nCycles)
	if(length(index)>0) PED$nCycles[index]=6
	

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
			tmp<-grep(PED$sire,pattern=id)
			if(length(tmp)>0){
				PED$sire[tmp]<-newID
			}
				
			tmp<-grep(PED$dam,pattern=id)
			if(length(tmp)>0){
				PED$dam[tmp]<-newID
			}
		
		}
		if(verbose){print(i)}
	}

	newPED<-editPed(sire=newPED$sire,dam=newPED$dam,label=newPED$label)
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
#########

