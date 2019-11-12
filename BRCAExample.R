#Data set information
DataSet="BRCA"
n=398
K=4
#Extract sample names
truedata=read.table(paste0("DataSets/",DataSet,"/TCGA Subtype Labels",n),stringsAsFactors=FALSE)
true=as.vector(truedata[,2],mode="numeric")
samples=as.vector(truedata[,1])
#Read multimodal data into a list
Data<-list()
Data[[1]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/mDNA",n),sep=" ",header=TRUE,row.names=1))
Data[[2]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/RNA",n),sep=" ",header=TRUE,row.names=1))
Data[[3]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/miRNA",n),sep=" ",header=TRUE,row.names=1))
Data[[4]] <- as.matrix(read.table(paste0("DataSets/",DataSet,"/RPPA",n),sep=" ",header=TRUE,row.names=1))
M=length(Data)
mod=c("mDNA","RNA","miRNA","RPPA")
#Log Transform of Sequence based Gene and miRNA modality
LogData=Data
LogData[[2]][LogData[[2]]==0]=1
LogData[[2]]=log(LogData[[2]],base=10)
LogData[[3]][LogData[[3]]==0]=1
LogData[[3]]=log(LogData[[3]],base=10)

#****************************************************************************** End of Data Import ***********************************************************************


cat("\nDataset=",DataSet,append=F)
cat("\n#Samples=",n,append=T)
cat("\nClusters=",K,append=T)

source("Normality.R")
out=Normality(LogData,mod)
Dmat=out$Dmat	
cat("\n\nFinal Joint Subspace written to file: JointSubspace.txt\n\n")
write.table(Dmat,col.names=FALSE,row.names=samples,quote=FALSE,file="JointSubspace.txt")
 
 
#Perform K-means clustering on joint subspace
JointSub=as.matrix(read.table("JointSubspace.txt",sep=" ",row.names=1,header=FALSE))
print("\n First few rows of joint subspace:\n")
print(JointSub[1:5,])
cat("\n Joint Subspace Dimension: ",dim(JointSub)[1]," rows",dim(JointSub)[2]," columns")
km=kmeans(JointSub,K)$cluster
df=data.frame(cbind(samples,km))
write.table(df,quote=FALSE,col.names=FALSE,row.names=FALSE,file=paste0(DataSet,"-ClusterAssignment.txt"))
cat("\n\nFinal cluster assignments written to file:",paste0(DataSet,"-ClusterAssignment.txt\n\n"))
