# Article:
A. Khan and P. Maji, "Low-Rank Joint Subspace Construction for Cancer Subtype Discovery," in *IEEE/ACM Transactions on Computational Biology and Bioinformatics*.

doi: 10.1109/TCBB.2019.2894635

URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8624406&isnumber=4359833


Two R pacakges are required for execution of the method: "nortest" and "e1071".
These packages can be installed by execution the following commands within the R environment in command line:

```r
install.pacakges("nortest")
install.pacakges("e1071")
```

For execution of the C program for the algorithm 
Input File format is as follows:
```c
Samples=  #the number of samples
Modalities=  #the number of modalities
Cluster= #the number of clusters
File1= #Filename for Modality1    Features= #Number of features in Modality1    logtransform= 1 for log transformation 0 for if not
File2= #Filename for Modality2    Features= #Number of features in Modality2    logtransform= 1 for log transformation 0 for if not
File3= #Filename for Modality3    Features= #Number of features in Modality3    logtransform= 1 for log transformation 0 for if not
File4= #Filename for Modality4    Features= #Number of features in Modality4    logtransform= 1 for log transformation 0 for if not
```
###### Example Input File (for BRCA dataset):
```c
Samples= 398
Modalities= 4
Cluster= 4
File1= DataSets/BRCA/mDNA		Features= 2000		logtransform= 0
File2= DataSets/BRCA/RNA		Features= 2000		logtransform= 1
File3= DataSets/BRCA/miRNA		Features= 278		logtransform= 1
File4= DataSets/BRCA/RPPA		Features= 212		logtransform= 0
```
After installing the R packages, execute C version of code by 
```shell
gcc Normality.c -lm -llapack
./a.out BRCA			#Input File name as command line argument
```

R code for the breast cancer (BRCA) data set is given in ``BRCAExample.R`` file.
To execute R code go to R environment and execute in command line:

>source("BRCAExample.R")

The R code for the cervical cancer (CESC) dataset is given in ``CESCExample.R`` file. To run this file within the R terminal execute
>source("CESCExample.R")

Joint subspace is written to file : ``JointSubspace.txt``
JointSubspace.txt contains a ``(n x r)`` where ``n`` is the number of samples in the data set and ``r`` is  the rank of the joint subspace.

The cluster assignments are written to the file ``BRCA-ClusterAssignment.txt`` for the BRCA data set and in file ``CESC-ClusterAssignment.txt`` for the CESC data set.

The file ``Normality.R`` contains the R implementation of the proposed method as a function `Normality`. 
Details of the fuctions is as follows:

Proposed Method Function Name: `Normality`

###### Usage 
`Normality(DataL, mod)`

Arguments

``DataL``:  A list object containing ``M`` data matrices representing ``M`` different omic data types measured in a set of ``n`` samples. For each matrix, the rows represent samples, and the columns represent genomic features.

``mod``: A string array of names of the modalities. Required for modality selection.
Example: ``mod=c("mDNA","RNA","miRNA","RPPA")``

Example call:

```r
Data<-list()
Data[[1]] <- as.matrix(read.table(paste0("DataSets/BRCA/mDNA",n),sep=" ",header=TRUE,row.names=1))
Data[[2]] <- as.matrix(read.table(paste0("DataSets/BRCA/RNA",n),sep=" ",header=TRUE,row.names=1))
Data[[3]] <- as.matrix(read.table(paste0("DataSets/BRCA/miRNA",n),sep=" ",header=TRUE,row.names=1))
Data[[4]] <- as.matrix(read.table(paste0("DataSets/BRCA/RPPA",n),sep=" ",header=TRUE,row.names=1))
mod=c("mDNA","RNA","miRNA","RPPA")
source("Normality.R")
out=Normality(Data,mod)
```


##### Contact Information

Aparajita Khan   
Senior Research Fellow   
Machine Intelligence Unit    
Indian Statistical Institute    
203, B. T. Road    
Kolkata- 700108    
E-mail: aparajitak_r@isical.ac.in,  aparajitakhan1107@gmail.com
