#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include "HeaderFiles/reqfuncAll.h"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define F_NAME_MAX 1000
#define SAMPLE_NAME_MAX 100
#define MAX_DIM 2000
#define MAX_RANK 11
struct EigenSpace
{
	double **U;
	double **S;
	int Rank;
};
int main(int argc, char** argv) 
{
	time_t t;
	srand((unsigned) time(&t));
	int m,i,j,k,r;
	int M, K, n, *Size=NULL, *Labels=NULL, *LogT=NULL;
	char **DataFnames=NULL, *LabelFile=NULL, **Samples=NULL;
	double ***XM=NULL, FrbP;
	char *File=argv[1];	
	printf("\n %-20s %s","Input File Name=",File);
	getKMn(&K, &M, &n, File);
	printf("\n %-20s %d","No. of Samples=",n);
	printf("\n %-20s %d","No. of Modalities=",M);
	printf("\n %-20s %d","No. of Clusters=",K);
	Size=dynint1d(Size,M);
	LogT=dynint1d(LogT,M);
	LabelFile=dynchar1d(LabelFile,F_NAME_MAX);
	DataFnames=dymemchar2d(DataFnames, M, F_NAME_MAX);
	Samples=dymemchar2d(Samples, n, SAMPLE_NAME_MAX);
	Labels=dynint1d(Labels,n);
	getfiles(File, M, n, DataFnames, Size, LogT, LabelFile, Samples, Labels);
	for(m=0;m<M;m++)
		printf("\n Modality %d File Name: %-80s #Features: %d  Log Transform: %d",m+1,DataFnames[m],Size[m],LogT[m]);
	printf("\n %-20s %s","Label File Name=",LabelFile);
	print2dchar(Samples, 5, "First 5 Samples");
	XM=dymemdouble3d(XM,M,n,MAX_DIM);
	getdata(XM, Size, DataFnames, n, M);	
	struct EigenSpace *ES=NULL;
	ES=(struct EigenSpace*)malloc(M*sizeof(struct EigenSpace));
	double **Uxi=NULL, alpha=0.05;
	double **Xm, **Xmcent, **U, **S, **Vt, **Umax, **Smax, **Um, **Sm, **UmStats=NULL;
	double *relv=NULL, *depend=NULL, *Phi=NULL, *Theta=NULL, Thetam, ThetaSum=0.0;
	UmStats=dymemdouble2d(UmStats, MAX_RANK, 5);
	int *order=NULL, *modmarker=NULL, *taken=NULL, modcount=0, col, Rank, iter=100, ntimes=20;
	relv=dyndouble1d(relv, M);
	depend=dyndouble1d(depend, M);
	Phi=dyndouble1d(Phi, M);
	Theta=dyndouble1d(Theta, M);
	order=dynint1d(order, M);
	modmarker=dynint1d(modmarker, M);
	taken=dynint1d(taken, M);
	for(m=0;m<M;m++)
	{
		Xm=NULL; U=NULL; S=NULL; Vt=NULL; Umax=NULL; Smax=NULL; Um=NULL; Sm=NULL; Xmcent=NULL;
		col=Size[m];
		Xm = dymemdouble2d(Xm, n, col);
		Xmcent = dymemdouble2d(Xmcent, n, col);
		U = dymemdouble2d(U, n, n);
		S = dymemdouble2d(S, n, col);
		Vt = dymemdouble2d(Vt, col, col);
		getXi(XM, Xm, m, n, col);
		if(LogT[m]==1)
			logtransform(Xm,n,col);
    		print2ddouble(Xm, 5, 5, "Xm");
		meancenter(Xm, n, col, Xmcent);
		svdlapack(Xmcent, n, col, U, S, Vt);
		Umax = dymemdouble2d(Umax, n, MAX_RANK);
		Smax = dymemdouble2d(Smax, MAX_RANK, MAX_RANK);
		copymat2ddouble(U, Umax, n, MAX_RANK);
		copymat2ddouble(S, Smax, MAX_RANK, MAX_RANK);
		freedouble2d(U, n, n);
		freedouble2d(S, n, col);
		freedouble2d(Vt, col, col);
		writeDouble2DMatrix(Umax,"Um",n,MAX_RANK);
		system("Rscript Royston.R");
		remove("Um");
		readDouble2DMatrix(UmStats,"UmStats",MAX_RANK,5);
		remove("UmStats");
		printf("\n\nModalities %d",m+1);
		for(i=0; i<MAX_RANK-1; i++)
			printf("\n Component %d H-stat= %3.5lf P-value= %1.10lf",i+1,UmStats[i][0],UmStats[i][1]);
		Rank=0;
		for(i=0; i<MAX_RANK-1; i++)
		{
			if(UmStats[i][1]>alpha && UmStats[i+1][1]>alpha)
				break;
			else
				Rank++;
		}
		printf("\n Rank of Modality %d= %d",m+1,Rank);	
		Um = dymemdouble2d(Um, n, Rank);
		Sm = dymemdouble2d(Sm, Rank, Rank);
		copymat2ddouble(Umax, Um, n, Rank);
		copymat2ddouble(Smax, Sm, Rank, Rank);
		freedouble2d(Umax, n, MAX_RANK);
		freedouble2d(Smax, MAX_RANK, MAX_RANK);		
		Phi[m]=0.5*(1+(UmStats[Rank][2]-UmStats[Rank][4])/MAX(UmStats[Rank][2],UmStats[Rank][4]));
		Thetam=0.0;
		for(r=0;r<Rank;r++)
			Thetam=Thetam+Sm[r][r];
		ThetaSum=ThetaSum+Thetam;
		Theta[m]=Thetam;
		(ES+m)->U=Um;
		(ES+m)->S=Sm;
		(ES+m)->Rank=Rank;
		if(Rank>0)
			modmarker[m]=1;
	}	
	freedouble2d(UmStats, MAX_RANK, 5);
	for(m=0;m<M;m++)
	{
		Theta[m]=Theta[m]/ThetaSum;
		relv[m]=Phi[m]*Theta[m];
		order[m]=m;
		printf("\n Modality %d Relevance= %lf Rank= %d",m+1,relv[m],(ES+m)->Rank);
	}
	isortwindex(relv, order, M, 1);
	struct EigenSpace *jointES=NULL;
	jointES=(struct EigenSpace *)malloc(sizeof(struct EigenSpace));
	double **jointU=NULL, **jointS=NULL, **P=NULL, **R=NULL, **jointUt=NULL, **ProjMat=NULL, *UselectStats=NULL, **takeU=NULL, **takeS=NULL, **tempU=NULL, **tempS=NULL;
	ProjMat=dymemdouble2d(ProjMat, n, n);
	int jointRank=(ES+order[0])->Rank, select, flag, countRes, r1,r2,  *takeRes=NULL;
	takeRes=dynint1d(takeRes, MAX_RANK);
	jointU=dymemdouble2d(jointU, n, jointRank);
	jointS=dymemdouble2d(jointS,  jointRank,  jointRank);
	copymat2ddouble((ES+order[0])->U, jointU, n, jointRank);
	copymat2ddouble((ES+order[0])->S, jointS, jointRank, jointRank);
	jointES->U=jointU;
	jointES->S=jointS;
	jointES->Rank=jointRank;
	print1dint(order, M, "ORDER");
	print1dint(modmarker, M, "modmarker");
	modmarker[order[0]]=0;
	taken[modcount]=order[0]+1;
	print1dint(modmarker, M, "modmarker");
	print1dint(taken, M, "taken");
	for(i=1;i<M;i++)
	{
	    jointRank=jointES->Rank;
	    constantarraydouble(depend, M, 0);
	    jointUt=dymemdouble2d(NULL, jointRank, n);
	    transpose(jointES->U, jointUt, n, jointRank);
	    matmul(jointES->U, jointUt, ProjMat, n, jointRank, jointRank, n);
	    freedouble2d(jointUt, jointRank, n);
	    flag=0;
	    for(j=0;j<M;j++)
	    {
	        if(modmarker[j]==1)
	        {
	           flag=1;
	           P=dymemdouble2d(NULL, n, (ES+j)->Rank);
	           matmul(ProjMat, (ES+j)->U, P, n, n, n, (ES+j)->Rank);
	           FrbP=forbnorm(P,n,(ES+j)->Rank);
	           depend[j]=FrbP/(ES+j)->Rank;
	           freedouble2d(P, n, (ES+j)->Rank);
	        }
	    }
	    print1ddouble(depend, M, "depend");
	    if(flag==1)
	    {
	        select=whichMax(depend, M);  
	        printf("Selected modality %d",select+1);
	        modmarker[select]=0;
            P=dymemdouble2d(NULL, n, (ES+select)->Rank);
            R=dymemdouble2d(NULL, n, (ES+select)->Rank);
            matmul(ProjMat, (ES+select)->U, P, n, n, n, (ES+select)->Rank);
            matdiff((ES+select)->U,P,R,n,(ES+select)->Rank);
            print2ddouble(R,5,(ES+select)->Rank,"Uselect");
    		writeDouble2DMatrix(R,"Uselect",n,(ES+select)->Rank);
    		freedouble2d(P, n, (ES+select)->Rank);
            freedouble2d(R, n, (ES+select)->Rank);
    		system("Rscript ShapiroWilk.R");
		    remove("Uselect");
		    UselectStats=dyndouble1d(NULL, (ES+select)->Rank);
		    readDouble1DArray(UselectStats,"UselectStats",(ES+select)->Rank);
		    remove("UselectStats");
		    countRes=0;
		    constantarrayint(takeRes, (ES+select)->Rank, 0);
		    for(j=0;j<(ES+select)->Rank;j++)
		    {
		        if(UselectStats[j]<alpha)
		        {
		            countRes++;
		            takeRes[j]=1;
		        }
		    }
		    free(UselectStats);
		    if(countRes>0)
		    {
		        takeU=dymemdouble2d(NULL, n, countRes);
		        takeS=dymemdouble2d(NULL, countRes, countRes);
		        r=-1;
		        for(j=0;j<(ES+select)->Rank;j++)
		        {
		            if(takeRes[j]==1)
		            {
		                printf("Taking Componnet %d",j);
		                r=r+1;
		                for(k=0;k<n;k++)
		                    takeU[k][r]=(ES+select)->U[k][j];
		                takeS[r][r]=(ES+select)->S[j][j];
		            }
		        }
			    r1=jointES->Rank;
			    r2=countRes;
			    tempU=jointES->U;
			    tempS=jointES->S;
			    jointES->Rank=r1+r2;
			    printf("\n No. Of Components to Take from Modality %d= %d",select+1,countRes);
			    printf("\n Joint rank changed to: %d",jointES->Rank);
			    jointES->U=dymemdouble2d(NULL, n, jointES->Rank);
			    jointES->S=dymemdouble2d(NULL, jointES->Rank, jointES->Rank);
			    appendU(jointES->U, tempU, takeU, n, r1, r2);
			    appendS(jointES->S, tempS, takeS, r1, r2);
			    modcount=modcount+1;
			    taken[modcount]=select+1;
			    freedouble2d(tempU, n, r1);
			    freedouble2d(tempS, r1, r1);
			    freedouble2d(takeU, n, r2);
			    freedouble2d(takeS, r2, r2);
		    }
		    else
			printf("\n Modality %d not updated as Residuals are all Normal", select);
	    }
	}
	double **FinalPC=NULL;
	int *Finalclust=NULL;
	FinalPC=dymemdouble2d(FinalPC, n, jointES->Rank);
	Finalclust=dynint1d(Finalclust, n);
	matmul(jointES->U, jointES->S, FinalPC, n, jointES->Rank, jointES->Rank, jointES->Rank);
	print2ddouble(FinalPC, 5, jointES->Rank, "Final Joint Subspace");
	writeDouble2DMatrix(FinalPC,"FinalJointSubspace.txt",n,jointES->Rank);
	double Rand, ARI, Dice, Jac, Purity, Fmeasure, NMI, withinss;
	HartiganKm(FinalPC, n, jointES->Rank, K, Finalclust, NULL, ntimes, iter, &withinss);
	confusion(Labels, Finalclust, n, &Rand, &ARI, &Dice, &Jac, &Purity, &Fmeasure, &NMI);
	printf("\n Fmeasure %lf, Rand %lf, Jac %lf, Dice %lf,  Purity %lf",Fmeasure, Rand, Jac, Dice, Purity);
	print1dint(modmarker, M, "ModalityMarker");
	print1dint(taken, M, "Taken");
	freedouble2d(ProjMat, n, n);
	free(relv);
	free(depend);
	free(taken);
	free(order);
	free(modmarker);
	free(Finalclust);
	freedouble2d(FinalPC, n, jointES->Rank);
	free(Size);
	free(LogT);
    printf("\n");
	return 0;
}
