#define CHAR_MAX 1000
double eddist(double *x, double *y, int M)
{
	double d=0.0;
	int i=0;
	for(i=0;i<M;i++)
	{
		d=d+(x[i]-y[i])*(x[i]-y[i]);
	}
	return(sqrt(d));
}
void isortwindex(double *val, int *index, int n, int flag)
{
		int i,j;
		double keyv;
		int keyin;
		for(j=1;j<n;j++)
		{
			keyv=val[j];
			keyin=index[j];
			i=j-1;
			if(flag==0)
			{
				while(i>=0 && keyv<val[i])
				{
					val[i+1]=val[i];
					index[i+1]=index[i];
					i=i-1;
				}
			}
			else
			{
				while(i>=0 && keyv>val[i])
				{
					val[i+1]=val[i];
					index[i+1]=index[i];
					i=i-1;
				}
			}
			val[i+1]=keyv;
			index[i+1]=keyin;
		}
}
int maxclust(int *arr, int n)
{
	int mx=0;
	int i;
	mx=arr[0];
	for(i=1;i<n;i++)
	{
		if(mx<arr[i])
			mx=arr[i];
	}
	return(mx);
}

double randind(int *org, int * clust, int N)
{
	double a,b,den;
	int i,j;
	double rind=0.0;
	den=(N*(N-1))/2;
	a=0;
	b=0;
	for(i=0;i<N;i++)
	{
		for(j=i+1;j<N;j++)
		{
			if(org[i]==org[j] && clust[i]==clust[j])
				a=a+1;
			else if(org[i]!=org[j] && clust[i]!=clust[j])
				b=b+1;
		}
	}
	rind=(a+b)/den;
	return(rind);
}
void centroid(double **data, double **cent, int *clust, int N, int M, int K) 
{
	int i,j,l,clustc;
	for(i=0;i<K;i++)
	{
		clustc=0;
		for(j=0;j<N;j++)
		{
			if(clust[j]==i+1)
			{
				clustc=clustc+1;
				for(l=0;l<M;l++)
					cent[i][l]=cent[i][l]+data[j][l];
			}
		}
		if(clustc!=0)
		{
			for(l=0;l<M;l++)
				cent[i][l]=cent[i][l]/clustc;
		}
	}
}
double maxm(double *arr, int n)
{
	double mx;
	int i;
	mx=arr[0];
	for(i=1;i<n;i++)
	{
		if(mx<arr[i])
			mx=arr[i];
	}
	return(mx);
}
double minm(double *arr, int n)
{
	double mn;
	int i;
	mn=arr[0];
	for(i=1;i<n;i++)
	{
		if(mn>arr[i])
			mn=arr[i];
	}
	return(mn);
}
int whichMax(double *arr, int n)
{
	double mx;
	int i, wh;
	mx=arr[0];
	wh=0;
	for(i=1;i<n;i++)
	{
		if(mx<arr[i])
		{
			mx=arr[i];
			wh=i;	
		}
	}
	return(wh);
}
int whichMin(double *arr, int n)
{
	double mn;
	int i, wh;
	mn=arr[0];
	wh=0;
	for(i=1;i<n;i++)
	{
		if(mn>arr[i])
		{
			mn=arr[i];
			wh=i;		
		}
	}
	return(wh);
}
void pointsincluster(int *clust, int *pc, int N, int K)
{
	int i,j,clustc;
	for(i=0;i<K;i++)
	{
		clustc=0;
		for(j=0;j<N;j++)
		{
			if(clust[j]==i+1)
				clustc=clustc+1;
		}
		pc[i]=clustc;
	}
}
double sigmaci(double **data, int *clust, double *cent,int ci, int N, int M)
{
	int i,j,count=0;
	double sum=0.0,sigmaci;
	for(i=0;i<N;i++)
	{
		if(clust[i]==ci)
		{
			count=count+1;
			sum=sum+eddist(data[i],cent,M);
		}
	}	
	return(sum/count);
}

void pairdist(double **data, int N, int M, double **datap)
{
	int i,j;
	for(i=0;i<N;i++)
	{
		for(j=0;j<=i;j++)
		{
			datap[i][j]=eddist(data[i],data[j],M);
			datap[j][i]=datap[i][j];
		}
	}
}
void isortval(double *val, int n, int flag)
{
	//flag=0 for ascending sort
		int i,j;
		double keyv;
		for(j=1;j<n;j++)
		{
			keyv=val[j];
			i=j-1;
			if(flag==0)
			{
				while(i>=0 && keyv<val[i])
				{
					val[i+1]=val[i];
					i=i-1;
				}
			}
			else
			{
				while(i>=0 && keyv>val[i])
				{
					val[i+1]=val[i];
					i=i-1;
				}
			}
			val[i+1]=keyv;
		}
}

double sumdoublearray(double *dArray, int len)
{
	double sum=0.0;
	int i;
	for(i=0;i<len;i++)
		sum=sum+dArray[i];
	return(sum);
}

int sumintarray(int *iArray, int len)
{
	int sum=0;
	int i;
	for(i=0;i<len;i++)
		sum=sum+iArray[i];
	return(sum);
}

char** dymemchar2d(char** cArray, int m, int n)
{
	int i;
	if(cArray==NULL)
	{
		cArray=(char**)malloc(m*sizeof(char*));		
		for(i=0;i<m;i++)
			cArray[i]=(char* )malloc(n*sizeof(char));
		cArray[i]='\0';
	}
	else
		printf("\n Variable already occupies memory");
	return(cArray);
}

int** dymemint2d(int** iArray, int m, int n)
{
	int i,j;
	if(iArray==NULL)
	{
		iArray=(int**)malloc(m*sizeof(int*));		
		for(i=0;i<m;i++)
		{
			iArray[i]=(int*)malloc(n*sizeof(int));
			for(j=0;j<n;j++)
				iArray[i][j]=0;
		}	
	}
	else
		printf("\n Variable already occupies memory");
	return(iArray);
}

double** dymemdouble2d(double** dArray, int m, int n)
{
	int i,j;
	if(dArray==NULL)
	{
		dArray=(double**)malloc(m*sizeof(double*));		
		for(i=0;i<m;i++)
		{
			dArray[i]=(double*)malloc(n*sizeof(double));
			for(j=0;j<n;j++)
				dArray[i][j]=0.0;
		}	
	}
	else
		printf("\n Variable already occupies memory");
	return(dArray);
}
double*** dymemdouble3d(double*** dArray, int M, int n, int maxd)
{
	int i,j,k;
	if(dArray==NULL)
	{
		dArray=(double***)malloc(M*sizeof(double**));		
		for(i=0;i<M;i++)
		{
			dArray[i]=(double**)malloc(n*sizeof(double*));		
			for(j=0;j<n;j++)
			{
				dArray[i][j]=(double*)malloc(maxd*sizeof(double));
				for(k=0;k<maxd;k++)
					dArray[i][j][k]=0.0;
			}
		}
	}
	else
		printf("\n Variable already occupies memory");
	return(dArray);
}

int* dynint1d(int *iArray, int l)
{
	int i;
	if(iArray==NULL)
	{
		iArray=(int*)malloc(l*sizeof(int));
	}
	else
		printf("\n Variable already occupies memory");
	for(i=0;i<l;i++)
		iArray[i]=0;
	return(iArray);
}

double* dyndouble1d(double *dArray, int l)
{
	int i;
	if(dArray==NULL)
	{
		dArray=(double*)malloc(l*sizeof(double));
	}
	else
		printf("\n Variable already occupies memory");
	for(i=0;i<l;i++)
		dArray[i]=0.0;
	return(dArray);
}

char* dynchar1d(char *cArray, int l)
{
	int i;
	if(cArray==NULL)
	{
		cArray=(char*)malloc(l*sizeof(char));
	}
	else
		printf("\n Variable already occupies memory");
	for(i=0;i<l;i++)
		cArray[i]=' ';
	return(cArray);
}
void getKMn(int* K, int* M, int* n, char* fileName)
{
	FILE *FP=NULL;
	FP=fopen(fileName,"r");
    char Header[CHAR_MAX];
    fscanf(FP,"%s %d",Header,n);
	fscanf(FP,"%s %d",Header,M);
	fscanf(FP,"%s %d",Header,K);
	printf("\n%d %d %d",*n, *M, *K);
	fclose(FP);
}
void getfiles(char* fileName, int M, int n, char** FArray, int* size, int* LogT)
{
	int i, tempI;
	FILE *FP=NULL;
	char tempC[CHAR_MAX];
	FP=fopen(fileName,"r");
	fscanf(FP,"%s %d %s %d %s %d", tempC, &tempI, tempC, &tempI, tempC, &tempI);
	for(i=0;i<M;i++)
	{
		fscanf(FP,"%s %s",tempC,FArray[i]);
		fscanf(FP,"%s %d",tempC,&size[i]);
		fscanf(FP,"%s %d",tempC,&LogT[i]);
	}
	fclose(FP);
}
void logtransform(double** logtrans, int r, int c)
{
	int i,j;
	double temp;
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
		{
			if(logtrans[i][j]==0)
				temp=1;
			else
				temp=logtrans[i][j];
			logtrans[i][j]=log10(temp);
		}
	}
}
void getdata(double*** XM, int* size, char** fnames, int n, int M)
{
	int i,j,k,dm;
	FILE *FP=NULL;
	char* fnameXi, buff[100000];
	for(i=0;i<M;i++)
	{
		fnameXi=fnames[i];
		FP=fopen(fnameXi,"r");
		fscanf(FP,"%[^\n]", buff);
		dm=size[i];
		for(j=0;j<n;j++)
		{
			fscanf(FP,"%s",buff);
			for(k=0;k<dm;k++)
			{
				fscanf(FP,"%lf",&XM[i][j][k]);
			}
		}
		fclose(FP);
	}
}
void print1dint(int *array, int len, char *name)
{
	int i;
	printf("\n %s length %d \n",name,len);
	for(i=0;i<len;i++)
		printf("%d ",array[i]);
}

void print1ddouble(double *array, int len, char *name)
{
	int i;
	printf("\n %s length %d \n",name,len);
	for(i=0;i<len;i++)
		printf("%lf ",array[i]);
}

void print2ddouble(double** X, int r, int c, char *name)
{
	int j,k;
	printf("\n %s : %d X %d \n",name,r,c);
	for(j=0;j<r;j++)
	{
		for(k=0;k<c;k++)
			printf("%lf ",X[j][k]);
		printf("\n");
	}
}

void print2ddoublepart(double** X, int startrow, int endrow, int startcol, int endcol, char *name)
{
	int j,k;
	printf("\n %s : \n",name);
	for(j=startrow;j<endrow;j++)
	{
		for(k=startcol;k<endcol;k++)
			printf("%lf ",X[j][k]);
		printf("\n");
	}
}

void print2dint(int** X, int r, int c, char *name)
{
	int j,k;
	printf("\n %s : %d X %d \n",name,r,c);
	for(j=0;j<r;j++)
	{
		for(k=0;k<c;k++)
			printf("%d ",X[j][k]);
		printf("\n");
	}
}

void print2dchar(char** X, int r, char *name)
{
	int j;
	printf("\n");
	puts(name);
	for(j=0;j<r;j++)
		puts(X[j]);
}

void print3d(double*** X, int m, int n, int d, char *name)
{
	printf("\n %s : %d X %d X %d ",name,m,n,d);
	int i,j,k;
	for(i=0;i<m;i++)
	{
		printf("\n Layer: %d \n",i);
		for(j=0;j<n;j++)
		{
			for(k=0;k<d;k++)
				printf("%lf ",X[i][j][k]);
			printf("\n");
		}
  }
}

void writeDouble2DMatrix(double **X, char *File, int r, int c)
{
	FILE *FP=NULL;
	FP=fopen(File,"w");
	int j,k;
	for(j=0;j<r;j++)
	{
		for(k=0;k<c;k++)
			fprintf(FP,"%lf ",X[j][k]);
		fprintf(FP,"\n");
	}
	fclose(FP);
}

void readDouble2DMatrix(double **X, char *File, int r, int c)
{
	FILE *FP=NULL;
	FP=fopen(File,"r");
	double temp;
	int j,k;
	for(j=0;j<r;j++)
	{
		for(k=0;k<c;k++)
			fscanf(FP,"%lf",&X[j][k]);			
	}
	fclose(FP);
}
void readDouble1DArray(double *X, char *File, int n)
{
	FILE *FP=NULL;
	FP=fopen(File,"r");
	double temp;
	int j;
	for(j=0;j<n;j++)
	    fscanf(FP,"%lf",&X[j]);			
	fclose(FP);
}


void freedouble2d(double** A, int m, int n)
{
	int i;
	for(i=0;i<m;i++)
		free(A[i]);
	free(A);
}

void freeint2d(int** A, int m, int n)
{
	int i;
	for(i=0;i<m;i++)
		free(A[i]);
	free(A);
}

void transpose(double** A, double **At, int m, int n)
{
	int i,j;
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			At[i][j]=A[j][i];
		}
	}
}

void **inverse(double **A, double **I, int n)
{
	int r,j,i;
	double temp;
	for(i=0;i<n;i++)								
		for(j=0;j<n;j++)
			if(i==j)
				I[i][j]=1;
			else
				I[i][j]=0;
	for(r=0;r<n;r++)				
	 {
		temp=A[r][r];										
		for(j=0;j<n;j++)
		{
			A[r][j]/=temp;
			I[r][j]/=temp;
		}	
		for(i=0;i<n;i++)
		{
			temp=A[i][r];
			for(j=0;j<n;j++)
			{
				if(i==r)
					break;
				A[i][j]-=A[r][j]*temp;
				I[i][j]-=I[r][j]*temp;
			}
		}
	 } 
}

void matmul(double **A,double **B,double **C,int m, int n, int r, int c)
{
	int i,j,k;
	double sum;
	if(n!=r)
		printf("\nMatrix Dimensions Mismatch");
	else
		for(i=0;i<m;i++)
	  {
			for(j=0;j<c;j++)
		  {
			 sum=0.0;
			 for(k=0;k<r;k++)
				sum=sum+A[i][k]*B[k][j];
			C[i][j]=sum;
		 }
	 }	
}

void getXi(double ***XM, double **Xi, int layer, int row, int col)
{
	int i,j,k;
	k=layer;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			Xi[i][j]=XM[k][i][j];
		}
	}
}

double getmean(double *array, int len)
{
	double mean=0.0;
	int i;
	for(i=0;i<len;i++)
		mean=mean+array[i];
	mean=mean/len;
	return(mean);
}

void meancenter(double **a, int row, int col, double **amcent)
{
	int i,j;
	double *matcol=NULL, *meanarr=NULL;
	matcol=dyndouble1d(matcol,row);
	meanarr=dyndouble1d(meanarr,col);
	for(j=0;j<col;j++)
	{
		for(i=0;i<row;i++)
			matcol[i]=a[i][j];
		meanarr[j]=getmean(matcol, row);
		for(i=0;i<row;i++)
			amcent[i][j]=a[i][j]-meanarr[j];
	}	
}

void copymat2ddouble(double **src, double **dest, int row, int col)
{
	int i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			dest[i][j]=src[i][j];
		}
	}
}
void copymat2ddoubletranspose(double **src, double **dest, int destrow, int destcol)
{
	int i,j;
	for(j=0;j<destcol;j++)
	{
		for(i=0;i<destrow;i++)
		{
			dest[i][j]=src[j][i];
		}
	}
}
void copymat1ddouble(double *src, double *dest, int len)
{
	int i;
	for(i=0;i<len;i++)
		dest[i]=src[i];
}

void copymat1dint(int *src, int *dest, int len)
{
	int i;
	for(i=0;i<len;i++)
		dest[i]=src[i];
}

void appendU(double **dest, double **src1, double **src2, int n, int r1, int r2)
{
	int i,j;
	for(i=0;i<n;i++)
	{
		for(j=0;j<r1;j++)
			dest[i][j]=src1[i][j];
	}
	for(i=0;i<n;i++)
	{
		for(j=0;j<r2;j++)
			dest[i][r1+j]=src2[i][j];
	}
}

void appendS(double **dest, double **src1, double **src2, int r1, int r2)
{
	int i;
	for(i=0;i<r1;i++)
		dest[i][i]=src1[i][i];
	for(i=0;i<r2;i++)
		dest[r1+i][r1+i]=src2[i][i];
}


void svdlapack(double **a, int row, int col, double **u, double **s, double **v)
{
	/* Locals  */
	int info, lwork;
	double *tmpa=NULL, *tmpu=NULL, *tmpv=NULL, *tmps=NULL;
	double wkopt;
	double* work;
	int i, j, ioff, iss;
 
	tmpa = dyndouble1d(tmpa, row*col);
	tmpu = dyndouble1d(tmpu, row*row);
	tmpv = dyndouble1d(tmpv, col*col);
	tmps = dyndouble1d(tmps, col);
 
	/* convert input to matrix 1D */
	ioff = 0;
	for(i=0; i<col; i++)
	{
		for(j=0; j<row; j++)
		{
			tmpa[ioff] = a[j][i];
			ioff = ioff+1;
		}
	}
 
	/* Query and allocate the optimal workspace */
	lwork = -1;
	dgesvd_( "All", "All", &row, &col, tmpa, &row, tmps, tmpu, &row, tmpv, &col, &wkopt,
		 &lwork, &info );
 
	lwork = (int)wkopt;
	work = (double*)malloc( lwork*sizeof(double) );
 
	/* Compute SVD */
	dgesvd_( "All", "All", &row, &col, tmpa, &row, tmps, tmpu, &row, tmpv, &col, work,
		 &lwork, &info );
 
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm computing SVD failed to converge.\n" );
		exit( 1 );
	}
 
	/* Convert from tmpu (matrix 1D) to u (matrix 2D) */
	ioff = 0;
	for(i=0; i<row; i++)
	{
		for(j=0; j<row; j++)
		{
			u[j][i] = tmpu[ioff];
			ioff = ioff+1;
		}
	}
 
	/* Convert from tmpv (matrix 1D) to v (matrix 2D) */
	ioff = 0;
	for(i=0; i<col; i++)
	{
		for(j=0; j<col; j++)
		{
			v[j][i] = tmpv[ioff];
			ioff = ioff+1;
		}
	}
 
	/* get minimum size from row and columns */
	if(row<col) iss = row;
	else iss=col;
 
	/* Convert from tmps (matrix 1D) to s (matrix 2D) */
	for(i=0; i<iss; i++)
	  s[i][i] = tmps[i];
 
	/* free allocated memory */
	free( (void*)work );
	free(tmpa);
	free(tmpu);
	free(tmpv);
	free(tmps);
}
// Simple K-means
void getinitcent(double **dmat, double **centset, double **initcent, int row, int col, int K)
{
	int i,j,rv,flag, *selc=NULL;
	selc=dynint1d(selc, K);
	if(centset==NULL)
	{
		for(i=0;i<K;i++)
		{
			while(1)
			{
				flag=0;
				rv=rand()%row;
				selc[i]=rv;
				for(j=0;j<=i-1;j++)
				{
					if(selc[i]==selc[j])
					{
						flag=1;
						break;
					}
				}
				if(flag==0)
					break;
			}
			rv=selc[i];
			for(j=0;j<col;j++)
				initcent[i][j]=dmat[rv][j];
		}
	}
	else
	{
		for(i=0;i<K;i++)
		{
			for(j=0;j<col;j++)
				initcent[i][j]=centset[i][j];
		}
	}
}
double** kmeans(double **data, int N, int M, int K, int *clust, int iter, double err, double **centset, int ntimes)
{
	int i,j,l,it,clustc,tm,flag,tc,ch;
	double **centP=NULL, **centN=NULL,objv,*objfn=NULL;
	double *temp=NULL,min,cdif,dfval;
	int *classN=NULL,**cass=NULL,cmin,rv,*index=NULL;;
	centP=(double **)malloc(K*sizeof(double *));
	centN=(double **)malloc(K*sizeof(double *));
	for(i=0;i<K;i++)
	{
		centP[i]=(double *)malloc(M*sizeof(double));
		centN[i]=(double *)malloc(M*sizeof(double));
	}
	cass=(int **)malloc(ntimes*sizeof(int *));
	for(i=0;i<ntimes;i++)
		cass[i]=(int *)malloc(N*sizeof(int));
	objfn=(double *)malloc(ntimes*sizeof(double));
	index=(int *)malloc(ntimes*sizeof(int));
	for(i=0;i<ntimes;i++)
		index[i]=i;
	classN=(int *)malloc(N*sizeof(int));
	temp=(double *)malloc(K*sizeof(double));
	for(tm=0;tm<ntimes;tm++)
	{
		getinitcent(data, NULL, centP, N, M, K);
		for(it=1;it<=iter;it++)
		{
			for(i=0;i<N;i++)
			{
				for(j=0;j<K;j++)
					temp[j]=eddist(data[i],centP[j],M);  
				min=temp[0];
				cmin=1;
				for(j=1;j<K;j++)
				{
					if(temp[j]<min)
					{
						min=temp[j];
						cmin=j+1;
					}
				}
				classN[i]=cmin;
			}
			for(i=0;i<K;i++)
			{
				for(j=0;j<M;j++)
					 centN[i][j]=0.0; 
			}
			for(i=0;i<K;i++)
			{
				clustc=0;
				for(j=0;j<N;j++)
				{
					if(classN[j]==i+1)
					{
						clustc=clustc+1;
						for(l=0;l<M;l++)
							centN[i][l]=centN[i][l]+data[j][l];
					}
				}
				if(clustc!=0)
				{
					for(l=0;l<M;l++)
						centN[i][l]=centN[i][l]/clustc;
				}
			}
			cdif=0.0;
			for(i=0;i<K;i++)
			{
				for(j=0;j<M;j++)
				{
					dfval=centN[i][j]-centP[i][j];
					if(dfval<0)
						dfval=-dfval;
					cdif=cdif+dfval;
				}
			}
			cdif=cdif/K;
			if(cdif<err)
				break;
			for(i=0;i<K;i++)
			{
				for(j=0;j<M;j++)
					centP[i][j]=centN[i][j];
			}
		}
		objv=0.0;
		for(i=0;i<N;i++)
		{
			cass[tm][i]=classN[i];
			objv=objv+eddist(data[i],centN[classN[i]-1],M);		
		}
		objfn[tm]=objv;
  }
 	isortwindex(objfn, index, ntimes, 0);
 	ch=index[0];
 	for(i=0;i<N;i++)
		clust[i]=cass[ch][i];
	for(i=0;i<K;i++)
		free(centP[i]);
	free(centP);
	for(i=0;i<ntimes;i++)
		free(cass[i]);
	free(cass);
	free(objfn);
	free(temp);
	free(classN);
/*	printf("\n Final Cluster Assignments : \n");*/
/*	for(i=0;i<N;i++)*/
/*		printf("%d ",clust[i]);*/
/*	printf("\n No of Data points in each Cluster : ");*/
/*	for(i=0;i<K;i++)*/
/*	{*/
/*		clustc=0;*/
/*		for(j=0;j<N;j++)*/
/*		{*/
/*			if(clust[j]==i+1)*/
/*				clustc=clustc+1;*/
/*		}*/
/*		printf("%d ",clustc);*/
/*	}*/
}
double SILHOUETTE(double **data, int *clust, int N, int M, int K, double **datap, int *pc)
{	
	int i,j,l,m,t;
	double a,b,sumn,sumd,num,den,sgval,sil, *denarr=NULL, *pointsil=NULL;
	denarr=(double *)malloc((K-1)*sizeof(double));
	pointsil=(double *)malloc(N*sizeof(double));
	int nk=0;
	double silsum=0.0;
	for(i=0;i<K;i++)
	{
		sgval=0.0;
		nk=0;
		for(j=0;j<N;j++)
		{
			if(clust[j]==i+1)
			{
				sumn=0.0;
				for(l=0;l<N;l++)
				{
					if(clust[l]==i+1)
						sumn=sumn+datap[j][l];
				}
				a=(double)sumn/(pc[i]-1);
				t=0;
				for(l=0;l<K;l++)
				{
					if(l!=i)
					{
					 sumd=0.0;
						for(m=0;m<N;m++)
						{
							if(clust[m]==l+1)
								sumd=sumd+datap[j][m];
						}
						denarr[t]=sumd/pc[l];
						t=t+1;
					}
				}
				b=minm(denarr,K-1);
				num=b-a;
				if(a>b)
					den=a;
				else
					den=b;
				pointsil[j]=num/den;
				sgval=sgval+(num/den);
				nk=nk+1;
			}
		}
		if(nk!=0)
		{
			silsum=silsum+sgval/nk;
		}
	}
	free(denarr);
	free(pointsil);
	sil=silsum/K;
	return(sil);
}
//Hartigan Wong K-means
#include "asa136.h"
void HartiganKm(double **dmat, int row, int col, int K, int *clust, double **centset, int ntimes, int iter, double* withinss)
{
	int i,j,flag,tm, *index=NULL, **cass=NULL;
	double **centroid=NULL, *wssntimes=NULL;
	centroid = dymemdouble2d(centroid, K, col);
	cass= dymemint2d(cass, ntimes, row); 
	wssntimes=dyndouble1d(wssntimes, ntimes);
	index=dynint1d(index, ntimes);
	for(i=0;i<ntimes;i++)
		index[i]=i;
	double *a;
   double *c;
   int *ic1;
   int ifault;
   int *nc;
   int nc_sum;
   double *wss;
   double wss_sum;
   a = ( double * ) malloc ( row * col * sizeof ( double ) );
   c = ( double * ) malloc ( K * col * sizeof ( double ) );
   ic1 = ( int * ) malloc ( row * sizeof ( int ) );
   nc = ( int * ) malloc ( K * sizeof ( int ) );
   wss = ( double * ) malloc ( K * sizeof ( double ) );
   for ( i = 1; i <= row; i++ )
	{
		for ( j = 1; j <= col; j++ )
		{
			a[i-1+(j-1)*row] = dmat[i-1][j-1];
		}
	}
	for(tm=0;tm<ntimes;tm++)
	{
/*		printf("\n Execution %d \n",tm);*/
		getinitcent(dmat, NULL, centroid, row, col, K);
		for ( i = 1; i <= K; i++ )
		{
			for ( j = 1; j <= col; j++ )
			{
				c[i-1+(j-1)*K] = centroid[i-1][j-1];
			}
		}
		kmns ( a, row, col, c, K, ic1, nc, iter, wss, &ifault );
	   if ( ifault != 0 )
	   {
		  printf ( "\n" );
		  printf ( "TEST02 - Fatal error!\n" );
		  printf ( "  KMNS returned IFAULT = %d\n", ifault );
		  return;
	   }
		nc_sum = 0;
  		wss_sum = 0.0;
	   for ( i = 1; i <= K; i++ )
	   {
		  nc_sum = nc_sum + nc[i-1];
		  wss_sum = wss_sum + wss[i-1];
	   }
	   wssntimes[tm]=wss_sum;
	   for(i=0;i<row;i++)
			cass[tm][i]=ic1[i];
	}
	isortwindex(wssntimes, index, ntimes, 0);
 	for(i=0;i<row;i++)
		clust[i]=cass[index[0]][i];
	*withinss=wssntimes[index[0]];
	freeint2d(cass, ntimes, row);	
	freedouble2d(centroid, K, col);
	free(wssntimes);
	free(index);
}

int maxarrint(int *arr, int n)
{
	int mx=0;
	int i;
	mx=arr[0];
	for(i=1;i<n;i++)
	{
		if(mx<arr[i])
			mx=arr[i];
	}
	return(mx);
}
int intersect(int *iArray1, int *iArray2, int n1, int n2)
{
	int count,i,j,check=0;
	count=0;
	for(i=0;i<n1;i++)
	{
		check=iArray1[i];
		for(j=0;j<n2;j++)
		{
			if(check==iArray2[j])
			{
				count=count+1;
				break;
			}	
		}
	}
	return(count);
}
int* findwhichint(int v, int *iArray, int n, int *len)
{
	int *found=NULL;
	int count=0,i,j;
	for(i=0;i<n;i++)
	{
		if(iArray[i]==v)
			count=count+1;
	}
	found=(int*)malloc(count*sizeof(int));
	count=0;	
	for(i=0;i<n;i++)
	{
		if(iArray[i]==v)
		{
			found[count]=i;
			count=count+1;
		}
	}
	*len=count;
	return(found);
}
void confusion(int *original, int *clust, int n, double *Rand, double *ARI, double *Dice, double *Jac, double *Purity, double *Fmeasure, double *NMI)
{
	int Kt, Kc, i ,j;
	Kt=maxarrint(original, n);
	Kc=maxarrint(clust, n);
	int n11,n10,n01,n00;
	n11=0; n10=0; n01=0; n00=0;
	for(i=0;i<n-1;i++)
	{
		for(j=i+1;j<n;j++)
		{
			if((original[i]==original[j])&(clust[i]==clust[j]))
				n11=n11+1;
			if((original[i]==original[j])&(clust[i]!=clust[j]))
				n10=n10+1;
			if((original[i]!=original[j])&(clust[i]==clust[j]))
				n01=n01+1;
			if((original[i]!=original[j])&(clust[i]!=clust[j]))
				n00=n00+1;
		}
	}
	*Rand=(double)(n11+n00)/(n11+n10+n01+n00);
	*ARI=(double)(2*(n00*n11-n01*n10))/((n00+n01)*(n01+n11)+(n00+n10)*(n10+n11));
	*Dice=(double)(2*n11)/(2*n11+n01+n10);
	*Jac=(double)n11/(n11+n10+n01);
	int **ContgCT=NULL, *ci=NULL, *tj=NULL, *tempr=NULL;
	double **Fij=NULL, *tempc=NULL, maxFi=0;
	ContgCT=dymemint2d(ContgCT, Kc, Kt);
	int ni,nj, nij, maxtj;
	for(i=0;i<Kc;i++)
	{
		for(j=0;j<Kt;j++)
		{
			ci=NULL;
			tj=NULL;
			ni=0; nj=0; nij=0;
			ci=findwhichint(i+1, clust, n, &ni);
			tj=findwhichint(j+1, original, n, &nj);
			nij=intersect(ci,tj, ni, nj);
			ContgCT[i][j]=nij;
		}
	}
	*Purity=0.0;
	tempr=dynint1d(tempr, Kt);
	for(i=0;i<Kc;i++)
	{
		for(j=0;j<Kt;j++)
			tempr[j]=ContgCT[i][j];
		maxtj=maxarrint(tempr, Kt);
		*Purity=(double)*Purity+maxtj;
	}
	*Purity=(double)(*Purity/n);
	*Fmeasure=0.0;
	Fij=dymemdouble2d(Fij, Kc, Kt);
	double *rowsumKc=NULL, *colsumKt=NULL, smv=0.0;
	rowsumKc=dyndouble1d(rowsumKc, Kc);
	colsumKt=dyndouble1d(colsumKt, Kt);
	for(i=0;i<Kc;i++)
	{
		smv=0.0;
		for(j=0;j<Kt;j++)
			smv=smv+ContgCT[i][j];
		rowsumKc[i]=smv;
	}
	for(j=0;j<Kt;j++)
	{
		smv=0.0;
		for(i=0;i<Kc;i++)
			smv=smv+ContgCT[i][j];
		colsumKt[j]=smv;
	}
	for(i=0;i<Kc;i++)
	{
		for(j=0;j<Kt;j++)			
			Fij[i][j]=(2*ContgCT[i][j])/(rowsumKc[i]+colsumKt[j]);
	}
	tempc=dyndouble1d(tempc, Kc);
	for(j=0;j<Kt;j++)		
	{
		nj=colsumKt[j];
		for(i=0;i<Kc;i++)
			tempc[i]=Fij[i][j];
		maxFi=maxm(tempc, Kc);
		*Fmeasure=(double)(*Fmeasure+nj*maxFi);
	}
	*Fmeasure=(double)(*Fmeasure/n);
	double dummy;
	double InfCT=0, HC=0, HT=0, wl, civ, tjv;
	for(i=0;i<Kc;i++)
	{
		for(j=0;j<Kt;j++)		
		{
		  wl=(n*ContgCT[i][j])/(rowsumKc[i]*colsumKt[j]);
		  if(wl!=0)
		  {
		  	dummy=((double)ContgCT[i][j])/((double)n);
		  	InfCT=InfCT+(double)dummy*log(wl);
		  }
				
		}
	}
	for(i=0;i<Kc;i++)
	{
		civ=rowsumKc[i]/n;
		if(civ!=0)
			HC=HC+civ*log(civ);
	}
	HC=-HC;
	for(j=0;j<Kt;j++)		
	{
		tjv=colsumKt[j]/n;
		if(tjv!=0)
			HT=HT+tjv*log(tjv);
	}
	HT=-HT;
	*NMI=2*InfCT/(HC+HT);
	//printf("\n Rand %lf, ARI %lf, Dice %lf, Jac %lf, Purity %lf, Fmeasure %lf, NMI %lf",*Rand, *ARI, *Dice, *Jac, *Purity, *Fmeasure, *NMI);
}
void constantarrayint(int *iArray, int l, int val)
{
	int i;
	for(i=0;i<l;i++)
		iArray[i]=val;
}
void constantarraydouble(double *dArray, int l, int val)
{
	int i;
	for(i=0;i<l;i++)
		dArray[i]=val;
}
double forbnorm(double** dmat, int r, int c)
{
	double sum=0.0;
	int i,j;
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
			sum=sum+(dmat[i][j]*dmat[i][j]);
	}
	return(sum);
}
double matdiff(double** amat, double** bmat, double** diff, int r, int c)
{
	int i,j;
	for(i=0;i<r;i++)
	{
		for(j=0;j<c;j++)
			diff[i][j]=amat[i][j]-bmat[i][j];
	}
}
void QRlapack(double **A, double **Q, int row, int col)
{
int i, j, ioff;
double *TAU=NULL;
int *JPVT=NULL, info, lwork, k=col;
double wkopt1, wkopt2;
double *work1=NULL, *work2=NULL;
JPVT=dynint1d(JPVT,col);
TAU=dyndouble1d(TAU,col);
double *tmpa=NULL;
tmpa = dyndouble1d(tmpa, row*col);
ioff = 0;
for(i=0; i<col; i++)
{
	for(j=0; j<row; j++)
	{
		tmpa[ioff] = A[j][i];
		ioff = ioff+1;
	}
}
lwork = -1;
dgeqp3_(&row, &col, tmpa, &row, JPVT, TAU, &wkopt1, &lwork, &info );
lwork = (int)wkopt1;
work1 = (double*)malloc( lwork*sizeof(double) );
dgeqp3_(&row, &col, tmpa, &row, JPVT, TAU, work1, &lwork, &info );

lwork = -1;
dorgqr_(&row, &col, &k, tmpa, &row, TAU, &wkopt2, &lwork, &info );
lwork = (int)wkopt2;
work2 = (double*)malloc( lwork*sizeof(double) );
dorgqr_(&row, &col, &k, tmpa, &row, TAU, work2, &lwork, &info );

ioff = 0;
for(i=0; i<col; i++)
{
	for(j=0; j<row; j++)
	{
			Q[j][i] = tmpa[ioff];
			ioff = ioff+1;
	}
}
free(tmpa);
free(TAU);
free(JPVT);
free(work1);
free(work2);
}

void cbind(double **A, double **B, double **AB, int row, int colA, int colB)
{
	int i,j,k;
	for(i=0;i<row;i++)
	{
		for(j=0;j<colA;j++)
		{
			AB[i][j]=A[i][j];
		}
	}
	for(i=0;i<row;i++)
	{
		for(j=0;j<colB;j++)
		{
			AB[i][j+colA]=B[i][j];
		}
	}
}
void rbind(double **A, double **B, double **AB, int rowA, int rowB, int col)
{
	int i,j,k;
	for(i=0;i<rowA;i++)
	{
		for(j=0;j<col;j++)
		{
			AB[i][j]=A[i][j];
		}
	}
	for(i=0;i<rowB;i++)
	{
		for(j=0;j<col;j++)
		{
			AB[i+rowA][j]=B[i][j];
		}
	}
}
void stackMatrixRowWise(double ***XM, double **stacked, int *size, int M, int n)
{
	double **Xi=NULL;
	int i,j,l, row, col, sumpi=0;
	for(i=0;i<M;i++)
	{
		Xi=NULL;
		row=n;
		col=size[i];
		Xi = dymemdouble2d(Xi, row, col);
		getXi(XM, Xi, i, row, col);
		if(i==1 || i==2)
			logtransform(Xi,row,col);
		for(j=0;j<size[i];j++)
		{
			for(l=0;l<n;l++)
			{
				stacked[j+sumpi][l]=Xi[l][j];
			}
		}
		freedouble2d(Xi, row, col); 
		sumpi=sumpi+size[i];
	}	
}
void constmat(double **Zmat, int row, int col, int val)
{
	int i,j;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
			Zmat[i][j]=val;
	}
}
void matColumnToVector(double **Mat, int n, int col, double *shapiroX)
{
	int i;
	for(i=0;i<n;i++)
		shapiroX[i]=Mat[i][col];
}
