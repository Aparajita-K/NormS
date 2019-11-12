Normality<-function(DataL,mod)
{

	alpha=0.05
	Thetas=list()
	M=length(Data)
	present=rep(0,M)
	rel=rep(0,M)
	red=rep(Inf,M)
	ranks=rep(0,M)
	sumD=0
	for(i in 1:M)
	{
		cat("\n\n",mod[i],"\n")
		Dmat=DataL[[i]]
		n=dim(Dmat)[1]
		d=dim(Dmat)[2]
		DmatSc=scale(Dmat,center=TRUE,scale=FALSE)
		dsv=svd(DmatSc)
		mxpc=11
		DmU=dsv$u[,1:mxpc]
		DmS=dsv$d[1:mxpc]
	    nrm=multinormHkd(DmU,mod=mod[i])
		npc=nrm$npc
		cat("\n #Components=",npc)
		if(length(npc)>0)
		{
			Thetas[[i]]=list(U=as.matrix(DmU[,npc]),D=DmS[npc],H=nrm$H,p.val=nrm$p.val,df=nrm$df,allnorm=nrm$allnorm)
			cat("\n Modality",i,"Hstat=",nrm$stat,"P-value=",nrm$p.val,"DF=",nrm$df)
			present[i]=1
			ranks[i]=length(npc)
		}
		else
		{
			present[i]=0
			Thetas[[i]]=list(U=as.matrix(DmU[,1]),D=DmS[1],H=nrm$H,p.val=nrm$p.val,df=nrm$df,allnorm=nrm$allnorm)
		}
		sumD=sumD+sum(Thetas[[i]]$D)
	}
	cat("\n Present Array=",present)
	modcount=sum(present)
	cat("\nModality Relevances=",rel)
	for(i in 1:M)
	{
			df=Thetas[[i]]$df
			mvnTemp<-Thetas[[i]]$H
			rel[i]=0.5*(1+(mvnTemp-qchisq((1-alpha),df))/max(mvnTemp,qchisq((1-alpha),df)))
			rel[i]=rel[i]*(sum(Thetas[[i]]$D)/sumD)
	}
	cat("\nIndividual Ranks=",ranks)
	cat("\nModality Relevances=",rel)
	for(m in 1:M)
	    cat("\n",mod[m],"Relevance=", rel[m], "Rank=",ranks[m])
	order=sort(rel,index.return=TRUE,decreasing=TRUE)$ix
	cat("\nModality Order:",mod[order],"\n\n")
	ThetaAs=list()
	ThetaAs=list(U=as.matrix(Thetas[[order[1]]]$U),D=Thetas[[order[1]]]$D,H=Thetas[[order[1]]]$H,df=Thetas[[order[1]]]$df)
	cat("\n Dim ThetaAs U-space dimension=",dim(ThetaAs$U))
	taken=order[1]
	present[taken]=0
	for(j in 2:modcount)
	{
		actualdf=rep(0,M)
		red=rep(0,M)
		redS=rep(0,M)
		redA=rep(0,M)
		for(i in 1:M)
		{
			if(present[i]==1)
			{
				G=as.matrix(t(ThetaAs$U)%*%Thetas[[i]]$U)
        		P=as.matrix(ThetaAs$U%*%G)
        		R=as.matrix(Thetas[[i]]$U-P)
        		redA[i]=sum(P^2)/sum(Thetas[[i]]$U^2)
			}
		}

        red=redA
		checkmod=which.max(red)
		present[checkmod]=0
		cat("\n Redundancy=",red)
		cat("\n\n Checking Residual of Module: ",checkmod,mod[checkmod])
		
		G=as.matrix(t(ThetaAs$U)%*%Thetas[[checkmod]]$U)
		P=as.matrix(ThetaAs$U%*%G)
		R=as.matrix(Thetas[[checkmod]]$U-P)
		nrmR=resNormTest(R,mod=paste0("Residual-",mod[checkmod]))
		whresnrm=nrmR$whresnrm
		cat("\n whresnrm=",whresnrm)
		if(whresnrm[1]!=0)
		{
			curJointMat=ThetaAs$U
			dimcur=dim(curJointMat)[2]
			mvncur<-roy(curJointMat)
			mvnpvcur=mvncur$p.value
			resJoinMat=cbind(ThetaAs$U,as.matrix(R[,whresnrm]))
			mvnresjoin<-roy(resJoinMat)
			mvnpvresjoin=mvnresjoin$p.value
			Uch=cbind(ThetaAs$U,Thetas[[checkmod]]$U[,whresnrm])
			Dch=c(ThetaAs$D,Thetas[[checkmod]]$D[whresnrm])
			ThetaAs=list(U=Uch,D=Dch,H=mvnresjoin$H,df=mvnresjoin$df)
			taken=c(taken,checkmod)			
			cat("\nCurrent Subspace: H-stat=",mvncur$H,"P-value=",mvncur$p.val,"DF=",mvncur$df)
			cat("\n Ki's=",mvncur$indvki)
			cat("\nUpdated Subspace: H-stat=",mvnresjoin$H,"P-value=",mvnresjoin$p.val,"DF=",mvnresjoin$df)
			cat("\n Ki's=",mvnresjoin$indvki)
			cat("\n Dimension ThetaAs:",length(ThetaAs$D))
			cat("\n\n Joint Module Updated with:",mod[checkmod])
		}
	    else
		    cat("\n\n Residuals are all Normal. \n Joint Module NOT Updated with:",mod[checkmod])
		
    }
	cat("\nModalities Considered for Joint Subspace:",mod[taken])
	if(length(ThetaAs$D)==1)
        Dmat=as.matrix(ThetaAs$U*(ThetaAs$D))
	else
		Dmat=as.matrix(ThetaAs$U%*%diag(ThetaAs$D))
	params=paste("Modalities=",paste(mod[taken],collapse=","),"DegreeFreedom=",ThetaAs$df,"Hstat=",ThetaAs$H)
	return(list(Dmat=Dmat,params=params))
}

resNormTest<-function(Dmat,mod="GENE",alpha=0.05)
{
	tc=dim(Dmat)[2]
	mvnt<-roy(Dmat)
	indvki=mvnt$indvki
	pvalki=pchisq(indvki, 1, lower.tail = FALSE)
	cat("\nUnivariate H-stat with 1-DF p-values of Residuals:\n",indvki)
	cat("\nUnivariate H-stat with 1-DF Statistic of Residuals:\n",pvalki)
	cat("\n Residual Space H-stat=",mvnt$H,"P-value=",mvnt$p.val,"Estimated DF=",mvnt$df)
	take=rep(0,tc)
	for(i in 1:tc)
	{
		if(pvalki[i]<alpha)
			take[i]=1
	}
	whresnrm=which(take==1)
	if(length(whresnrm)==0)
		whresnrm=0		
	return(list(whresnrm=whresnrm))
}

multinormHkd<-function(Dmat,mod="GENE",alpha=0.05)
{	
	n=dim(Dmat)[1]
	d=dim(Dmat)[2]-1
	npc=NULL
	H=0
	p.val=1
	mvnt<-roy(Dmat)
	indvki=mvnt$indvki
	pvalki=pchisq(indvki, 1, lower.tail = FALSE)
	df=0
	count=NULL
	for(i in 1:d)
	{
		cat("\n Component=",i,"ki=",indvki[i],"P-value for ki=",pvalki[i])
		if(pvalki[i]>alpha && pvalki[i+1]>alpha)
			break
        else
		{	    
		        count=c(count,i)
		}	
	}
	allnorm=0
	for(i in 1:d)
	    if(pvalki[i]<alpha && indvki[i]<30)
	        allnorm=allnorm+1
	cat("\nallnorm=",allnorm)
	#allnorm=1
	if(length(count)>0)
	{
		npc=count
		df=length(count)
		H=sum(indvki[count])
		p.val=pchisq(H, count, lower.tail = FALSE)
	}
	if(length(count)==0)
	{
	    H=indvki[1]
	    df=1
	}
	cat("\n #components by Ki=",npc)		
	return(list(npc=npc,H=H,p.val=p.val,df=df,allnorm=allnorm))
}

roy<-function (data, qqplot = FALSE) 
{
    library(nortest)
	library(e1071)
    if (!is.data.frame(data) && !is.matrix(data)) 
        stop("Input must be one of classes \"data frame\" or \"matrix\"")
    dataframe = as.data.frame(data)
    dname <- deparse(substitute(data))
    data <- data[complete.cases(data), ]
    data <- as.matrix(data)
    p <- dim(data)[2]
    n <- dim(data)[1]
    z <- matrix(nrow <- p, ncol = 1)
    z <- as.data.frame(z)
    w <- matrix(nrow <- p, ncol = 1)
    w <- as.data.frame(w)
    data.org <- data
    if (n <= 3) {
        stop("n must be greater than 3")
    }
    else if (n >= 4 || n <= 11) {
        x <- n
        g <- -2.273 + 0.459 * x
        m <- 0.544 - 0.39978 * x + 0.025054 * x^2 - 0.0006714 * 
            x^3
        s <- exp(1.3822 - 0.77857 * x + 0.062767 * x^2 - 0.0020322 * 
            x^3)
        for (i in 1:p) {
            a2 <- data[, i]
            {
                if (kurtosis(a2) > 3) {
                  w <- sf.test(a2)$statistic
                }
                else {
                  w <- shapiro.test(a2)$statistic
                }
            }
            z[i, 1] <- (-log(g - (log(1 - w))) - m)/s
        }
    }
    if (n > 2000) {
        stop("n must be less than 2000")
    }
    else if (n >= 12 || n <= 2000) {
        x <- log(n)
        g <- 0
        m <- -1.5861 - 0.31082 * x - 0.083751 * x^2 + 0.0038915 * 
            x^3
        s <- exp(-0.4803 - 0.082676 * x + 0.0030302 * x^2)
        for (i in 1:p) {
            a2 <- data[, i]
            {
                if (kurtosis(a2) > 3) {
                  w <- sf.test(a2)$statistic
                }
                else {
                  w <- shapiro.test(a2)$statistic
                }
            }
            z[i, 1] <- ((log(1 - w)) + g - m)/s
        }
    }
    else {
        stop("n is not in the proper range")
    }
    if(p==1)
    {
    	 RH<-(qnorm((pnorm(-z[p,1]))/2))^2
    	 pv <- pchisq(RH, 1, lower.tail = FALSE)
    	 indvki=RH
    	 edf=p
    }
    else
    {
       u <- 0.715
		 v <- 0.21364 + 0.015124 * (log(n))^2 - 0.0018034 * (log(n))^3
		 l <- 5
		 C <- cor(data)
		 NC <- (C^l) * (1 - (u * (1 - C)^u)/v)
		 T <- sum(sum(NC)) - p
    	 mC <- T/(p^2 - p)
    	 edf <- p/(1 + (p - 1) * mC)
    	 Res <- matrix(nrow = p, ncol = 1)
		 Res <- as.data.frame(Res)
		 for (i in 1:p) {
		     Res <- (qnorm((pnorm(-z[, ]))/2))^2
		 }
		 RH <- (edf * (sum(Res)))/p
		 pv <- pchisq(RH, edf, lower.tail = FALSE)
		 indvki=as.vector(Res)
    }
    
    result <- list(H = RH, p.value = pv, dname = dname, indvki=indvki, df=edf)
    return(result)
}
