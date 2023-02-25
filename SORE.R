
#---------------------------------------------------------------------------------------------------------------------------------
#---------------------- Functions to Fit SORE models  ------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------
#---------------------- Fit SORE model ------------------------------------------------
library(Formula)
SORE <- function(formula, W=NULL, xstar=NULL, digit=NULL, RR.GC=F, stdR=F, stdE=F, data, parstart=NULL, hess=F, prof=F,  trace=0, tolprof=1e-8){
  ## formula: a two-part formula with the first part describing the model for y and the second part giving names of endogenous regressors X or its normal score transformations.
  ## W: a one-sided formula to specify exogenous regressors
  ## xstar: the column positions of the matrix x that undergo normal transformtion to model GC regressor-error dependence. 
  ## RR.GC: A flag: =T if regressor-regressor follows a GC relationship.  
  ## digits: a vector denoting the decimal place after rounding continuous endogenous regressors in SORE modeling. No rounding if digits=NULL
  ## stdR: a flag: =T is standardize regressors, =F if not. 
  ## stdE: a flag: =T is standardize error term by its SD estimate, =F if not
  ## data: name of data frame containing all data elements
  ## parstart: starting values of parameters. Default to OLS estimates. 
  ## hess: =T if computing hessian matrix
  ## prof: =T if using profile likelihood
  ## trace: the amount of printed information at intermediate iteration steps
  ## tolprof: the threshold value used in the profile likelihood to determine the convergence of  baseline parameters. 

  F.formula <- as.Formula(formula)
  ymodel = formula(F.formula, lhs = 1, rhs = 1)
  x <- as.matrix(model.matrix(F.formula, data = data, rhs = 1))
  y <- model.part(F.formula, data = data, lhs = 1)[, 1]
  if (ncol(x) < 1)  stop("No predictor variables specified for the outcome")
  endox <- as.matrix(model.matrix(F.formula, data = data, rhs = length(F.formula)[2]))
  endox <- endox[,colnames(endox)!= "(Intercept)", drop=F] ## remove the intercept to obtain endogenous regressors only.
  nendox <- ncol(endox); 
  if (!is.null(digit)) for (i in 1:ncol(endox)) endox[,i]= round(endox[,i], digit[i])
  if (stdR) endox<-scale(endox)
  x4w <- endox
  
  ns=nrow(endox)/(nrow(endox)+1)
  if (!is.null(xstar))  {
    for (k in 1:length(xstar)) {
      loc=xstar[k]
      endox[,loc]=ecdf(endox[,loc])(endox[,loc]); endox[endox[,loc]==1,loc]=ns
      endox[,loc]= qnorm(endox[,loc])
    }
  }

  if (is.null(W)) wmat <- NULL else
  {
    wmat <- as.matrix(model.matrix(W, data = data))
    wmat <- wmat[,colnames(wmat)!= "(Intercept)", drop=F]
    if (stdR) wmat=scale(wmat)
    if (RR.GC==T) {
        for (k in 1:ncol(x4w)) {x4w[,k]=ecdf(x4w[,k])(x4w[,k]); x4w[x4w[,k]==1,k]=ns;x4w[,k]=qnorm(x4w[,k])}
        for (k in 1:ncol(wmat)) {wmat[,k]=ecdf(wmat[,k])(wmat[,k]); wmat[wmat[,k]==1,k]=ns;wmat[,k]=qnorm(wmat[,k])}
      }
      
  }
  
  if (is.null(wmat)) nw<-0 else 
    nw<- ncol(wmat)
  
  ## set the reference point
  for (i in 1:ncol(endox)) endox[,i]=endox[,i]-min(endox[,i])
  for (i in 1:ncol(x4w)) x4w[,i]=x4w[,i]- min(x4w[,i])  

  
  ## create parameters 
  data.ols <- lm(ymodel, data=data)
  beta <- coef(data.ols); lsigma <- log(summary(data.ols)$sigma)
  parx <- ux<- ux4w<-ux4wstar<- uxcount<- nparx<- vector("list", ncol(endox))
  for (i in 1:nendox) {
    endox.table <-table(endox[,i])
    ux[[i]] = as.numeric(names(table(endox[,i])))
    uxcount[[i]]<-as.data.frame(endox.table)$Freq
    ux4w[[i]] = as.numeric(names(table(x4w[,i])))
    parx[[i]] = list(gamma=rep(0,i+nw), lambda= log(uxcount[[i]])[-length(uxcount[[i]])]-log(uxcount[[i]])[length(uxcount[[i]])])
    nparx[[i]] =list(ngamma=length(parx[[i]]$gamma), nlambda=length(parx[[i]]$lambda))
  }
  
  paryx <- c(beta, lsigma, unlist(parx))
  if (prof==T) {paryx<-c(beta, lsigma); for (i in 1:nendox) paryx<-c(paryx,parx[[i]]$gamma)}
  parnames<-names(paryx)
  if (!is.null(parstart)) {paryx=parstart;} 


  if (prof==T)  fitSORE<-nlminb(paryx, obj=fun.profnllSORE,  y=y, X=x, W=wmat, endox=endox, ux=ux, uxcount=uxcount, ux4z=ux4w, nparx=nparx, xstar=xstar,stdE=stdE,  tolprof=tolprof,
                    control=list(trace=trace, eval.max=10000,iter.max=100000)) else
                      fitSORE<-nlminb(paryx, obj=fun.nllSORE,  y=y, X=x, W=wmat, endox=endox, ux=ux, ux4z=ux4w,nparx=nparx, xstar=xstar,stdE=stdE,
                                      control=list(trace=trace, eval.max=10000,iter.max=100000))

  if (hess==T) {
    print("Compute Hessian Matrix now...")
    
    if (prof==T) {
      fitHess=optimHess(par=fitSORE$par, fn=fun.profnllSORE,  y=y, X=x, W=wmat, endox=endox, ux=ux, uxcount=uxcount,ux4z=ux4w, nparx=nparx, xstar=xstar,stdE=stdE, tolprof=tolprof,
                        control=list(trace=1, maxit=100000))
    }else
    {
      fitHess=optimHess(par=fitSORE$par, fn=fun.nllSORE,  y=y, X=x, W=wmat, endox=endox, ux=ux, ux4z=ux4w, nparx=nparx, xstar=xstar,stdE=stdE,
                        control=list(trace=1, maxit=100000))
    }
    
    return(list(par=fitSORE$par, se=sqrt(diag(solve(fitHess))), objective=fitSORE$obj, convergence=fitSORE$converg, iterations=fitSORE$iterations, evaluations=fitSORE$evaluations))
  }
  else return(fitSORE) 

}


fun.profnllSORE <- function(paryx, y, X, W, endox, ux, uxcount, ux4z,  nparx, xstar,stdE, tolprof, lambda=NULL){
  res<- fun.findlambda(paryx, y, X, W, endox, ux, uxcount, ux4z, nparx, xstar,stdE,tolprof,lambda=lambda)$res
  res
}

fun.findlambda <- function(paryx, y, X, W, endox, ux, uxcount,ux4z, nparx, xstar,stdE,tolprof, lambda=NULL){
  
  beta=paryx[1:ncol(X)]
  sig <-exp(paryx[length(beta)+1])
  parx <-  vector("list", ncol(endox))
  
  eps<- y-X%*%beta
  eps <- eps[,1]
  ## Obtain the lambda parameters 
  lambdaest = vector("list", ncol(endox))
  
  k=length(beta)+1; 
  for (i in 1:ncol(endox)) {
    parx[[i]] = paryx[(k+1):(k+nparx[[i]]$ngamma)]
    k= k+nparx[[i]]$ngamma  
  }
  
  
  res <-  sum(dnorm(eps, sd=sig, log=T))
  
  if(stdE) eps<-eps/sig
  for (i in 1:ncol(endox)) {
    if (is.null(lambda))  temp= profllk.SORE(eps=eps, X=endox[,1:i],Z=W, ux=ux[[i]], uxcount=uxcount[[i]], ux4z=ux4z[[i]],  gamma=parx[[i]], tolprof=tolprof) else
      temp= profllk.SORE(eps=eps, X=endox[,1:i],Z=W, ux=ux[[i]], uxcount=uxcount[[i]], ux4z=ux4z[[i]], gamma=parx[[i]], lambda=lambda[[i]], tolprof=tolprof)
    
    res= res+ temp$res
    lambdaest[[i]]= temp$lambdak 
  }
  return(list(res=-res, lambda=lambdaest))
}



profllk.SORE<- function(eps, X, Z, ux,uxcount, ux4z, gamma,  lambda=NULL, tolprof) {

  X<-as.matrix(X); 
  if (is.null(Z)) {
    if (ncol(X)==1) emax=as.matrix(eps)%*%as.matrix(gamma)%*%t(as.matrix(ux)) else 
      if (ncol(X)>1) emax=as.matrix(cbind(X[,1:(ncol(X)-1)],eps))%*%as.matrix(gamma)%*%t(as.matrix(ux)) else
        stop("The dimension of X in SOR model is less than 1")
  } else {
    if (ncol(X)==1) {
      emax= Z %*%as.matrix(gamma[-length(gamma)])%*%t(as.matrix(ux4z))+ (as.matrix(eps)*gamma[length(gamma)])%*%t(as.matrix(ux))
    } else 
      if (ncol(X)>1) emax=Z %*%as.matrix(gamma[1:ncol(Z)])%*%t(as.matrix(ux4z))+ as.matrix(cbind(X[,1:(ncol(X)-1)],eps))%*%as.matrix(gamma[-(1:ncol(Z))])%*%t(as.matrix(ux)) else
        stop("The dimension of X in SOR model is less than 1")
  }
  emax=exp(emax)
  
  
  lambda1 <- log(uxcount[-length(uxcount)]/uxcount[length(uxcount)])
  lambda0 <- lambda1+1  
  if (is.null(lambda)){
    emax[is.infinite(emax)]=exp(700) ## avoid infinite
 
    itermax=0
    while (max(abs(lambda1-lambda0)) > tolprof & itermax<1000 ) {
      
      lambda0=lambda1
      p0= exp(c(lambda0,0))
      p0= p0/sum(p0); 
      emax1=emax*matrix(rep(as.matrix(p0),each=nrow(X)),nrow=nrow(X))
      emax2=apply(emax1,1,sum);   
      emax2= t(matrix(rep(emax2,each=length(p0)),nrow=length(p0)))
      p1= emax/emax2;
      p1= apply(p1, 2, sum)
      p1= uxcount/p1
      p1=p1/sum(p1)    
      lambda1= log(p1[-length(p1)]/p1[length(p1)])
      itermax=itermax+1  
    }
    
    emax= emax1/emax2
  } else
  {
    emax=emax+ matrix(rep(as.matrix(c(lambda,0)),each=nrow(X)),nrow=nrow(X))
    emax[is.infinite(emax)]=exp(700) ## avoid infinite
    emax=exp(emax); 
    emax=t(apply(emax, 1, function(i) i/sum(i))) 
  }
  
  res=0
  for (i in 1:nrow(X)) res=res+log(emax[i,round(X[i,ncol(X)],10)==round(ux,10)])
  return(list(res=res, lambdak=lambda1))  
}

fun.nllSORE <- function(paryx, y, X, W, endox, ux, ux4z, nparx, xstar,stdE){
  res=fun.nllSOREvec(paryx, y, X, W, endox, ux, ux4z, nparx,xstar,stdE)
  return(sum(res))
}

fun.nllSOREvec <- function(paryx, y, X, W, endox, ux, ux4z,nparx, xstar,stdE){
  
  beta=paryx[1:ncol(X)]
  sig <-exp(paryx[length(beta)+1])
  parx <-  vector("list", ncol(endox))
  
  ## convert parameter vectors back to lists
  k=length(beta)+1; 
  for (i in 1:ncol(endox)) {
    parx[[i]] = list(gamma=paryx[(k+1):(k+nparx[[i]]$ngamma)],
                     lambda=paryx[(k+nparx[[i]]$ngamma+1):(k+nparx[[i]]$ngamma+nparx[[i]]$nlambda)] )
    k= k+nparx[[i]]$ngamma+nparx[[i]]$nlambda
  }
  
  eps<- y-X%*%beta;
  eps <- eps[,1]
  res <-  dnorm(eps, sd=sig, log=T)
  if(stdE) eps<-eps/sig; 
  for (i in 1:ncol(endox)) {
    res= res+ llk.SORE(eps=eps, X=endox[,1:i],Z=W,ux=ux[[i]],ux4z=ux4z[[i]], lambda=parx[[i]]$lambda, gamma=parx[[i]]$gamma)
  }
  -res
}

llk.SORE<- function(eps, X, Z, ux, ux4z,  lambda, gamma) {
  
  X<-as.matrix(X)
  if (is.null(Z)) {
    if (ncol(X)==1) emax=as.matrix(eps)%*%as.matrix(gamma)%*%t(as.matrix(ux)) else 
      if (ncol(X)>1) emax=as.matrix(cbind(X[,1:(ncol(X)-1)],eps))%*%as.matrix(gamma)%*%t(as.matrix(ux)) else
        stop("The dimension of X in SOR model is less than 1")
  } else {
    if (ncol(X)==1) {
         emax= Z %*%as.matrix(gamma[-length(gamma)])%*%t(as.matrix(ux4z))+ (as.matrix(eps)*gamma[length(gamma)])%*%t(as.matrix(ux)) ##else
    } else 
      if (ncol(X)>1) emax=Z %*%as.matrix(gamma[1:ncol(Z)])%*%t(as.matrix(ux4z))+ as.matrix(cbind(X[,1:(ncol(X)-1)],eps))%*%as.matrix(gamma[-(1:ncol(Z))])%*%t(as.matrix(ux)) else
        stop("The dimension of X in SOR model is less than 1")
  }
  
  emax=emax+ matrix(rep(as.matrix(c(lambda,0)),each=nrow(X)),nrow=nrow(X)); 
  emax=exp(emax); 
  emax=t(apply(emax, 1, function(i) i/sum(i)))
  res=numeric(nrow(X))
  for (i in 1:nrow(X)) res[i]=log(emax[i,round(X[i,ncol(X)],10)==round(ux,10)])
  return(res)

}


#---------------------------------------------------------------------------------------------------------
#---------------------- Use SORE to fit sample data. ------------------------------------------------
#----------------------------------------------------------------------------------------------------------

#################### Example 1: Small data set ###############################
## Read in the sample data
data<-read.table(file="example1.dat", header=T)

## SORE estimation
results<-SORE(formula=y~x|x, xstar=1, stdE=T, data=data, prof=T)

results$par[3] <- exp(results$par[3])
results$par
-results$obj

## JCM estimation
summary(lm(y~x+xstar, data=data))


############################# Example 2: Big data ####################################
data<-read.table(file="example2.dat", header=T)

## Rounding to reduce the unique values of the endogenous regressors to be modeled in SORE, original x values are used in the outcome model
## Rounding  not only can yield more stabilized estimates in small to moderate samples,
## but also substantially reduce computational time when the endogenous regressor has many unique values. 
## Recommend to have at least 30 unique values remaining after rounding. 
data$xstar = round(data$xstar,1) 

fitSORE<-SORE(formula=y~x|xstar,  data=data, prof=T, tolprof=1e-3)

fitSORE$par[3] <- exp(fitSORE$par[3])
fitSORE$par

############################# Example 3: Correlated x and w ###############################
data<-read.table(file="example3.dat", header=T)

## Rounding to reduce the unique values of the endogenous regressors to reduce computational time.
## no need to round w as W is not modeled
data$xstar = round(data$xstar,1) 


fitSORE<-SORE(formula=y~x+w|xstar, W=~wstar, data=data, prof=T, tolprof=1e-3)

fitSORE$par[4] <- exp(fitSORE$par[4])
fitSORE$par

