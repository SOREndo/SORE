
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
  ## hess: =T if computing standard error by inverting hessian matrix
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
    
    return(list(par=fitSORE$par, se=sqrt(diag(solve(fitHess))), loglik=-fitSORE$obj, convergence=fitSORE$converg, iterations=fitSORE$iterations, evaluations=fitSORE$evaluations))
  }
  else {fitSORE$loglik=-fitSORE$obj; return(fitSORE)} 

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
## (1) formula=y~x|x below specifies a structural model y=intercept+ a*x+Error, 
##     where x is endogenous and so is specified after the vertical bar "|".
## (2) the function argument xstar specifies a numeric vector denoting the positions of endogenous regressors that need normal score transformation for GC 
##            regressor-error dependence structure. In this example, xstar=1 means the first endogenous regressor follows GC
##            dependence structure. 
##  (3) stdE=T means the structural error term is standardized by its standard deviation. The two possible values (T and F) of stdE resut in 
##            reparameterization of the same SORE model, and simply scale  the log-odds-ratio parameter for regressor-error dependence 
##            differently in relative to the structural error standard deviation.
##  (4) prof=T means the estimation is done via maximizing the profile likelihood that eliminates the nuisance parameters in baseline functions. 
##          =F means that estimation is done via maximizing the regular likelihood. 
results<-SORE(formula=y~x|x, xstar=1, stdE=T, data=data, prof=T)

## Because the estimation routine minimizes the negative likelihood function over unconstrained parameter space, we used the log
## transformation for the error standard deviation. The step below transform this log error standard deviation back to the original
## error standard deviation. 
results$par[3] <- exp(results$par[3])
## The step below prints out the parameter estimate for intercept, coefficient for x, error standard deviation and the log odds ratio parameter \gamma
results$par  
## This step prints out the maximized log likelihood for the SORE model
results$loglik

## JCM estimation as the benchmark analysis
summary(lm(y~x+xstar, data=data))


############################# Example 2: Big data ####################################
data<-read.table(file="example2.dat", header=T)

## Rounding to reduce the unique values of the endogenous regressors to be modeled in SORE, original x values are used in the outcome model
## Rounding  not only can yield more stabilized estimates in small to moderate samples,
## but also substantially reduce computational time when the endogenous regressor has many unique values. 
## Recommend to have at least 30 unique values remaining after rounding. 
data$xstar = round(data$xstar,1) 

## (1) data contains three variables: y (the outcome), x (the endogenous regressor) and xstar (the GC copula transformation of x)
## (2) formula=y~x|xstar below specifies the structural model as y=intercept+a*x+Error, 
##     where the regressor x is endogenous, and its latent copula data, xstar, is specified after the vertical bar "|".
##     Including the GC copula (normal score) transformation of an endogenous regressor after the vertical bar "|" permits modeling GC regressor-error dependence
##     and requires creating this variable and including it in data. 
##     Alternatively, one can specify the GC regressor-error dependence without creating and including the normal score transformation variable in the data.
##     See Example 1 above for an example of this alternative way to specify GC regressor-error dependence, which makes use of the "xstar" argument in the SORE() function. 
##     See Web Appendix Table W3 for more details.
## (3) tolprof specifies the threshold value used to determine convergence  of nuisance baseline parameters within the profile likelihood algorithm (Step 2.c in Web Appendix A.2). 
##     Default=1e-8. We find a larger threhold value is often sufficient and results in essentially the same converged value.         
fitSORE<-SORE(formula=y~x|xstar,  data=data, prof=T, tolprof=1e-3)

fitSORE$par[3] <- exp(fitSORE$par[3])
fitSORE$par

############################# Example 3: Correlated x and w ###############################
data<-read.table(file="example3.dat", header=T)


## (1) formula=y~x+w|x below specifies the structural model as y=intercept+a*x+b*w+Error,  where the regressor x is endogenous
## (2) W=~w specifies the exogenous variable w is used to model the endogenous regressor.
## (3) xstar=1 specifies that the first endogenous regressor and the error term jointly follow a GC dependence. xstar=NULL (default) yields LB regressor-error dependence. 
## (4) RR.GC=T specifies that the endogenous regressor and the exogenous regressors follows a GC dependence. RR.GC=F yields LB regressor-regressor dependence. 
fitSORE<-SORE(formula=y~x+w|x, W=~w, xstar=1, RR.GC=T, data=data, prof=T, tolprof=1e-3)

## print out model estimates
fitSORE$par[4] <- exp(fitSORE$par[4])
fitSORE$par
## print out maximized model likelihood
fitSORE$loglik

############################# Example 4: Binary Endogenous Regressors #####################################
data<-read.table(file="example4.dat", header=T)

## OLS analysis
ols<- lm(y~x, data=data); summary(ols)

## SORE estimation
## To safeguard against the possibility of multiple local maxima in likelihood function and ensure locating the global maximum, a common practice is to run the estimation algorithm
## from different starting values and compare the maximized likelihood(s) obtained from these starting values. The code below sets starting values of structural model
## parameters as the OLS estimates and a grid of values for the endogeneity log-OR parameter gamma0
gamma0<- 10 ## the gamma0 value is the starting value value for the log OR endogeneity parameer gamma and   can be varied in a grid of  values, e.g., from 10 to -10. 
par0<-c(coef(summary(ols))[,1], log(summary(ols)$sigma), gamma0)
## (1) formula=y~x|x below specifies the structural model as y=intercept+a*x+Error, 
##     where the regressor x is endogenous, and  is specified after the vertical bar "|". This specifies an LB regressor-error dependence since xstar=NULL by default
## (2) parstart=par0 specifies the starting values of estimation algorithm as the par0 vector. The OLS estimate and gamma=0 are set as the default starting values. 
fitSORE<-SORE(formula=y~x|x, data=data, prof=T, tolprof=1e-3, hess=T, parstart=par0)

## print the parameter estimates for intercept, coefficient for x, error standard deviation and the log odds ratio parameter \gamma
fitSORE$par[3] <- exp(fitSORE$par[3])
fitSORE$par
## Print the likelihood 
fitSORE$loglik

