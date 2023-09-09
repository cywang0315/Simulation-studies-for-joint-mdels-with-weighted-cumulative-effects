##########Data is generated under the setting sigma_1=10,sigma_2=0.5
##########Fit a joint model with B^{\omega}(t_0, t; \sigma_1)) replaced by  B^U (t_0, t),denote by case 4* in our paper
##########That is, we approximate the flat weighted biomarker level with the uniformly weighted biomarker level
library(snowfall)
library(parallel)

################
sfInit(parallel = TRUE, cpus = 16)
sfLibrary(survival)
sfLibrary(JM)
sfLibrary(joineR)
sfLibrary(dplyr) # near
sfLibrary(statmod)
sfLibrary(progress)
sfLibrary(MASS)
sfLibrary(mvtnorm)
sfLibrary(tensor)
sfSource(here::here("source_case4+.R"))
####################

SIMULATE=function(s){
  Q.2=matrix(0,nrow=q,ncol=q-2)
  for(l in 1:(q-2)){
    Q.2[l,l]=6/((knots[l+4]-knots[l+2])*(knots[l+4]-knots[l+1]))
    Q.2[l+2,l]=6/((knots[l+4]-knots[l+2])*(knots[l+5]-knots[l+2]))
    Q.2[l+1,l]=-(Q.2[l,l]+Q.2[l+2,l])
  }
  
  data=try(gendat(s,obstim=seq(0,20/2,by=1/2),obsmax=10,gammatrue=-2,alpha1true=0.3,alpha2true=0.5,betatrue=c(6,3,7,1,8,5,4),D0=diag(c(3,4,4,5,4,3,4)),sigmatrue=1,
              knots,Q.2,sig1true=10,sig2true=0.5))

  if ('try-error' %in% class(data)){
    return(numeric(39)+1)
  }
  M=as.data.frame(mycubicbs(data$time))
  names(M)=paste0("time",c(1:7))
  data=cbind(data,M)
  data.id=data[!duplicated(data$id),]
  
  initialvalue=try(inival(data,data.id,ctl=lmeControl (msMaxIter=100),knots,Q.2,low=0.01,up=15))
  if ('try-error' %in% class(initialvalue)) {
    return(numeric(39)+2)
  }
  
  beta=initialvalue$beta
  sigma2=initialvalue$sigma2
  D=initialvalue$D
  gamma=initialvalue$gamma
  cumbase=initialvalue$cumbase
  alpha1=initialvalue$alpha1
  alpha2=initialvalue$alpha2
  coxts=initialvalue$coxts
  res.ts=c(coxts,beta,sigma2,diag(D))
  sig2=initialvalue$sig2
  
  loglik=try(logLik(data,data.id,gamma,alpha1,alpha2,beta,sigma2,D,cumbase,knots,Q.2,sig2,L=2000))
  if ('try-error' %in% class(loglik)) {
    return(numeric(39)+3)
  }
  
  for (k in 1:10){
    ### UPDATE PARAMETERS VIA EM ALGORITHM
    res.EM=try(est.EM(loglik,data,data.id,gamma,alpha1,alpha2,cumbase,beta,sigma2,D,knots,Q.2,sig2))
    if ('try-error' %in% class(res.EM)) {
      return(numeric(39)+4)
    }
    if(is.null(res.EM)){
      res.jm=c(gamma,alpha1,alpha2,beta,sigma2,diag(D),sig2)
      break
    }
    
    gamma=res.EM$gamma
    alpha1=res.EM$alpha1
    alpha2=res.EM$alpha2
    cumbase=res.EM$cumbase
    beta=res.EM$beta
    sigma2=res.EM$sigma2
    D=res.EM$D
    loglik=res.EM$logLik
    
    ### UPDATE  SIG2 
    res.sig2=try(est.sig2(data,data.id,gamma,alpha1,alpha2,cumbase,beta,sigma2,D,knots,Q.2,sig2))
    if ('try-error' %in% class(res.sig2)) {
      return(numeric(39)+5)
    }

    sig2new=res.sig2$sig2
    flag=res.sig2$flag
    if(flag){
      res.jm=c(gamma,alpha1,alpha2,beta,sigma2,diag(D),sig2new)
      break
    }
    
    sig2=sig2new
  }
  
  if(k==10){
    return(numeric(39)+6)
  }else{
    return(c(k,c(res.jm),c(res.ts)))
  }
  
}

RES=sfLapply(1:20,SIMULATE)
sfStop()
