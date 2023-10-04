m=1000
q=7
internal_knots=c(5,10,15)/2
boundary_knots=c(0,20)/2
out_knots=c(seq(from=min(boundary_knots),by=-0.1,length=4),seq(from=max(boundary_knots),by=0.1,length=4))
knots=sort(c(internal_knots,out_knots))

####### B-spline basis #########
b_sp=function(s,k,l){
  if(l==1){
    return(ifelse((s<knots[k+1])&(s>=knots[k]),1,0))
  }else{
    (s-knots[k])/(knots[k+l-1]-knots[k])*b_sp(s,k,l-1)+(knots[k+l]-s)/(knots[k+l]-knots[k+1])*b_sp(s,k+1,l-1)
  }
  
}

## my cubic B-spline basis functions ####
mycubicbs=function(s){
  mycubicbs1_t=function(s){
    return(sapply(1:(length(internal_knots)+3+1),b_sp,s=s,l=4))
  }
  return(t(sapply(s,mycubicbs1_t)))
}



#### truncated [t1,t2] normal density as weight function
Nor_wei=function(s,sigma,t1,t2){
  dnorm(s,sd=sigma)/(pnorm(t2,sd=sigma)-pnorm(t1,sd=sigma))
}



#################### Calucate  \int_{t_0}^{t} w(t-s)B_{k,4}(s) ds for k=1,...,q 
cumwei_cubicbs=function(t0,tt,sigma){ # tt can be a vector
  if(any(t0>tt)){
    return("error: upper limit smaller than lower limit")
  }else{
    cumwei_cubicbs_sing=function(k,l,t){# t is a scalar, k is a scalar
      if(l==1){
        ifelse((knots[k+1]<t0)|(knots[k]>t),0,(pnorm(t-max(t0,knots[k]),sd=sigma)-pnorm(t-min(t,knots[k+1]),sd=sigma))/(pnorm(t-t0,sd=sigma)-pnorm(0,sd=sigma)))
      }else{
        if(l==2){
          r1=ifelse((knots[k+1]<t0)|(knots[k]>t),0,
                    (t-knots[k])/(knots[k+1]-knots[k])*(pnorm(t-max(t0,knots[k]),sd=sigma)-pnorm(t-min(t,knots[k+1]),sd=sigma))/(pnorm(t-t0,sd=sigma)-pnorm(0,sd=sigma))+
                      sigma^2/(knots[k+1]-knots[k])*(Nor_wei(t-max(t0,knots[k]),sigma,0,t-t0)-Nor_wei(t-min(t,knots[k+1]),sigma,0,t-t0)))
          r2=ifelse((knots[k+2]<t0)|(knots[k+1]>t),0,
                    (knots[k+2]-t)/(knots[k+2]-knots[k+1])*(pnorm(t-max(t0,knots[k+1]),sd=sigma)-pnorm(t-min(t,knots[k+2]),sd=sigma))/(pnorm(t-t0,sd=sigma)-pnorm(0,sd=sigma))-
                      sigma^2/(knots[k+2]-knots[k+1])*(Nor_wei(t-max(t0,knots[k+1]),sigma,0,t-t0)-Nor_wei(t-min(t,knots[k+2]),sigma,0,t-t0)))
          return(r1+r2)
          
        }else{
          (t-knots[k])/(knots[k+l-1]-knots[k])*cumwei_cubicbs_sing(k,l-1,t)+
            (knots[k+l]-t)/(knots[k+l]-knots[k+1])*cumwei_cubicbs_sing(k+1,l-1,t)-
            sigma^2*(Nor_wei(0,sigma,0,t-t0)*b_sp(t,k,l-1)-Nor_wei(t-t0,sigma,0,t-t0)*b_sp(t0,k,l-1))/(knots[k+l-1]-knots[k])+
            sigma^2*(Nor_wei(0,sigma,0,t-t0)*b_sp(t,k+1,l-1)-Nor_wei(t-t0,sigma,0,t-t0)*b_sp(t0,k+1,l-1))/(knots[k+l]-knots[k+1])+
            sigma^2*(l-2)/(knots[k+l-1]-knots[k])*(cumwei_cubicbs_sing(k,l-2,t)/(knots[k+l-2]-knots[k])-cumwei_cubicbs_sing(k+1,l-2,t)/(knots[k+l-1]-knots[k+1]))-
            sigma^2*(l-2)/(knots[k+l]-knots[k+1])*(cumwei_cubicbs_sing(k+1,l-2,t)/(knots[k+l-1]-knots[k+1])-cumwei_cubicbs_sing(k+2,l-2,t)/(knots[k+l]-knots[k+2]))
          
        }
      }
      
    }
    cumwei_cubicbs_t=function(t){
      return(sapply(1:(length(internal_knots)+3+1),cumwei_cubicbs_sing,l=4,t=t))
    }
    
    return(t(sapply(tt,cumwei_cubicbs_t)))
  }
  
}



cumwei_R.2=function(t0,t,sigma){
  cumwei_R.2=matrix(0,ncol=q-2,nrow=q-2)
  for (l in 1:(q-2)){
    rr1=(pnorm(t-max(t0,knots[l+2]),sd=sigma)-pnorm(t-min(t,knots[l+3]),sd=sigma))/(pnorm(t-t0,sd=sigma)-pnorm(0,sd=sigma))
    
    r1=ifelse(t<knots[l+2]|t0>knots[l+3],0,1/(knots[l+3]-knots[l+2])^2*
                (sigma^2*(Nor_wei(t-min(t,knots[l+3]),sigma,0,t-t0)*(t-min(t,knots[l+3]))-Nor_wei(t-max(t0,knots[l+2]),sigma,0,t-t0)*(t-max(t0,knots[l+2])))+
                   (sigma^2+(t-knots[l+2])^2)*rr1-
                   2*sigma^2*(t-knots[l+2])*(Nor_wei(t-min(t,knots[l+3]),sigma,0,t-t0)-Nor_wei(t-max(t0,knots[l+2]),sigma,0,t-t0))))
    
    rr2=(pnorm(t-max(t0,knots[l+3]),sd=sigma)-pnorm(t-min(t,knots[l+4]),sd=sigma))/(pnorm(t-t0,sd=sigma)-pnorm(0,sd=sigma))
    r2=ifelse(t<knots[l+3]|t0>knots[l+4],0,1/(knots[l+4]-knots[l+3])^2*
                (sigma^2*(Nor_wei(t-min(t,knots[l+4]),sigma,0,t-t0)*(t-min(t,knots[l+4]))-Nor_wei(t-max(t0,knots[l+3]),sigma,0,t-t0)*(t-max(t0,knots[l+3])))+
                   (sigma^2+(t-knots[l+4])^2)*rr2-
                   2*sigma^2*(t-knots[l+4])*(Nor_wei(t-min(t,knots[l+4]),sigma,0,t-t0)-Nor_wei(t-max(t0,knots[l+3]),sigma,0,t-t0))))
    cumwei_R.2[l,l]=r1+r2
    if(l<(q-2)){
      cumwei_R.2[l,l+1]=ifelse(t<knots[l+3]|t0>knots[l+4],0,1/(knots[l+4]-knots[l+3])^2*
                                 (-sigma^2*(Nor_wei(t-min(t,knots[l+4]),sigma,0,t-t0)*(t-min(t,knots[l+4]))-Nor_wei(t-max(t0,knots[l+3]),sigma,0,t-t0)*(t-max(t0,knots[l+3])))+
                                    ((knots[l+4]-t)*(t-knots[l+3])-sigma^2)*rr2+
                                    sigma^2*(2*t-knots[l+3]-knots[l+4])*(Nor_wei(t-min(t,knots[l+4]),sigma,0,t-t0)-Nor_wei(t-max(t0,knots[l+3]),sigma,0,t-t0))))
      
      cumwei_R.2[l+1,l]=cumwei_R.2[l,l+1]
      
    }
  }
  return(cumwei_R.2)
}



R.2=function(t){
  R.2=matrix(0,ncol=q-2,nrow=q-2)
  for (l in 1:(q-2)){
    if((t>knots[l+2])&(t<=knots[l+3])){
      R.2[l,l]=(t-knots[l+2])^3/(3*(knots[l+3]-knots[l+2])^2)
    }else{
      if((t>knots[l+3])&(t<knots[l+4])){
        R.2[l,l]=(knots[l+4]-knots[l+2])/3-(knots[l+4]-t)^3/(3*(knots[l+4]-knots[l+3])^2)
        R.2[l,l+1]=(-t^3/3+(knots[l+4]+knots[l+3])/2*t^2-knots[l+4]*knots[l+3]*t-(knots[l+3])^3/6+
                      (knots[l+3])^2*knots[l+4]/2)/(knots[l+4]-knots[l+3])^2
        R.2[l+1,l]= R.2[l,l+1]
      }else{
        if(t>=knots[l+4]){
          R.2[l,l]=(knots[l+4]-knots[l+2])/3
          R.2[l,l+1]=(knots[l+4]-knots[l+3])/6
          R.2[l+1,l]=(knots[l+4]-knots[l+3])/6
          
        }
      }
    }
  }
  R.2[1,1]=ifelse(t>=knots[5],(knots[5]-knots[4])/3,(knots[5]-knots[4])/3-(knots[5]-t)^3/(3*(knots[5]-knots[4])^2))
  return(R.2)
}

# true value of longitudinal biomarker (no measurement error)
longit.true=function(t,fix_eff, ran_eff){
  desmat=mycubicbs(t)# design matrix
  desmat%*%(fix_eff+ran_eff)
}




adv=function(x){
  grad=diag(1,length(x))
  grad2=matrix(0,ncol=length(x),nrow = length(x))
  attr(x,"grad")=grad
  attr(x,"grad2")=grad2
  class(x)="adv"
  x
}


quad.adv=function(A,b){
  grad.b=attr(b,"grad")
  grad2.b=attr(b,"grad2")
  b=as.numeric(b)
  d=c(t(b)%*%A%*%b) # scalar
  attr(d,"grad")=grad.b%*%(2*A%*%b)
  attr(d,"grad2")=grad.b%*%(2*A)
  class(d)="adv"
  d
}


"^.adv"=function(b,alpha){ # the input is a scalar function of b
  grad.b=attr(b,"grad")
  grad2.b=attr(b,"grad2")
  b=as.numeric(b)
  d=b^alpha
  attr(d,"grad")=alpha*b^(alpha-1)*grad.b
  attr(d,"grad2")=alpha*(alpha-1)*b^(alpha-2)*tcrossprod(grad.b)+alpha*b^(alpha-1)*grad2.b
  class(d)="adv"
  d
}


linear_sca.adv=function(A,b,j){ ### derivative of <Ab>j with respect to b, A is matrix independent of b
  grad.b=attr(b,"grad")
  b=as.numeric(b) ## transform to a vector
  d=c(A%*%b)[j]
  K=numeric(nrow(A))
  K[j]=1
  attr(d,"grad")=t(A)%*%K
  attr(d,"grad2")=matrix(0,ncol=length(b),nrow=length(b))
  class(d)="adv"
  d
}

linear.adv=function(A,b){## derivative of Ab w.r.t b, A is a matrix or a column vector
  grad.b=attr(b,"grad")
  b=as.numeric(b)
  check=length(A)==length(b)
  if(check){
    d=t(A)%*%b
    attr(d,"grad")=A
    attr(d,"grad2")=matrix(0,ncol=length(b),nrow = length(b))
  }else{
    d=A%*%b
    attr(d,"grad")=A
  }
  
  class(d)="adv"
  d
}

exp.adv=function(b){ ## the input is a scalar function of b
  grad.b=attr(b,"grad")
  grad2.b=attr(b,"grad2")
  b=as.numeric(b)
  d=exp(b)
  attr(d,"grad")=exp(b)*grad.b
  attr(d,"grad2")=exp(b)*tcrossprod(grad.b)+exp(b)*grad2.b
  class(d)="adv"
  d
}



"*.adv"=function(a1,a2){#derivative of a1(b)*a2(b) w.r.t. b, where a1 is a scalar function of b and a2 is a scalar or vector function of b.
  
  if(is.null(attr(a1,"grad"))){
    grad.a2=attr(a2,"grad")
    grad2.a2=attr(a2,"grad2")
    a2=as.numeric(a2)
    d=a1*a2
    attr(d,"grad")=a1*grad.a2
    attr(d,"grad2")=a1*grad2.a2
    
  }else{
    grad.a1=attr(a1,"grad")
    grad2.a1=attr(a1,"grad2")
    grad.a2=attr(a2,"grad")
    grad2.a2=attr(a2,"grad2")
    a1=as.numeric(a1)
    a2=as.numeric(a2)
    d=a1*a2
    check=length(a2)==1
    if(check){
      attr(d,"grad")=grad.a1*a2+a1*grad.a2
      attr(d,"grad2")=grad.a1%*%t(grad.a2)+grad.a2%*%t(grad.a1)+a2*grad2.a1+a1*grad2.a2
    }else{
      attr(d,"grad")=a2%*%t(grad.a1)+a1*grad.a2
    }
  }
  
  class(d)="adv"
  d
}


"+.adv"=function(a,b){ #derivative of a(b)+b(b) where a and b are vector or scalar functions of b, and a(.) can be not function of b
  if(is.null(attr(a,"grad"))){
    grad.b=attr(b,"grad")
    grad2.b=attr(b,"grad2")
    b=as.numeric(b)
    d=a+b
    attr(d,"grad")=grad.b
    attr(d,"grad2")=grad2.b
    
  }else{
    grad.a=attr(a,"grad")
    grad2.a=attr(a,"grad2")
    
    grad.b=attr(b,"grad")
    grad2.b=attr(b,"grad2")
    
    a=as.numeric(a)
    b=as.numeric(b)
    d=a+b
    attr(d,"grad")=grad.a+grad.b
    attr(d,"grad2")=grad2.a+grad2.b
    
  }
  
  class(d)="adv"
  d
}

##################Generate data 

gendat=function(s,obstim,obsmax,gammatrue,alpha1true,alpha2true,betatrue,D0,sigmatrue,knots,Q.2,sig1true,sig2true){
  set.seed(s)
  basetrue=function(t){ifelse(t<10/2,0,exp(-2))} ## 0 when t<5
  b=mvrnorm(m,mu=rep(0,q),Sigma=D0) # random effects in longitudinal submodel
  W=sample(c(0,1),size=m,replace = TRUE)
  desmat=cumwei_cubicbs(1,seq(10/2,obsmax,by=0.005),sig1true)# design matrix  
  #set.seed(1)
  TD.fun=function(i,t){# i is a scale, t can be a vector
    tt=function(t){
      c(basetrue(t)*exp(gammatrue*W[i]+alpha1true*desmat[match(t,seq(10/2,obsmax,by=0.005)),]%*%(betatrue+b[i,])+
                          alpha2true*sqrt(t(betatrue+b[i,])%*%Q.2%*%cumwei_R.2(1,t,sig2true)%*%t(Q.2)%*%(betatrue+b[i,]))))
    }
    #hazard=sapply(t,tt)
    #p=1-exp(-hazard*0.01)
    r=rbinom(n=length(t),size=1,prob=1-exp(-sapply(t,tt)*0.005)) ##prob: event occurs
    Time1=obsmax
    delta1=0
    if(max(r)>0){
      Time1=t[min(which(r==1))]
      delta1=1
    }
    return(list(Time1=Time1,delta1=delta1))
  }
  
  TD=sapply(c(1:m),TD.fun,t=seq(10/2,obsmax,by=0.005))
  Time1=as.numeric(TD[1,])
  delta1=as.numeric(TD[2,])
  ####generate censoring 
  #censtim=sample(seq(12,20,by=0.001),size=m,replace = TRUE)
  censtim=sample(seq(12/2,30/2,by=0.001),size=m,replace = TRUE)
  delta=ifelse(Time1<=censtim,delta1,0)
  sum(delta)
  Time=ifelse(Time1<=censtim,Time1,censtim)
  #############generate longitudinal data
  #set.seed(1)
  Y=c()
  l=numeric(m)
  time=c()
  for (i in 1:m){
    l[i]=floor(2*Time[i])+1
    for(j in seq(0,Time[i],by=0.5)){
      YY=longit.true(j,fix_eff=betatrue,ran_eff=b[i,])+rnorm(1,sd=sigmatrue)
      Y=rbind(Y,YY)
      time=rbind(time,j)
    }
  }
  
  Y=as.vector(Y)
  time=as.vector(time)
  ### construct data frame. Note here we use "id" not "sub"
  data=data.frame(id=rep(c(1:m),l),Time=rep(Time,l),delta=rep(delta,l),Y=Y,W=rep(W,l),
                  time=time, start=time,stop=time+1/2,event=numeric(length(Y))) ##should adjust "stop" and "event"
  
  data$stop=ifelse(data$stop<=data$Time,data$stop,data$Time)
  data$event=ifelse((data$stop==data$Time)&(data$delta==1),1,data$event)
  return(data)
  
}

############################# Get initial value ##################
inival=function(data,data.id,ctl,knots,Q.2,low,up){
  lme.data=lme(fixed=Y~-1+time1+time2+time3+time4+time5+time6+time7,
               random=list(id=pdDiag(form=~-1+time1+time2+time3+time4+time5+time6+time7)),data = data,control=ctl)
  cox.data=coxph(Surv(Time,delta)~W,data=data.id,x=TRUE)
  beta=lme.data$coefficients[[1]]
  sigma2=lme.data$sigma^2
  D=diag(c(as.numeric(VarCorr(lme.data)[1]),as.numeric(VarCorr(lme.data)[2]),as.numeric(VarCorr(lme.data)[3]),
           as.numeric(VarCorr(lme.data)[4]),as.numeric(VarCorr(lme.data)[5]),as.numeric(VarCorr(lme.data)[6]),as.numeric(VarCorr(lme.data)[7])))
  

  conloglik_3=function(par){
    gamma=par[1]
    alpha1=par[2]
    alpha2=par[3]
    conloglik=0
    for(i in data.id$id[data.id$delta==1]){
      Ti=data.id$Time[data.id$id==i]
      rs.id=data.id$id[data.id$Time>=Ti]
      conloglik=conloglik+gamma*data.id$W[i]+alpha1*sum(c(cumwei_cubicbs(1,Ti,sig1))*(beta+t(random.effects(lme.data)[i,])))+
        alpha2*sqrt(1/(Ti-1)*sum((beta+t(random.effects(lme.data)[i,]))*(Q.2%*%(R.2(Ti)-R.2(1))%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[i,])))))-
        log(sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(cumwei_cubicbs(1,Ti,sig1))*(beta+t(random.effects(lme.data)[rs.id,])))+
                      alpha2*sqrt(1/(Ti-1)*colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%(R.2(Ti)-R.2(1))%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))))
    }
    return(-conloglik)
    
  }
  
  
  gg_conloglik_3=function(par){
    gamma=par[1]
    alpha1=par[2]
    alpha2=par[3]
    gg_gamma=0
    gg_alpha1=0
    gg_alpha2=0
    i_alpha1=0
    i_alpha2=0
    for(i in data.id$id[data.id$delta==1]){
      Ti=data.id$Time[data.id$id==i]
      rs.id=data.id$id[data.id$Time>=Ti]
      s=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(cumwei_cubicbs(1,Ti,sig1))*(beta+t(random.effects(lme.data)[rs.id,])))+
                  alpha2*sqrt(1/(Ti-1)*colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%(R.2(Ti)-R.2(1))%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,])))))))
      s1=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(cumwei_cubicbs(1,Ti,sig1))*(beta+t(random.effects(lme.data)[rs.id,])))+
                   alpha2*sqrt(1/(Ti-1)*colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%(R.2(Ti)-R.2(1))%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))*data.id$W[rs.id])
      
      s2=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(cumwei_cubicbs(1,Ti,sig1))*(beta+t(random.effects(lme.data)[rs.id,])))+
                   alpha2*sqrt(1/(Ti-1)*colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%(R.2(Ti)-R.2(1))%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))*colSums(c(cumwei_cubicbs(1,Ti,sig1))*(beta+t(random.effects(lme.data)[rs.id,]))))
      
      
      s3=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(cumwei_cubicbs(1,Ti,sig1))*(beta+t(random.effects(lme.data)[rs.id,])))+
                   alpha2*sqrt(1/(Ti-1)*colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%(R.2(Ti)-R.2(1))%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))*sqrt(1/(Ti-1)*colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%(R.2(Ti)-R.2(1))%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))
      
      
      gg_gamma=gg_gamma+data.id$W[i]-s1/s
      gg_alpha1=gg_alpha1+sum(c(cumwei_cubicbs(1,Ti,sig1))*(beta+t(random.effects(lme.data)[i,])))-s2/s

      gg_alpha2=gg_alpha2+sqrt(1/(Ti-1)*sum((beta+t(random.effects(lme.data)[i,]))*(Q.2%*%(R.2(Ti)-R.2(1))%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[i,])))))-s3/s

      
    }
    return(c(-gg_gamma,-gg_alpha1,-gg_alpha2))
    #return(c(-gg_alpha1,-gg_alpha2))
    
  }
  
  conloglik_sig1=function(sig1){
    sig1=sig1
    conloglik=0
    for(i in data.id$id[data.id$delta==1]){
      Ti=data.id$Time[data.id$id==i]
      rs.id=data.id$id[data.id$Time>=Ti]
      conloglik=conloglik+gamma*data.id$W[i]+alpha1*sum(c(cumwei_cubicbs(1,Ti,sig1))*(beta+t(random.effects(lme.data)[i,])))+
        alpha2*sqrt(1/(Ti-1)*sum((beta+t(random.effects(lme.data)[i,]))*(Q.2%*%(R.2(Ti)-R.2(1))%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[i,])))))-
        log(sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(cumwei_cubicbs(1,Ti,sig1))*(beta+t(random.effects(lme.data)[rs.id,])))+
                      alpha2*sqrt(1/(Ti-1)*colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%(R.2(Ti)-R.2(1))%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,]))))))))
    }
    return(-conloglik)
    
  }
  
  
  gamma=0
  alpha1=0
  alpha2=0
  sig1_set=data.frame(sig1=c(0.1,2,5,10))

  for(k in 1:dim(sig1_set)[1]){
    sig1=sig1_set$sig1[k]
    res=optim(c(gamma,alpha1,alpha2),conloglik_3,gg_conloglik_3,method = "BFGS")
    sig1_set$value[k]=res$value
    sig1_set$gamma[k]=res$par[1]
    sig1_set$alpha1[k]=res$par[2]
    sig1_set$alpha2[k]=res$par[3]
  }
  
  sig1=sig1_set[which.min(sig1_set$value),]$sig1
  #gamma=sig12_set[which.min(sig12_set$value),]$gamma
  #alpha1=sig12_set[which.min(sig12_set$value),]$alpha1
  #alpha2=sig12_set[which.min(sig12_set$value),]$alpha2
  
   for(k in 1:10){
    
    res=optim(c(gamma,alpha1,alpha2),conloglik_3,gg_conloglik_3,method = "BFGS")
    gamma=res$par[1]
    alpha1=res$par[2]
    alpha2=res$par[3]
    sig1_res=optim(c(sig1),conloglik_sig1,lower=c(0.01),upper = c(15),method = "L-BFGS-B")
    sig1=sig1_res$par
    if(res$value<sig1_res$value){
      break
    }else{
      if(abs(res$value-sig1_res$value)/res$value<10^(-6)){
        break
      }
    }
  }

  coxts=c(gamma,alpha1,alpha2,sig1)
  
  ### CALCULATE CUMULATIVE HAZARD
  riskset=function(t){# individuals in the risk set
    unique(data$id[(data$Time>=t)])## "=" is important
  }
  
  cc=c()
  for(t in sort(unique(data$Time[data$delta==1]))){
    rs.id=riskset(t)
    s=sum(exp(gamma*data.id$W[rs.id]+alpha1*colSums(c(cumwei_cubicbs(1,t,sig1))*(beta+t(random.effects(lme.data)[rs.id,])))+
                alpha2*sqrt(1/(t-1)*colSums((beta+t(random.effects(lme.data)[rs.id,]))*(Q.2%*%(R.2(t)-R.2(1))%*%t(Q.2)%*%(beta+t(random.effects(lme.data)[rs.id,])))))))
    
    cc=c(cc,length(data.id$id[(data.id$Time==t)&(data.id$delta==1)])/s)
  }
  
  cumbase=data.frame(hazard=0,time=c(0,sort(unique(data$Time[data$delta==1]))))
  cumbase[,1]=cumsum(c(0,cc))
  
  
  out=list(beta=beta,sigma2=sigma2,D=D,gamma=gamma,alpha1=alpha1,alpha2=alpha2,cumbase=cumbase,coxts=coxts,sig1=sig1)
  
}

###calculate log likelihood 

logLik=function(data,data.id,gamma,alpha1,alpha2,beta,sigma2,D,cumbase,knots,Q.2,sig1,L){
  
  lambda0=function(t){
    if(t %in% unique(data$Time[data$delta==1])){
      x=cumbase[which(near(cumbase[,2],t)),1]-cumbase[which(near(cumbase[,2],t))-1,1] ##warn: options(digits = 22) default 7
      if(length(x)>1){
        x=x[1]
      }
      
    }else{
      x=0}
    return(x)
  }
  
  
  des.Y=mycubicbs(data$time)
  des.T=cumwei_cubicbs(1,data.id$Time[data.id$delta==1],sig1)
  Time=data.id$Time
  delta=data.id$delta
  W=data.id$W
  K.2=array(c(sapply(data.id$Time[data.id$delta==1],function(t) 1/(t-1)*Q.2%*%(R.2(t)-R.2(1))%*%t(Q.2))),
            dim=c(q,q,length(data.id$Time[data.id$delta==1]))) 
  l=as.data.frame(table(data$id))$Freq
  m=dim(data.id)[1]
  
  #set.seed(1)
  log.p.Y=numeric(m)
  log.s_set=matrix(0,nrow=m,ncol=L)
  ZTimeb=matrix(0,nrow=sum(delta==1),ncol=L) ##Z_Time%*%b
  bKb=matrix(0,nrow=sum(delta==1),ncol=L)
  rmc=rmvnorm(L,mean=rep(0,q))
  for(i in data.id$id){
    Zi=des.Y[data$id==i,]
    Xi=Zi
    Yi=data$Y[data$id==i]
    
    chol_invSig=chol(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
    Q=solve(chol_invSig)
    
    mu=c(Q%*%t(Q)%*%(t(Zi)/c(sigma2))%*%(Yi-Xi%*%c(beta))) #E(bi|Yi,theta),initial value of bi
    
    Sigma=solve(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
    
    log.p.Y[data.id$id==i]=dmvnorm(Yi,mean=c(Xi%*%c(beta)),sigma=sigma2*diag(length(Yi))+Zi%*%D%*%t(Zi),log=TRUE)
    rmci=t(mu+t(chol(Sigma))%*%t(rmc))## matrix with ncol=q,nrow=L, MC sample from bi|yi
    if(delta[data.id$id==i]){
      ZTimeb[match(i,data.id$id[data.id$delta==1]),]=des.T[match(i,data.id$id[data.id$delta==1]),]%*%t(rmci)
      bKb[match(i,data.id$id[data.id$delta==1]),]=colSums((t(rmci)+beta)*(K.2[,,match(Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(t(rmci)+beta)))
    }
    
    a=unique(data.id$Time[(data.id$Time<=Time[data.id$id==i])&(data.id$delta==1)])
    if(length(a)>0){
      if(length(a)==1){
        btKb=colSums((t(rmci)+beta)*(K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(t(rmci)+beta)))
      }else{
        btKb=apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) colSums((t(rmci)+beta)*(y%*%(t(rmci)+beta))))
      }
      
      log.s_set[data.id$id==i,]=apply(sapply(a,lambda0)*exp(c(gamma*W[data.id$id==i]))*
                                        exp(alpha1*(des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(t(rmci)+beta))+
                                              alpha2*t(sqrt(btKb))),2,sum)
      
    }
    
  }
  
  log.hazard=matrix(0,ncol=L,nrow=m)
  log.hazard[which(delta==1),]=log(sapply(Time[delta==1],lambda0))+c(W[delta==1]*gamma)+alpha1*(c(des.T%*%beta)+ZTimeb)+
    alpha2*sqrt(bKb)
  log.p.Tb=log.hazard-log.s_set
  p.Tb=exp(log.p.Tb)
  logLikmc=sum(log.p.Y+log(rowMeans(p.Tb)))
  return(logLikmc)              
  
}


###############Estimation##########

est.EM=function(loglik,data,data.id,gamma,alpha1,alpha2,cumbase,beta,sigma2,D,knots,Q.2,sig1){
  
  lambda0=function(t){
    if(t %in% unique(data$Time[data$delta==1])){
      x=cumbase[which(near(cumbase[,2],t)),1]-cumbase[which(near(cumbase[,2],t))-1,1] ##warn: options(digits = 22) default 7
      if(length(x)>1){
        x=x[1]
      }
      
    }else{
      x=0}
    return(x)
  }
  
  riskset=function(t){# individuals in the risk set
    unique(data$id[(data$Time>=t)])## "=" is important
  }
  
  Time=data$Time[!duplicated(data$id)]
  delta=data$delta[!duplicated(data$id)]
  W=data$W[!duplicated(data$id)]
  des.Y=as.matrix(data[,c(paste0("time",c(1:q)))])
  des.T=cumwei_cubicbs(1,data.id$Time[data.id$delta==1],sig1)
  N=nrow(des.Y)
  l=as.data.frame(table(data$id))$Freq
  ### collection of K(t)=Q.2%*%R.2(t)%*%t(Q.2) for t=data.id$Time[data.id$delta==1]
  K.2=array(c(sapply(data.id$Time[data.id$delta==1],function(t) 1/(t-1)*Q.2%*%(R.2(t)-R.2(1))%*%t(Q.2))),
            dim=c(q,q,length(data.id$Time[data.id$delta==1]))) 
  ### EM algorithm  ###
  K=10
  gamma_set=numeric(K)
  alpha1_set=numeric(K)
  alpha2_set=numeric(K)
  beta_set=matrix(0,ncol=q,nrow=K)
  sigma2_set=numeric(K)
  D_set=array(0,dim=c(q,q,K))
  Q.fun_set=numeric(K)
  logLik_set=numeric(K)
  b_set=matrix(0,ncol=q,nrow=m)
  inv.Fish_set=array(0,dim=c(q,q,m)) ## V(bi)~=inv.Fishi
  ## Estimation
  for (k in 1:K){
    
    for(i in 1:m){
      Zi=des.Y[data$id==i,]
      Xi=Zi
      Yi=data$Y[data$id==i]
      chol_invSig=chol(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
      Q=solve(chol_invSig)
      bi=c(Q%*%t(Q)%*%(t(Zi)/c(sigma2))%*%(Yi-Xi%*%beta)) #E(bi|Yi,theta),initial value of bi
      a=unique(data.id$Time[(data.id$Time<=Time[i])&(data.id$delta==1)])
      if(length(a)>0){
        if(length(a)>1){
          ### initial value of p(Ti,deltai,Yi,bi;theta)
          log.p.b=dmvnorm(bi,mean=rep(0,q),sigma=D,log=TRUE)
          log.p.Yb=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%bi),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
          log.s=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                            alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(beta+bi)%*%y%*%(beta+bi)))))
          
          log.hazard=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                              alpha2*sqrt(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)))
          
          log.p.YTb=log.hazard-log.s+log.p.Yb+log.p.b
          
          for(kk in 1:10){
            wei=c(exp(gamma*W[i])*sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                          alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(beta+bi)%*%y%*%(beta+bi)))))
            ZK=alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]+alpha2*t(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-1/2)*y%*%(beta+bi)))
            Sbi=-colSums(wei*ZK)+t(Zi)%*%(Yi-Xi%*%beta-Zi%*%bi)/sigma2-solve(D)%*%bi
            
            Fishi=tensor(array(apply(ZK,1,tcrossprod),dim=c(q,q,length(a)))+
                           alpha2*array(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-1/2)*y),dim=c(q,q,length(a)))-
                           alpha2*array(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-3/2)*tcrossprod(y%*%(beta+bi))),dim=c(q,q,length(a))),
                         wei,3,1)+t(Zi)%*%Zi/sigma2+solve(D)
            if(delta[i]==1){
              Sbi=alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)+Sbi
              Fishi=-alpha2*(-c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-3/2)*tcrossprod(K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))+
                               c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])])+Fishi
            }
            ## update bi using Newton-Raphson: bi=bi+solve(I)S
            binew=c(bi+solve(Fishi)%*%Sbi)
            if(all(abs((binew-bi)/bi)[which(is.nan(abs((binew-bi)/bi))!=TRUE)]<10^(-3))){
              break
            }
            ## Check
            log.p.b.new=dmvnorm(binew,mean=rep(0,q),sigma=D,log=TRUE)
            log.p.Yb.new=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%binew),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
            log.s.new=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                                                  alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(beta+binew)%*%y%*%(beta+binew)))))
            log.hazard.new=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                    alpha2*sqrt(t(beta+binew)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+binew)))
            
            log.p.YTb.new=log.hazard.new-log.s.new+log.p.Yb.new+log.p.b.new
            if(log.p.YTb.new<=log.p.YTb){
              #print("oh")
              break
            }
            bi=binew
            log.p.YTb=log.p.YTb.new
            #print(log.p.YTb)
            
          }
        }else{
          ### initial value of p(Ti,deltai,Yi,bi;theta)
          log.p.b=dmvnorm(bi,mean=rep(0,q),sigma=D,log=TRUE)
          log.p.Yb=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%bi),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
          
          log.s=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                            alpha2*sqrt(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))))
          
          log.hazard=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                              alpha2*sqrt(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)))
          
          log.p.YTb=log.hazard-log.s+log.p.Yb+log.p.b
          
          for(kk in 1:10){
            wei=c(exp(gamma*W[i])*sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi)+
                                                          alpha2*sqrt(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))))
            ZK=alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi)
            Sbi=-colSums(wei*ZK)+t(Zi)%*%(Yi-Xi%*%beta-Zi%*%bi)/sigma2-solve(D)%*%bi
            
            Fishi=wei*(tcrossprod(ZK)-alpha2*c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-3/2)*tcrossprod(K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))+
                         alpha2*c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(a,data.id$Time[data.id$delta==1])])+t(Zi)%*%Zi/sigma2+solve(D)
            if(delta[i]==1){
              Sbi=alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi)+Sbi
              Fishi=-alpha2*(-c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-3/2)*tcrossprod(K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))+
                               c(t(beta+bi)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(Time[i],data.id$Time[data.id$delta==1])])+Fishi
            }
            ## update bi using Newton-Raphson: bi=bi+solve(I)S
            binew=c(bi+solve(Fishi)%*%Sbi)
            if(all(abs((binew-bi)/bi)[which(is.nan(abs((binew-bi)/bi))!=TRUE)]<10^(-3))){
              break
            }
            ## Check
            log.p.b.new=dmvnorm(binew,mean=rep(0,q),sigma=D,log=TRUE)
            log.p.Yb.new=sum(dnorm(Yi,mean=c(Xi%*%beta+Zi%*%binew),sd=sqrt(sigma2),log=TRUE)) # p(Y|b)
            log.s.new=exp(gamma*W[i])*sum(sapply(a,lambda0)*exp(alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                                                  alpha2*sqrt(t(beta+binew)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+binew))))
            log.hazard.new=ifelse(delta[i]==0,0,log(lambda0(Time[i]))+gamma*W[i]+alpha1*des.T[match(Time[i],data.id$Time[data.id$delta==1]),]%*%(beta+binew)+
                                    alpha2*sqrt(t(beta+binew)%*%K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]%*%(beta+binew)))
            
            log.p.YTb.new=log.hazard.new-log.s.new+log.p.Yb.new+log.p.b.new
            if(log.p.YTb.new<=log.p.YTb){
              #print("oh")
              break
            }
            bi=binew
            log.p.YTb=log.p.YTb.new
            #print(log.p.YTb)
            
          }
          
        }
        
      }else{
        Fishi=t(Zi)%*%Zi/sigma2+solve(D)
      }
      
      b_set[i,]=bi
      # Fish_set[,,i]=Fishi
      inv.Fish_set[,,i]=solve(Fishi)
    }
    
    tr=0 #sum of trace
    for(i in 1:m){
      tra=sum(diag(t(des.Y[data$id==i,])%*%des.Y[data$id==i,]
                   %*%inv.Fish_set[,,i]))
      tr=tr+tra
    }
    
    ###################### Updata parameters using b_set
    sgamma=0
    salpha1=0
    salpha2=0
    sbeta=numeric(q)
    Igamma=0
    Ialpha1=0
    Ialpha2=0
    Ibeta=matrix(0,ncol=q,nrow=q)
    Igamalp1=0
    Igamalp2=0
    Ialp12=0
    
    
    for(i in data.id$id[data.id$delta==1]){
      
      rs=riskset(Time[i])
      Exp.f=numeric(length(rs))
      Exp.fBb=numeric(length(rs))
      Exp.fBb2=numeric(length(rs))
      Exp.fbetaKb=numeric(length(rs))
      sqrbetaKb=numeric(length(rs))
      Exp.sqrfbetaKb=numeric(length(rs))
      Exp.fBbbetaKb=numeric(length(rs))
      Exp.scorebeta=matrix(0,ncol=q,nrow=length(rs))
      
      
      K2=K.2[,,match(Time[i],data.id$Time[data.id$delta==1])]
      desT=des.T[match(Time[i],data.id$Time[data.id$delta==1]),]
      
      betabi=adv(c(beta+b_set[i,]))
      misqrbetaKbKb=c(c(t(beta+b_set[i,])%*%K2%*%(beta+b_set[i,]))^(-1/2)*K2%*%(beta+b_set[i,]))+
        sapply(c(1:q),function(s) sum(diag(inv.Fish_set[,,i]%*%attr((quad.adv(K2,betabi))^(-1/2)*linear_sca.adv(K2,betabi,s),"grad2")))/2)
      
      
      
      for(j in rs){
        f=c(exp(alpha1*desT%*%b_set[j,]+
                  alpha2*sqrt(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))))
        inv.Fishj=inv.Fish_set[,,j]
        
        
        fg1=c(f*(alpha1*desT+alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,])))
        
        fg2=f*(tcrossprod(alpha1*desT+alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,]))-
                 alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-3/2)*tcrossprod(K2%*%(beta+b_set[j,]))+
                 alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2)
        
        h=c(sqrt(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,])))
        
        hg1=c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,])
        
        hg2=-c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-3/2)*tcrossprod(K2%*%(beta+b_set[j,]))+
          c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2
        
        Exp.f[which(rs==j)]=f+sum(diag(inv.Fishj%*%fg2))/2
        
        Bb=c(t(desT)%*%b_set[j,])
        
        Exp.fBb[which(rs==j)]=f*Bb+
          sum(diag(inv.Fishj%*%(fg2*Bb+fg1%*%t(desT)+desT%*%t(fg1))))/2
        
        Exp.fBb2[which(rs==j)]=f*Bb^2+
          sum(diag(inv.Fishj%*%(Bb^2*fg2+2*Bb*(fg1%*%t(desT)+desT%*%t(fg1))+2*f*desT%*%t(desT))))/2
        
        betaKb=c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))
        
        Exp.fbetaKb[which(rs==j)]=f*betaKb+
          sum(diag(inv.Fishj%*%(betaKb*fg2+2*fg1%*%t(K2%*%(beta+b_set[j,]))+2*(K2%*%(beta+b_set[j,]))%*%t(fg1)+2*f*K2)))/2
        
        sqrbetaKb[which(rs==j)]=h+sum(diag(inv.Fishj%*%hg2))/2
        
        Exp.sqrfbetaKb[which(rs==j)]=f*h+sum(diag(inv.Fishj%*%(h*fg2+fg1%*%t(hg1)+hg1%*%t(fg1)+f*hg2)))/2
        
        Exp.fBbbetaKb[which(rs==j)]=f*h*Bb+sum(diag(inv.Fishj%*%(Bb*(h*fg2+fg1%*%t(hg1)+hg1%*%t(fg1)+f*hg2)+(h*fg1+f*hg1)%*%t(desT)+desT%*%t((h*fg1+f*hg1)))))/2
        
        betabj=adv(c(beta+b_set[j,]))
        Exp.scorebeta[which(rs==j),]=c(f*alpha2*c(t(beta+b_set[j,])%*%K2%*%(beta+b_set[j,]))^(-1/2)*K2%*%(beta+b_set[j,]))+
          sapply(c(1:q),function(s) sum(diag(inv.Fishj%*%attr(exp(alpha1*linear.adv(desT,adv(b_set[j,])))*exp(alpha2*(quad.adv(K2,betabj))^(1/2))*(alpha2*(quad.adv(K2,betabj))^(-1/2))*linear_sca.adv(K2,betabj,s),
                                                              "grad2")))/2)
        
      }
      
      ### score of gamma; Fisher of gamma
      d=exp(gamma*W[rs])#*(exp(alpha*xbeta))# a vector, can be cancelled 
      Deno=c(t(d)%*%Exp.f)
      ssgamma=W[i]-t(W[rs]*d)%*%Exp.f/Deno
      sgamma=sgamma+ssgamma
      igamma=t(W[rs]^2*d)%*%Exp.f/Deno-(t(W[rs]*d)%*%Exp.f/Deno)^2
      Igamma=Igamma+igamma
      
      ### score of alpha1; Fisher of alpha1
      xbeta=matrix(rep(c(t(desT)%*%beta),length(rs)),ncol=1)#collection of Xj(Ti)%*%beta
      ssalpha1=desT%*%(beta+b_set[i,])-t(d)%*%(xbeta*Exp.f+Exp.fBb)/Deno
      salpha1=salpha1+ssalpha1
      ialpha1=t(d)%*%(xbeta*Exp.fBb+Exp.fBb2)/Deno-(t(d)%*%(xbeta*Exp.f+Exp.fBb))*(t(d)%*%Exp.fBb)/(Deno^2)
      Ialpha1=Ialpha1+ialpha1
      
      ### score of alpha2; Fisher of alpha2
      ssalpha2=sqrbetaKb[which(rs==i)]-t(d)%*%Exp.sqrfbetaKb/Deno
      salpha2=salpha2+ssalpha2
      ialpha2=t(d)%*%Exp.fbetaKb/Deno-(t(d)%*%Exp.sqrfbetaKb/Deno)^2
      Ialpha2=Ialpha2+ialpha2
      
      ##############
      igamalp1=t(d)%*%(Exp.fBb*W[rs])/Deno-(t(d)%*%Exp.fBb)*(t(W[rs]*d)%*%Exp.f)/(Deno^2)
      Igamalp1=Igamalp1+igamalp1
      
      igamalp2=t(d)%*%(Exp.sqrfbetaKb*W[rs])/Deno-(t(d)%*%Exp.sqrfbetaKb)*(t(W[rs]*d)%*%Exp.f)/(Deno^2)
      Igamalp2=Igamalp2+igamalp2
      
      ialp12=t(d)%*%(Exp.fBbbetaKb+c(desT%*%beta)*Exp.sqrfbetaKb)/Deno-(t(d)%*%Exp.sqrfbetaKb)*(t(d)%*%(Exp.fBb+c(desT%*%beta)*Exp.f))/(Deno^2)
      Ialp12=Ialp12+ialp12
      
      ### part of score of beta; part of Fisher of beta
      
      ssbeta=alpha1*c(desT)+alpha2*misqrbetaKbKb-(Deno*alpha1*desT+colSums(d*Exp.scorebeta))/Deno
      sbeta=sbeta+ssbeta
      ibeta=tcrossprod(ssbeta)
      Ibeta=Ibeta+ibeta
      
      
    }
    
    res=c(data$Y-apply(des.Y*t(beta+t(b_set[rep(c(1:m),l),])),1,sum))
    
    sigma.2=1/N*((data$Y-c(des.Y%*%beta))%*%(data$Y-c(des.Y%*%beta)-2*apply(des.Y*b_set[rep(c(1:m),l),],1,sum))+
                   tr+sum(apply(des.Y*b_set[rep(c(1:m),l),],1,sum)^2))
    sigma.2=c(sigma.2)
    sbeta=c(sbeta)+1/sigma.2*t(des.Y)%*%res
    
    Ibeta=Ibeta+(-(sigma.2)^(-2)*2/N)*t(des.Y)%*%res%*%t(res)%*%des.Y+(sigma.2)^(-1)*t(des.Y)%*%des.Y
    
    
    
    ### UPDATE PARAMETERS
    Igamalp12=matrix(0,ncol=3,nrow=3)
    Igamalp12[upper.tri(Igamalp12)]=c(Igamalp1,Igamalp2,Ialp12)
    Igamalp12=Igamalp12+t(Igamalp12)
    diag(Igamalp12)=c(Igamma,Ialpha1,Ialpha2)
    if(k>1){
      loglik=logLik_set[k-1]
    }
    
    
    for(v in 0:10){
      stepsize=2^(-v)
      par.sur=c(gamma,alpha1,alpha2)+stepsize*solve(Igamalp12,c(sgamma,salpha1,salpha2))
      gammanew=par.sur[1]
      alpha1new=par.sur[2]
      alpha2new=par.sur[3]
      betanew=c(beta+stepsize*solve(Ibeta)%*%sbeta)
      
      sigma2new=1/N*(t(data$Y-des.Y%*%betanew)%*%(data$Y-des.Y%*%betanew-2*apply(des.Y*b_set[rep(c(1:m),l),],1,sum))+tr+sum(apply(des.Y*b_set[rep(c(1:m),l),],1,sum)^2))
      sigma2new=c(sigma2new)
      Dnew=1/m*(apply(inv.Fish_set,c(1,2),sum)+t(b_set)%*%b_set)
      Dnew=diag(diag(Dnew))
      
      
      cc=c()
      for(t in sort(unique(data$Time[data$delta==1]))){
        rs=riskset(t)
        Exp.f=numeric(length(rs))
        for(j in rs){
          f=c(exp(alpha1new*des.T[match(t,data.id$Time[data.id$delta==1]),]%*%b_set[j,]+
                    alpha2new*sqrt(t(betanew+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(betanew+b_set[j,]))))
          
          inv.Fishj=inv.Fish_set[,,j]
          bb=c(t(betanew+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(betanew+b_set[j,]))

          fg2=f*(tcrossprod(alpha1new*des.T[match(t,data.id$Time[data.id$delta==1]),]+alpha2new*bb^(-1/2)*K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(betanew+b_set[j,]))-
                   alpha2new*bb^(-3/2)*tcrossprod(K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(betanew+b_set[j,]))+
                   alpha2new*bb^(-1/2)*K.2[,,match(t,data.id$Time[data.id$delta==1])])
          
          
          Exp.f[which(rs==j)]=f+sum(diag(inv.Fishj%*%fg2))/2
          if(Exp.f[which(rs==j)]<0){
            Exp.f[which(rs==j)]=f
          }
        }
        
        xbeta=matrix(rep(des.T[match(t,data.id$Time[data.id$delta==1]),],length(rs)),ncol=q,byrow=TRUE)%*%betanew #collection of Xj(t)%*%beta1
        cc=c(cc,length(data.id$id[(data.id$Time==t)&(data.id$delta==1)])/(c(t(exp(gammanew*W[rs]+alpha1new*xbeta))%*%Exp.f)))
      }
      
      cumbase=data.frame(hazard=0,time=c(0,sort(unique(data$Time[data$delta==1]))))
      cumbase[,1]=cumsum(c(0,cc))
      
      
      logLikmcnew=logLik(data,data.id,gammanew,alpha1new,alpha2new,betanew,sigma2new,Dnew,cumbase,knots,Q.2,sig1,L=2000)
      
      cat("UPDATED LOGLIKELIHOOD=",logLikmcnew,"\n")
      if(logLikmcnew>loglik){
        print(v)
        break
      }
      
    }
    
    if(k==1){
      if(logLikmcnew<loglik){
        return(NULL)
      }else{
        if(((abs(logLikmcnew-loglik)/abs(loglik))<2*10^(-6))){
          out=list(iter=k,logLik=logLikmcnew,beta=betanew,sigma2=sigma2new,D=Dnew,gamma=gammanew,alpha1=alpha1new,alpha2=alpha2new,cumbase=cumbase)
          return(out)
        }
      }
    }else{
if(logLikmcnew<loglik){
        cc=c()
        for(t in sort(unique(data$Time[data$delta==1]))){
          rs=riskset(t)
          Exp.f=numeric(length(rs))
          for(j in rs){
            f=c(exp(alpha1*des.T[match(t,data.id$Time[data.id$delta==1]),]%*%b_set[j,]+
                      alpha2*sqrt(t(beta+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))))
            
            inv.Fishj=inv.Fish_set[,,j]
            
            fg2=f*(tcrossprod(alpha1*des.T[match(t,data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))^(-1/2)*K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))-
                     alpha2*c(t(beta+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))^(-3/2)*tcrossprod(K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))+
                     alpha2*c(t(beta+b_set[j,])%*%K.2[,,match(t,data.id$Time[data.id$delta==1])]%*%(beta+b_set[j,]))^(-1/2)*K.2[,,match(t,data.id$Time[data.id$delta==1])])
            
            
            Exp.f[which(rs==j)]=f+sum(diag(inv.Fishj%*%fg2))/2
            if(Exp.f[which(rs==j)]<0){
              Exp.f[which(rs==j)]=f
            }
          }
          xbeta=matrix(rep(des.T[match(t,data.id$Time[data.id$delta==1]),],length(rs)),ncol=q,byrow=TRUE)%*%beta #collection of Xj(t)%*%beta1
          cc=c(cc,length(data.id$id[(data.id$Time==t)&(data.id$delta==1)])/(c(t(exp(gamma*W[rs]+alpha1*xbeta))%*%Exp.f)))
        }
        
        cumbase=data.frame(hazard=0,time=c(0,sort(unique(data$Time[data$delta==1]))))
        cumbase[,1]=cumsum(c(0,cc))
        out=list(iter=k-1,logLik=loglik,beta=beta,sigma2=sigma2,D=D,gamma=gamma,alpha1=alpha1,alpha2=alpha2,cumbase=cumbase)
        return(out)
      }else{
        if(((abs(logLikmcnew-loglik)/abs(loglik))<2*10^(-6))){
          out=list(iter=k,logLik=logLikmcnew,beta=betanew,sigma2=sigma2new,D=Dnew,gamma=gammanew,alpha1=alpha1new,alpha2=alpha2new,cumbase=cumbase)
          return(out)
        }
        
      }}
    
    
    
    gamma=gammanew
    alpha1=alpha1new
    alpha2=alpha2new
    beta=betanew
    D=Dnew
    sigma2=sigma2new
    logLikmc=logLikmcnew
    
    logLik_set[k]=logLikmc
    gamma_set[k]=gamma
    alpha1_set[k]=alpha1
    alpha2_set[k]=alpha2
    beta_set[k,]=beta
    sigma2_set[k]=sigma2
    D_set[,,k]=D
    cat(k,"th element in LOGLIK_set is",logLikmc,"\n")
    
  }
  
}


########## METHOD 2: UPDADE  SIG1 AND SIG2 BY OPTIMIZING log of E^{bi|Yi}[p(Ti|bi)]
est.sig1=function(data,data.id,gamma,alpha1,alpha2,cumbase,beta,sigma2,D,knots,Q.2,sig1){
  lambda0=function(t){
    if(t %in% unique(data$Time[data$delta==1])){
      x=cumbase[which(near(cumbase[,2],t)),1]-cumbase[which(near(cumbase[,2],t))-1,1] ##warn: options(digits = 22) default 7
      if(length(x)>1){
        x=x[1]
      }
      
    }else{
      x=0}
    return(x)
  }
  
  flag=0
  des.Y=as.matrix(data[,c(paste0("time",c(1:q)))])
  bnew_set=matrix(0,ncol=q,nrow=m)
  for(i in 1:m){
    Zi=des.Y[data$id==i,]
    Xi=Zi
    Yi=data$Y[data$id==i]
    chol_invSig=chol(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
    Q=solve(chol_invSig)
    bnew_set[i,]=c(Q%*%t(Q)%*%(t(Zi)/c(sigma2))%*%(Yi-Xi%*%beta)) #E(bi|Yi,theta),initial value of bi
  }
  

  conlogliknew_sig1=function(sig1){
    des.T=cumwei_cubicbs(1,data.id$Time[data.id$delta==1],sig1)

    K.2=array(c(sapply(data.id$Time[data.id$delta==1],function(t) 1/(t-1)*Q.2%*%(R.2(t)-R.2(1))%*%t(Q.2))),
              dim=c(q,q,length(data.id$Time[data.id$delta==1])))
    Q.fun.up.new=0
    for(i in data.id$id){
      Zi=des.Y[data$id==i,]
      Xi=Zi
      Yi=data$Y[data$id==i]
      bi=bnew_set[i,]
      chol_invSig=chol(solve(D)+c(sigma2^(-1))*t(Zi)%*%Zi)
      Q=solve(chol_invSig)
      inv.Fishi=Q%*%t(Q)
      a=unique(data.id$Time[(data.id$Time<=data.id$Time[data.id$id==i])&(data.id$delta==1)])
      log.hazard=0
      log.s=0
      Fishinew=matrix(0,ncol=q,nrow = q) ##length(a)=0 means only longitudinal part left
      if(length(a)>0){
        log.hazard=ifelse(data.id$delta[data.id$id==i]==0,0,log(lambda0(data.id$Time[data.id$id==i]))+c(gamma*data.id$W[data.id$id==i])+alpha1*(des.T[match(data.id$Time[data.id$id==i],data.id$Time[data.id$delta==1]),]%*%(beta+bi))+
                            alpha2*sqrt(t(beta+bi)%*%K.2[,,match(data.id$Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(beta+bi)))
        if(length(a)>1){
          weinew=c(exp(c(gamma*data.id$W[data.id$id==i]))*sapply(a,lambda0)*exp(alpha1*(des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi))+
                                                                                  alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) sqrt(t(beta+bi)%*%y%*%(beta+bi)))))
          
          ZKnew=alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]+alpha2*apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-1/2))*
            t(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) t(y)%*%beta)+apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) y%*%bi))
          
          Sbi=-colSums(weinew*ZKnew)
          
          Fishinew=tensor(array(apply(ZKnew,1,tcrossprod),dim=c(q,q,length(a)))+
                            alpha2*array(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-1/2)*y),dim=c(q,q,length(a)))-
                            alpha2*array(apply(K.2[,,match(a,data.id$Time[data.id$delta==1])],3,function(y) (c(t(beta+bi)%*%y%*%(beta+bi)))^(-3/2)*tcrossprod(t(y)%*%beta+y%*%bi)),dim=c(q,q,length(a))),
                          weinew,3,1)+Fishinew
          
        }else{
          weinew=c(exp(c(gamma*data.id$W[data.id$id==i]))*sapply(a,lambda0)*exp(alpha1*(des.T[match(a,data.id$Time[data.id$delta==1]),]%*%(beta+bi))+
                                                                                  alpha2*sqrt(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi))))
          
          ZKnew=alpha1*des.T[match(a,data.id$Time[data.id$delta==1]),]+
            alpha2*((c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi)))^(-1/2)*(t(K.2[,,match(a,data.id$Time[data.id$delta==1])])%*%beta+K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%bi))
          
          Sbi=-weinew*ZKnew
          
          Fishinew=weinew*(tcrossprod(ZKnew)-alpha2*(c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi)))^(-3/2)*tcrossprod(t(K.2[,,match(a,data.id$Time[data.id$delta==1])])%*%beta+K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%bi)+
                             alpha2*(c(t(beta+bi)%*%K.2[,,match(a,data.id$Time[data.id$delta==1])]%*%(beta+bi)))^(-1/2)*K.2[,,match(a,data.id$Time[data.id$delta==1])])+Fishinew
          
          
        }
        
        log.s=sum(weinew)
      }
      
      
      if(data.id$delta[data.id$id==i]==1){
        Sbi=alpha1*des.T[match(data.id$Time[data.id$id==i],data.id$Time[data.id$delta==1]),]+alpha2*c(t(beta+bi)%*%K.2[,,match(data.id$Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(beta+bi))^(-1/2)*K.2[,,match(data.id$Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(beta+bi)+Sbi
        
        Fishinew=-alpha2*(-(c(t(beta+bi)%*%K.2[,,match(data.id$Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(beta+bi)))^(-3/2)*tcrossprod(t(K.2[,,match(data.id$Time[data.id$id==i],data.id$Time[data.id$delta==1])])%*%beta+K.2[,,match(data.id$Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%bi)+
                            +(c(t(beta+bi)%*%K.2[,,match(data.id$Time[data.id$id==i],data.id$Time[data.id$delta==1])]%*%(beta+bi)))^(-1/2)*K.2[,,match(data.id$Time[data.id$id==i],data.id$Time[data.id$delta==1])])+Fishinew
      }
      
      aa=exp(log.hazard-log.s)+sum(diag(inv.Fishi%*%(exp(log.hazard-log.s)*(-Fishinew+tcrossprod(Sbi)))))/2
      Q.fun.up.new=Q.fun.up.new+log(ifelse(aa<0,exp(log.hazard-log.s),aa))
      
    }
    return(-Q.fun.up.new)
  }
  
  conloglik=conlogliknew_sig1(sig1)
  sig1_opt=optim(c(sig1),conlogliknew_sig1,lower=0.01,upper=15,method = "L-BFGS-B")
  
  conloglik_up=sig1_opt$value
  if(conloglik_up<conloglik){
    sig1=sig1_opt$par
    if(abs(conloglik-conloglik_up)/conloglik<10^(-6)){
      flag=1
      return(list(flag=flag,sig1=sig1_opt$par))
    }else{
      return(list(flag=flag,sig1=sig1_opt$par))
    }
  }else{
    flag=1
    return(list(flag=flag,sig1=sig1))
    
  }
  
}


