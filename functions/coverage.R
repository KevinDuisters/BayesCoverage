#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Feb 2020
#--------------------------------------------------------------------------------------------------------------------#
# update Feb 20: exclude inf/sup region, code is no longer used

#--------------------------------------------------------------------------------------------------------------------#
# Coverage function
coverage <- function(alpha,lambda,w,dist,thetamax,h=0.01,plot.cov=F,plot.cov.title=T,plot.neg=F){
  
  excl <- T # update Feb 20: exclude inf/sup region, code is no longer used
  
  thetaseq <- seq(lambda,thetamax,h)
  xgrid <- seq(min(thetaseq)-30,max(thetaseq)+30,0.005)
  if(dist=="Cauchy"| dist=="t1"){xgrid <- seq(min(thetaseq)-250,max(thetaseq)+250,0.005)}
  
  Ugrid <- Lgrid <- regimeU <- regimeL <- numeric(length(xgrid))
  
  for(i in 1:length(xgrid)){
    x <- xgrid[i]
    out.U <-getU(x,alpha,lambda,w,dist)
    Ugrid[i] <- out.U$val
    regimeU[i] <- out.U$regime
    
    out.L <- getU(-x,alpha,lambda,w,dist)
    Lgrid[i] <- -(out.L$val)  
    regimeL[i] <- out.L$regime
  }
  
  
  if(dist=="Lap"){
    g <- function(x){dlaplace(x,m=0,s=1)}
    G <- function(x){plaplace(x,m=0,s=1)}
    distname<-"Laplace(0,1)"
    }
  if(dist=="Normal"){
    g <- function(x){dnorm(x,0,1)} 
    G <- function(x){pnorm(x,0,1)} 
    distname<-"N(0,1)"
    }
  if(dist=="t3"){
    g <- function(x){dt(x,3)} 
    G <- function(x){pt(x,3)} 
  distname<-"t(3)"
  }
  if(dist=="t2"){
    g <- function(x){dt(x,2)} 
    G <- function(x){pt(x,2)} 
    distname<-"t(2)"
  }
  if(dist=="t1"){
    g <- function(x){dt(x,1)} 
    G <- function(x){pt(x,1)} 
    distname<-"t(1)"
    }
  if(dist=="Cauchy"){
    g <- function(x){dcauchy(x,0,1)} 
    G <- function(x){pcauchy(x,0,1)} 
    distname <- "Cauchy"
    }
  
  # determine t_alpha
  ta <- -1e8
  if(2*G(-lambda) <= ((1-w)/w)*(alpha/(1-alpha))*g(0)){
    ta <-max(0,xgrid[((w/(1-w))*(G(xgrid - lambda) + G(-xgrid - lambda))/g(xgrid)) <= (alpha/(1-alpha))])
  }
  
  C.inf <- C.sup <- NULL
  if(excl==F){
    XU.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Ugrid>theta0))]})
    XU.sup <- sapply(thetaseq,function(theta0){xgrid[max(which(Ugrid<theta0))]  })
    XL.sup <-sapply(thetaseq,function(theta0){xgrid[max(which(Lgrid<theta0))]})
    XL.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Lgrid>theta0))]})
  
    C.inf <- G(XL.inf-thetaseq)-G(XU.sup-thetaseq)
    C.sup <- G(XL.sup-thetaseq)-G(XU.inf-thetaseq)
  
  
    for(t in 1:length(thetaseq)){
      theta0 <- thetaseq[t]
      xstar <- XU.sup[t]
      xtilde <- XU.inf[t]
      xcheck <- XL.inf[t]
      xhat <- XL.sup[t]
    
      # C.inf
      if(abs(xstar) <= ta){
        if(xcheck > ta){C.inf[t] <- G(xcheck - theta0) - G(ta-theta0)}
        if(xcheck < (-ta)){C.inf[t] <- G(xcheck - theta0) - G(-ta-theta0)}
        if(abs(xcheck) < ta){C.inf[t] <- 0}
      }else{
        if(xstar > ta){C.inf[t] <- G(xcheck - theta0) - G(xstar-theta0)}else{
           if(xstar < (-ta)){
              if(xcheck > ta){C.inf[t] <- G(-ta - theta0) - G(xstar-theta0) + G(xcheck - theta0) - G(ta - theta0)}
              if(xcheck < (-ta)){C.inf[t] <- G(xcheck - theta0) - G(xstar-theta0)}
              if(abs(xcheck) < ta){C.inf[t] <- G(-ta - theta0) - G(xstar-theta0)}
      }
      }
    }
      # C.sup
      if(abs(xtilde) <= ta){
        if(xhat > ta){C.sup[t] <- G(xhat - theta0) - G(ta-theta0)}
        if(xhat < (-ta)){C.sup[t] <- G(xhat - theta0) - G(-ta-theta0)}
        if(abs(xhat) < ta){C.sup[t] <- 0}
      }else{
        if(xtilde > ta){C.sup[t] <- G(xhat - theta0) - G(xtilde-theta0)}else{
          if(xtilde < (-ta)){
              if(xhat > ta){C.sup[t] <- G(-ta - theta0) - G(xtilde-theta0) + G(xhat - theta0) - G(ta - theta0)}
              if(xhat < (-ta)){C.sup[t] <- G(xhat - theta0) - G(xtilde-theta0)}
              if(abs(xhat) < ta){C.sup[t] <- G(-ta - theta0) - G(xtilde-theta0)}
        }
        }
      }
  }
  ## 
  }
  
  ## numeric proxy for frequentist coverage
  C.num <- sapply(thetaseq,function(theta0){
              Xcov <- xgrid[which(Ugrid >= theta0 & Lgrid <= theta0 & abs(xgrid) > ta)]
              return(sum(g(Xcov - theta0))/sum(g(xgrid-theta0)))
              })  
  # color numeric proxy
  col.tvec <- sapply(thetaseq,function(theta0){
    reg.t <- regimeU[which(Ugrid >= theta0 & Lgrid <= theta0 & abs(xgrid) > ta)]
    if(is.na(mean(reg.t))){col.t <- "white"}else{
      col.t <- rgb(red=sum(reg.t%in%c(2,4))/length(reg.t),green=sum(reg.t%in%c(1,5))/length(reg.t),blue=sum(reg.t==3)/length(reg.t),alpha=1)
    }
    return(col.t)
  })
  
  
  regimeCinf<-regimeCsup <- NULL
  if(excl==F){
  # coverage regime
  regimeCinf <- regimeCsup <- numeric(length(thetaseq))
  for(r in 5:0){
    reg <- (sapply(XU.sup,function(x) regimeU[which(xgrid==x)])==r) # inf
    if(sum(reg)>0){regimeCinf[reg] <- r}
    
    reg <- (sapply(XU.inf,function(x) regimeU[which(xgrid==x)])==r)# Sup
    if(sum(reg)>0){regimeCsup[reg] <- r}
  }
  }
  
  ## output results
  if(plot.cov == F){
    return(list(C.inf=C.inf, C.sup=C.sup,C.num= C.num, xgrid=xgrid,regimeU=regimeU,regimeL=regimeL,regimeCinf=regimeCinf,regimeCsup=regimeCsup,colvec=col.tvec))
    }else{
      if(plot.neg==F){
      plot(thetaseq,C.sup,type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])),ylim=c(1-2*alpha,1),xlim=c(0,thetamax))
      if(plot.cov.title==T){title(bquote(paste(.(distname),", ",lambda==.(lambda),", ",w==.(w))))}
      abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),lty=rep(3,4),col=rep("grey",4))
      text(x=rep(max(thetaseq)-2,4),y=0.005+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=1,adj=0)
      abline(v=lambda,lty=3,col="darkgrey")
      text(expression(lambda),x=lambda+0.4,y=1-2*alpha+0.005)
      for(i in 2:length(thetaseq)){
        segments(x0=thetaseq[i-1],x1=thetaseq[i],y0=C.num[i-1],y1=C.num[i],col=col.tvec[i])
      }
      }else{
        plot(thetaseq,C.sup,type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])),ylim=c(1-2*alpha,1),xlim=c(-thetamax,thetamax))
        if(plot.cov.title==T){title(bquote(paste(.(distname),", ",lambda==.(lambda),", ",w==.(w))))}
        abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),lty=rep(3,4),col=rep("grey",4))
        text(x=rep(max(thetaseq)-2,4),y=0.005+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=1,adj=0)
        abline(v=c(-lambda,lambda),lty=3,col="darkgrey")
        text(c(expression(-lambda),expression(lambda)),x=c(-lambda-0.4,lambda+0.4),y=rep(1-2*alpha+0.005,2))
        for(i in 2:length(thetaseq)){
          segments(x0=thetaseq[i-1],x1=thetaseq[i],y0=C.num[i-1],y1=C.num[i],col=col.tvec[i])
          segments(x0=-thetaseq[i-1],x1=-thetaseq[i],y0=C.num[i-1],y1=C.num[i],col=col.tvec[i])
        } 
      }
      
    }
}




