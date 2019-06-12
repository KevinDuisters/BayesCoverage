#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Oct 2018
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
# Coverage function

coverage <- function(thetaseq,alpha,lambda,w,dist,plot.cov=F,cols=rep("black",5)){
  
  #xgrid <- seq(min(thetaseq)-20,max(thetaseq)+15,0.005)
  xgrid <- seq(min(thetaseq)-30,max(thetaseq)+30,0.005)
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
  
  XU.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Ugrid>theta0))]})
  XU.sup <- sapply(thetaseq,function(theta0){xgrid[max(which(Ugrid<theta0))]  })
  XL.sup <-sapply(thetaseq,function(theta0){xgrid[max(which(Lgrid<theta0))]})
  XL.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Lgrid>theta0))]})
  
  if(dist=="Lap"){G <- function(x,theta0=0){plaplace(x,m=theta0,s=1)} }
  if(dist=="Normal"){G <- function(x,theta0=0){pnorm(x,theta0,1)} }
  if(dist=="t3"){G <- function(x){pt(x,3)} }
  if(dist=="t5"){G <- function(x){pt(x,5)} }
  if(dist=="Cauchy"){G <- function(x,theta0=0){pcauchy(x,theta0,1)} }
  
  #plot(thetaseq,(1-alpha)/2*G(lambda-XL.inf),ylim=c(0,alpha/2))
  #lines(thetaseq,alpha*G(lambda-XU.sup),col="red")
  #abline(h=(1-alpha)/2*G(Ginv(alpha/(1+alpha))),lty=2,col="green")
  #abline(h=(1-alpha)/2*G(3*Ginv(alpha/(1+alpha))),lty=2,col="blue")
  #abline(h=(1-alpha)/2*G(Ginv(alpha/(1+alpha)) + Ginv(2*alpha/(2+alpha))),lty=2)
  #abline(h=alpha^2/(2+alpha),lty=2,col="red")
  
  
  C.inf <- G(XL.inf-thetaseq)-G(XU.sup-thetaseq)
  C.sup <- G(XL.sup-thetaseq)-G(XU.inf-thetaseq)
  
  if(plot.cov == F){
    return(list(C.inf=C.inf, C.sup=C.sup, xgrid=xgrid,regimeU=regimeU,regimeL=regimeL))
    }else{
      plot(thetaseq,C.sup,xlim=range(thetaseq),type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])),ylim=c(1-2*alpha,1))
      abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),lty=rep(3,4),col=rep("grey",4))
      polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),col="white",border="white")
      polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),col="grey90",border="white",density=20,angle=45)
      polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),col="grey90",border="white",density=20,angle=-45)
      if(lambda==0.75 & w==1 & dist=="Normal"){segments(x0=3.41,x1=3.41,y0=0.952,y1=0.956,col="grey90",lwd=1.5)}# manual tweak to make top-left plot Normal look nicer, be careful generalizing this
      text(x=rep(lambda+9,4),y=0.005+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=0.8,adj=0)
      
      
      for(r in 5:1){
          reg <- (sapply(XU.sup,function(x) regimeU[which(xgrid==x)])==r) # inf
          if(sum(reg)>0){lines(thetaseq[reg],C.inf[reg],col=cols[r],lty=2)}
          reg <- (sapply(XU.inf,function(x) regimeU[which(xgrid==x)])==r)# Sup
          if(sum(reg)>0){lines(thetaseq[reg],C.sup[reg],col= cols[r])}
      }
    }
}




