#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Oct 2018
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
# Coverage function

coverage <- function(thetaseq,alpha,lambda,w,dist,plot.cov=F,cols=rep("black",5+1)){
  
  #xgrid <- seq(min(thetaseq)-20,max(thetaseq)+15,0.005)
  xgrid <- seq(min(thetaseq)-25,max(thetaseq)+25,0.01)
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
  
  if(dist=="Lap"){
    g <- function(x){dlaplace(x,m=0,s=1)}
    G <- function(x){plaplace(x,m=0,s=1)} 
    }
  if(dist=="Normal"){
    g <- function(x){dnorm(x,0,1)} 
    G <- function(x){pnorm(x,0,1)} }
  if(dist=="t3"){
    g <- function(x){dt(x,3)} 
    G <- function(x){pt(x,3)} }
  if(dist=="t5"){
    g <- function(x){dt(x,5)} 
    G <- function(x){pt(x,5)} }
  if(dist=="Cauchy"){
    g <- function(x){dcauchy(x,0,1)} 
    G <- function(x){pcauchy(x,0,1)} }
  
  C.inf <- G(XL.inf-thetaseq)-G(XU.sup-thetaseq)
  C.sup <- G(XL.sup-thetaseq)-G(XU.inf-thetaseq)
  
  # numeric proxy for frequentist coverage
  C.num <- sapply(thetaseq,function(theta0){
              Xcov <- xgrid[which(Ugrid >= theta0 & Lgrid <= theta0)]
              return(sum(g(Xcov - theta0))/sum(g(xgrid-theta0)))
              })  
  
  
  # coverage regime
  regimeCinf <- regimeCsup <- numeric(length(thetaseq))
  for(r in 5:0){
    reg <- (sapply(XU.sup,function(x) regimeU[which(xgrid==x)])==r) # inf
    if(sum(reg)>0){regimeCinf[reg] <- r}
    
    reg <- (sapply(XU.inf,function(x) regimeU[which(xgrid==x)])==r)# Sup
    if(sum(reg)>0){regimeCsup[reg] <- r}
  }
  
  if(plot.cov == F){
    return(list(C.inf=C.inf, C.sup=C.sup, C.num, xgrid=xgrid,regimeU=regimeU,regimeL=regimeL,regimeCinf=regimeCinf,regimeCsup=regimeCsup))
    }else{
      plot(thetaseq,C.sup,xlim=range(thetaseq),type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])),ylim=c(1-2*alpha,1))
      abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),lty=rep(3,4),col=rep("grey",4))
      polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),col="white",border="white")
      polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),col="grey90",border="white",density=20,angle=45)
      polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),col="grey90",border="white",density=20,angle=-45)
      if(lambda==0.75 & w==1 & dist=="Normal"){segments(x0=3.41,x1=3.41,y0=0.952,y1=0.956,col="grey90",lwd=1.5)}# manual tweak to make top-left plot Normal look nicer, be careful generalizing this
      text(x=rep(lambda+16,4),y=0.005+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=0.8,adj=0)
      
      lines(thetaseq,C.num,col='green',lwd=2)
      #points(thetaseq,C.num,col='green')
      
      for(r in 5:0){
        reg <- (regimeCinf==r) # inf
        if(sum(reg)>0){lines(thetaseq[reg],C.inf[reg],col=cols[r+1],lty=2)}
        reg <- (regimeCsup==r) # sup
        if(sum(reg)>0){lines(thetaseq[reg],C.sup[reg],col= cols[r+1])}
      }
      
    }
}




