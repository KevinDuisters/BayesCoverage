#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Feb 2020
#--------------------------------------------------------------------------------------------------------------------#
# update Feb 20: exclude inf/sup region, code is no longer used

#--------------------------------------------------------------------------------------------------------------------#
# Coverage function
coverage <- function(alpha,lambda,w,dist,thetamax,h=0.01,plot.cov=F,plot.cov.title=T,plot.neg=F,line.col=NULL){
  
  # localize functions
  if(dist=="Lap"){
    g <- function(x){dlaplace(x,m=0,s=1)}
    G <- function(x){plaplace(x,m=0,s=1)}
    sampleg<-function(n){rlaplace(n,m=0,s=1)}
    distname<-"Laplace(0,1)"
  }
  if(dist=="Normal"){
    g <- function(x){dnorm(x,0,1)} 
    G <- function(x){pnorm(x,0,1)} 
    sampleg<-function(n){rnorm(n,0,1)}
    distname<-"N(0,1)"
  }
  if(dist=="t3"){
    g <- function(x){dt(x,3)} 
    G <- function(x){pt(x,3)} 
    sampleg<-function(n){rt(n,3)}
    distname<-"t(3)"
  }
  if(dist=="t2"){
    g <- function(x){dt(x,2)} 
    G <- function(x){pt(x,2)} 
    sampleg<-function(n){rt(n,2)}
    distname<-"t(2)"
  }
  if(dist=="t1"){
    g <- function(x){dt(x,1)} 
    G <- function(x){pt(x,1)} 
    sampleg<-function(n){rt(n,1)}
    distname<-"t(1)"
  }
  if(dist=="Cauchy"){
    g <- function(x){dcauchy(x,0,1)} 
    G <- function(x){pcauchy(x,0,1)} 
    sampleg<-function(n){rcauchy(n,0,1)}
    distname <- "Cauchy"
  }
  
  # generate thetas
  thetaseq <- seq(lambda,thetamax,h)
  
  
  # One strategy would be to sample grid for each theta. However, this is computationally slow.
  # instead, fix a (fine) grid of x and compute L,U, and regime once
  #fix xgrid
  xgrid <- seq(min(thetaseq)-30,max(thetaseq)+30,0.005)
  Ugrid <- Lgrid <- numeric(length(xgrid))
  regime <- numeric(length(xgrid))
  for(i in 1:length(xgrid)){
    x <- xgrid[i]
    out.U <-getU(x,alpha,lambda,w,dist)
    Ugrid[i] <- out.U$val
    regime[i] <- out.U$regime
    
    out.L <- getU(-x,alpha,lambda,w,dist)
    Lgrid[i] <- -(out.L$val)  
  }
  
  
  # loop over thetaseq to get coverage
  
  ## numeric proxy for frequentist coverage
  
  C.num <- sapply(thetaseq,function(theta0){
    Xcov <- xgrid[which(Ugrid >= theta0 & Lgrid <= theta0 )]
    return(sum(g(Xcov - theta0))/sum(g(xgrid-theta0)))
  })
  # color gradient for C.num
  if(is.null(line.col)){
       col.tvec <- sapply(thetaseq,function(theta0){
                    reg.t <- regime[which(Ugrid >= theta0 & Lgrid <= theta0)]
                    if(is.na(mean(reg.t))){
                        col.t <- "white"
                        }else{
                        col.t <- rgb(red=sum(reg.t%in%c(2,4))/length(reg.t),green=sum(reg.t%in%c(1,5))/length(reg.t),blue=sum(reg.t==3)/length(reg.t),alpha=1)          
                        }
                  return(col.t)
                  })
  }else{
      col.tvec <- rep(line.col,length(thetaseq))
  }
  
  
  ## output results
  if(plot.cov == F){
    return(list(C.num= C.num, colvec=col.tvec))
    }else{
      if(plot.neg==F){
      plot(thetaseq,C.num,type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])),ylim=c(1-2*alpha,1),xlim=c(0,thetamax))
      if(plot.cov.title==T){title(bquote(paste(.(distname),", ",lambda==.(lambda),", ",w==.(w))))}
      abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),lty=rep(3,4),col=rep("grey",4))
      text(x=rep(max(thetaseq)-2,4),y=0.005+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=1,adj=0)
      abline(v=lambda,lty=3,col="darkgrey")
      text(expression(lambda),x=lambda+0.4,y=1-2*alpha+0.005)
      for(i in 2:length(thetaseq)){
        segments(x0=thetaseq[i-1],x1=thetaseq[i],y0=C.num[i-1],y1=C.num[i],col=col.tvec[i])
      }
      }else{
        plot(thetaseq,C.num,type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])),ylim=c(1-2*alpha,1),xlim=c(-thetamax,thetamax))
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
      
    } # end plot
} # end function




