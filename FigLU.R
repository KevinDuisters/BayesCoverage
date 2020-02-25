#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Oct 2019
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
# source functions
library(rmutil) #dlaplace in densities
source("functions/getU.R")
#--------------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------------#
## Fig: Visualize L,U 
## Manual input Laplace density
#--------------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------------#
# set universal parameters
alpha <- 0.05
#--------------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------------#
# Universal code chunk for plotting for some w and lambda
code.chunk <- function(lambda,w,alpha,h=0.01){
  
  xgrid <- seq(-20-lambda,20+lambda,h)
  ranges <- c(-10,10)
  Ugrid <- Lgrid <- regimeU <- regimeL <- numeric(length(xgrid))
  
  for(i in 1:length(xgrid)){
    x <- xgrid[i]
    out.U <-getU(x,alpha,lambda,w,dist="Lap")
    Ugrid[i] <- out.U$val
    regimeU[i] <- out.U$regime
    
    out.L <- getU(-x,alpha,lambda,w,dist="Lap")
    Lgrid[i] <- -(out.L$val)
    regimeL[i] <- out.L$regime
  }
  
  
  theta0seq <- c(seq(-lambda-20,-lambda-h,h),seq(lambda+h,lambda+20,h))
  xU.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Ugrid>=theta0))]})
  xU.sup <- sapply(theta0seq,function(theta0){xgrid[max(which(Ugrid<=theta0))]})
  xL.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Lgrid>=theta0))]})
  xL.sup <-sapply(theta0seq,function(theta0){xgrid[max(which(Lgrid<=theta0))]})
  
  
  plot(xgrid,xgrid,ylab=expression(paste(L[alpha](x)," and ",U[alpha](x))),xlab=expression(x),type="n",xlim=ranges,ylim=1.5*ranges)
  abline(h=c(-lambda,lambda),lty=3,col="grey")
  ta <-max(0,xgrid[((w/(1-w))*(plaplace(xgrid - lambda) + plaplace(-xgrid - lambda))/dlaplace(xgrid)) <= (alpha/(1-alpha))])
  
  if(ta==0){ # this is a manual trick, not exact
    abline(v=c(max(xgrid[regimeU==5]),max(xgrid[regimeU==4]),min(xgrid[regimeU==2]),min(xgrid[regimeU==1])),lty=3,col="grey")
    mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])+max(xgrid[regimeU==4]))/2,(max(xgrid[regimeU==4])+min(xgrid[regimeU==2]))/2,(min(xgrid[regimeU==2])+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2))
  }else{
    abline(v=c(max(xgrid[regimeU==5]),-ta,ta,min(xgrid[regimeU==1])),lty=3,col="grey")
    mtext(c("I","IV","","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])-ta)/2,(-ta+ta)/2,(ta+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2))
    mtext(c(expression(-t[alpha]),expression(t[alpha])),side=1,at=c(-ta,ta),line=1)
  }
  
  
  for(r in 0:5){
    reg <- which(regimeU==r)
    lines(xgrid[reg],Ugrid[reg],col="black")
    
    reg <- which(regimeL==r)
    lines(xgrid[reg],Lgrid[reg],col="blue")
  }
  
  mtext(c(expression(-lambda),expression(lambda)),line=-1.5,side=2,at=c(-lambda,lambda))
  legend("topleft",c(expression(U[alpha]),expression(L[alpha])),lty=c(1,1),col=c("black","blue"),horiz=F,bg="white",box.col="white",seg.len=2)
  box()
  
}

#--------------------------------------------------------------------------------------------------------------------#
# Visualize
pdf("Figures/figLU.pdf",width=12,height=4.5)

par(mfrow=c(1,3))
code.chunk(lambda=1,w=0.25,alpha) # left
code.chunk(lambda=5,w=0.25,alpha) # mid
code.chunk(lambda=5,w=1,alpha) # right

dev.off()
#------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
# Extend LU code with figsize

pdf("Figures/figsize.pdf",width=9)
#pdf("Figures/figsizeSI.pdf",width=9) #w=0.2

par(mfrow=c(1,length(lambdas)),pty="s")
for(lambda in lambdas){
  plot(xgrid,xgrid,type="n",xlab="x",ylab=expression(nu(x)),ylim=c(0,10))
  #polygon(x=c(-100,100,100,-100,-100),y=c(-100,-100,0,0,-100),col="lightgrey")
  #box()
  
  
  ta <-max(0,xgrid[((w/(1-w))*(G(xgrid - lambda) + G(-xgrid - lambda))/g(xgrid)) <= (alpha/(1-alpha))])
  
  x1 <- xgrid[(xgrid > lambda + R1f(xgrid,alpha,w,lambda))]
  x2 <- xgrid[(xgrid > 0) & (-lambda + R3f(xgrid,alpha,w,lambda) < xgrid) & (xgrid <= lambda + R1f(xgrid,alpha,w,lambda))  ]
  x3 <- xgrid[(abs(xgrid) <= -lambda + R3f(xgrid,alpha,w,lambda))]
  x4 <- (-x2)
  x5 <- xgrid[(-xgrid > lambda + R1f(xgrid,alpha,w,lambda))]
  
  lines(x1[abs(x1)>ta],2*R1f(x1[abs(x1)>ta],alpha,w,lambda),col="blue") # black
  lines(x2[abs(x2)>ta],x2[x2>ta]+R2f(x2[abs(x2)>ta],alpha,w,lambda)-lambda,col="blue") #green
  lines(x3[abs(x3)>ta],2*R3f(x3[abs(x3)>ta],alpha,w,lambda)-2*lambda,col="blue") # blue
  lines(x4[abs(x4)>ta],-x4[abs(x4)>ta]+R2f(-x4[abs(x4)>ta],alpha,w,lambda)-lambda,col="blue") # red
  lines(x5[abs(x5)>ta],2*R1f(x5[abs(x5)>ta],alpha,w,lambda),col="blue") # black
  
  lapply(list(x1,x2,x3,x4,x5),function(x)lines(x[abs(x)<=ta],rep(0,sum(abs(x)<=ta)),col="blue")) # draw zero size for HPD(x)={0}.
  
  if(ta==0){ # this is a manual trick, not exact
    abline(v=c(max(x5),max(x4),min(x2),min(x1)),lty=3)
    mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(xgrid)+max(x5))/2,(max(x5)+max(x4))/2,(max(x4)+min(x2))/2,(min(x2)+min(x1))/2,(min(x1)+max(xgrid))/2))
  }else{
    abline(v=c(max(x5),-ta,ta,min(x1)),lty=3,col="grey")
    mtext(c("I","IV","-","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(x5))/2,(max(x5)-ta)/2,(-ta+ta)/2,(ta+min(x1))/2,(min(x1)+max(ranges))/2))
    mtext(c(expression(-t[alpha]),expression(t[alpha])),side=1,at=c(-ta,ta),line=1)
  }
  
  
}

dev.off()


