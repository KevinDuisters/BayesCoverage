#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Oct 2018
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
# source functions
library(rmutil) #dlaplace

#--------------------------------------------------------------------------------------------------------------------#
# Plots
#--------------------------------------------------------------------------------------------------------------------#
# Fig 2: Visualize L,U and inverse xL, xU 
dist <- "Lap"

# localize g and G functions
if(dist=="Norm"){
  g <- function(x){dnorm(x)}
  G <- function(u){pnorm(u)}
  Ginv <- function(u){qnorm(pmin(u,1))}
  Delta <- function(x,lambda){G(lambda-x) - G(-lambda-x)}
}
if(dist=="Lap"){
  g <- function(x){dlaplace(x)}
  G <- function(u){plaplace(u)}
  Ginv <- function(u){qlaplace(pmin(u,1))}
  Delta <- function(x,lambda){G(lambda-x) - G(-lambda-x)}
}
if(dist=="t3"){
  g <- function(x){dt(x,3)}
  G <- function(u){pt(u,3)}
  Ginv <- function(u){qt(pmin(u,1),3)}
  Delta <- function(x,lambda){G(lambda-x) - G(-lambda-x)}
}


R1 <- function(x,alpha,lambda,w){Ginv(1 - alpha/2 - ((1-w)/(2*w))*alpha*g(x) - ((1-alpha)/2)*Delta(x,lambda))}
R2 <- function(x,alpha,lambda,w){Ginv(1 - alpha - ((1-w)/(w))*alpha*g(x) + alpha*Delta(x,lambda) + G(-lambda-x))}
R3 <- function(x,alpha,lambda,w){Ginv(1 - alpha/2 - ((1-w)/(2*w))*alpha*g(x) + (alpha/2)*Delta(x,lambda))}


alpha <- 0.05
lambdas <- c(0.75,7.5)
w <- 1
#w<-0.2

pdf("Figures/figR.pdf",width=9)

xgrid <- seq(-15,15,0.01)
par(mfrow=c(1,length(lambdas)),pty="s")
for(lambda in lambdas){
  plot(xgrid,xgrid,type="n",xlab="x",ylab="R",ylim=c(-15,15),asp=1) #main=bquote(lambda==.(lambda))
  
  x1 <- xgrid[xgrid > lambda + R1(xgrid,alpha,lambda,w)]
  x2 <- xgrid[(xgrid > 0) & (-lambda + R3(xgrid,alpha,lambda,w) < xgrid) & (xgrid <= lambda + R1(xgrid,alpha,lambda,w))  ]
  x3 <- xgrid[abs(xgrid) <= -lambda + R3(xgrid,alpha,lambda,w)]
  x4 <- (-x2)
  x5 <- xgrid[-xgrid > lambda + R1(xgrid,alpha,lambda,w)]
  
  abline(v=c(max(x5),min(x3),max(x3),min(x1)),lty=3)
  mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(xgrid)+max(x5))/2,(max(x5)+min(x3))/2,(min(x3)+max(x3))/2,(max(x3)+min(x1))/2,(min(x1)+max(xgrid))/2))
  
  lines(xgrid,R1(xgrid,alpha,lambda,w),lty=2)
  lines(xgrid,R2(xgrid,alpha,lambda,w),col="green",lty=2)
  lines(xgrid,R2(-xgrid,alpha,lambda,w),col="red",lty=2)
  lines(xgrid,R3(xgrid,alpha,lambda,w),col="blue",lty=2)
  
  # exclude |x| <= t_alpha
  ta <-xgrid[((w/(1-w))*(G(xgrid - lambda) + G(-xgrid - lambda))/g(xgrid)) <= (alpha/(1-alpha))]
  if(length(ta)>0){
    polygon(x=c(-ta,ta,ta,-ta,-ta),y=c(-100,-100,100,100,-100),col="grey")
  }
  
  
  lines(x1,R1(x1,alpha,lambda,w),lty=1)
  lines(x2,R2(x2,alpha,lambda,w),col="green",lty=1)
  lines(x3,R3(x3,alpha,lambda,w),col="blue",lty=1)
  lines(x4,R2(-x4,alpha,lambda,w),col="red",lty=1)
  lines(x5,R1(x5,alpha,lambda,w),lty=1)
  
  legend("topleft",c(expression(R[1](x)),expression(R[2](-x)),expression(R[2](x)),expression(R[3](x))),
         lty=c(1,1,1,1),col=c(1,"red","green","blue"),bty="n",seg.len=2)
}

dev.off()

# check, should give the same as FigLU :
# plot U
par(mfrow=c(1,length(lambdas)))
for(lambda in lambdas){
  plot(xgrid,xgrid,type="n",xlab="x",ylab="U",ylim=c(-15,15),main=bquote(lambda==.(lambda)))
  
  lines(xgrid,xgrid+R1(xgrid,alpha,lambda,w),lty=2)
  lines(xgrid,xgrid+R2(xgrid,alpha,lambda,w),col="red",lty=2)
  lines(xgrid,xgrid+R3(xgrid,alpha,lambda,w),col="blue",lty=2)
  
  abline(v=lambda,lty=3)
  abline(h=lambda,lty=3)
  
  x1 <- x2 <- x3 <- x5 <-  numeric(0)
  x1 <- xgrid[xgrid > lambda + R1(xgrid,alpha,lambda,w)]
  x5 <- xgrid[-xgrid > lambda + R1(xgrid,alpha,lambda,w)]
  x2 <- xgrid[(xgrid > 0) & (-lambda + R3(xgrid,alpha,lambda,w) < xgrid) & (xgrid <= lambda + R1(xgrid,alpha,lambda,w))]
  x3 <- xgrid[abs(xgrid) <= -lambda + R3(xgrid,alpha,lambda,w)]
  
  if(length(x1)>0){lines(x1,x1+R1(x1,alpha,lambda,w),lty=1)}
  if(length(x5)>0){lines(x5,x5+R1(x5,alpha,lambda,w),lty=1)}
  if(length(x2)>0){lines(x2,x2+R2(x2,alpha,lambda,w),col="red",lty=1)}
  if(length(x3)>0){lines(x3,x3+R3(x3,alpha,lambda,w),col="blue",lty=1)}
  
}

# ---------------------------------------------------------------------------------------- #
# Lebesgue measure HPD credible set


alpha <- 0.05
lambdas <- c(0.75,7.5)
w <- 1
#w<-0.2

pdf("Figures/figsize.pdf",width=9)

xgrid <- seq(-15,15,0.01)
par(mfrow=c(1,length(lambdas)),pty="s")
for(lambda in lambdas){
  plot(xgrid,xgrid,type="n",xlab="x",ylab=expression(paste("Lebesgue measure of ",HPD[alpha])),ylim=c(-15,15),asp=1)
  
  x1 <- xgrid[xgrid > lambda + R1(xgrid,alpha,lambda,w)]
  x2 <- xgrid[(xgrid > 0) & (-lambda + R3(xgrid,alpha,lambda,w) < xgrid) & (xgrid <= lambda + R1(xgrid,alpha,lambda,w))  ]
  x3 <- xgrid[abs(xgrid) <= -lambda + R3(xgrid,alpha,lambda,w)]
  x4 <- (-x2)
  x5 <- xgrid[-xgrid > lambda + R1(xgrid,alpha,lambda,w)]
  
  lines(x1,2*R1(x1,alpha,lambda,w)) # black
  lines(x2,x2+R2(x2,alpha,lambda,w)-lambda,col="black") #green
  lines(x3,2*R3(x3,alpha,lambda,w)-2*lambda,col="black") # blue
  lines(x4,-x4+R2(-x4,alpha,lambda,w)-lambda,col="black") # red
  lines(x5,2*R1(x5,alpha,lambda,w)) # black
  
  abline(v=c(max(x5),min(x3),max(x3),min(x1)),lty=3)
  mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(xgrid)+max(x5))/2,(max(x5)+min(x3))/2,(min(x3)+max(x3))/2,(max(x3)+min(x1))/2,(min(x1)+max(xgrid))/2))
  
}

dev.off()

