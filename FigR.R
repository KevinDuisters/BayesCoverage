#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Feb 2020
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
# source functions
library(rmutil) #dlaplace

#--------------------------------------------------------------------------------------------------------------------#
# Plots
#--------------------------------------------------------------------------------------------------------------------#
# Fig: Visualize R functions
dist <- "Lap"

# localize g and G functions
if(dist=="Lap"){
  g <- function(x,theta0=0){dlaplace(x,m=theta0,s=1)}
  G <- function(x,theta0=0){plaplace(x,m=theta0,s=1)}  
  Ginv <- function(p){qlaplace(p)}
  distname <- "Laplace(0,1)"
}

#--------------------------------------------------------------------------------------------------------------------#
# R functions
R1f <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha/2 - (1-w)/(2*w)*alpha*g(x)-((1-alpha)/2)*(G(lambda-x)-G(-lambda-x)))))}
R2f <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha - (1-w)/w*alpha*g(x)+alpha*(G(lambda-x)-G(-lambda-x)) + G(-lambda-x))))}
R3f <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha/2 - (1-w)/(2*w)*alpha*g(x)+(alpha/2)*(G(lambda-x)-G(-lambda-x)))))}

#--------------------------------------------------------------------------------------------------------------------#
# Settings for plot
alpha <- 0.05
lambda <- 5
w <- 1


pdf("Figures/figR.pdf",width=4,height=4.5)

xgrid <- seq(0,15,0.0001)
ranges <- c(0,15)

par(mfrow=c(1,1))

# get R
R1 <- R1f(xgrid,alpha,w,lambda)
R2 <- R2f(xgrid,alpha,w,lambda)
R3 <- R3f(xgrid,alpha,w,lambda)
  
plot(xgrid,xgrid,type="n",xlab="x",ylab="R",xlim=ranges,ylim=c(0,10)) 
title(bquote(paste(.(distname),", ",lambda==.(lambda),", ",w==.(w))))

  
X1 <- xgrid[xgrid > lambda + R1]
X2 <- xgrid[(xgrid > 0) & (-lambda + R3 < xgrid) & (xgrid <= lambda + R1)  ]
X3 <- xgrid[abs(xgrid) <= -lambda + R3]
X4 <- -X2
X5 <- -X1

  # exclude |x| <= t_alpha
  ta <-max(0,xgrid[((w/(1-w))*(G(xgrid - lambda) + G(-xgrid - lambda))/g(xgrid)) <= (alpha/(1-alpha))])
  
  if(ta==0){ # this is a manual trick, not exact
  abline(v=c(max(X5),max(X4),min(X2),min(X1)),lty=3,col="darkgrey")
  mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(X5))/2,(max(X5)+max(X4))/2,(max(X4)+min(X2))/2,(min(X2)+min(X1))/2,(min(X1)+max(ranges))/2))
  }else{
    abline(v=c(max(X5),-ta,ta,min(X1)),lty=3,col="darkgrey")  
    mtext(c("I","IV","-","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(X5))/2,(max(X5)-ta)/2,(-ta+ta)/2,(ta+min(X1))/2,(min(X1)+max(ranges))/2))
    mtext(c(expression(-t[alpha]),expression(t[alpha])),side=1,at=c(-ta,ta),line=1)
  }
  
  lines(xgrid,R1,col="blue",lty=2)
  lines(xgrid,R2,col="red",lty=2)
  lines(xgrid,R3,col="green",lty=2)
  
  # Active x in regime
  lines(X1[abs(X1)>ta],R1f(X1[abs(X1)>ta],alpha,w,lambda),col="green",lty=1)
  lines(X2[abs(X2)>ta],R2f(X2[abs(X2)>ta],alpha,w,lambda),col="red",lty=1)
  lines(X3[abs(X3)>ta],R3f(X3[abs(X3)>ta],alpha,w,lambda),col="blue",lty=1)
  

dev.off()

# note, could use this code to generate L,U figures as well (replace getU function). Gives the same results