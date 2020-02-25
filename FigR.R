#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# June 2019
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
# source functions
library(rmutil) #dlaplace

#--------------------------------------------------------------------------------------------------------------------#
# Plots
#--------------------------------------------------------------------------------------------------------------------#
# Fig 3: Visualize R functions
dist <- "Lap"

# localize g and G functions
if(dist=="Lap"){
  g <- function(x,theta0=0){dlaplace(x,m=theta0,s=1)}
  G <- function(x,theta0=0){plaplace(x,m=theta0,s=1)}  
  Ginv <- function(p){qlaplace(p)}
}
if(dist=="Normal"){
  g <- function(x,theta0=0){dnorm(x,theta0,1)}
  G <- function(x,theta0=0){pnorm(x,theta0,1)}
  Ginv <- function(p){qnorm(p)}  # make precise by always considering small p as opposed to large  
}
if(dist=="t3"){
  g <- function(x,theta0=0){dt(x,3)}
  G <- function(x,theta0=0){pt(x,3)}  
  Ginv <- function(p){qt(p,3)}
}
if(dist=="t5"){
  g <- function(x,theta0=0){dt(x,5)}
  G <- function(x,theta0=0){pt(x,5)}  
  Ginv <- function(p){qt(p,5)}
}
if(dist=="Cauchy"){
  g <- function(x,theta0=0){dcauchy(x,theta0,1)}
  G <- function(x,theta0=0){pcauchy(x,theta0,1)}  
  Ginv <- function(p){qcauchy(p)}
}

# R functions
R1f <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha/2 - (1-w)/(2*w)*alpha*g(x)-((1-alpha)/2)*(G(lambda-x)-G(-lambda-x)))))}
R2f <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha - (1-w)/w*alpha*g(x)+alpha*(G(lambda-x)-G(-lambda-x)) + G(-lambda-x))))}
R3f <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha/2 - (1-w)/(2*w)*alpha*g(x)+(alpha/2)*(G(lambda-x)-G(-lambda-x)))))}


alpha <- 0.05
lambdas <- c(0.5,5)
w <- 1
#w<-0.2


pdf("Figures/figR.pdf",width=9)

xgrid <- seq(-15,15,0.01)
ranges <- c(-15,15)

par(mfrow=c(1,length(lambdas)),pty="s")
for(lambda in lambdas){
  
  # get R
  R1 <- R1f(xgrid,alpha,w,lambda)
  R2 <- R2f(xgrid,alpha,w,lambda)
  R3 <- R3f(xgrid,alpha,w,lambda)
  
  plot(xgrid,xgrid,type="n",xlab="x",ylab="R",xlim=ranges,asp=1) #main=bquote(lambda==.(lambda))
  
  x1 <- xgrid[xgrid > lambda + R1]
  x2 <- xgrid[(xgrid > 0) & (-lambda + R3 < xgrid) & (xgrid <= lambda + R1)  ]
  x3 <- xgrid[abs(xgrid) <= -lambda + R3]
  x4 <- (-x2)
  x5 <- xgrid[-xgrid > lambda + R1]
  
  
  # exclude |x| <= t_alpha
  ta <-max(0,xgrid[((w/(1-w))*(G(xgrid - lambda) + G(-xgrid - lambda))/g(xgrid)) <= (alpha/(1-alpha))])
  #if(ta!=0){polygon(x=c(-ta,ta,ta,-ta,-ta),y=c(-100,-100,100,100,-100),col="grey")}
  
  if(ta==0){ # this is a manual trick, not exact
  abline(v=c(max(x5),max(x4),min(x2),min(x1)),lty=3)
  mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(x5))/2,(max(x5)+max(x4))/2,(max(x4)+min(x2))/2,(min(x2)+min(x1))/2,(min(x1)+max(ranges))/2))
  }else{
    abline(v=c(max(x5),-ta,ta,min(x1)),lty=3)  
    mtext(c("I","IV","-","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(x5))/2,(max(x5)-ta)/2,(-ta+ta)/2,(ta+min(x1))/2,(min(x1)+max(ranges))/2))
    mtext(c(expression(-t[alpha]),expression(t[alpha])),side=1,at=c(-ta,ta),line=1)
  }
  
  lines(xgrid,R1,lty=2)
  lines(xgrid,R2,col="red",lty=2)
  lines(xgrid,R2f(-xgrid,alpha,w,lambda),col="red",lty=2)
  lines(xgrid,R3,col="blue",lty=2)
  
  # Active x in regime
  lines(x1[abs(x1)>ta],R1f(x1[abs(x1)>ta],alpha,w,lambda),lty=1)
  lines(x2[abs(x2)>ta],R2f(x2[abs(x2)>ta],alpha,w,lambda),col="red",lty=1)
  lines(x3[abs(x3)>ta],R3f(x3[abs(x3)>ta],alpha,w,lambda),col="blue",lty=1)
  lines(x4[abs(x4)>ta],R2f(-x4[abs(x4)>ta],alpha,w,lambda),col="red",lty=1)
  lines(x5[abs(x5)>ta],R1f(x5[abs(x5)>ta],alpha,w,lambda),lty=1)
  
  #legend("topleft",c(expression(R[1](x)),expression(paste(R[2](-x),",",R[2](x))),expression(R[3](x))),lty=c(1,1,1),col=c(1,"red","blue"),bg="white",box.col="white",seg.len=2)
  
}

dev.off()

# check, should give the same as FigLU :
# will replace getU and figLU
# plot U
par(mfrow=c(1,length(lambdas)))
for(lambda in lambdas){
  plot(xgrid,xgrid,type="n",xlab="x",ylab="U",ylim=c(-15,15),main=bquote(lambda==.(lambda)))
  
  lines(xgrid,xgrid+R1f(xgrid,alpha,w,lambda),lty=2)
  lines(xgrid,xgrid+R2f(xgrid,alpha,w,lambda),col="red",lty=2)
  lines(xgrid,xgrid+R3f(xgrid,alpha,w,lambda),col="blue",lty=2)

  x1 <- x2 <- x3 <- x5 <-  numeric(0)
  x1 <- xgrid[xgrid > lambda + R1f(xgrid,alpha,w,lambda)]
  x5 <- xgrid[-xgrid > lambda + R1f(xgrid,alpha,w,lambda)]
  x2 <- xgrid[(xgrid > 0) & (-lambda + R3f(xgrid,alpha,w,lambda) < xgrid) & (xgrid <= lambda + R1f(xgrid,alpha,w,lambda))]
  x3 <- xgrid[abs(xgrid) <= -lambda + R3f(xgrid,alpha,w,lambda)]
  
  abline(h=c(-lambda,lambda),lty=3)
  #abline(v=c(max(x5),min(x3),max(x3),min(x1)),lty=3)
  #mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(xgrid)+max(x5))/2,(max(x5)+min(x3))/2,(min(x3)+max(x3))/2,(max(x3)+min(x1))/2,(min(x1)+max(xgrid))/2))
  abline(v=c(max(x5),max(x4),min(x2),min(x1)),lty=3)
  mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(xgrid)+max(x5))/2,(max(x5)+max(x4))/2,(max(x4)+min(x2))/2,(min(x2)+min(x1))/2,(min(x1)+max(xgrid))/2))
  
  
  if(length(x1)>0){lines(x1,x1+R1f(x1,alpha,w,lambda),lty=1)}
  if(length(x5)>0){lines(x5,x5+R1f(x5,alpha,w,lambda),lty=1)}
  if(length(x2)>0){lines(x2,x2+R2f(x2,alpha,w,lambda),col="red",lty=1)}
  if(length(x2)>0){lines(-x2,rep(-lambda,length(x2)),col="red",lty=1)}
  if(length(x3)>0){lines(x3,x3+R3f(x3,alpha,w,lambda),col="blue",lty=1)}
  
}



# ---------------------------------------------------------------------------------------- #

