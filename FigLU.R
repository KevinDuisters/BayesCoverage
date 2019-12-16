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
# Plots
#--------------------------------------------------------------------------------------------------------------------#
# Fig: Visualize L,U and inverse xL, xU 
#--------------------------------------------------------------------------------------------------------------------#
# set density of interest
dist <- "Lap"
#dist <- "t3"
#dist <- "Normal"
#--------------------------------------------------------------------------------------------------------------------#
# code chunk to localize g, G, Ginv to memory
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
#--------------------------------------------------------------------------------------------------------------------#
# set parameters
alpha <- 0.05
#w <- 1
w<-0.2
##--------------------------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------------------------#
# Visualize

if(w==1){pdf("Figures/figLU.pdf",width=9)}else{pdf("Figures/figSI1.pdf",width=9)}  # w=0.2}

right <- F # excluded inversion figures (December '19)
#par(mfrow=c(1,1))
par(mfrow=c(1,2))

ranges <- c(-10,10)

#--------------------------------------------------------------------------------------------------------------------#
# New Left (lambda = 0.5)
lambda <- 0.5
h <- 0.01
xgrid <- seq(-20-lambda,20+lambda,h)
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


theta0seq <- c(seq(-lambda-20,-lambda-h,h),seq(lambda+h,lambda+20,h))
xU.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Ugrid>=theta0))]})
xU.sup <- sapply(theta0seq,function(theta0){xgrid[max(which(Ugrid<=theta0))]})
xL.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Lgrid>=theta0))]})
xL.sup <-sapply(theta0seq,function(theta0){xgrid[max(which(Lgrid<=theta0))]})


plot(xgrid,xgrid,ylab=expression(paste(L[alpha](x)," and ",U[alpha](x))),xlab=expression(x),type="n",xlim=ranges,ylim=ranges,asp=1)
abline(h=c(-lambda,lambda),lty=3,col="grey")
ta <-max(0,xgrid[((w/(1-w))*(G(xgrid - lambda) + G(-xgrid - lambda))/g(xgrid)) <= (alpha/(1-alpha))])

if(ta==0){ # this is a manual trick, not exact
  abline(v=c(max(xgrid[regimeU==5]),max(xgrid[regimeU==4]),min(xgrid[regimeU==2]),min(xgrid[regimeU==1])),lty=3,col="grey")
  mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])+max(xgrid[regimeU==4]))/2,(max(xgrid[regimeU==4])+min(xgrid[regimeU==2]))/2,(min(xgrid[regimeU==2])+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2))
}else{
  abline(v=c(max(xgrid[regimeU==5]),-ta,ta,min(xgrid[regimeU==1])),lty=3,col="grey")
  mtext(c("I","IV","-","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])-ta)/2,(-ta+ta)/2,(ta+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2))
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

# manual things to make plot look nice (be careful in generalizing this code)
box()
if(w==1){
  points(y=c(-lambda,lambda,-lambda,lambda),x=c(max(xgrid[regimeU==4]),min(xgrid[regimeU==3]),max(xgrid[regimeL==3]),min(xgrid[regimeL==4])),
         col=c("black","black","blue","blue"),pch=rep(21,4),bg=c("black","white","blue","white")) # discontinuity jump at lambda
}
if(w < 1){
  if(ta==0){points(y=c(-lambda,lambda,-lambda,lambda),x=c(min(xgrid[regimeU==3]),min(xgrid[regimeU==3]),max(xgrid[regimeL==3]),max(xgrid[regimeL==3])),
                   col=c("black","black","blue","blue"),pch=rep(21,4),bg=c("black","white","blue","white")) # discontinuity jump at lambda
  }else{
    points(y=c(-lambda,0,0,lambda),x=c(-ta,-ta,ta,ta),
           col=c("black","blue","black","blue"),pch=rep(21,4),bg=c("black","white","black","white")) # discontinuity jump at lambda
  }
}

#------------------------------------------------------------------------------------------------------------------------#

# Left
lambda <- 5
h <- 0.01
xgrid <- seq(-20-lambda,20+lambda,h)
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


theta0seq <- c(seq(-lambda-20,-lambda-h,h),seq(lambda+h,lambda+20,h))
xU.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Ugrid>=theta0))]})
xU.sup <- sapply(theta0seq,function(theta0){xgrid[max(which(Ugrid<=theta0))]})
xL.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Lgrid>=theta0))]})
xL.sup <-sapply(theta0seq,function(theta0){xgrid[max(which(Lgrid<=theta0))]})


plot(xgrid,xgrid,ylab=expression(paste(L[alpha](x)," and ",U[alpha](x))),xlab=expression(x),type="n",xlim=ranges,ylim=ranges,asp=1)
abline(h=c(-lambda,lambda),lty=3,col="grey")
ta <-max(0,xgrid[((w/(1-w))*(G(xgrid - lambda) + G(-xgrid - lambda))/g(xgrid)) <= (alpha/(1-alpha))])

if(ta==0){ # this is a manual trick, not exact
  abline(v=c(max(xgrid[regimeU==5]),max(xgrid[regimeU==4]),min(xgrid[regimeU==2]),min(xgrid[regimeU==1])),lty=3,col="grey")
  mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])+max(xgrid[regimeU==4]))/2,(max(xgrid[regimeU==4])+min(xgrid[regimeU==2]))/2,(min(xgrid[regimeU==2])+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2))
}else{
  abline(v=c(max(xgrid[regimeU==5]),-ta,ta,min(xgrid[regimeU==1])),lty=3,col="grey")
  mtext(c("I","IV","-","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])-ta)/2,(-ta+ta)/2,(ta+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2))
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

# manual things to make plot look nice (be careful in generalizing this code)
box()
if(w==1){
points(y=c(-lambda,lambda,-lambda,lambda),x=c(max(xgrid[regimeU==4]),min(xgrid[regimeU==3]),max(xgrid[regimeL==3]),min(xgrid[regimeL==4])),
       col=c("black","black","blue","blue"),pch=rep(21,4),bg=c("black","white","blue","white")) # discontinuity jump at lambda
}
if(w < 1){
  if(ta==0){points(y=c(-lambda,lambda,-lambda,lambda),x=c(min(xgrid[regimeU==3]),min(xgrid[regimeU==3]),max(xgrid[regimeL==3]),max(xgrid[regimeL==3])),
                   col=c("black","black","blue","blue"),pch=rep(21,4),bg=c("black","white","blue","white")) # discontinuity jump at lambda
  }else{
    points(y=c(-lambda,0,0,lambda),x=c(-ta,-ta,ta,ta),
           col=c("black","blue","black","blue"),pch=rep(21,4),bg=c("black","white","black","white")) # discontinuity jump at lambda
  }
}

#------------------------------------------------------------------------------------------------------------------------#
if(right==T){

# Right
plot(theta0seq,xU.sup,type="n",xlim=ranges,ylim=ranges,ylab=expression(paste(X[U],(theta[0])," and ",{X[L]}(theta[0]))),col="black",xlab=expression(theta[0]),asp=1)
if(ta==0){ # this is a manual trick, not exact
  abline(h=c(max(xgrid[regimeU==5]),max(xgrid[regimeU==4]),min(xgrid[regimeU==2]),min(xgrid[regimeU==1])),lty=3,col="grey")
  mtext(c("I","IV","III","II","I"),side=2,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])+max(xgrid[regimeU==4]))/2,(max(xgrid[regimeU==4])+min(xgrid[regimeU==2]))/2,(min(xgrid[regimeU==2])+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2))
}else{
  abline(h=c(max(xgrid[regimeU==5]),-ta,ta,min(xgrid[regimeU==1])),lty=3,col="grey")
  mtext(c("I","IV","-","II","I"),side=2,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])-ta)/2,(-ta+ta)/2,(ta+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2))
  mtext(c(expression(-t[alpha]),expression(t[alpha])),side=2,at=c(-ta,ta),line=1)
}
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=-1.5,side=1,at=c(-lambda,lambda))

polygon(x=c(theta0seq[theta0seq<(-lambda)],sort(theta0seq[theta0seq<(-lambda)],decreasing=T)),y=c(xL.sup[theta0seq<(-lambda)],(xL.inf[theta0seq<(-lambda)])[order(theta0seq[theta0seq<(-lambda)],decreasing=T)]),col="blue",border=NA,density=20)
polygon(x=c(theta0seq[theta0seq>(lambda)],sort(theta0seq[theta0seq>(lambda)],decreasing=T)),y=c(xU.sup[theta0seq>(lambda)],(xU.inf[theta0seq>(lambda)])[order(theta0seq[theta0seq>(lambda)],decreasing=T)]),col="black",border=NA,density=20)





#plot(theta0seq[theta0seq>(lambda)],xU.sup[theta0seq>(lambda)],type="l")
#lines(theta0seq[theta0seq>(lambda)],xU.inf[theta0seq>(lambda)],type="l",col="red")
#polygon(c(theta0seq[theta0seq>lambda],rev(theta0seq[theta0seq>lambda])),c(xU.inf[theta0seq>lambda],rev(xU.sup[theta0seq>lambda])))

lines(theta0seq[theta0seq<(-lambda)],xU.sup[theta0seq<(-lambda)])
lines(theta0seq[theta0seq>(lambda)],xU.sup[theta0seq>(lambda)])
lines(theta0seq[theta0seq<(-lambda)],xU.inf[theta0seq<(-lambda)],lty=2)
lines(theta0seq[theta0seq>(lambda)],xU.inf[theta0seq>(lambda)],lty=2)
lines(theta0seq[theta0seq<(-lambda)],xL.sup[theta0seq<(-lambda)],col="blue")
lines(theta0seq[theta0seq>(lambda)],xL.sup[theta0seq>(lambda)],col="blue")
lines(theta0seq[theta0seq<(-lambda)],xL.inf[theta0seq<(-lambda)],lty=2,col="blue")
lines(theta0seq[theta0seq>(lambda)],xL.inf[theta0seq>(lambda)],lty=2,col="blue")

legend("topleft",c(expression(paste("sup ",X[U])),expression(paste("inf ",X[U])),expression(paste("sup ",X[L])),expression(paste("inf ",X[L]))),
       lty=c(1,2,1,2),col=c("black","black","blue","blue"),ncol=2,bg="white",box.col="white",seg.len=2)


# manual things to make plot look nice (be careful in generalizing this code)
box()
if(w==1){
points(x=c(-lambda,lambda,-lambda,lambda),y=c(min(xgrid[regimeU==4]),min(xgrid[regimeU==3]),max(xgrid[regimeL==3]),max(xgrid[regimeL==4])),
       col=c("black","black","blue","blue"),pch=rep(21,4),bg=c("black","white","blue","white")) 
}
if(w < 1){
  if(ta==0){
  points(x=c(-lambda,lambda,-lambda,lambda),y=c(min(xgrid[regimeU==4]),max(xgrid[regimeU==3]),min(xgrid[regimeL==3]),max(xgrid[regimeL==4])),
         col=c("black","black","blue","blue"),pch=rep(21,4),bg=c("black","white","blue","white"))
  segments(x0=0,x1=0,y0=min(xgrid[Ugrid==0]),y1=max(xgrid[Ugrid==0]),lwd=1.5,col="blue")
  segments(x0=0,x1=0,y0=min(xgrid[Lgrid==0]),y1=max(xgrid[Lgrid==0]),lwd=1.5,col="blue")
  points(x=c(0,0,0,0),y=c(min(xgrid[Ugrid==0]),min(xgrid[Ugrid==0]),min(xgrid[Lgrid==0]),max(xgrid[Lgrid==0])),col=c(1,"blue",1,"blue"),pch=rep(21,4),bg=c("white","black","white","blue"))
  }else{
    points(x=c(-lambda,lambda,-lambda,lambda),y=c(min(xgrid[regimeU==4]),ta,-ta,max(xgrid[regimeL==4])),
           col=c("black","black","blue","blue"),pch=rep(21,4),bg=c("black","white","blue","white"))
    segments(x0=0,x1=0,y0=min(xgrid[Ugrid==0]),y1=max(xgrid[Ugrid==0]),lwd=1.5,col="blue")
    segments(x0=0,x1=0,y0=min(xgrid[Lgrid==0]),y1=max(xgrid[Lgrid==0]),lwd=1.5,col="blue")
    points(x=c(0,0,0,0),y=c(min(xgrid[Ugrid==0]),min(xgrid[Ugrid==0]),min(xgrid[Lgrid==0]),max(xgrid[Lgrid==0])),col=c(1,"blue",1,"blue"),pch=rep(21,4),bg=c("white","black","white","blue"))
  }
}

} # end if right==T

dev.off()


