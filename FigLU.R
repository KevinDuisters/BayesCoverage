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
cx <- 1 # label scale
#--------------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------------#
# Universal code chunk for plotting for some w and lambda
code.chunk <- function(lambda,w,alpha,output,h=0.01){
  
  xgrid <- seq(-20-lambda,20+lambda,h)
  ranges <- c(-10,10)
  Ugrid <- Lgrid <- regimeU <- regimeL <- numeric(length(xgrid))
  #cols <- c("green","red","blue","red","green")
  cols <- rep("black",5)
  
  for(i in 1:length(xgrid)){
    x <- xgrid[i]
    out.U <-getU(x,alpha,lambda,w,dist="Lap")
    Ugrid[i] <- out.U$val
    regimeU[i] <- out.U$regime
    
    out.L <- getU(-x,alpha,lambda,w,dist="Lap")
    Lgrid[i] <- -(out.L$val)
    regimeL[i] <- out.L$regime
  }
  
  
  if(output=="LU"){
  theta0seq <- c(seq(-lambda-20,-lambda-h,h),seq(lambda+h,lambda+20,h))
  xU.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Ugrid>=theta0))]})
  xU.sup <- sapply(theta0seq,function(theta0){xgrid[max(which(Ugrid<=theta0))]})
  xL.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Lgrid>=theta0))]})
  xL.sup <-sapply(theta0seq,function(theta0){xgrid[max(which(Lgrid<=theta0))]})
  
  
  plot(xgrid,xgrid,ylab=expression(paste(L[alpha](x)," and ",U[alpha](x))),xlab=expression(x),type="n",xlim=ranges,ylim=1.5*ranges,cex.lab=cx)
  abline(h=c(-lambda,lambda),lty=3,col="darkgrey")
  ta <-max(0,xgrid[((w/(1-w))*(plaplace(xgrid - lambda) + plaplace(-xgrid - lambda))/dlaplace(xgrid)) <= (alpha/(1-alpha))])
  
  if(ta==0){ # this is a manual trick, not exact
    abline(v=c(max(xgrid[regimeU==5]),max(xgrid[regimeU==4]),min(xgrid[regimeU==2]),min(xgrid[regimeU==1])),lty=3,col="darkgrey")
    mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])+max(xgrid[regimeU==4]))/2,(max(xgrid[regimeU==4])+min(xgrid[regimeU==2]))/2,(min(xgrid[regimeU==2])+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2),cex=cx)
  }else{
    abline(v=c(max(xgrid[regimeU==5]),-ta,ta,min(xgrid[regimeU==1])),lty=3,col="darkgrey")
    mtext(c("I","IV","","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])-ta)/2,(-ta+ta)/2,(ta+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2),cex=cx)
    mtext(c(expression(-t[alpha]),expression(t[alpha])),side=1,at=c(-ta,ta),line=1,cex=cx)
  }
  
  
  for(r in 0:5){
    reg <- which(regimeU==r)
    lines(xgrid[reg],Ugrid[reg],col=cols[r])
    
    reg <- which(regimeL==r)
    lines(xgrid[reg],Lgrid[reg],col=cols[r],lty=2)
  }
  
  mtext(c(expression(-lambda),expression(lambda)),line=-1.5,side=2,at=c(-lambda,lambda),cex=cx)
  legend("topleft",c(expression(U[alpha]),expression(L[alpha])),lty=c(1,2),col=c("black","black"),horiz=F,bg="white",box.col="white",seg.len=2,cex=cx)
  box()
  title(bquote(paste("Laplace(0,1), ",lambda==.(lambda),", ",w==.(w))))
  }
  
  if(output=="size"){
    plot(xgrid,xgrid,type="n",xlab="x",ylab=expression(paste("length of ",HPD[alpha](x))),xlim=ranges,ylim=c(0,2*qlaplace(1-alpha/2)),cex.lab=cx)
    abline(h=2*qlaplace(1-alpha/2),lty=2,col="red")
    ta <-max(0,xgrid[((w/(1-w))*(plaplace(xgrid - lambda) + plaplace(-xgrid - lambda))/dlaplace(xgrid)) <= (alpha/(1-alpha))])
    
    if(ta==0){ # this is a manual trick, not exact
      abline(v=c(max(xgrid[regimeU==5]),max(xgrid[regimeU==4]),min(xgrid[regimeU==2]),min(xgrid[regimeU==1])),lty=3,col="darkgrey")
      mtext(c("I","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])+max(xgrid[regimeU==4]))/2,(max(xgrid[regimeU==4])+min(xgrid[regimeU==2]))/2,(min(xgrid[regimeU==2])+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2),cex=cx)
    }else{
      abline(v=c(max(xgrid[regimeU==5]),-ta,ta,min(xgrid[regimeU==1])),lty=3,col="darkgrey")
      mtext(c("I","IV","","II","I"),side=1,adj=0.5,line=-1.5,at=c((min(ranges)+max(xgrid[regimeU==5]))/2,(max(xgrid[regimeU==5])-ta)/2,(-ta+ta)/2,(ta+min(xgrid[regimeU==1]))/2,(min(xgrid[regimeU==1])+max(ranges))/2),cex=cx)
      mtext(c(expression(-t[alpha]),expression(t[alpha])),side=1,at=c(-ta,ta),line=1,cex=cx)
    }
    
    sizes <- numeric(length(xgrid))
    for(i in 1:length(xgrid)){
      Ui <- Ugrid[i]
      Li <- Lgrid[i]
      if(Ui > lambda & Li < (-lambda) ){sizes[i] <- Ui - Li - 2*lambda}else{
        sizes[i] <- Ui - Li
      }
    }
    lines(xgrid,sizes)
    title(bquote(paste("Laplace(0,1), ",lambda==.(lambda),", ",w==.(w))))
    
    
  }
}

#--------------------------------------------------------------------------------------------------------------------#
# Visualize LU
#pdf("Figures/figLU.pdf",width=12,height=4.5)
png("Figures/figLU.png",width=12,height=4.5,unit="in",res=300)

par(mfrow=c(1,3))
code.chunk(lambda=0.5,w=0.25,alpha,output="LU") # left
code.chunk(lambda=5,w=0.25,alpha,output="LU") # mid
code.chunk(lambda=5,w=1,alpha,output="LU") # right

dev.off()
#------------------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------------------------#
# Visualize LU
#pdf("Figures/figsize.pdf",width=12,height=4.5)
png("Figures/figsize.png",width=12,height=4.5,unit="in",res=300)

par(mfrow=c(1,3))
code.chunk(lambda=0.5,w=0.25,alpha,output="size") # left
code.chunk(lambda=5,w=0.25,alpha,output="size") # mid
code.chunk(lambda=5,w=1,alpha,output="size") # right

dev.off()
#------------------------------------------------------------------------------------------------------------------------#




