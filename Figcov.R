#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Nov 2019
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
# libraries
library(rmutil) #dlaplace

# source functions
source("functions/getU.R")
source("functions/coverage.R")

#--------------------------------------------------------------------------------------------------------------------#
# Code chunk for plotting multiple w's in one figure
code.chunk <- function(lambda,alpha,wseq,dist,thetamax,h){
  
  
  thetaseq <- seq(lambda,thetamax,h) # same as inside coverage function
  if(dist=="Lap"){distname <- "Laplace(0,1)"}
  if(dist=="Normal"){distname <- "N(0,1)"}
  if(dist=="t3"){distname <- "t(3)"}
  if(dist=="t2"){distname <- "t(2)"}
  if(dist=="t1"){distname <- "t(1)"}
  if(dist=="Cauchy"){distname <- "Cauchy"}
  
  yr <- c(1-3*alpha,1)
  
  plot(thetaseq,thetaseq,type="n",ylim=range(yr),xlab=expression(theta[0]),ylab=expression(C(theta[0])),xlim=c(0,thetamax))
  title(bquote(paste(.(distname),", ",lambda==.(lambda),", ",alpha==.(alpha))))
  abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha,1-5*alpha/2),lty=rep(3,3),col=rep("grey",3))
  text(x=rep(max(thetaseq)-3,5),y=0.03*(max(yr)-min(yr))+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha,1-5*alpha/2),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha),expression(1-5*alpha/2)),cex=1,adj=0)
  abline(v=lambda,lty=3,col="darkgrey")
  text(expression(lambda),x=lambda+0.4,y=0.03*(max(yr)-min(yr))+ min(yr))
  
  for(w in wseq){
    obj <- coverage(alpha,lambda,w,dist,thetamax,h,plot.cov=F)
    for(i in 2:length(thetaseq)){
      segments(x0=thetaseq[i-1],x1=thetaseq[i],y0=obj$C.num[i-1],y1=obj$C.num[i],col=obj$colvec[i])
    }
  }
}
#--------------------------------------------------------------------------------------------------------------------#
# Universal parameters
#h <- 0.001 # theta grid stepsize (change to 0.05 for speed when testing)
h <- 0.1 # theta grid stepsize (change to 0.05 for speed when testing)
thetamax <- 15 # thetaseq endpoint
wseq <- c(0.1,0.25,0.5,0.75,0.9,1)

#--------------------------------------------------------------------------------------------------------------------#
# Visualze panel
pdf("Figures/figcov.pdf",width=12,height=9)

par(mfcol=c(2,3))
#code.chunk(alpha=0.05,lambda=0.5,wseq=0.25,thetamax,h,dist="Lap")
#code.chunk(alpha=0.05,lambda=0.5,wseq=0.25,thetamax,h,dist="t3")
#code.chunk(alpha=0.05,lambda=5,wseq=0.25,thetamax,h,dist="Lap")
#code.chunk(alpha=0.05,lambda=5,wseq=1,thetamax,h,dist="Normal")
#code.chunk(alpha=0.05,lambda=5,wseq=1,thetamax,h,dist="Lap")
#code.chunk(alpha=0.05,lambda=0.5,wseq=1,thetamax,h,dist="Lap")

code.chunk(alpha=0.05,lambda=0.5,wseq,thetamax,h,dist="Normal")
code.chunk(alpha=0.05,lambda=5,wseq,thetamax,h,dist="Normal")

code.chunk(alpha=0.05,lambda=0.5,wseq,thetamax,h,dist="Lap")
code.chunk(alpha=0.05,lambda=5,wseq,thetamax,h,dist="Lap")

code.chunk(alpha=0.05,lambda=0.5,wseq,thetamax,h,dist="t2")
code.chunk(alpha=0.05,lambda=5,wseq,thetamax,h,dist="t2")



dev.off()
