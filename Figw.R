#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Feb 2020
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
# Coverage for changing w

#--------------------------------------------------------------------------------------------------------------------#
library(rmutil) #dlaplace
# source functions
source("functions/getU.R")
source("functions/coverage.R")
#--------------------------------------------------------------------------------------------------------------------#

# universal settings
dist <- "Lap"


#h <- 0.001 # theta grid stepsize (change to 0.05 for speed when testing)
h <- 0.05
thetamax <- 15 # thetaseq endpoint


#wseq <- c(1e-5,1e-4,1e-3,1e-2,0.025,0.05,0.075,0.1,0.2,0.5,1)
wseq <- c(0.01,0.1,0.25,0.5,0.75,0.9,1)

#--------------------------------------------------------------------------------------------------------------------#
# Code to plot
code.chunk <- function(lambda,alpha,wseq,dist,thetamax,h){

  
thetaseq <- seq(lambda,thetamax,h) # same as inside coverage function
if(dist=="Lap"){distname <- "Laplace(0,1)"}
if(dist=="Normal"){distname <- "N(0,1)"}

plot(thetaseq,thetaseq,type="n",ylim=c(0.84,1),xlab=expression(theta[0]),ylab=expression(C(theta[0])),xlim=c(0,thetamax))
title(bquote(paste(.(distname),", ",lambda==.(lambda),", ",alpha==.(alpha))))
abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2),lty=rep(3,3),col=rep("grey",3))
text(x=rep(max(thetaseq)-2,4),y=0.005+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=1,adj=0)
abline(v=lambda,lty=3,col="darkgrey")
text(expression(lambda),x=lambda+0.4,y=1-2*alpha+0.005)

for(w in wseq){
  obj <- coverage(alpha,lambda,w,dist,thetamax,h,plot.cov=F)
  for(i in 2:length(thetaseq)){
    segments(x0=thetaseq[i-1],x1=thetaseq[i],y0=obj$C.num[i-1],y1=obj$C.num[i],col=obj$colvec[i])
  }
}
}
#--------------------------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------------------------#
# Visualize coverage for different w
pdf("Figures/figw.pdf",width=12,height=4.5)

par(mfrow=c(1,3))
code.chunk(lambda=0.5,alpha=0.05,wseq,dist,thetamax,h)
code.chunk(lambda=5,alpha=0.05,wseq,dist,thetamax,h)
code.chunk(lambda=5,alpha=0.01,wseq,dist,thetamax,h)


dev.off()