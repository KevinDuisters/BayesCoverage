#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Oct 2019
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
# Coverage for changing w

#--------------------------------------------------------------------------------------------------------------------#
library(rmutil) #dlaplace
# source functions
source("functions/getU.R")
source("functions/coverage.R")
#--------------------------------------------------------------------------------------------------------------------#

dist <- "Lap"
alpha <- 0.05
cols <- c("grey","black","red","blue","orange","green")
h <- 0.005

wseq <- c(1e-5,1e-4,1e-3,1e-2,0.025,0.05,0.075,0.1,0.2,0.5,1)

pdf("Figures/figw.pdf",width=9)
par(mfrow=c(1,2))

lambdas <- c(0.75,7.5)
qmat <- matrix(NA,nrow=length(wseq),ncol=length(lambdas),dimnames=list(rows=wseq,cols=lambdas))
for(lambda in lambdas){
  thetaseq <- seq(lambda+h,lambda+18,h)
  plot(thetaseq,thetaseq,type="n",ylim=c(0,1),xlab=expression(theta[0]),ylab=expression(C(theta[0])))
  abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),lty=rep(3,4),col=rep("grey",4))
  text(x=rep(lambda+16,4),y=0.005+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=0.8,adj=0)
  for(w in wseq){
  qmat[which(wseq==w),which(lambdas==lambda)] <- qlaplace(min(1-1e-8,alpha/(1+alpha)*(1 + ((1-w)/w)*dlaplace(lambda))))  
  obj <- coverage(thetaseq,alpha,lambda,w,dist,plot.cov=F) 
  lines(thetaseq,obj$C.inf,col=cols[0+1]) # detach this for better visual appearance
  for(r in 5:1){
  lines(thetaseq[obj$regimeCinf==r],obj$C.inf[obj$regimeCinf==r],col=cols[r+1])
  }
  
  }
}

dev.off()