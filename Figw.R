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
h <- 0.1

wseq <- c(0.00001,0.001,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5,0.75,1)

#pdf("Figures/figw.pdf",width=9)
par(mfrow=c(1,2))

lambdas <- c(0.75,7.5)
for(lambda in lambdas){
  thetaseq <- seq(lambda+h,lambda+10,h)
  plot(thetaseq,thetaseq,type="n",ylim=c(0,1),xlab=expression(theta[0]),ylab=expression(C(theta[0])))
  for(w in wseq){
  obj <- coverage(thetaseq,alpha,lambda,w,dist,plot.cov=F) 
  for(r in 5:1){
  lines(thetaseq[obj$regimeCinf==r],obj$C.inf[obj$regimeCinf==r],col=c("black","red","blue","orange","green")[r])
  }
  
  }
}

#dev.off()