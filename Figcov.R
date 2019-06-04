#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# May 2019
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
# libraries
library(rmutil) #dlaplace

# source functions
source("functions/getU.R")
source("functions/coverage.R")

#--------------------------------------------------------------------------------------------------------------------#
# Plot Coverage
alpha <- 0.05
w <- 1
#w=0.2 # Supplemental

if(w==1){pdf("Figures/figcov.pdf",width=9)}else{pdf("Figures/SFcov.pdf",width=9)}
par(mfcol=c(2,2),xpd=F,mar=c(4,4,2,2))

lambdas <- c(0.75,7.5)
for(lambda in lambdas){
  thetaseq <- seq(lambda+0.005,lambda+10,0.005)
  coverage(thetaseq,alpha,lambda,w,dist="Normal",plot.cov=T,cols=c("black","red","blue","orange","green")) 
  coverage(thetaseq,alpha,lambda,w,dist="Lap",plot.cov=T,cols=c("black","red","blue","orange","green"))
}

dev.off()

# bottom
if(w==1){pdf("Figures/figcovbottom.pdf",width=9)}else{pdf("Figures/SFcovbottom.pdf",width=9)}
par(mfrow=c(2,2),xpd=F,mar=c(4,4,2,2))

lambdas <- c(0.75,7.5)
for(lambda in lambdas){
  thetaseq <- seq(lambda+0.005,lambda+10,0.005)
  coverage(thetaseq,alpha,lambda,w,dist="t5",plot.cov=T,cols=c("black","red","blue","orange","green")) 
  #coverage(thetaseq,alpha,lambda,w,dist="t3",plot.cov=T,cols=c("black","red","blue","orange","green")) 
  #coverage(thetaseq,alpha,lambda,w,dist="Cauchy",plot.cov=T,cols=c("black","red","blue","orange","green"))
}

dev.off()

