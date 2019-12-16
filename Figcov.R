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
# Plot Coverage
alpha <- 0.05
cols <- c("grey","green","red","blue","orange","magenta")
w <- 1
#w <- 0.2 # Supplemental

if(w==1){pdf("Figures/figcov.pdf",width=9)}else{pdf("Figures/SFcov.pdf",width=9)}
par(mfcol=c(2,2),xpd=F,mar=c(4,4,2,2))

lambdas <- c(0.5,5)
for(lambda in lambdas){
  thetaseq <- seq(lambda+0.005,lambda+14,0.005)
  coverage(thetaseq,alpha,lambda,w,dist="Normal",plot.cov=T,cols) 
  coverage(thetaseq,alpha,lambda,w,dist="Lap",plot.cov=T,cols) 
}

dev.off()

# bottom
if(w==1){pdf("Figures/figcovbottom.pdf",width=9)}else{pdf("Figures/SFcovbottom.pdf",width=9)}
par(mfrow=c(2,2),xpd=F,mar=c(4,4,2,2))

lambdas <- c(0.5,5)
for(lambda in lambdas){
  thetaseq <- seq(lambda+0.005,lambda+14,0.005)
  coverage(thetaseq,alpha,lambda,w,dist="t3",plot.cov=T,cols) 
}

dev.off()



