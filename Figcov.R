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
# Universal parameters
alpha <- 0.05
h <- 0.01 # theta grid stepsize (slows down code!) for testing
#h <- 0.0001 # final version of figures (30 mins)
thetamax <- 15 # thetaseq endpoint
#--------------------------------------------------------------------------------------------------------------------#


pdf("Figures/figcov.pdf",width=12,height=9)

par(mfcol=c(2,3))
coverage(alpha,lambda=0.5,w=0.25,thetamax,h,dist="Lap",plot.cov=T)
coverage(alpha,lambda=0.5,w=0.25,thetamax,h,dist="t3",plot.cov=T)
coverage(alpha,lambda=5,w=0.25,thetamax,h,dist="Lap",plot.cov=T)
coverage(alpha,lambda=5,w=1,thetamax,h,dist="Normal",plot.cov=T)
coverage(alpha,lambda=5,w=1,thetamax,h,dist="Lap",plot.cov=T)
coverage(alpha,lambda=0.5,w=1,thetamax,h,dist="Lap",plot.cov=T)


dev.off()
