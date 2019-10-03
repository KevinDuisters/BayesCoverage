#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# June 2019
#--------------------------------------------------------------------------------------------------------------------#
# Credible set efficiency: C(\theta_0)/(E_{\theta_0}(nu(X)))

#--------------------------------------------------------------------------------------------------------------------#
# source
source("functions/getU.R")
source("functions/coverage.R")

#--------------------------------------------------------------------------------------------------------------------#
# settings
dist <- "Lap"
alpha <- 0.05
lambda <- 0.75
w <- 1

h <- 0.01
thetaseq <- seq(lambda+h,lambda+18,h)
#--------------------------------------------------------------------------------------------------------------------#
# Coverage
Cseq <- coverage(thetaseq,alpha,lambda,w,dist="Lap",plot.cov=F,cols=c("black","red","blue","orange","green"))
}
