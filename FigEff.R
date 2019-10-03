#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# June 2019
#--------------------------------------------------------------------------------------------------------------------#
# Credible set efficiency: C(\theta_0)/(E_{\theta_0}(nu(X)))

#--------------------------------------------------------------------------------------------------------------------#
# source
library(rmutil) # laplace
source("functions/getU.R")
source("functions/coverage.R")

#--------------------------------------------------------------------------------------------------------------------#
# settings
dist <- "Lap"
draw <- function(n,theta){rlaplace(n,m=theta)}
n <- 1e3
  
alpha <- 0.05
lambda <- 0.75
w <- 1

h <- 0.05
thetaseq <- seq(lambda+h,lambda+10,h)
xgrid <-seq(min(thetaseq)-5,max(thetaseq)+5,h)

#--------------------------------------------------------------------------------------------------------------------#
# Coverage
Cseq <- coverage(thetaseq,alpha,lambda,w,dist,plot.cov=F,cols=c("black","red","blue","orange","green"))

# size
E.nu.seq <- sapply(thetaseq,function(theta){
  xvec <- draw(n,theta)
  Uvec <- sapply(xvec,function(x)getU(x,alpha,lambda,w,dist)$val)
  Lvec <- -sapply(xvec,function(x)getU(-x,alpha,lambda,w,dist)$val)
  ex <- sapply(1:n,function(i){
                    out <- 0
                    if(Uvec[i] >= lambda & Lvec[i] <= - lambda){out <- 2*lambda}
                    return(out)
        })
  return(mean(Uvec-Lvec - ex))}
  )

# plot efficiency
plot(smooth.spline(thetaseq,Cseq$C.sup/E.nu.seq),type="l",ylab="Efficiency",xlab=expression(theta[0]),
     ylim=c(0.1,0.25),xlim=c(0,10))

#--------------------------------------------------------------------------------------------------------------------#
# inverse: observe x. Expected coverage?
# suppose w=1. Then E(theta_0) = x for |x| > lambda; lambda for x in (0,lambda) and - lambda for x in (-lambda,0).
# For w < 1, this needs to be finetuned a bit with weights.
set.seed(1234)
xmat <- sapply(thetaseq,function(theta)draw(n,theta))
colnames(xmat) <- thetaseq

thetamat <- matrix(NA,nrow(xmat),ncol(xmat))
for(j in 1:length(thetaseq)){
  thetamat[,j] <- thetaseq[j]
}

# expected coverage given x
s <- 0.05
E.cov.grid <- sapply(xgrid,function(x){
  thetavec <- as.numeric(thetamat[(xmat > x - s) & (xmat < x+s)])
  print(100*(which(x==xgrid)/length(xgrid)))
  return(table(thetavec)%*%coverage(unique(thetavec),alpha,lambda,w,dist,plot.cov=F,cols=c("black","red","blue","orange","green"))$C.sup/length(thetavec))
})

# size
nu.x <- sapply(xgrid,function(x){
  U <- getU(x,alpha,lambda,w,dist)$val
  L <- -getU(-x,alpha,lambda,w,dist)$val
  return(U-L-(U>=lambda & L <= - lambda)*2*lambda)
})

# plot expected efficiency
plot(smooth.spline(xgrid,E.cov.grid/nu.x),type="l",ylab="Efficiency",xlab="x",ylim=c(0.1,0.25),col="blue")
plot(smooth.spline(xgrid,nu.x),type="l")





