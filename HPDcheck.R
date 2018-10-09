# ------------------------------------------------------------------------------------------------------------ 
# Duisters & Schmidt-Hieber, 2018
# Source functions Coverage paper
# ------------------------------------------------------------------------------------------------------------ 

# ------------------------------------------------------------------------------------------------------------ 
# source functions
library(rmutil) #dlaplace
source("getU.R")
source("post.R")

# ------------------------------------------------------------------------------------------------------------ 
# HPD procedure indeed gives same results as analytical expressions
# ------------------------------------------------------------------------------------------------------------ 
dist <- "Normal"
lambda<-0.75
alpha <- 0.05
w <- 1

xgrid <- seq(-lambda-15,lambda+15,0.01)
Ugrid <- Lgrid <- regimeU <- regimeL <- numeric(length(xgrid))

for(i in 1:length(xgrid)){
  x <- xgrid[i]
  out.U <-getU(x,alpha,lambda,w,dist)
  Ugrid[i] <- out.U$val
  regimeU[i] <- out.U$regime
  
  out.L <- getU(-x,alpha,lambda,w,dist)
  Lgrid[i] <- -(out.L$val)
  regimeL[i] <- out.L$regime
}


# ------------------------------------------------------------------------------------------------------------ 
# Check (changepoints) U, L with HPD procedure

Ut <- Lt <- numeric(length(xgrid))
step <- 0.05
for(i in 1:length(xgrid)){
  x <- xgrid[i]
  hpd.mass <-0
  iter <- 0
  Us <- Ls <- x # alternative HPD procedure (not top down), equal for unimodal g (verify!)
  while(hpd.mass<1-alpha & abs(Us)<30){
    
    hpd.mass <- post(0,x,lambda,w,dist)
    
    Ls <- Ls - step
    Us <- Us + step
    if(Ls > -lambda & Ls < lambda ){L <- lambda}else{L <- Ls}
    if(Us < lambda & Us > -lambda){U <- -lambda}else{U <- Us}
    
    cat('\r',paste(round(i/length(xgrid)*100,0),"U=",U," L=",L))
    if(((L==lambda) & (U==(-lambda))) | (U-L)<0.1){hpd.mass <- 0}else{
      hpd.mass <- integrate(post,lower=L,upper=U,x=x,lambda=lambda,w=w,dist=dist)$value + post(0,x,lambda,w,dist)}
  }
  Lt[i] <- L
  Ut[i] <- U
}
par(mfrow=c(1,1))
plot(xgrid,Ut,pch=16,ylim=c(-5,5),xlim=c(-5,5))
points(xgrid,Lt,pch=16,col="blue")
abline(v=xgrid[sapply(5:1,function(r) max(which(regimeU==r)))],lty=3)
abline(v=0,col="red")

plot(xgrid,Ugrid,pch=16,ylim=c(-5,5),xlim=c(-5,5))
points(xgrid,Lgrid,pch=16,col="blue")


# ------------------------------------------------------------------------------------------------------------ 
# check credible set mass produced by getU
# Ugrid, Lgrid
mass.check <- numeric(length(Ugrid))
for(i in 1:length(Ugrid)){
 mass.check[i]<- integrate(post,lower=Lgrid[i],upper=Ugrid[i],x=xgrid[i],lambda=lambda,w=w,dist=dist)$value + post(0,x=xgrid[i],lambda,w,dist)
}
plot(xgrid,mass.check,ylim=c(0,1))
abline(v=xgrid[sapply(5:1,function(r) max(which(regimeU==r)))],lty=3)
abline(h=1-alpha,lty=2)

# Ut,Lt
mass.check <- numeric(length(Ut))
for(i in 1:length(Ut)){
  mass.check[i]<- integrate(post,lower=Lt[i],upper=Ut[i],x=xgrid[i],lambda=lambda,w=w,dist=dist)$value + post(0,x=xgrid[i],lambda,w,dist)
}
plot(xgrid,mass.check,ylim=c(0,1))
abline(v=xgrid[sapply(5:1,function(r) max(which(regimeU==r)))],lty=3)
abline(h=1-alpha,lty=2)
# ------------------------------------------------------------------------------------------------------------ 
