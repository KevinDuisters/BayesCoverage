# ------------------------------------------------------------------------------------------------------------ 
# Duisters & Schmidt-Hieber, 2018
# Source functions Coverage paper
# ------------------------------------------------------------------------------------------------------------ 

library(rmutil) #dlaplace

# ------------------------------------------------------------------------------------------------------------ 
# get U function
source("getU.R")

# ------------------------------------------------------------------------------------------------------------ 
# Plots
# ------------------------------------------------------------------------------------------------------------ 
# Fig 2: Visualize L,U and inverse xtilde, xstar 
dist <- "Normal"
#lambda<-7.5
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

theta0seq <- c(seq(-lambda-10,-lambda-0.01,0.01),seq(lambda+0.01,lambda+10,0.01))
xstar.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Ugrid>theta0))]})
xstar.sup <- sapply(theta0seq,function(theta0){xgrid[max(which(Ugrid<theta0))]})
xtilde.sup <-sapply(theta0seq,function(theta0){xgrid[max(which(Lgrid<theta0))]})
xtilde.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Lgrid>theta0))]})
#xtilde.inf <- xtilde.sup

# for "Fig2.R"
loc <- which(theta0seq==1.75)
c(xstar.inf[loc],xtilde.sup[loc])
loc2 <- which(theta0seq==4)
c(xstar.inf[loc2],xtilde.sup[loc2])


# ----- 

output <- T
if(output==T){
#pdf("fig2.pdf",width=9)
  pdf("figSI1.pdf",width=9)  # w=0.2

ranges <- c(-10-lambda,10+lambda)
par(mfrow=c(1,2))
plot(xgrid,xgrid,ylab=expression(paste(L[alpha](x)," and ",U[alpha](x))),xlab=expression(x),type="n",xlim=ranges,ylim=ranges,asp=1)
polygon(y=c(-lambda,lambda,lambda,-lambda,-lambda),x=c(1.04*min(ranges),1.04*min(ranges),1.04*max(ranges),1.04*max(ranges),1.04*min(ranges)),col="grey90",border=NA)
abline(h=c(-lambda,lambda),lty=3,col="grey")
abline(v=xgrid[sapply(5:1,function(r) max(which(regimeU==r)))],lty=3)
mtext(c("V","IV","III","II","I"),side=1,adj=0.5,line=-2,at=c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU==r)))][-5])+diff(c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU==r)))][-5],ranges[2]))/2   )

for(r in 1:5){
  reg <- which(regimeU==r)
  lines(xgrid[reg],Ugrid[reg],col="black")
  
  reg <- which(regimeL==r)
  lines(xgrid[reg],Lgrid[reg],col="blue")
}

mtext(c(expression(-lambda),expression(lambda)),line=-2,side=2,at=c(-lambda,lambda))
legend("topleft",c(expression(U[alpha]),expression(L[alpha])),lty=c(1,1),col=c("black","blue"),horiz=F,bty="n",seg.len=2)


plot(theta0seq,xstar.sup,type="n",xlim=ranges,ylim=ranges,ylab=expression(paste(x,"*",(theta[0])," and ",tilde(x)(theta[0]))),col="black",xlab=expression(theta[0]),asp=1)
polygon(x=c(-lambda,lambda,lambda,-lambda,-lambda),y=c(1.04*min(ranges),1.04*min(ranges),1.04*max(ranges),1.04*max(ranges),1.04*min(ranges)),col="grey90",border=NA)
abline(h=xgrid[sapply(5:1,function(r) max(which(regimeU==r)))],lty=3)
mtext(c("V","IV","III","II","I"),side=2, adj=0.5, line=-2,at=c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU==r)))][-5])+diff(c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU==r)))][-5],ranges[2]))/2   )

polygon(x=c(theta0seq[theta0seq<(-lambda)],sort(theta0seq[theta0seq<(-lambda)],decreasing=T)),y=c(xtilde.sup[theta0seq<(-lambda)],(xtilde.inf[theta0seq<(-lambda)])[order(theta0seq[theta0seq<(-lambda)],decreasing=T)]),col="blue",border=NA,density=20)
polygon(x=c(theta0seq[theta0seq>(lambda)],sort(theta0seq[theta0seq>(lambda)],decreasing=T)),y=c(xstar.sup[theta0seq>(lambda)],(xstar.inf[theta0seq>(lambda)])[order(theta0seq[theta0seq>(lambda)],decreasing=T)]),col="black",border=NA,density=20)

lines(theta0seq[theta0seq<(-lambda)],xstar.sup[theta0seq<(-lambda)])
lines(theta0seq[theta0seq>(lambda)],xstar.sup[theta0seq>(lambda)])
lines(theta0seq[theta0seq<(-lambda)],xstar.inf[theta0seq<(-lambda)],lty=2)
lines(theta0seq[theta0seq>(lambda)],xstar.inf[theta0seq>(lambda)],lty=2)
lines(theta0seq[theta0seq<(-lambda)],xtilde.sup[theta0seq<(-lambda)],col="blue")
lines(theta0seq[theta0seq>(lambda)],xtilde.sup[theta0seq>(lambda)],col="blue")
lines(theta0seq[theta0seq<(-lambda)],xtilde.inf[theta0seq<(-lambda)],lty=2,col="blue")
lines(theta0seq[theta0seq>(lambda)],xtilde.inf[theta0seq>(lambda)],lty=2,col="blue")
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=-2,side=1,at=c(-lambda,lambda))
#abline(h=c(-x1,-x2,x2,x1),lty=3,col="grey")
#mtext(c(expression(-x[1]),expression(-x[2]),expression(x[2]),expression(x[1])),line=-2,side=2,at=c(-x1,-x2,x2,x1))

legend("topleft",c(expression(paste("sup ",x,"*")),expression(paste("inf ",x,"*")),expression(paste("sup ",tilde(x))),expression(paste("inf ",tilde(x)))),
       lty=c(1,2,1,2),col=c("black","black","blue","blue"),ncol=2,bty="n",
       seg.len=2)


dev.off()
}

# ------------------------------------------------------------------------------------------------------------ 
# Check (changhepoints) U, L with HPD procedure
source("post.R")

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
plot(xgrid,Ut,pch=16)
points(xgrid,Lt,pch=16,col="blue")
abline(v=xgrid[sapply(5:1,function(r) max(which(regimeU==r)))],lty=3)
abline(v=0,col="red")

# ------------------------------------------------------------------------------------------------------------ 
# check credible set mass produced by getU
# Ugrid, Lgrid
masscheck <- numeric(length(Ugrid))
for(i in 1:length(Ugrid)){
 masscheck[i]<- integrate(post,lower=Lgrid[i],upper=Ugrid[i],x=xgrid[i],lambda=lambda,w=w,dist=dist)$value + post(0,x=xgrid[i],lambda,w,dist)
}
plot(xgrid,masscheck,ylim=c(0,1))
abline(v=xgrid[sapply(5:1,function(r) max(which(regimeU==r)))],lty=3)
abline(h=1-alpha,lty=2)

# Ut,Lt
masscheck <- numeric(length(Ut))
for(i in 1:length(Ugrid)){
  masscheck[i]<- integrate(post,lower=Lt[i],upper=Ut[i],x=xgrid[i],lambda=lambda,w=w,dist=dist)$value + post(0,x=xgrid[i],lambda,w,dist)
}
plot(xgrid,masscheck,ylim=c(0,1))
abline(v=xgrid[sapply(5:1,function(r) max(which(regimeU==r)))],lty=3)
abline(h=1-alpha,lty=2)



# ------------------------------------------------------------------------------------------------------------ 
# Coverage

# not the most elegant code ever in terms of environment handling ,but it works
codechunk <- function(cols,w,dist,new){

xgrid <- seq(-lambda-15,lambda+15,0.005)
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

theta0seq <- c(seq(lambda+0.01,lambda+10,0.005))
xstar.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Ugrid>theta0))]})
xstar.sup <- sapply(theta0seq,function(theta0){xgrid[max(which(Ugrid<theta0))]  })
xtilde.sup <-sapply(theta0seq,function(theta0){xgrid[max(which(Lgrid<theta0))]})
xtilde.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Lgrid>theta0))]})


if(dist=="Lap"){G <- function(x,theta0=0){plaplace(x,m=theta0,s=1)} }
if(dist=="Normal"){G <- function(x,theta0=0){pnorm(x,theta0,1)} }
C.inf <- G(xtilde.inf-theta0seq)-G(xstar.sup-theta0seq)
C.sup <- G(xtilde.sup-theta0seq)-G(xstar.inf-theta0seq)

if(new==T){
  plot(theta0seq,C.sup,xlim=range(theta0seq),type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])),ylim=c(1-2*alpha,1))
  abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),lty=rep(3,4),col=rep("grey",4))
  abline(h=1-w+w*(w-alpha),lty=3,col="grey")
  text(x=rep(lambda+9,4),y=0.01+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=0.8,adj=0)
}

polygon(x=c(theta0seq,sort(theta0seq,decreasing=T)),y=c(C.inf,C.sup[order(theta0seq,decreasing=T)]),
        col="grey90",border="white",density=20,angle=45)

polygon(x=c(theta0seq,sort(theta0seq,decreasing=T)),y=c(C.inf,C.sup[order(theta0seq,decreasing=T)]),
        col="grey90",border="white",density=20,angle=-45)

for(r in 5:1){
  
  reg <- (sapply(xstar.sup,function(x) regimeL[which(xgrid==x)])==r) # inf
  if(sum(reg)>0){lines(theta0seq[reg],C.inf[reg],col=cols,lty=2)}
  
  reg <- (sapply(xstar.inf,function(x) regimeL[which(xgrid==x)])==r)# Sup
  if(sum(reg)>0){lines(theta0seq[reg],C.sup[reg],col=cols)}

  }

} # end codechunk


pdf("figcov.pdf",width=9)
#pdf("SFcov.pdf",width=9)

par(mfrow=c(2,2))
alpha <- 0.05

## top Normal
lambda <- 0.75 # left
codechunk(cols="black",w=1,dist="Normal",new=T)
#codechunk(cols="black",w=0.2,dist="Normal",new=T)

lambda <- 7.5 # right
codechunk(cols="black",w=1,dist="Normal",new=T)
#codechunk(cols="black",w=0.2,dist="Normal",new=T)


## Bottom Laplace
lambda <- 0.75 # left
codechunk(cols="blue",w=1,dist="Lap",new=T)
#codechunk(cols="blue",w=0.2,dist="Lap",new=T)

lambda <- 7.5 # right
codechunk(cols="blue",w=1,dist="Lap",new=T)
#codechunk(cols="blue",w=0.2,dist="Lap",new=T)

dev.off()



# ------------------------------------------------------------------------------------------------------------ 
# Fig: Interval length
cardinal <- function(xgrid,lambda,w,dist){
Ugrid <- Lgrid <- regimeU <- regimeL <- UL.length <- numeric(length(xgrid))

for(i in 1:length(xgrid)){
  x <- xgrid[i]
  out.U <-getU(x,alpha,lambda,w,dist)
  Ugrid[i] <- out.U$val
  regimeU[i] <- out.U$regime
  
  out.L <- getU(-x,alpha,lambda,w,dist)
  Lgrid[i] <- -(out.L$val)
  regimeL[i] <- out.L$regime
  
  
  if(is.na(Ugrid[i]) | is.na(Lgrid[i]) ){UL.length[i]=0}else{
    UL.length[i] <- Ugrid[i] - Lgrid[i]
    if(Ugrid[i] >= lambda & Lgrid[i] <= (-lambda)){UL.length[i] <- UL.length[i] - 2*lambda}
  }
  
}
return(UL.length)
}


# Figure length
xgrid <- seq(-lambda-15,lambda+15,0.005)

pdf("figdiam.pdf",width=9)
#pdf("SFdiam.pdf",width=9)

par(mfrow=c(1,2))
ylims <- c(0,6)

# w=1, lambda = 0.75
plot(xgrid,cardinal(xgrid,lambda=0.75,w=1,dist="Normal"),col="black",ylim=ylims,xlab="x",ylab=expression(paste("diam(",U[alpha](x),",",L[alpha](x),")")),type="l")
lines(xgrid,cardinal(xgrid,lambda=0.75,w=1,dist="Lap"),col="blue",lty=2)

## w=1, lmabda=7.5
plot(xgrid,cardinal(xgrid,lambda=7.5,w=1,dist="Normal"),col="black",ylim=ylims,xlab="x",ylab=expression(paste("diam(",U[alpha](x),",",L[alpha](x),")")),type="l")
lines(xgrid,cardinal(xgrid,lambda=7.5,w=1,dist="Lap"),col="blue",lty=2)

## w=0.2, lambda=7.5
#plot(xgrid,cardinal(xgrid,lambda=0.75,w=0.2,dist="Normal"),col="black",ylim=ylims,xlab="x",ylab=expression(paste("diam(",U[alpha](x),",",L[alpha](x),")")),type="l")
#lines(xgrid,cardinal(xgrid,lambda=0.75,w=0.2,dist="Lap"),col="blue",lty=2)

#plot(xgrid,cardinal(xgrid,lambda=7.5,w=0.2,dist="Normal"),col="black",ylim=ylims,xlab="x",ylab=expression(paste("diam(",U[alpha](x),",",L[alpha](x),")")),type="l")
#lines(xgrid,cardinal(xgrid,lambda=7.5,w=0.2,dist="Lap"),col="blue",lty=2)



dev.off()

# ------------------------------------------------------------------------------------------------------------ 
# "efficiency"
# coverage / length
# not the most elegant code ever in terms of environment handling ,but it works
effchunk <- function(cols,w,dist,new){
  
  xgrid <- seq(-lambda-15,lambda+15,0.005)
  Ugrid <- Lgrid <- regimeU <- regimeL  <- UL.length <- numeric(length(xgrid))
  
  for(i in 1:length(xgrid)){
    x <- xgrid[i]
    out.U <-getU(x,alpha,lambda,w,dist)
    Ugrid[i] <- out.U$val
    regimeU[i] <- out.U$regime
    
    out.L <- getU(-x,alpha,lambda,w,dist)
    Lgrid[i] <- -(out.L$val)  
    regimeL[i] <- out.L$regime
  }
  
  # Cov
  theta0seq <- c(seq(lambda+0.01,lambda+10,0.005))
  xstar.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Ugrid>theta0))]})
  xstar.sup <- sapply(theta0seq,function(theta0){xgrid[max(which(Ugrid<theta0))]  })
  xtilde.sup <-sapply(theta0seq,function(theta0){xgrid[max(which(Lgrid<theta0))]})
  xtilde.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Lgrid>theta0))]})
  
  
  if(dist=="Lap"){G <- function(x,theta0=0){plaplace(x,m=theta0,s=1)} }
  if(dist=="Normal"){G <- function(x,theta0=0){pnorm(x,theta0,1)} }
  C.inf <- G(xtilde.inf-theta0seq)-G(xstar.sup-theta0seq)
  C.sup <- G(xtilde.sup-theta0seq)-G(xstar.inf-theta0seq)
  
  # Diam
  for(i in 1:length(xgrid)){
  if(is.na(Ugrid[i]) | is.na(Lgrid[i]) ){UL.length[i]=0}else{
    UL.length[i] <- Ugrid[i] - Lgrid[i]
    if(Ugrid[i] >= lambda & Lgrid[i] <= (-lambda)){UL.length[i] <- UL.length[i] - 2*lambda}
  }
  }
  
  # Expected diam (function of theta0) to match with C(theta0)
  Ediam <- sapply(theta0seq,function(t){
    
    if(dist=="Lap"){evec <- rlaplace(1000)}
    if(dist=="Normal"){evec <- rnorm(1000)}
    
    return(mean(sapply(evec,function(e){
      x <- t + e
      return(UL.length[which.min(abs(x - xgrid))])
    }))
    )
  })
  
  if(new==T){
    plot(theta0seq,C.sup/Ediam,xlim=range(theta0seq),type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])/E[theta[0]](diam)),ylim=c(0,0.5))
  }
  
  polygon(x=c(theta0seq,sort(theta0seq,decreasing=T)),y=c(C.inf/Ediam,(C.sup/Ediam)[order(theta0seq,decreasing=T)]),
          col="grey90",border="white",density=20,angle=45)
  
  polygon(x=c(theta0seq,sort(theta0seq,decreasing=T)),y=c(C.inf/Ediam,(C.sup/Ediam)[order(theta0seq,decreasing=T)]),
          col="grey90",border="white",density=20,angle=-45)
  
  for(r in 5:1){
    
    reg <- (sapply(xstar.sup,function(x) regimeL[which(xgrid==x)])==r) # inf
    if(sum(reg)>0){lines(theta0seq[reg],C.inf[reg]/Ediam[reg],col=cols,lty=2)}
    
    reg <- (sapply(xstar.inf,function(x) regimeL[which(xgrid==x)])==r)# Sup
    if(sum(reg)>0){lines(theta0seq[reg],C.sup[reg]/Ediam[reg],col=cols)}
    
  }
  
} # end effchunk



alpha <- 0.05

#pdf("figeff.pdf",width=9)
pdf("SFeff.pdf",width=9)

par(mfrow=c(2,2))
## top Normal
lambda <- 0.75 # left
#effchunk(cols="black",w=1,dist="Normal",new=T)
effchunk(cols="black",w=0.2,dist="Normal",new=T)

lambda <- 7.5 # right
#effchunk(cols="black",w=1,dist="Normal",new=T)
effchunk(cols="black",w=0.2,dist="Normal",new=T)


## Bottom Laplace
lambda <- 0.75 # left
#effchunk(cols="blue",w=1,dist="Lap",new=T)
effchunk(cols="blue",w=0.2,dist="Lap",new=T)

lambda <- 7.5 # right
#effchunk(cols="blue",w=1,dist="Lap",new=T)
effchunk(cols="blue",w=0.2,dist="Lap",new=T)

dev.off()
