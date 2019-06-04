#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# June 2019
#--------------------------------------------------------------------------------------------------------------------#
# Source functions
source("functions/post.R")
source("functions/getU.R")
#--------------------------------------------------------------------------------------------------------------------#
# Settings
lambda <- 0.75
w <- 0.2
alpha <- 0.05
xl<-7
ranges <- c(-xl,xl)

thetaseq <- c(seq(-xl,-lambda-0.005,0.005),seq(lambda+0.005,xl,0.005))
sub1 <- thetaseq<(-lambda)
sub2 <- thetaseq>=(-lambda) & thetaseq <= lambda
sub3 <- thetaseq>lambda

# prep
dist <- "Normal"
G <- function(x){pnorm(x)}
xgrid <- seq(-2*xl,2*xl,0.001)
Ugrid <- Lgrid <- regimeU <- matrix(NA,length(xgrid),2)

for(k in 1:2){
    if(k==1){w<-1}else{w<-0.2}
       for(i in 1:length(xgrid)){
       x <- xgrid[i]
       out.U <-getU(x,alpha,lambda,w,dist)
       Ugrid[i,k] <- out.U$val
       regimeU[i,k] <- out.U$regime
       out.L <- getU(-x,alpha,lambda,w,dist)
       Lgrid[i,k] <- -(out.L$val)
}
}
 
#--------------------------------------------------------------------------------------------------------------------#
# plot
pdf("Figures/figIntrotop.pdf",width=9)
par(mfrow=c(2,2),xpd=F,mar=c(4,4,2,2))

x <- 1.25 

# a top
plot(0,0,type="n",ylim=c(0,0.7),xlim=c(-xl,xl),xaxt="s",yaxt="n",xlab=expression(theta),ylab="Density")
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=1,side=1,at=c(-lambda,lambda))
arrows(-xl,1/(6-2*lambda),-lambda,1/(6-2*lambda),lwd=1.5,angle=15,length=0.10,code=1)
arrows(lambda,1/(6-2*lambda),xl,1/(6-2*lambda),lwd=1.5,angle=15,length=0.10,code=2)
segments(-lambda,0,lambda,0,lwd=1.5)
points(c(-lambda,lambda),rep(0,2),pch=21,bg=1)
points(c(-lambda,lambda),rep(1/(6-2*lambda),2),pch=21,bg="white")
text(x=3,y=0.17,expression(paste(pi,"(",theta,")")),adj=0)


# b top (# the 0.65 for the spike is just a choice for visual illustration, has nothing to do with 1-w)
plot(0,0,type="n",ylim=c(0,0.7),xlim=c(-xl,xl),xaxt="s",yaxt="n",xlab=expression(theta),ylab="Density*")
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=1,side=1,at=c(-lambda,lambda))
arrows(-xl,w*1/(6-2*lambda),-lambda,w*1/(6-2*lambda),lwd=1.5,angle=15,length=0.10,code=1)
arrows(lambda,w*1/(6-2*lambda),xl,w*1/(6-2*lambda),lwd=1.5,angle=15,length=0.10,code=2)
segments(-lambda,0,lambda,0,lwd=1.5)
points(c(-lambda,lambda),rep(0,2),pch=21,bg=1)
points(c(-lambda,lambda),rep(w*1/(6-2*lambda),2),pch=21,bg="white")
arrows(x0=0,x1=0,y0=0,y1=0.65,lwd=3,angle=15,length=0.10)
text(x=0.5,y=0.65,expression(paste(P,"(",theta,"=0) = 1-w")),adj=0)
text(x=3,y=0.105,expression(paste(pi,"(",theta,")")),adj=0)

# a middle
plot(0,0,type="n",ylim=c(0,0.7),xlim=c(-xl,xl),xaxt="s",yaxt="n",xlab=expression(theta),ylab="Density")
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=1,side=1,at=c(-lambda,lambda))
polygon(c(-lambda,seq(-lambda,lambda,0.01),lambda,-lambda),c(0,dnorm(seq(-lambda,lambda,0.01)-x),0,0),col="grey",density=20,border=NA)
segments(x0=c(-lambda,lambda),x1=c(-lambda,lambda),y0=c(0,0),y1=c(dnorm(-lambda-x),dnorm(lambda-x)),col="grey")
lines(seq(-xl,xl,0.01),dnorm(x-seq(-xl,xl,0.01)),lty=2)
segments(-lambda,x1=lambda,y0=0,lwd=1.5,col="blue")
lines(thetaseq[sub1],sapply(thetaseq[sub1],function(t)post(t,x=x,lambda,w=1,"Normal")),col="blue",lwd=1.5)
lines(thetaseq[sub2],sapply(thetaseq[sub2],function(t)post(t,x=x,lambda,w=1,"Normal")),col="blue",lwd=1.5)
lines(thetaseq[sub3],sapply(thetaseq[sub3],function(t)post(t,x=x,lambda,w=1,"Normal")),col="blue",lwd=1.5)
points(c(-lambda,lambda),rep(0,2),pch=21,bg="blue",col="blue")
points(c(-lambda,lambda),c(post(-lambda-1e-6,x,lambda,w=1,"Normal"),post(lambda+1e-6,x,lambda,w=1,"Normal")),pch=21,bg="white",col="blue")
text(x=2.45,y=0.4,expression(paste(pi,"(",theta,"|X=x)")),adj=0)
text(x=1.3,y=0.12,expression(paste("g(x-",theta,")")),adj=0)
polygon(c(-4,-3.5,-3.5,-4,-4),c(0.275,0.275,0.325,0.325,0.275),col="grey",density=20)
text(x=-3.25,y=0.3,expression(Delta[lambda](x)),adj=0)

# b middle
plot(0,0,type="n",ylim=c(0,0.7),xlim=c(-xl,xl),xaxt="s",yaxt="n",xlab=expression(theta),ylab="Density*")
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=1,side=1,at=c(-lambda,lambda))
polygon(c(-lambda,seq(-lambda,lambda,0.01),lambda,-lambda),c(0,dnorm(seq(-lambda,lambda,0.01)-x),0,0),col="grey",density=20,border=NA)
segments(x0=c(-lambda,lambda),x1=c(-lambda,lambda),y0=c(0,0),y1=c(dnorm(-lambda-x),dnorm(lambda-x)),col="grey")
lines(seq(-xl,xl,0.01),dnorm(x-seq(-xl,xl,0.01)),lty=2)
segments(-lambda,x1=lambda,y0=0,lwd=1.5,col="blue")
lines(thetaseq[sub1],sapply(thetaseq[sub1],function(t)post(t,x=x,lambda,w,"Normal")),col="blue",lwd=1.5)
lines(thetaseq[sub2],sapply(thetaseq[sub2],function(t)post(t,x=x,lambda,w,"Normal")),col="blue",lwd=1.5)
lines(thetaseq[sub3],sapply(thetaseq[sub3],function(t)post(t,x=x,lambda,w,"Normal")),col="blue",lwd=1.5)
points(c(-lambda,lambda),rep(0,2),pch=21,bg="blue",col="blue")
points(c(-lambda,lambda),c(post(-lambda-1e-6,x,lambda,w,"Normal"),post(lambda+1e-6,x,lambda,w,"Normal")),pch=21,bg="white",col="blue")
arrows(x0=0,x1=0,y0=0,y1=0.65,lwd=3,angle=15,length=0.10,col="blue")
text(x=0.5,y=0.65,expression(paste(P,"(",theta,"=0|X=x)")),adj=0)
text(x=1.9,y=0.4,expression(paste("g(x-",theta,")")),adj=0)
text(x=3.2,y=0.12,expression(paste(pi,"(",theta,"|X=x)")),adj=0)
segments(x0=2.65, x1=3, y0=0.1, y1=0.12,col="blue")
polygon(c(-4,-3.5,-3.5,-4,-4),c(0.275,0.275,0.325,0.325,0.275),col="grey",density=20)
text(x=-3.25,y=0.3,expression(Delta[lambda](x)),adj=0)


dev.off()
pdf("Figures/figIntrobottom.pdf",width=9)
par(mfrow=c(2,2),xpd=F,mar=c(4,4,2,2))

# a bottom
xU.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Ugrid[,1]>theta0))]})
xU.sup <- sapply(thetaseq,function(theta0){xgrid[max(which(Ugrid[,1]<theta0))]})
xL.sup <-sapply(thetaseq,function(theta0){xgrid[max(which(Lgrid[,1]<theta0))]})
xL.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Lgrid[,1]>theta0))]})

C.inf <- G(xL.inf-thetaseq)-G(xU.sup-thetaseq)
C.sup <- G(xL.sup-thetaseq)-G(xU.inf-thetaseq)

plot(thetaseq,C.sup,xlim=range(thetaseq),type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])),ylim=c(1-2*alpha,1))
abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),lty=rep(3,4),col=rep("grey",4))
text(x=rep(lambda+4.75,4),y=0.005+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=0.8,adj=0)


polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),col="grey90",border="white",density=20,angle=45)
polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),col="grey90",border="white",density=20,angle=-45)



# translate to theta regimes
starsupgrid <- sapply(xU.sup,function(x) regimeU[which(xgrid==x),1])
starinfgrid <- sapply(xU.inf,function(x) regimeU[which(xgrid==x),1])
tildesupgrid <- sapply(xL.sup,function(x) regimeU[which(xgrid==x),1])
tildeinfgrid <- sapply(xL.inf,function(x) regimeU[which(xgrid==x),1])

for(r in 5:1){
  
  # sup cov  (sup xL, inf xU); for theta0>lambda regime change in xU, for theta0 < - lambda regime change in xL
  regmin <- (thetaseq < (- lambda)) & (tildesupgrid==r)
  regplus <- (thetaseq > (lambda))  & (starinfgrid==r)
  lines(thetaseq[regmin], C.sup[regmin],col="black")
  lines(thetaseq[regplus], C.sup[regplus],col="black")
  
  
  # inf cov  (inf xL, sup xU); for theta0>lambda regime change in xU, for theta0 < - lambda regime change in xL
  #regmin <- (thetaseq < (- lambda)) & (tildeinfgrid==r)
  #regplus <- (thetaseq > (lambda))  & (starsupgrid==r)
  #lines(thetaseq[regmin], C.inf[regmin],lty=2,col="blue",lwd=1.5)
  #lines(thetaseq[regplus], C.inf[regplus],lty=2,col="blue",lwd=1.5)
}

abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=1,side=1,at=c(-lambda,lambda))

# some manual adjustments for visualization finetuning
points(c(-lambda,lambda),rep(C.sup[which.min((thetaseq-lambda)^2)],2),pch=16) # discontinuity points

# discontinuity between regimes
bool <- (((thetaseq > lambda) & (starinfgrid==3)))
points(c(-max(thetaseq[bool]),max(thetaseq[bool])),rep(C.sup[min(which(thetaseq==max(thetaseq[bool])))],2),pch=21,bg="white")
bool <- (((thetaseq > lambda) & (starinfgrid==2)))
points(c(-min(thetaseq[bool]),min(thetaseq[bool])),rep(C.sup[min(which(thetaseq==min(thetaseq[bool])))],2),pch=21,bg=1)



# b bottom
xU.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Ugrid[,2]>theta0))]})
xU.sup <- sapply(thetaseq,function(theta0){xgrid[max(which(Ugrid[,2]<theta0))]})
xL.sup <-sapply(thetaseq,function(theta0){xgrid[max(which(Lgrid[,2]<theta0))]})
xL.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Lgrid[,2]>theta0))]})

C.inf <- G(xL.inf-thetaseq)-G(xU.sup-thetaseq)
C.sup <- G(xL.sup-thetaseq)-G(xU.inf-thetaseq)

plot(thetaseq,C.sup,xlim=range(thetaseq),type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])),ylim=c(1-2*alpha,1))
abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),lty=rep(3,4),col=rep("grey",4))
text(x=rep(lambda+4.7,4),y=0.005+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=0.8,adj=0)

polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),
        col="grey90",border="white",density=20,angle=45)
polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),
        col="grey90",border="white",density=20,angle=-45)

# translate to theta regimes
starsupgrid <- sapply(xU.sup,function(x) regimeU[which(xgrid==x),2])
starinfgrid <- sapply(xU.inf,function(x) regimeU[which(xgrid==x),2])
tildesupgrid <- sapply(xL.sup,function(x) regimeU[which(xgrid==x),2])
tildeinfgrid <- sapply(xL.inf,function(x) regimeU[which(xgrid==x),2])

for(r in 5:1){
  
  # sup cov  (sup xL, inf xU); for theta0>lambda regime change in xU, for theta0 < - lambda regime change in xL
  regmin <- (thetaseq < (- lambda)) & (tildesupgrid==r)
  regplus <- (thetaseq > (lambda))  & (starinfgrid==r)
  lines(thetaseq[regmin], C.sup[regmin],col="black")
  lines(thetaseq[regplus], C.sup[regplus],col="black")
  
  # inf cov  (inf xL, sup xU); for theta0>lambda regime change in xU, for theta0 < - lambda regime change in xL
  #regmin <- (thetaseq < (- lambda)) & (tildeinfgrid==r)
  #regplus <- (thetaseq > (lambda))  & (starsupgrid==r)
  #lines(thetaseq[regmin], C.inf[regmin],lty=2,col="blue",lwd=1.5)
  #lines(thetaseq[regplus], C.inf[regplus],lty=2,col="blue",lwd=1.5)
  
}

abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=1,side=1,at=c(-lambda,lambda))


# some manual adjustments for visualization finetuning
points(0,1,pch=21,bg=1) # theta_0 = 0 always covered since 0 always in HPD for w<1
points(c(-lambda,lambda),rep(C.sup[which.min((thetaseq-lambda)^2)],2),pch=21,bg=1) # discontinuity points

dev.off()

