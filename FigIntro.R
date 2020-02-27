#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Feb 2020
#--------------------------------------------------------------------------------------------------------------------#
# Source functions
source("functions/post.R")
source("functions/getU.R")
source("functions/coverage.R")
#--------------------------------------------------------------------------------------------------------------------#
# Settings
lambda <- 0.5
w <- 0.2
dist <- "Normal"
distname <- "N(0,1)"
alpha <- 0.05
xl<-7
ranges <- c(-xl,xl)


x <- 1.25 


thetaseq <- c(seq(-xl,-lambda-0.005,0.005),seq(lambda+0.005,xl,0.005))
sub1 <- thetaseq<(-lambda)
sub2 <- thetaseq>=(-lambda) & thetaseq <= lambda
sub3 <- thetaseq>lambda

#--------------------------------------------------------------------------------------------------------------------#
# plot
#pdf("Figures/figIntro.pdf",width=12,height=9)
png("Figures/figIntro.png",width=12,height=9,units="in",res=300)
par(mfcol=c(2,3))


for(w in c(1,0.25)){
if(w==1){ylw <- "Density"}else{ylw <- "Density*"}
plot(0,0,type="n",ylim=c(0,0.7),xlim=c(-xl,xl),xaxt="s",yaxt="n",xlab=expression(theta),ylab=ylw)
title(bquote(paste(.(distname),", ",lambda==.(lambda),", ",w==.(w))))
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=1,side=1,at=c(-lambda,lambda))
arrows(-xl,w*1/(6-2*lambda),-lambda,w*1/(6-2*lambda),lwd=1,angle=15,length=0.10,code=1)
arrows(lambda,w*1/(6-2*lambda),xl,w*1/(6-2*lambda),lwd=1,angle=15,length=0.10,code=2)
segments(-lambda,0,lambda,0,lwd=1.5)
#points(c(-lambda,lambda),rep(0,2),pch=21,bg=1)
#points(c(-lambda,lambda),rep(w*1/(6-2*lambda),2),pch=21,bg="white")

if(w==0.25){
arrows(x0=0,x1=0,y0=0,y1=0.65,lwd=3,angle=15,length=0.10) # this 0.65 is arbitrary
text(x=0.5,y=0.65,expression(paste(P,"(",theta,"=0) = 1-w")),adj=0)
text(x=3,y=0.08,expression(paste(pi,"(",theta,")")),adj=0)
}
if(w==1){text(x=3,y=0.17,expression(paste(pi,"(",theta,")")),adj=0)}
}


# middle
for(w in c(1,0.25)){
if(w==1){ylw <- "Density"}else{ylw <- "Density*"}
plot(0,0,type="n",ylim=c(0,0.7),xlim=c(-xl,xl),xaxt="s",yaxt="n",xlab=expression(theta),ylab=ylw)
title(bquote(paste(.(distname),", ",lambda==.(lambda),", ",w==.(w))))
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=1,side=1,at=c(-lambda,lambda))
polygon(c(-lambda,seq(-lambda,lambda,0.01),lambda,-lambda),c(0,dnorm(seq(-lambda,lambda,0.01)-x),0,0),col="grey",density=20,border=NA)
segments(x0=c(-lambda,lambda),x1=c(-lambda,lambda),y0=c(0,0),y1=c(dnorm(-lambda-x),dnorm(lambda-x)),col="grey")
lines(seq(-xl,xl,0.01),dnorm(x-seq(-xl,xl,0.01)),lty=2)
segments(-lambda,x1=lambda,y0=0,lwd=1.5,col="black")
lines(thetaseq[sub1],sapply(thetaseq[sub1],function(t)post(t,x=x,lambda,w,dist)),col="black",lwd=1)
lines(thetaseq[sub2],sapply(thetaseq[sub2],function(t)post(t,x=x,lambda,w,dist)),col="black",lwd=1)
lines(thetaseq[sub3],sapply(thetaseq[sub3],function(t)post(t,x=x,lambda,w,dist)),col="black",lwd=1)
#points(c(-lambda,lambda),rep(0,2),pch=21,bg="black",col="black")
#points(c(-lambda,lambda),c(post(-lambda-1e-6,x,lambda,w,dist),post(lambda+1e-6,x,lambda,w,dist)),pch=21,bg="white",col="black")

if(w==1){
text(x=2.45,y=0.4,expression(paste(pi,"(",theta,"|X=x)")),adj=0,col="black")
text(x=1.3,y=0.12,expression(paste("g(x-",theta,")")),adj=0)
}
if(w==0.25){
  arrows(x0=0,x1=0,y0=0,y1=0.65,lwd=3,angle=15,length=0.10,col="black")
  text(x=0.5,y=0.65,expression(paste(P,"(",theta,"=0|X=x)")),adj=0,col="black")
  text(x=1.9,y=0.4,expression(paste("g(x-",theta,")")),adj=0)
  text(x=1,y=0.07,expression(paste(pi,"(",theta,"|X=x)")),adj=0,col="black")
}
polygon(c(-4,-3.5,-3.5,-4,-4),c(0.275,0.275,0.325,0.325,0.275),col="grey",density=20)
text(x=-3.25,y=0.3,expression(Delta[lambda](x)),adj=0)

}

# bottom
for(w in c(1,0.25)){
coverage(alpha,lambda,w,dist="Normal",thetamax=xl,h=0.001,plot.cov=T,plot.cov.title=T,plot.neg=T,line.col="black")
if(w<1){points(0,1,pch=16) # theta_0 = 0 always covered since 0 always in HPD for w<1}
}
}

dev.off()



