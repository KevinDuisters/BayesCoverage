#--------------------------------------------------------------------------------------------------------------------#
# Introduction figure
# coverage paper
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Sep 2018
#--------------------------------------------------------------------------------------------------------------------#
# Source functions
source("post.R")
source("getU.R")
#--------------------------------------------------------------------------------------------------------------------#
# Settings
lambda <- 0.75
w <- 0.2
alpha <- 0.05
xl<-7
ranges <- c(-xl,xl)

thetaseq <- c(seq(-xl,-lambda-0.001,0.001),seq(lambda+0.001,xl,0.001))
sub1 <- thetaseq<(-lambda)
sub2 <- thetaseq>=(-lambda) & thetaseq <= lambda
sub3 <- thetaseq>lambda

# prep
dist <- "Normal"
G <- function(x){pnorm(x)}
#xgrid <- seq(-lambda-15,lambda+15,0.01)
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
 



# plot

pdf("Figures/figIntroA.pdf",width=9)
par(mfrow=c(2,2),xpd=F,mar=c(4,4,2,2))
#layout(matrix(c(1,2,3,4,5,6),2,3,byrow=T),widths = c(1,2,2),heights=rep(1,3) )

#plot(-10:10,-10:10,bty="n",xaxt="n",yaxt="n",type="n",xlab="",ylab="")
#text(0,0,"Prior")


x <- 1.25 

# a top
plot(0,0,type="n",ylim=c(0,0.7),xlim=c(-xl,xl),xaxt="s",yaxt="n",xlab=expression(theta),ylab="Density")
#segments(x0=-xl,x1=xl,y0=0)
#mtext(side=1,at=c(-lambda,0,lambda),c(expression(-lambda),"0",expression(lambda)))

abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=-2,side=1,at=c(-lambda,lambda))

arrows(-xl,1/(6-2*lambda),-lambda,1/(6-2*lambda),lwd=1.5,angle=15,length=0.10,code=1)
arrows(lambda,1/(6-2*lambda),xl,1/(6-2*lambda),lwd=1.5,angle=15,length=0.10,code=2)
segments(-lambda,0,lambda,0,lwd=1.5)
points(c(-lambda,lambda),rep(0,2),pch=21,bg=1)
points(c(-lambda,lambda),rep(1/(6-2*lambda),2),pch=21,bg="white")
text(x=3,y=0.17,expression(paste(pi,"(",theta,")")),adj=0)


# b top (# the 0.65 for the spike is just a choice for visual illustration, has nothing to do with 1-w)
plot(0,0,type="n",ylim=c(0,0.7),xlim=c(-xl,xl),xaxt="s",yaxt="n",xlab=expression(theta),ylab="Density*")
#segments(x0=-xl,x1=xl,y0=0)
#mtext(side=1,at=c(-lambda,0,lambda),c(expression(-lambda),"0",expression(lambda)))

abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=-2,side=1,at=c(-lambda,lambda))

arrows(-xl,w*1/(6-2*lambda),-lambda,w*1/(6-2*lambda),lwd=1.5,angle=15,length=0.10,code=1)
arrows(lambda,w*1/(6-2*lambda),xl,w*1/(6-2*lambda),lwd=1.5,angle=15,length=0.10,code=2)
segments(-lambda,0,lambda,0,lwd=1.5)
points(c(-lambda,lambda),rep(0,2),pch=21,bg=1)
points(c(-lambda,lambda),rep(w*1/(6-2*lambda),2),pch=21,bg="white")
arrows(x0=0,x1=0,y0=0,y1=0.65,lwd=3,angle=15,length=0.10)
text(x=0.5,y=0.65,expression(paste(P,"(",theta,"=0) = 1-w")),adj=0)
text(x=3,y=0.105,expression(paste(pi,"(",theta,")")),adj=0)
#segments(x0=-8.75,x1=xl,y0=1/(6-2*lambda),lty=3,col="grey")


#plot(-10:10,-10:10,bty="n",xaxt="n",yaxt="n",type="n",xlab="",ylab="")
#text(0,0,"Posterior")

# a bottom
plot(0,0,type="n",ylim=c(0,0.7),xlim=c(-xl,xl),xaxt="s",yaxt="n",xlab=expression(theta),ylab="Density")
#mtext(side=1,at=c(-lambda,0,lambda),c(expression(-lambda),"0",expression(lambda)))

abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=-2,side=1,at=c(-lambda,lambda))


polygon(c(-lambda,seq(-lambda,lambda,0.01),lambda,-lambda),c(0,dnorm(seq(-lambda,lambda,0.01)-x),0,0),col="grey",density=20,border=NA)
segments(x0=c(-lambda,lambda),x1=c(-lambda,lambda),y0=c(0,0),y1=c(dnorm(-lambda-x),dnorm(lambda-x)),col="grey")
lines(seq(-xl,xl,0.01),dnorm(x-seq(-xl,xl,0.01)),lty=2)
#segments(x0=-xl,x1=xl,y0=0)
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

# b bottom
plot(0,0,type="n",ylim=c(0,0.7),xlim=c(-xl,xl),xaxt="s",yaxt="n",xlab=expression(theta),ylab="Density*")
#mtext(side=1,at=c(-lambda,0,lambda),c(expression(-lambda),"0",expression(lambda)))

abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=-2,side=1,at=c(-lambda,lambda))

polygon(c(-lambda,seq(-lambda,lambda,0.01),lambda,-lambda),c(0,dnorm(seq(-lambda,lambda,0.01)-x),0,0),col="grey",density=20,border=NA)
segments(x0=c(-lambda,lambda),x1=c(-lambda,lambda),y0=c(0,0),y1=c(dnorm(-lambda-x),dnorm(lambda-x)),col="grey")
lines(seq(-xl,xl,0.01),dnorm(x-seq(-xl,xl,0.01)),lty=2)
#segments(x0=-xl,x1=xl,y0=0)

segments(-lambda,x1=lambda,y0=0,lwd=1.5,col="blue")
lines(thetaseq[sub1],sapply(thetaseq[sub1],function(t)post(t,x=x,lambda,w,"Normal")),col="blue",lwd=1.5)
lines(thetaseq[sub2],sapply(thetaseq[sub2],function(t)post(t,x=x,lambda,w,"Normal")),col="blue",lwd=1.5)
lines(thetaseq[sub3],sapply(thetaseq[sub3],function(t)post(t,x=x,lambda,w,"Normal")),col="blue",lwd=1.5)
points(c(-lambda,lambda),rep(0,2),pch=21,bg="blue",col="blue")
points(c(-lambda,lambda),c(post(-lambda-1e-6,x,lambda,w,"Normal"),post(lambda+1e-6,x,lambda,w,"Normal")),pch=21,bg="white",col="blue")

#segments(x0=-19+xl+x+0.25,x1=xl,y0=post(x,x,lambda,w=1,"Normal"),lty=3,col="grey")


arrows(x0=0,x1=0,y0=0,y1=0.65,lwd=3,angle=15,length=0.10,col="blue")
text(x=0.5,y=0.65,expression(paste(P,"(",theta,"=0|X=x)")),adj=0)

text(x=1.9,y=0.4,expression(paste("g(x-",theta,")")),adj=0)
text(x=3.2,y=0.12,expression(paste(pi,"(",theta,"|X=x)")),adj=0)
segments(x0=2.65, x1=3, y0=0.1, y1=0.12,col="blue")


polygon(c(-4,-3.5,-3.5,-4,-4),c(0.275,0.275,0.325,0.325,0.275),col="grey",density=20)
text(x=-3.25,y=0.3,expression(Delta[lambda](x)),adj=0)

dev.off()

pdf("figures/figIntroB.pdf",width=9)
par(mfcol=c(2,2),xpd=F,mar=c(4,4,2,2))
#layout(matrix(c(1,2,3,4,5,6),2,3,byrow=F),widths = c(1,2,2),heights=rep(1,3) )


#plot(-10:10,-10:10,bty="n",xaxt="n",yaxt="n",type="n",xlab="",ylab="")
#text(0,0,"Inverse HPD")

#plot(-10:10,-10:10,bty="n",xaxt="n",yaxt="n",type="n",xlab="",ylab="")
#text(0,0,"Freq. Cov.")


# (3,1)
xstar.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Ugrid[,1]>theta0))]})
xstar.sup <- sapply(thetaseq,function(theta0){xgrid[max(which(Ugrid[,1]<theta0))]})
xtilde.sup <-sapply(thetaseq,function(theta0){xgrid[max(which(Lgrid[,1]<theta0))]})
xtilde.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Lgrid[,1]>theta0))]})

#plot(thetaseq,xstar.sup,type="n",xlim=ranges,ylim=ranges,ylab=expression(paste(x,"*",(theta[0])," and ",tilde(x)(theta[0]))),col="black",xlab=expression(theta[0]),asp=1)
plot(thetaseq,xstar.sup,type="n",xlim=ranges,ylim=ranges,ylab=expression(paste(x,"*",(theta[0])," and ",tilde(x)(theta[0]))),col="black",xlab=expression(theta[0]))
#polygon(x=c(-lambda,lambda,lambda,-lambda,-lambda),y=c(1.04*min(ranges),1.04*min(ranges),1.04*max(ranges),1.04*max(ranges),1.04*min(ranges)),col="grey90",border=NA)
abline(h=xgrid[sapply(5:2,function(r) max(which(regimeU[,1]==r)))],lty=3)
mtext(c("V","IV","III","II","I"),side=4, adj=0.5, line=-2,at=c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU[,1]==r)))][-5])+diff(c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU[,1]==r)))][-5],ranges[2]))/2   )

polygon(x=c(thetaseq[thetaseq<(-lambda)],sort(thetaseq[thetaseq<(-lambda)],decreasing=T)),y=c(xtilde.sup[thetaseq<(-lambda)],(xtilde.inf[thetaseq<(-lambda)])[order(thetaseq[thetaseq<(-lambda)],decreasing=T)]),col="blue",border=NA,density=20)
polygon(x=c(thetaseq[thetaseq>(lambda)],sort(thetaseq[thetaseq>(lambda)],decreasing=T)),y=c(xstar.sup[thetaseq>(lambda)],(xstar.inf[thetaseq>(lambda)])[order(thetaseq[thetaseq>(lambda)],decreasing=T)]),col="black",border=NA,density=20)

# technically the following are not lines (only piecewise continuous), but this is corrected by polygon
lines(thetaseq[thetaseq<(-lambda)],xstar.sup[thetaseq<(-lambda)])
lines(thetaseq[thetaseq>(lambda)],xstar.sup[thetaseq>(lambda)])
lines(thetaseq[thetaseq<(-lambda)],xstar.inf[thetaseq<(-lambda)],lty=2)
lines(thetaseq[thetaseq>(lambda)],xstar.inf[thetaseq>(lambda)],lty=2)
lines(thetaseq[thetaseq<(-lambda)],xtilde.sup[thetaseq<(-lambda)],col="blue")
lines(thetaseq[thetaseq>(lambda)],xtilde.sup[thetaseq>(lambda)],col="blue")
lines(thetaseq[thetaseq<(-lambda)],xtilde.inf[thetaseq<(-lambda)],lty=2,col="blue")
lines(thetaseq[thetaseq>(lambda)],xtilde.inf[thetaseq>(lambda)],lty=2,col="blue")
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=-2,side=1,at=c(-lambda,lambda))


legend("topleft",c(expression(paste("sup ",x,"*")),expression(paste("inf ",x,"*")),expression(paste("sup ",tilde(x))),expression(paste("inf ",tilde(x)))),
       lty=c(1,2,1,2),col=c("black","black","blue","blue"),ncol=2,bty="n",
       seg.len=2)

# (4,1)
C.inf <- G(xtilde.inf-thetaseq)-G(xstar.sup-thetaseq)
C.sup <- G(xtilde.sup-thetaseq)-G(xstar.inf-thetaseq)

plot(thetaseq,C.sup,xlim=range(thetaseq),type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])),ylim=c(1-2*alpha,1))
abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),lty=rep(3,4),col=rep("grey",4))
text(x=rep(lambda+4.75,4),y=0.005+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=0.8,adj=0)


polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),
        col="grey90",border="white",density=20,angle=45)
polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),
        col="grey90",border="white",density=20,angle=-45)



# translate to theta regimes
starsupgrid <- sapply(xstar.sup,function(x) regimeU[which(xgrid==x),1])
starinfgrid <- sapply(xstar.inf,function(x) regimeU[which(xgrid==x),1])
tildesupgrid <- sapply(xtilde.sup,function(x) regimeU[which(xgrid==x),1])
tildeinfgrid <- sapply(xtilde.inf,function(x) regimeU[which(xgrid==x),1])

for(r in 5:1){
  
  # sup cov  (sup xtilde, inf xstar); for theta0>lambda regime change in xstar, for theta0 < - lambda regime change in xtilde
  regmin <- (thetaseq < (- lambda)) & (tildesupgrid==r)
  regplus <- (thetaseq > (lambda))  & (starinfgrid==r)
  lines(thetaseq[regmin], C.sup[regmin],col="black")
  lines(thetaseq[regplus], C.sup[regplus],col="black")
  
  # ablines (sup)
  #if(r<5){abline(v=min(thetaseq[regmin | regplus]),lty=3,col="red")}
  
  # inf cov  (inf xtilde, sup xstar); for theta0>lambda regime change in xstar, for theta0 < - lambda regime change in xtilde
  regmin <- (thetaseq < (- lambda)) & (tildeinfgrid==r)
  regplus <- (thetaseq > (lambda))  & (starsupgrid==r)
  lines(thetaseq[regmin], C.inf[regmin],lty=2,col="blue")
  lines(thetaseq[regplus], C.inf[regplus],lty=2,col="blue")
  
  # ablines (inf)
  if(r==5){mini <- min(thetaseq)}
  maxi <- max(thetaseq[regmin | regplus])
  if(r>1){abline(v=maxi,lty=3,col="grey")}
  if(r==1){maxi <- max(thetaseq)}
  if(r!=3){mtext(as.roman(r),line=-2,side=3,adj=0.5,at=(maxi+mini)/2)}else{mtext(c("III","III"),line=-2,side=3,adj=0.5,at=c((maxi+lambda)/2,(-lambda-maxi)/2))}
  mini <- maxi        
}




#polygon(x=c(-lambda,lambda,lambda,-lambda,-lambda),y=c(1.04*min(ranges),1.04*min(ranges),1.04*max(ranges),1.04*max(ranges),1.04*min(ranges)),col="grey90",border=NA)
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=-2,side=1,at=c(-lambda,lambda))


# (3,2)
xstar.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Ugrid[,2]>theta0))]})
xstar.sup <- sapply(thetaseq,function(theta0){xgrid[max(which(Ugrid[,2]<theta0))]})
xtilde.sup <-sapply(thetaseq,function(theta0){xgrid[max(which(Lgrid[,2]<theta0))]})
xtilde.inf <- sapply(thetaseq,function(theta0){xgrid[min(which(Lgrid[,2]>theta0))]})

#plot(thetaseq,xstar.sup,type="n",xlim=ranges,ylim=ranges,ylab=expression(paste(x,"*",(theta[0])," and ",tilde(x)(theta[0]))),col="black",xlab=expression(theta[0]),asp=1)
plot(thetaseq,xstar.sup,type="n",xlim=ranges,ylim=ranges,ylab=expression(paste(x,"*",(theta[0])," and ",tilde(x)(theta[0]))),col="black",xlab=expression(theta[0]))

#polygon(x=c(-lambda,lambda,lambda,-lambda,-lambda),y=c(1.04*min(ranges),1.04*min(ranges),1.04*max(ranges),1.04*max(ranges),1.04*min(ranges)),col="grey90",border=NA)
abline(h=xgrid[sapply(5:2,function(r) max(which(regimeU[,2]==r)))],lty=3)
mtext(c("V","IV","III","II","I"),side=4, adj=0.5, line=-2,at=c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU[,2]==r)))][-5])+diff(c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU[,2]==r)))][-5],ranges[2]))/2   )

polygon(x=c(thetaseq[thetaseq<(-lambda)],sort(thetaseq[thetaseq<(-lambda)],decreasing=T)),y=c(xtilde.sup[thetaseq<(-lambda)],(xtilde.inf[thetaseq<(-lambda)])[order(thetaseq[thetaseq<(-lambda)],decreasing=T)]),col="blue",border=NA,density=20)
polygon(x=c(thetaseq[thetaseq>(lambda)],sort(thetaseq[thetaseq>(lambda)],decreasing=T)),y=c(xstar.sup[thetaseq>(lambda)],(xstar.inf[thetaseq>(lambda)])[order(thetaseq[thetaseq>(lambda)],decreasing=T)]),col="black",border=NA,density=20)

lines(thetaseq[thetaseq<(-lambda)],xstar.sup[thetaseq<(-lambda)])
lines(thetaseq[thetaseq>(lambda)],xstar.sup[thetaseq>(lambda)])
lines(thetaseq[thetaseq<(-lambda)],xstar.inf[thetaseq<(-lambda)],lty=2)
lines(thetaseq[thetaseq>(lambda)],xstar.inf[thetaseq>(lambda)],lty=2)
lines(thetaseq[thetaseq<(-lambda)],xtilde.sup[thetaseq<(-lambda)],col="blue")
lines(thetaseq[thetaseq>(lambda)],xtilde.sup[thetaseq>(lambda)],col="blue")
lines(thetaseq[thetaseq<(-lambda)],xtilde.inf[thetaseq<(-lambda)],lty=2,col="blue")
lines(thetaseq[thetaseq>(lambda)],xtilde.inf[thetaseq>(lambda)],lty=2,col="blue")
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=-2,side=1,at=c(-lambda,lambda))
#abline(h=c(-x1,-x2,x2,x1),lty=3,col="grey")
#mtext(c(expression(-x[1]),expression(-x[2]),expression(x[2]),expression(x[1])),line=-2,side=2,at=c(-x1,-x2,x2,x1))

legend("topleft",c(expression(paste("sup ",x,"*")),expression(paste("inf ",x,"*")),expression(paste("sup ",tilde(x))),expression(paste("inf ",tilde(x)))),
       lty=c(1,2,1,2),col=c("black","black","blue","blue"),ncol=2,bty="n",
       seg.len=2)

# (4,2)
C.inf <- G(xtilde.inf-thetaseq)-G(xstar.sup-thetaseq)
C.sup <- G(xtilde.sup-thetaseq)-G(xstar.inf-thetaseq)

plot(thetaseq,C.sup,xlim=range(thetaseq),type="n",xlab=expression(theta[0]),ylab=expression(C(theta[0])),ylim=c(1-2*alpha,1))
abline(h=c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),lty=rep(3,4),col=rep("grey",4))
text(x=rep(lambda+4.7,4),y=0.005+c(1-alpha/2,1-alpha,1-3*alpha/2,1-2*alpha),labels=c(expression(1-alpha/2),expression(1-alpha),expression(1-3*alpha/2),expression(1-2*alpha)),cex=0.8,adj=0)

polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),
        col="grey90",border="white",density=20,angle=45)
polygon(x=c(thetaseq,sort(thetaseq,decreasing=T)),y=c(C.inf,C.sup[order(thetaseq,decreasing=T)]),
        col="grey90",border="white",density=20,angle=-45)

# translate to theta regimes
starsupgrid <- sapply(xstar.sup,function(x) regimeU[which(xgrid==x),2])
starinfgrid <- sapply(xstar.inf,function(x) regimeU[which(xgrid==x),2])
tildesupgrid <- sapply(xtilde.sup,function(x) regimeU[which(xgrid==x),2])
tildeinfgrid <- sapply(xtilde.inf,function(x) regimeU[which(xgrid==x),2])

for(r in 5:1){
  
  # sup cov  (sup xtilde, inf xstar); for theta0>lambda regime change in xstar, for theta0 < - lambda regime change in xtilde
  regmin <- (thetaseq < (- lambda)) & (tildesupgrid==r)
  regplus <- (thetaseq > (lambda))  & (starinfgrid==r)
  lines(thetaseq[regmin], C.sup[regmin],col="black")
  lines(thetaseq[regplus], C.sup[regplus],col="black")
  
  # inf cov  (inf xtilde, sup xstar); for theta0>lambda regime change in xstar, for theta0 < - lambda regime change in xtilde
  regmin <- (thetaseq < (- lambda)) & (tildeinfgrid==r)
  regplus <- (thetaseq > (lambda))  & (starsupgrid==r)
  lines(thetaseq[regmin], C.inf[regmin],lty=2,col="blue")
  lines(thetaseq[regplus], C.inf[regplus],lty=2,col="blue")
  
  # ablines (inf)
  if(r==5){mini <- min(thetaseq)}
  maxi <- max(thetaseq[regmin | regplus])
  if(r>1){abline(v=maxi,lty=3,col="grey")}
  if(r==1){maxi <- max(thetaseq)}
  if(r!=3){mtext(as.roman(r),line=-2,side=3,adj=0.5,at=(maxi+mini)/2)}else{mtext(c("III","III"),line=-2,side=3,adj=0.5,at=c((maxi+lambda)/2,(-lambda-maxi)/2))}
  mini <- maxi    
  
}


#polygon(x=c(-lambda,lambda,lambda,-lambda,-lambda),y=c(1.04*min(ranges),1.04*min(ranges),1.04*max(ranges),1.04*max(ranges),1.04*min(ranges)),col="grey90",border=NA)
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=-2,side=1,at=c(-lambda,lambda))

dev.off()