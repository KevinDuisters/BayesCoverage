#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Oct 2018
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#
# source functions
library(rmutil) #dlaplace in densities
source("functions/getU.R")
#--------------------------------------------------------------------------------------------------------------------#
# Plots
#--------------------------------------------------------------------------------------------------------------------#
# Fig 2: Visualize L,U and inverse xL, xU 
#dist <- "Lap"
dist <- "t3"
lambda<-7.5
alpha <- 0.05
w <- 1
#w<-0.9


#xgrid <- seq(-lambda-15,lambda+15,0.01)
xgrid <- seq(-lambda-25,lambda+25,0.01)
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

#theta0seq <- c(seq(-lambda-10,-lambda-0.01,0.01),seq(lambda+0.01,lambda+10,0.01))
theta0seq <- c(seq(-lambda-20,-lambda-0.01,0.01),seq(lambda+0.01,lambda+20,0.01))
xU.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Ugrid>theta0))]})
xU.sup <- sapply(theta0seq,function(theta0){xgrid[max(which(Ugrid<theta0))]})
xL.sup <-sapply(theta0seq,function(theta0){xgrid[max(which(Lgrid<theta0))]})
xL.inf <- sapply(theta0seq,function(theta0){xgrid[min(which(Lgrid>theta0))]})


#pdf("Figures/figLU.pdf",width=9)
pdf("Figures/figSI1.pdf",width=9)  # w=0.2

par(mfrow=c(1,2))
#ranges <- c(-15,15)
ranges <- c(-25,25)

# Left
plot(xgrid,xgrid,ylab=expression(paste(L[alpha](x)," and ",U[alpha](x))),xlab=expression(x),type="n",xlim=ranges,ylim=ranges,asp=1)
abline(h=c(-lambda,lambda),lty=3,col="grey")
abline(v=xgrid[sapply(5:1,function(r) max(which(regimeU==r)))],lty=3,col="grey")
mtext(c("V","IV","III","II","I"),side=1,adj=0.5,line=-1.5,at=c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU==r)))][-5])+diff(c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU==r)))][-5],ranges[2]))/2   )

for(r in 1:5){
  reg <- which(regimeU==r)
  lines(xgrid[reg],Ugrid[reg],col="black")
  
  reg <- which(regimeL==r)
  lines(xgrid[reg],Lgrid[reg],col="blue")
}

mtext(c(expression(-lambda),expression(lambda)),line=-1.5,side=2,at=c(-lambda,lambda))
legend("topleft",c(expression(U[alpha]),expression(L[alpha])),lty=c(1,1),col=c("black","blue"),horiz=F,bg="white",box.col="white",seg.len=2)

# manual things to make plot look nice (be careful in generalizing this code)
box()
if(w==1){
points(y=c(-lambda,lambda,-lambda,lambda),x=c(max(xgrid[regimeU==4]),min(xgrid[regimeU==3]),max(xgrid[regimeL==3]),min(xgrid[regimeL==4])),
       col=c("black","black","blue","blue"),pch=rep(21,4),bg=c("black","white","blue","white")) # discontinuity jump at lambda
}
if(w < 1){points(y=c(-lambda,0,0,lambda),x=c(min(xgrid[regimeU==3]),min(xgrid[regimeU==3]),max(xgrid[regimeL==3]),max(xgrid[regimeL==3])),
        col=c("black","blue","black","blue"),pch=rep(21,4),bg=c("black","white","black","white")) # discontinuity jump at lambda
}

# Right
plot(theta0seq,xU.sup,type="n",xlim=ranges,ylim=ranges,ylab=expression(paste(X^U,(theta[0])," and ",{X^L}(theta[0]))),col="black",xlab=expression(theta[0]),asp=1)
abline(h=xgrid[sapply(5:1,function(r) max(which(regimeU==r)))],lty=3,col="grey")
mtext(c("V","IV","III","II","I"),side=2, adj=0.5, line=-1.5,at=c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU==r)))][-5])+diff(c(ranges[1],xgrid[sapply(5:1,function(r) max(which(regimeU==r)))][-5],ranges[2]))/2   )

polygon(x=c(theta0seq[theta0seq<(-lambda)],sort(theta0seq[theta0seq<(-lambda)],decreasing=T)),y=c(xL.sup[theta0seq<(-lambda)],(xL.inf[theta0seq<(-lambda)])[order(theta0seq[theta0seq<(-lambda)],decreasing=T)]),col="blue",border=NA,density=20)
polygon(x=c(theta0seq[theta0seq>(lambda)],sort(theta0seq[theta0seq>(lambda)],decreasing=T)),y=c(xU.sup[theta0seq>(lambda)],(xU.inf[theta0seq>(lambda)])[order(theta0seq[theta0seq>(lambda)],decreasing=T)]),col="black",border=NA,density=20)

lines(theta0seq[theta0seq<(-lambda)],xU.sup[theta0seq<(-lambda)])
lines(theta0seq[theta0seq>(lambda)],xU.sup[theta0seq>(lambda)])
lines(theta0seq[theta0seq<(-lambda)],xU.inf[theta0seq<(-lambda)],lty=2)
lines(theta0seq[theta0seq>(lambda)],xU.inf[theta0seq>(lambda)],lty=2)
lines(theta0seq[theta0seq<(-lambda)],xL.sup[theta0seq<(-lambda)],col="blue")
lines(theta0seq[theta0seq>(lambda)],xL.sup[theta0seq>(lambda)],col="blue")
lines(theta0seq[theta0seq<(-lambda)],xL.inf[theta0seq<(-lambda)],lty=2,col="blue")
lines(theta0seq[theta0seq>(lambda)],xL.inf[theta0seq>(lambda)],lty=2,col="blue")
abline(v=c(-lambda,lambda),lty=3,col="grey")
mtext(c(expression(-lambda),expression(lambda)),line=-1.5,side=1,at=c(-lambda,lambda))

legend("topleft",c(expression(paste("sup ",X^U)),expression(paste("inf ",X^U)),expression(paste("sup ",X^L)),expression(paste("inf ",X^L))),
       lty=c(1,2,1,2),col=c("black","black","blue","blue"),ncol=2,bg="white",box.col="white",seg.len=2)


# manual things to make plot look nice (be careful in generalizing this code)
box()
if(w==1){
points(x=c(-lambda,lambda,-lambda,lambda),y=c(min(xgrid[regimeU==4]),min(xgrid[regimeU==3]),max(xgrid[regimeL==3]),max(xgrid[regimeL==4])),
       col=c("black","black","blue","blue"),pch=rep(21,4),bg=c("black","white","blue","white")) 
}
if(w < 1){
  points(x=c(-lambda,lambda,-lambda,lambda),y=c(min(xgrid[regimeU==4]),max(xgrid[regimeU==3]),min(xgrid[regimeL==3]),max(xgrid[regimeL==4])),
         col=c("black","black","blue","blue"),pch=rep(21,4),bg=c("black","white","blue","white"))
  segments(x0=0,x1=0,y0=min(xgrid[Ugrid==0]),y1=max(xgrid[Ugrid==0]),lwd=1.5,col="blue")
  segments(x0=0,x1=0,y0=min(xgrid[Lgrid==0]),y1=max(xgrid[Lgrid==0]),lwd=1.5,col="blue")
  points(x=c(0,0,0,0),y=c(min(xgrid[Ugrid==0]),min(xgrid[Ugrid==0]),min(xgrid[Lgrid==0]),max(xgrid[Lgrid==0])),col=c(1,"blue",1,"blue"),pch=rep(21,4),bg=c("white","black","white","blue"))
}



dev.off()


