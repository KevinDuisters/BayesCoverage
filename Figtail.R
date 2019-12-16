#----------------------------------------------------------------------------------------------------------------------------#
# Duisters, Schmidt-Hieber
# Sep 2019
#----------------------------------------------------------------------------------------------------------------------------#
# Figure tail decay
#----------------------------------------------------------------------------------------------------------------------------#
# load library
library(rmutil)
#library(RColorBrewer)
#----------------------------------------------------------------------------------------------------------------------------#
# settings
h <- 0.001
alphaseq <- seq(h,1-h,h)
#wseq <- c(alphaseq,1)
wseq <- seq(0,1,h)

lambda1 <- 0.5
densities <- c("N(0,1)","Lap(0,1)","t10","t5","t3","t2","Cauchy(0,1)")
cols <- c("black","blue","darkgreen","chartreuse4","aquamarine3","green","red")
ltys <- c(5,1,3,3,3,3,4)
ld <- length(densities)

# lambda sensitivity
lambdaseq <- c(0.25,0.5,0.75,1,2.5,5,7.5)
ll <- length(lambdaseq)
lambdaltys <- c(2,1,3,4,5,1,6)
lambdacols <- rep(cols[which(densities=="Lap(0,1)")],ll)


#----------------------------------------------------------------------------------------------------------------------------#
# result
mat <- matrix(NA,length(alphaseq),length(wseq))
rownames(mat) <- alphaseq
colnames(mat) <- wseq
tail<- lapply(1:ld,function(l) mat)
Laptails <- lapply(1:ll,function(l)mat)

#----------------------------------------------------------------------------------------------------------------------------#
# loop
for(d in 1:ld){
dens <- densities[d]
  
if(dens=="N(0,1)"){
    g <- function(x){dnorm(x,0,1)}
    G <- function(x){pnorm(x,0,1)}
    Ginv <- function(p){qnorm(p,0,1)}
  }
  if(dens=="t10"){
    g <- function(x){dt(x,10)}
    G <- function(x){pt(x,10)}
    Ginv <- function(p){qt(p,10)}
  }
  if(dens=="t5"){
    g <- function(x){dt(x,5)}
    G <- function(x){pt(x,5)}
    Ginv <- function(p){qt(p,5)}
  }
  if(dens=="t3"){
    g <- function(x){dt(x,3)}
    G <- function(x){pt(x,3)}
    Ginv <- function(p){qt(p,3)}
  }
  if(dens=="t2"){
    g <- function(x){dt(x,2)}
    G <- function(x){pt(x,2)}
    Ginv <- function(p){qt(p,2)}
  }
  if(dens=="Lap(0,1)"){
    g <- function(x){dlaplace(x)}
    G <- function(x){plaplace(x)}
    Ginv <- function(p){qlaplace(p)}
  }
  if(dens=="Cauchy(0,1)"){
    g <- function(x){dcauchy(x)}
    G <- function(x){pcauchy(x)}
    Ginv <- function(p){qcauchy(p)}
  }
  
for(i in 1:(length(alphaseq)-1)){  
  alpha <- alphaseq[i]
  for(j in (i+1):length(wseq)){
      w <- wseq[j]
      q <- Ginv((alpha/(1+alpha))*(1 + ((1-w)/w)*g(lambda1)) )
      if(q < 0){tail[[d]][i,j] <- (G(2*q)<(2*alpha^2)/(1-alpha^2))}else{tail[[d]][i,j] <- NA}
      
      # lambda plot
      if(dens=="Lap(0,1)"){
      for(l in 1:ll){
        ql <- Ginv((alpha/(1+alpha))*(1 + ((1-w)/w)*g(lambdaseq[l])))
        if(ql < 0){Laptails[[l]][i,j] <- (G(2*ql)<(2*alpha^2)/(1-alpha^2))}else{Laptails[[l]][i,j] <- NA}
      }
      }
    }
  }
}

# tail decay
pdf("Figures/figtail.pdf",width=9,height=9)
plot(0,0,xlim=range(wseq),ylim=range(alphaseq),ylab=expression(alpha),xlab="w",type="n",asp=1,cex.axis=2,cex.lab=2)
#polygon(x=c(-2,2,-2,-2),y=c(-2,2,2,-2),col="lightgrey",border="white")
abline(h=seq(0,1,0.1),lty=2,col="lightgrey")
abline(v=seq(0,1,0.1),lty=2,col="lightgrey")
for(a in seq(0,1,0.1)){segments(x0=-1,x1=a,y0=a,y1=a,lty=2,col="white")}
box()
for(d in 1:7){
  alphaline <- apply(tail[[d]],2,function(c)alphaseq[min(which(c==T))])
  lines(wseq,alphaline,lty=ltys[d],col=cols[d],lwd=2)
}
legend(x=0,y=1,densities[ld:1],lty=ltys[ld:1],col=cols[ld:1],lwd=rep(2,ld),bg="white",cex=2)
dev.off()



# lambda sensitivity
pdf("Figures/figtaillambdaLap.pdf",width=9,height=9)
plot(0,0,xlim=range(wseq),ylim=range(alphaseq),ylab=expression(alpha),xlab="w",type="n",asp=1,cex.axis=2,cex.lab=2)
#polygon(x=c(-2,2,-2,-2),y=c(-2,2,2,-2),col="lightgrey",border="white")
abline(h=seq(0,1,0.1),lty=2,col="lightgrey")
abline(v=seq(0,1,0.1),lty=2,col="lightgrey")
for(a in seq(0,1,0.1)){segments(x0=-1,x1=a,y0=a,y1=a,lty=2,col="white")}
box()
for(l in 1:ll){
  alphaline <- apply(Laptails[[l]],2,function(c)alphaseq[min(which(c==T))])
  lines(wseq,alphaline,lty=lambdaltys[l],col=lambdacols[l],lwd=2)
}
mylabels <- sapply(1:ll,function(l)as.expression(bquote(lambda==.(lambdaseq[l]))))
legend(x=0,y=1,mylabels,lty=lambdaltys,col=lambdacols,lwd=rep(2,ll),bg="white",cex=2)
dev.off()

