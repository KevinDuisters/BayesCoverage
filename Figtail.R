#----------------------------------------------------------------------------------------------------------------------------#
# Duisters, Schmidt-Hieber
# June 2019
#----------------------------------------------------------------------------------------------------------------------------#
# Figure tail decay
#----------------------------------------------------------------------------------------------------------------------------#
# load library
library(rmutil)
#----------------------------------------------------------------------------------------------------------------------------#
# settings
lambda <-  0.75
alphaseq <- seq(0.01,0.5,0.01)
x <- seq(lambda,lambda+20,0.05)    # X_{II} \cap (\lambda,\infty)
#wseq <- seq(0.01,1,0.005)
wseq <- seq(0.5,1,0.005)

#----------------------------------------------------------------------------------------------------------------------------#
# result
mat <- matrix(NA,length(alphaseq),length(wseq))
rownames(mat) <- alphaseq
colnames(mat) <- wseq
result <-tail<- lapply(1:7,function(l) mat)

#----------------------------------------------------------------------------------------------------------------------------#
# loop
for(d in 1:7){

  if(d==1){
    g <- function(x){dnorm(x,0,1)}
    G <- function(x){pnorm(x,0,1)}
    Ginv <- function(p){qnorm(p,0,1)}
  }
 
  if(d==2){
    g <- function(x){dlaplace(x)}
    G <- function(x){plaplace(x)}
    Ginv <- function(p){qlaplace(p)}
  }
  if(d==3){
    g <- function(x){dt(x,10)}
    G <- function(x){pt(x,10)}
    Ginv <- function(p){qt(p,10)}
  }
  if(d==4){
    g <- function(x){dt(x,5)}
    G <- function(x){pt(x,5)}
    Ginv <- function(p){qt(p,5)}
  }
  if(d==5){
    g <- function(x){dt(x,3)}
    G <- function(x){pt(x,3)}
    Ginv <- function(p){qt(p,3)}
  }
  if(d==6){
    g <- function(x){dt(x,2)}
    G <- function(x){pt(x,2)}
    Ginv <- function(p){qt(p,2)}
  }
  if(d==7){
    g <- function(x){dcauchy(x)}
    G <- function(x){pcauchy(x)}
    Ginv <- function(p){qcauchy(p)}
  }
  
  R1 <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha/2 - (1-w)/(2*w)*alpha*g(x)-((1-alpha)/2)*(G(lambda-x)-G(-lambda-x)))))}
  R2 <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha - (1-w)/w*alpha*g(x)+alpha*(G(lambda-x)-G(-lambda-x)) + G(-lambda-x))))}
  R3 <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha/2 - (1-w)/(2*w)*alpha*g(x)+(alpha/2)*(G(lambda-x)-G(-lambda-x)))))}
  
  
  for(j in 1:length(wseq)){
  w <- wseq[j]
  for(i in 1:length(alphaseq)){
    alpha <- alphaseq[i]
    
    R1vec <- R2vec <- R3vec <- targetvec <-numeric(length(x))
    target <- (w/alpha)*((1-alpha)/(1-w))
    R1vec <- R1(x,alpha,w,lambda)
    R2vec <- R2(x,alpha,w,lambda)
    R3vec <- R3(x,alpha,w,lambda)
  
    # check if x in X_{II} \cap (\lambda,\infty) (otherwise no target has to be satisfied)
    region <- ((-lambda + R3vec < x)  & (x < lambda + R1vec) & (x > lambda))
    result[[d]][i,j]<-(sum((g(x[region]))/G(-lambda-x[region]) < target) == sum(region))
    #result[[d]][i,j]<-(sum((g(x[region]))/(G(lambda-x[region])-G(-lambda-x[region])) < (w/(1-w))   ) == sum(region))
    
    
    # tail decay condition 
    #q <- Ginv(alpha/(1+alpha)) + Ginv(alpha/((2-w)*alpha+w  ))
    #tail[[d]][i,j]<-((G(q)/G(q/3)) < (2*alpha)/(1-alpha))
    
    # new tail decay
    #xstar <- max(x[region])
    #q <- min(1,((2*alpha)/(2+alpha))*(1 + (1-w)/w*g(xstar)))
    #tail[[d]][i,j]<-(G(2*Ginv(q)) < alpha/(1-alpha)*q)
    #tail[[d]][i,j] <- (G(3*Ginv((alpha/(1+alpha))*(1+(1-w)/w*g(xstar)))) < alpha/(1-alpha)*q)
    tail[[d]][i,j] <- (G(2*Ginv(alpha/w))<(2*alpha^2)/(1-alpha^2))
  }
}
}

# assumption
plot(0,0,xlim=range(wseq),ylim=range(alphaseq),ylab="alpha",xlab="w",type="n")
title("g(x)/G(-lambda-x) < (w/alpha)*(1-alpha)/(1-w) for x in X_II")
mtext(paste("lambda =",lambda))

for(d in 1:7){
alphaline <- apply(result[[d]],2,function(c)alphaseq[min(which(c==F))])
lines(smooth.spline(wseq[is.na(alphaline)==F],alphaline[is.na(alphaline)==F]),col=d)
}
legend("topleft",c("N(0,1)","Lap(0,1)","t10","t5","t3","t2","Cauchy"),lty=rep(1,7),col=1:7,bty="n")
abline(h=c(0,0.05,0.1),lty=3)


# tail decay
plot(0,0,xlim=range(wseq),ylim=range(alphaseq),ylab="alpha",xlab="w",type="n")
title("tail decay")
mtext("does not depend on lambda")
for(d in 1:7){
  alphaline <- apply(tail[[d]],2,function(c)alphaseq[min(which(c==T))])
  #lines(wseq,alphaline,col=d)
  lines(smooth.spline(wseq[is.na(alphaline)==F],alphaline[is.na(alphaline)==F]),col=d)
}
legend("topleft",c("N(0,1)","Lap(0,1)","t10","t5","t3","t2","Cauchy"),lty=rep(1,7),col=1:7,bty="n")
abline(h=c(0,0.05,0.1),lty=3)
