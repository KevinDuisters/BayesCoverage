post <- function(thetaseq,x,lambda,w,dist){
  
  if(dist=="Lap"){
    g <- function(x,theta0=0){dlaplace(x,m=theta0,s=1)}
    G <- function(x,theta0=0){plaplace(x,m=theta0,s=1)}  
    Ginv <- function(p){qlaplace(p)}
  }
  if(dist=="Normal"){
    g <- function(x,theta0=0){dnorm(x,theta0,1)}
    G <- function(x,theta0=0){pnorm(x,theta0,1)}
    Ginv <- function(p){qnorm(p)}  # make precise by always considering small p as opposed to large  
  }
  Delta <- function(x,lambda){G(lambda-x) - G(-lambda-x)  }
  Deltabar <- function(x,lambda){G(x-lambda)+G(-lambda-x)}  # 1-Delta is much more precise
  
  return(sapply(thetaseq,function(theta){
    if(theta <= lambda & theta >= -lambda & theta!=0){return(0)}else{
      z <- (1-w)*g(x) + w*Deltabar(x,lambda)
      if(theta==0){(1-w)*g(x)/z}else{
        w*g(x-theta)/z
      }
    }
  }))
} # vectorized function for integrate()