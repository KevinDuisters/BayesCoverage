#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# June 2019
#--------------------------------------------------------------------------------------------------------------------#

getU <- function(x,alpha,lambda,w,dist){  
  
  # code chunk to localize g, G, Ginv to memory
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
  if(dist=="t10"){
    g <- function(x,theta0=0){dt(x,10)}
    G <- function(x,theta0=0){pt(x,10)}  
    Ginv <- function(p){qt(p,10)}
  }
  if(dist=="t5"){
    g <- function(x,theta0=0){dt(x,5)}
    G <- function(x,theta0=0){pt(x,5)}  
    Ginv <- function(p){qt(p,5)}
  }
  if(dist=="t3"){
    g <- function(x,theta0=0){dt(x,3)}
    G <- function(x,theta0=0){pt(x,3)}  
    Ginv <- function(p){qt(p,3)}
  }
  if(dist=="t2"){
    g <- function(x,theta0=0){dt(x,2)}
    G <- function(x,theta0=0){pt(x,2)}  
    Ginv <- function(p){qt(p,2)}
  }
  if(dist=="Cauchy"){
    g <- function(x,theta0=0){dcauchy(x,theta0,1)}
    G <- function(x,theta0=0){pcauchy(x,theta0,1)}  
    Ginv <- function(p){qcauchy(p)}
  }
  #
  
  # R functions
  R1f <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha/2 - (1-w)/(2*w)*alpha*g(x)-((1-alpha)/2)*(G(lambda-x)-G(-lambda-x)))))}
  R2f <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha - (1-w)/w*alpha*g(x)+alpha*(G(lambda-x)-G(-lambda-x)) + G(-lambda-x))))}
  R3f <- function(x,alpha,w,lambda){Ginv(pmax(0,pmin(1,1 - alpha/2 - (1-w)/(2*w)*alpha*g(x)+(alpha/2)*(G(lambda-x)-G(-lambda-x)))))}
  D <- function(x,w,lambda){((1-w)/w)*g(x) + 1-(G(lambda - x) - G(-lambda - x))}
  
  
  #aux1 <- 1 - (1/2)*(alpha*(z/w)+Delta(x,lambda))
  #aux2 <- 1 + G(-lambda-x) - alpha*(z/w)
  #aux3 <- 1 - (1/2)*(alpha*(z/w))
  
  # re-used values
  Pi0x <- (1-w)*g(x)/(w*D(x,w,lambda)) 
  R1 <- R1f(x,alpha,w,lambda)
  R2 <- R2f(x,alpha,w,lambda)
  R3 <- R3f(x,alpha,w,lambda)
  
  
  if(Pi0x >= 1 - alpha){# point mass assumes entire credible set already 
    U <- 0
    regime <- 3 # proxy
  }else{
    if(x + R3 > lambda){ # U > lambda; 
      if(x - R1 > lambda){# L(x) > lambda
        U <- x + R1
        regime <- 1
      }else{
        if(x - R3 < (-lambda)){ # L(x) <  -lambda
          U <- x + R3
          regime <- 3
        }else{# L(x) =  lambda
          regime <- 2
          U <- x + R2
        }
      }
    }else{# U <= -lambda
      if(x+R1 < (-lambda)){#U < -lambda
        U <- x + R1
        regime <- 5
      }else{#U = -lambda
        regime <- 4
        U <- (-lambda)
      }
    }
  }
  return(list(val=U,regime=regime)) 
}
