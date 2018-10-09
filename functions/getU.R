getU <- function(x,alpha,lambda,w,dist){  
  
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
  
  z <- (1-w)*g(x) +w*Deltabar(x,lambda)
  aux1 <- 1 - (1/2)*(alpha*(z/w)+Delta(x,lambda))
  aux2 <- 1 + G(-lambda-x) - alpha*(z/w)
  #aux3 <- 1 - (1/2)*(alpha*(z/w) + ((1-w)/w)*g(x))
  aux3 <- 1 - (1/2)*(alpha*(z/w))
  
  if((1-w)*(g(x)/z)>= 1 - alpha){# point mass assumes entire credible set already 
    
    #U <- lambda # or put it at zero? Better for interval length plot later?
    U <- NA
    regime <- 3 # proxy
    
  }else{
    if((G(lambda-x) < aux3)){ # U > lambda; 
    #if((x>0)|(G(lambda-x) < aux3)){ # U > lambda; 
      if(G(x-lambda) > aux1){# L(x) > lambda
        U <- x + Ginv(aux1)
        regime <- 1
      }else{
        if(G(x+lambda)<aux3){ # L(x) <  -lambda
          U <- x + Ginv(aux3)
          regime <- 3
        }else{# L(x) =  lambda
          regime <- 2
          U <- x + Ginv(aux2)
        }
      }
    }else{# U <= -lambda
      if(G(-lambda-x) > aux1){#U < -lambda
        U <- x + Ginv(aux1)
        regime <- 5
      }else{#U = -lambda
        regime <- 4
        U <- -lambda
      }
    }
  }
  return(list(val=U,regime=regime)) 
}
