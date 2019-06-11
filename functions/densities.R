#--------------------------------------------------------------------------------------------------------------------#
# On frequentist coverage of Bayesian credible sets for estimation of the mean under constraints
# Duisters & Schmidt-Hieber
# Math Institute, Leiden University
# Spring 2019
#--------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------#

# 'function' code chunk to localize g,G,Ginv in memory

densities <- function(dist){
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
  if(dist=="t3"){
    g <- function(x,theta0=0){dt(x,3)}
    G <- function(x,theta0=0){pt(x,3)}  
    Ginv <- function(p){qt(p,3)}
  }
  if(dist=="t5"){
    g <- function(x,theta0=0){dt(x,5)}
    G <- function(x,theta0=0){pt(x,5)}  
    Ginv <- function(p){qt(p,5)}
  }
  if(dist=="Cauchy"){
    g <- function(x,theta0=0){dcauchy(x,theta0,1)}
    G <- function(x,theta0=0){pcauchy(x,theta0,1)}  
    Ginv <- function(p){qcauchy(p)}
  }
}