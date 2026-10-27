#rm(list=ls())

#### Fun??es para o processo Gama
library(maxLik)
library(invgamma)
library(nortest)
library(tseries)
library(forecast)
#library(ggplot2)
#library(reshape2)
library(purrr)
library(robustbase)
library(lamW)

# <log-likelihood function for K-Bessel regression >
loglink_KBR = function(Y,par,L){
  # L: variavel global (numero de looks)
  # Xobs: variavel global (vetor de dados)
  # matrix_model: variavel global (matriz modelo)
  # Xobs: variavel global (vetor de respostas)
  
  alpha = par[1]; 
  delta = par[2];
  phi = par[3];
  
  mu[1] = 
  
  res <- N * log(2) +
    N * 0.5 * (L + alpha) * log(L) -
    N * lgamma(alpha) - N * lgamma(L) +
    N * 0.5 * (L + alpha) * log(alpha) -
    N * 0.5 * (L + alpha) * bet0 -
    0.5 * (L + alpha) * bet1 * sum(matrix_model[,2]) +
    0.5 * (L + alpha-2) * sum(log(Xobs)) +
    sum(
      log(
        besselK(
          2 * sqrt( L * alpha * Xobs / exp(auxi) ),
          alpha - L)
      ),na.rm=T
    )
  
  return(res)  
}

Chute.inicial <- function(Y,x,L){
  X = cbind(rep(1,length(x)),x)
  
  
  alpha = seq(1,100,10)
  b1 = solve( t(X)%*%X   )%*%t(X)%*%log(Y)
  if(b1[2] > 20){
    b1[2] = 20
  }else if( b1[2] < -20 ){
    b1[2] = -20
  }  
  
  if(b1[1] > 30){
    b1[1] = 30
  }else if( b1[1] < -30 ){
    b1[1] = -30
  }  
  
  par = cbind( alpha, rep(  b1[1],length(alpha)  ),rep(  b1[2],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = loglink_KBR( Y, x, par[i,], L=4  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, b1=b1  ) )
  
}
cml.K <- function(z, x, par, L, method){
  # L: variavel global (numero de looks)
  # x: variavel global (matriz modelo)
  # z: variavel global (vetor de respostas)
  # par: chute inicial
  #method: metodo de otimizacao: BFGS ou SMCO
  #A = matrix(c(1,0),1,2, byrow = T)
  #B = c(-0.001)
  
  A = matrix(c(1,0,0,-1,0,0,0,1,0,0,-1,0,0,0,1,0,0,-1),6,3, byrow = T)
  B = c(-0.01, 100,100,100,30,30)
  matrix_model = cbind(rep(1,length(x)),x)
  
  lk.new = function(par)  loglink_KBR(z,x, par,L)
  
  if( method == "BFGS" ){
    EST = maxLik(lk.new,start=par,method="BFGS",constraints=list(ineqA = A, ineqB = B))
    Estimate = EST$estimate; Code = EST$code;eta = (matrix_model)%*%c(Estimate[2],Estimate[3]) ; mu.hat = exp( eta )
    
  }else if( method == "smco"){
    EST=smco(par = par,
             LB = c(0.1,-20,-20), UB = c(50,20,20),
             fn = -lk.new,
             maxiter = 10000, trc = FALSE)
    Estimate = EST$par; Code = 0
  }
  
  return(list( estimate = Estimate, code = Code, mu.hat=mu.hat  ))
  
}
    DistK = function( alpha, L, mu, y ){
  X =  ( 1/gamma(L)*gamma(alpha) ) * meijerg( alpha, L, mu, y  )[,1]
  return(X)
}

#### Geradores ####

# Gera??o distribui??o K


rK_AR1_log = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3]
  mu=c();y=c();lambda=c()
  y0 = 1
  
  mu[1] = exp(
    delta + phin*( log(y0) )
  )
  lambda[1] = alpha/mu[1]
  y[1] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1]  )
  
  for(t in 2:(N+1)){
    mu[t] = exp(
      delta + phin*( log(y[t-1]) ) 
    )
    lambda[t] = alpha/mu[t]
    y[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )
  }
  
  return( y[2:(N+1)] )
  
}
rK_AR2_log = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p = 2
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  z0 = 1
  
  mu[1] = exp(
    delta + phin[1]*( log(z0) ) + phin[2]*( log(z0) )
  )
  z[1] = rinvgamma( 1, shape=-alpha, rate=mu[1]*(-alpha-1) ) 
  
  mu[2] = exp(
    delta + phin[1]*( log(z[1]) ) + phin[2]*( log(z0) ) 
  )
  z[2] = rinvgamma( 1, shape=-alpha, rate=mu[2]*(-alpha-1) ) 
  
  for(t in 3:(N+1)){
    mu[t] = exp(
      delta + phin[1]*( log(z[t-1]) ) + phin[2]*( log(z[t-2]) ) 
    )
    z[t] = rinvgamma( 1, shape=-alpha, rate=mu[t]*(-alpha-1) )
  }
  
  au = rgamma(N, shape=look, scale=1/look); 
  x = au*z
  return(x[3:(N+1)])
}
rK_MA1_log = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  alpha = vP[1]; look = L; delta = vP[2]; theta = vP[3]
  mu=c();y=c();lambda=c()
  y0 = 1; mu0 = 1
  
  mu[1] = exp(
    delta + theta*( log(y0) - log(mu0) )
  )
  lambda[1] = alpha/mu[1]
  y[1] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1]  )
  
  for(t in 2:(N+1)){
    mu[t] = exp(
      delta + theta*( log(y[t-1]) - log(mu[t-1]) )
    )
    lambda[t] = alpha/mu[t]
    y[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )
  }
  
  return( y[2:(N+1)] )
  
}
rK_ARMA11_log = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  alpha = vP[1]; look = L; delta = vP[2]; phi = vP[3];theta = vP[4]
  mu=c();y=c();lambda=c()
  y0 = 1; mu0 = 1
  
  mu[1] = exp(
    delta + phi*log(y0) + theta*( log(y0) - log(mu0) )
  )
  lambda[1] = alpha/mu[1]
  y[1] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1]  )
  
  for(t in 2:(N+1)){
    mu[t] = exp(
      delta + phi*log(y[t-1]) + theta*( log(y[t-1]) - log(mu[t-1]) )
    )
    lambda[t] = alpha/mu[t]
    y[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )
  }
  
  return( y[2:(N+1)] )
  
}

rK_AR1_W = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3]
  mu=c();z=c()
  z0 = 1
  
  require(lamW)
  # to find optimal lambda
  eta=c()
  lambda = c()

  eta[1] = delta + phin*( log(lambertW0( z0 )) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )  
  lambda[1] = alpha/mu[1]
  z[1] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1]  )
  
  for(t in 2:(N+1)){
    eta[t] = delta + phin*( log(lambertW0( z[t-1] )) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    lambda[t] = alpha/mu[t]
    z[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )

  }
  
  return( z[2:(N+1)] )
  
}
rK_AR2_W = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p = 2
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  z0 = 1
  lambda=c()
  require(lamW)
  eta=c()
  eta[1] = delta + phin[1]*( log(lambertW0(1)) ) + phin[2]*( log(lambertW0( 1 )) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  lambda[1] = alpha/mu[1]
  z[1] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1]  )
  
  eta[2] = delta + phin[1]*( log(lambertW0(z[1])) ) + phin[2]*( log(lambertW0( 1 )) ) 
  mu[2] = exp( eta[2] + exp(eta[2])  )  
  lambda[2] = alpha/mu[2]
  z[2] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[2]  )
  

  for(t in 3:(N+1)){
    eta[t] = delta + phin[1]*( log(lambertW0( z[t-1] )) ) + phin[2]*( log(lambertW0( z[t-2] )) )
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    lambda[t] = alpha/mu[t]
    z[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )
    }
  
  return(z[3:(N+1)])
}
rK_MA1_W = function(N,par,L){
  # N: number of outcomes
  # par: vector of parameters  
  
  q=1; look = L
  alpha = par[1]; delta = par[2]; theta = par[3]
  mu=c();z=c();
  lambda=c()
  
  require(lamW)
  eta=c()
  eta[1] = delta 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  lambda[t] = alpha/mu[1]
  z[1] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1]  )
  
  for(t in 2:(N+1)){
    eta[t] = delta + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    lambda[t] = alpha/mu[t]
    z[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )
    }
  
  return(z[2:(N+1)])
}
rK_ARMA11_W = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3];theta = vP[4]
  mu=c();z=c()
  z0 = 1;mu0=1
  eta=c()
  lambda = c()
  require(lamW)
  
  eta[1] = delta + phin*( log(lambertW0(1)) ) + theta*( log(lambertW0(1)) - log(lambertW0(1)) ) 
  mu[1] = 1
  lambda[1] = alpha/mu[1]
  z[1] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1]  )  
  
  for(t in 2:(N+1)){
    eta[t] = delta + phin*( log(lambertW0( z[t-1] )) ) + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    lambda[t] = alpha/mu[t]
    z[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )
    }
  
  return(z[2:(N+1)])
}

rK_AR1_sqrt = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3]
  mu=c();z=c()
  z0 = 1
  eta=c();lambda=c()
  
  eta[1] = delta + phin*( sqrt(z0) )
  mu[1] = eta[1]^2
  lambda[1] = alpha/mu[1]
  z[1] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1]  )  
  
  for(t in 2:(N+1)){
    eta[t] = delta + phin*( sqrt(z[t-1]) ) 
    mu[t] = eta[t]^2
    lambda[t] = alpha/mu[t]
    z[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )  
    }
  
  return(z[2:(N+1)])
}
rK_AR2_sqrt = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p = 2
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  z0 = 1
  eta=c();lambda=c()
  
  eta[1] = delta + phin[1]*( sqrt(z0) ) + phin[2]*( sqrt(z0) )
  mu[1] = eta[1]^2
  lambda[1] = alpha/mu[1]
  z[1] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1]  )
  
  eta[2] = delta + phin[1]*( sqrt(z[1]) ) + phin[2]*( sqrt(z0) ) 
  mu[2] = eta[2]^2
  lambda[2] = alpha/mu[2]
  z[2] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[2]  )
  
  for(t in 3:(N+2)){
    eta[t] = delta + phin[1]*( sqrt(z[t-1]) ) + phin[2]*( sqrt(z[t-2]) ) 
    mu[t] = eta[t]^2
    lambda[t] = alpha/mu[t]
    z[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )
  }
  
  return(z[3:(N+2)])
}
rK_MA1_sqrt = function(N,par,L){
  # N: number of outcomes
  # par: vector of parameters  
  
  q=1; look = L
  alpha = par[1]; delta = par[2]; theta = par[3]
  mu=c();z=c();
  eta=c();lambda=c()
  
  eta[1] = delta + theta*(  sqrt( 1 ) - sqrt( 1 )  )
  mu[1] = 1
  lambda[1] = alpha/mu[1]
  z[1] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1]  )
  
  for(t in 2:(N+1)){
    eta[t] = delta + theta*( sqrt( z[t-1] ) - sqrt( mu[t-1] )  )
    mu[t] = eta[t]^2
    lambda[t] = alpha/mu[t]
    z[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )  }
  
  return(z[2:(N+1)])
}
rK_MA2_sqrt = function(N,par,L){
  # N: number of outcomes
  # par: vector of parameters  
  
  q=1; look = L
  alpha = par[1]; delta = par[2]; theta = par[3:4]
  mu=c();z=c();
  eta=c();lambda=c()
  
  eta[1:2] = delta + theta[1]*(  sqrt( 1 ) - sqrt( 1 )  ) + theta[2]*(  sqrt( 1 ) - sqrt( 1 )  )
  mu[1:2] = 1
  lambda[1:2] = alpha/mu[1:2]
  z[1:2] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1:2]  )
  
  for( t in 3:(N+1) ){
    eta[t] = delta + theta[1]*( sqrt( z[t-1] ) - sqrt( mu[t-1] )  ) + theta[2]*( sqrt( z[t-2] ) - sqrt( mu[t-2] )  )
    mu[t] = eta[t]^2
    lambda[t] = alpha/mu[t]
    z[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )  }
  
  return(z[2:(N+1)])
}
rK_ARMA11_sqrt = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3];theta = vP[4]
  mu=c();z=c()
  z0 = 1;mu0=1
  eta=c();lambda=c()
  
  eta[1] = delta + phin*( sqrt(z0) ) + theta*( sqrt(z0) - sqrt(mu0) )
  mu[1] =eta[1]^2
  lambda[1] = alpha/mu[1]
  z[1] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[1]  )  
  
  for(t in 2:(N+1)){
    eta[t] = delta + phin*( sqrt(z[t-1]) ) + theta*( sqrt(z[t-1]) - sqrt(mu[t-1]) ) 
    mu[t] = eta[t]^2
    lambda[t] = alpha/mu[t]
    z[t] = rgamma(1, shape=L, scale=1/L) * rgamma(1, shape=alpha, scale=1/lambda[t]  )  }
  
  return(z[2:(N+1)])
}

## Gera??o da distribui??o gamma


####    Log verossimilhan?as     #####

# Distribui??o K

#AR(1)
lk.K.AR1_log = function(y,par,L){
  
  ## ---> Input
  ## y: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  alpha = par[1]; delta = par[2]; phin=par[3]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  y0 = rep(mean(y),p)
  
  mu[1] = exp(
    delta + phin*( log(y0) )
  )
  
  lk[1] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(y[1]) - log(mu[1])   )  -
    log(y[1]) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (y[1]/mu[1]) ), alpha - L) )
  
  for(t in (p+1):n){
    mu[t] = exp(
      delta + phin*( log(y[t-1]) ) 
    )
  }
  
    lk = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(y) - log(mu)   )  -
       log(y) - lgamma(L) - lgamma(alpha) + 
       log( besselK( 2 * sqrt( L * alpha * (y/mu) ), alpha - L) )
  
  return( sum( lk[2:n], na.rm = T) )
  
}
lk.K.AR2_log = function(y,par,L){
  ## ---> Input
  ## y: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## <-- Output
  ## log.verossimilhan?a 
  
  p = 2
  n = length(y)
  alpha = par[1]; delta = par[2]; phin=par[3:4]
  mu = c();v=c();lk=c()
  eta=c()
  y0 = mean(y);mu0=1
  
  mu[1] = exp(
    delta + phin[1]*( log(y0) ) + phin[2]*( log(y0) ) 
  )
  mu[2] = exp(
    delta + phin[1]*( log(y[1]) ) + phin[2]*( log(y0) ) 
  )
  
  lk[1:2] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(y[1:2]) - log(mu[1:2])   )  -
    log(y[1:2]) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (y[1:2]/mu[1:2]) ), alpha - L) )
  
  for(t in (p+1):n){
    mu[t] = exp(
      delta + phin[1]*( log(y[t-1]) ) + phin[2]*( log(y[t-2]) ) 
    )
  }
  
  lk = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(y) - log(mu)   )  -
    log(y) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (y/mu) ), alpha - L) )
  
  return( sum( lk[3:n], na.rm = T) )
  
}
lk.K.MA1_log = function(y,par,L){
  
  ## ---> Input
  ## y: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  alpha = par[1]; delta = par[2]; theta=par[3]
  mu=c();lk=c()
  
  mu0 = 1
  y0 = mean(y)
  
  mu[1] = exp(
    delta + theta*( log(y0) - log(mu0) )
  )
  
  lk[1] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(y[1]) - log(mu[1])   )  -
    log(y[1]) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (y[1]/mu[1]) ), alpha - L) )
  
  for(t in (p+1):n){
    mu[t] = exp(
      delta + theta*( log(y[t-1]) - log(mu[t-1]) ) 
    )
  }
  
  lk = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(y) - log(mu)   )  -
    log(y) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (y/mu) ), alpha - L) )
  
  return( sum( lk[2:n], na.rm = T) )
  
}
lk.K.MA2_log = function(y,par,L){
  
  ## ---> Input
  ## y: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=2
  n=length(y)
  alpha = par[1]; delta = par[2]; theta=par[3:4]
  mu=c();lk=c()
  
  mu0 = 1
  y0 = mean(y)
  
  mu[1] = exp(
    delta + theta[1]*( log(y0) - log(mu0) ) + theta[2]*( log(y0) - log(mu0) )
  )
  mu[2] = exp(
    delta + theta[1]*( log(y[1]) - log(mu[1]) ) + theta[2]*( log(y0) - log(mu0) )
  )
  
  lk[1:2] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(y[1:2]) - log(mu[1:2])   )  -
    log(y[1:2]) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (y[1:2]/mu[1:2]) ), alpha - L) )
  
  for(t in (p+1):n){
    mu[t] = exp(
      delta + theta[1]*( log(y[t-1]) - log(mu[t-1]) ) + theta[2]*( log(y[t-2]) - log(mu[t-2]) )
    )
  }
  
  lk = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(y) - log(mu)   )  -
    log(y) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (y/mu) ), alpha - L) )
  
  return( sum( lk[3:n], na.rm = T) )
  
}
lk.K.ARMA11_log = function(y,par,L){
  
  ## ---> Input
  ## y: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  alpha = par[1]; delta = par[2]; phi=par[3]; theta = par[4]
  mu=c();v=c();lk=c()
  
  mu0 = 1
  y0 = mean(y)
  
  mu[1] = exp(
    delta + phi*( log(y0) ) + theta*( log(y0) - log(mu0) )
  )
  
  lk[1] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(y[1]) - log(mu[1])   )  -
    log(y[1]) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (y[1]/mu[1]) ), alpha - L) )
  
  for(t in (p+1):n){
    mu[t] = exp(
      delta + phi*( log(y[t-1]) ) + theta*( log(y[t-1]) - log(mu[t-1]) )
    )
  }
  
  lk = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(y) - log(mu)   )  -
    log(y) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (y/mu) ), alpha - L) )
  
  return( sum( lk[2:n], na.rm = T) )
  
}


lk.K.AR1_W = function(z,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(z)
  alpha = par[1]; delta = par[2]; phin=par[3]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  z0 = rep(mean(z),p)
  
  require(lamW)
  # to find optimal lambda
  eta=c()
  
  eta[1] = delta + phin*( log(lambertW0(1)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )  
  
  lk[1] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(z[1]) - log(mu[1])   )  -
    log(z[1]) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (z[1]/mu[1]) ), alpha - L) )  
  
  for(t in (p+1):n){
    eta[t] = delta + phin*( log(lambertW0(z[t-1])) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    lk[t] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(z[t]) - log(mu[t])   )  -
      log(z[t]) - lgamma(L) - lgamma(alpha) + 
      log( besselK( 2 * sqrt( L * alpha * (z[t]/mu[t]) ), alpha - L) )    
  }
  
  return(sum(lk[2:n]))
  
}
lk.K.AR2_W = function(z,par,L){
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## <-- Output
  ## log.verossimilhan?a 
  
  p = 2
  n = length(z)
  alpha = par[1]; delta = par[2]; phin=par[3:length(par)]
  mu = c();v=c();lk=c()
  z0 =mean(z)
  
  require(lamW)
  eta=c()
  eta[1] = delta + phin[1]*( log(lambertW0(z0)) ) + phin[2]*( log(lambertW0( z0 )) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  eta[2] = delta + phin[1]*( log(lambertW0(z[1])) ) + phin[2]*( log(lambertW0( z0 )) ) 
  mu[2] = exp( eta[2] + exp(eta[2])  )  
  
  
  lk[1:2] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(z[1:2]) - log(mu[1:2])   )  -
    log(z[1:2]) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (z[1:2]/mu[1:2]) ), alpha - L) )
  
  for(t in (p+1):n){
    eta[t] = delta + phin[1]*( log(lambertW0(z[t-1])) ) + phin[2]*( log(lambertW0( z[t-2] )) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    lk[t] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(z[t]) - log(mu[t])   )  -
      log(z[t]) - lgamma(L) - lgamma(alpha) + 
      log( besselK( 2 * sqrt( L * alpha * (z[t]/mu[t]) ), alpha - L) )    
  }
  
  return( sum(lk[3:n]) )
  
}
lk.K.MA1_W <- function(z, par,L){
  
  n=length(z)
  alpha = par[1]
  delta<- par[2]
  theta <- par[3]
  m=1;q=1
  
  require(lamW)
  eta <- mu <- v <- c()
  
  eta[1] = delta + theta*( log(lambertW0(1)) - log(lambertW0(1)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  
  for(t in (m+1):n){
    eta[t] <- delta + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) )
    mu[t] <- exp( eta[t] + exp(eta[t])  )   
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  lv = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(zn) - log(mun)   )  -
    log(zn) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (zn/mun) ), alpha - L) )
  
    return(sum(lv))
} 
lk.K.MA2_W <- function(z, par,L){
  
  n=length(z)
  alpha = par[1]
  delta<- par[2]
  theta <- par[3:4]
  m=2;z0=mean(z);mu0=1
  
  require(lamW)
  eta <- mu <- v <- c()
  
  eta[1] <- delta + theta[1]*( log(lambertW0(z0)) - log(lambertW0(mu0)) )+ theta[2]*( log(lambertW0(z0)) - log(lambertW0(mu0)) )
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  eta[2] <- delta + theta[1]*( log(lambertW0(z[1])) - log(lambertW0(mu[1])) ) + theta[2]*( log(lambertW0(z0)) - log(lambertW0(mu0)) )
  mu[2] = exp( eta[2] + exp(eta[2])  )   
  
  t=3
  for(t in (m+1):n){
    eta[t] <- delta + theta[1]*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) + theta[2]*( log(lambertW0(z[t-2])) - log(lambertW0(mu[t-2])) )
    mu[t] <- exp( eta[t] + exp(eta[t])  )   
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]

  lv = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(zn) - log(mun)   )  -
    log(zn) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (zn/mun) ), alpha - L) )
  
  return(sum(lv))
} 
lk.K.ARMA11_W <- function(z, par,L){
  
  n=length(z)
  alpha <- par[1];delta <- par[2];phi <- par[3];theta <- par[4]
  
  m=1;q=1;p=1
  
  z0 = mean(z)
  
  eta <- mu <- v <- c()
  eta[1] = delta + phi*( log(lambertW0(z0)) ) + theta*( log(lambertW0(z0)) - log(lambertW0(z0)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )   

  for(t in (m+1):n){
    eta[t] <- delta + phi*( log(lambertW0(z[t-1])) ) + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) 
    mu[t]   <- exp( eta[t] + exp(eta[t])  )   
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]

  lv = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(zn) - log(mun)   )  -
    log(zn) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (zn/mun) ), alpha - L) )
  
  return(sum(lv))
} 


lk.K.AR1_sqrt = function(z,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(z)
  alpha = par[1]; delta = par[2]; phin=par[3]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  z0 = rep(mean(z),p)
  eta=c()
  eta[1] = delta + phin*( sqrt(z0) )
  mu[1] = eta[1]^2
  
  lk[1] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(z[1]) - log(mu[1])   )  -
    log(z[1]) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (z[1]/mu[1]) ), alpha - L) )    
  
  for(t in (p+1):n){
    eta[t] = delta + phin*( sqrt(z[t-1]) ) 
    mu[t] = eta[t]^2
    lk[t] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(z[t]) - log(mu[t])   )  -
      log(z[t]) - lgamma(L) - lgamma(alpha) + 
      log( besselK( 2 * sqrt( L * alpha * (z[t]/mu[t]) ), alpha - L) )  
  }
  
  return(sum(lk[2:n]))
  
}
lk.K.AR2_sqrt = function(z,par,L){
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## <-- Output
  ## log.verossimilhan?a 
  
  p = 2
  n = length(z)
  alpha = par[1]; delta = par[2]; phin=par[3:4]
  mu = c();v=c();lk=c()
  eta=c()
  z0 = mean(z)
  
  eta[1] = delta + phin[1]*( sqrt(z0) ) + phin[2]*( sqrt(z0) ) 
  mu[1] = eta[1]^2
  eta[2] = delta + phin[1]*( sqrt(z[1]) ) + phin[2]*( sqrt(z0) ) 
  mu[2] = eta[2]^2
  
  lk[1] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(z[1]) - log(mu[1])   )  -
    log(z[1]) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (z[1]/mu[1]) ), alpha - L) )    
  lk[2] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(z[2]) - log(mu[2])   )  -
    log(z[2]) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (z[2]/mu[2]) ), alpha - L) )  
  
  
    for(t in (p+1):n){
    eta[t] = delta + phin[1]*( sqrt(z[t-1]) ) + phin[2]*( sqrt(z[t-2]) ) 
    mu[t] = eta[t]^2
    lk[t] = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(z[t]) - log(mu[t])   )  -
      log(z[t]) - lgamma(L) - lgamma(alpha) + 
      log( besselK( 2 * sqrt( L * alpha * (z[t]/mu[t]) ), alpha - L) )      
  }
  
  return(sum(lk[3:n]))
  
}
lk.K.MA1_sqrt <- function(z, par,L){
  
  n=length(z)
  alpha <- par[1]
  delta<- par[2]
  theta <- par[3]
  
  m=1;q=1
  
  eta <- mu <- v <- c()
  eta[1] = delta + theta*(  sqrt(mean(z)) - sqrt(1)    )
  mu[1] = eta[1]^2

  for(t in (m+1):n){
    eta[t] <- delta + theta*( sqrt(z[t-1]) - sqrt(mu[t-1])    )
    mu[t] <- eta[t]^2
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]

  
  lv = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(zn) - log(mun)   )  -
    log(zn) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (zn/mun) ), alpha - L) )
  
  return(sum(lv))
} 
lk.K.MA2_sqrt <- function(z, par,L){
  
  n=length(z)
  alpha = par[1]
  delta<- par[2]
  theta <- par[3:4]
  z0=mean(z);mu0=1
  v=c()
  m=2;q=2
  eta <- mu <- v <- c()
  eta[1] <- delta + theta[1]*( sqrt(z0) - sqrt(mu0)  ) + theta[2]*( sqrt(z0) - sqrt(mu0)  ) 
  mu[1]  <- eta[1]^2
  eta[2] <- delta + theta[1]*( sqrt(z[1]) - sqrt(mu[1])  ) + theta[2]*( sqrt(z0) - sqrt(mu0)  ) 
  mu[2]   <- eta[2]^2
  
  for(t in (m+1):n){
    eta[t] <- delta + theta[1]*( sqrt(z[t-1]) - sqrt(mu[t-1])  ) + theta[2]*( sqrt(z[t-2]) - sqrt(mu[t-2])  ) 
    mu[t]  <- eta[t]^2
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]

  
  lv = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(zn) - log(mun)   )  -
    log(zn) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (zn/mun) ), alpha - L) )
  
  return(sum(lv))
} 
lk.K.ARMA11_sqrt <- function(z, par,L){
  
  n=length(z)
  alpha <- par[1]
  delta <- par[2]
  phi <- par[3]
  theta <- par[4]
  
  m=1;q=1;p=1
  
  eta <- mu <- v <- c()
  eta[1] = delta + phi*sqrt(mean(z)) + theta*(  sqrt(mean(z)) - sqrt(1)   ) 
  mu[1] = eta[1]^2
  
  for(t in (m+1):n){
    eta[t] <- delta + phi*sqrt(z[t-p]) + theta*(  sqrt(z[t-1]) - sqrt(mu[t-1])   )
    mu[t]   <- eta[t]^2
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  
  lv = log(2) + ( (L +alpha)/2  )*(  log(L) + log(alpha) + log(zn) - log(mun)   )  -
    log(zn) - lgamma(L) - lgamma(alpha) + 
    log( besselK( 2 * sqrt( L * alpha * (zn/mun) ), alpha - L) )
  
  return(sum(lv))
} 

### chute inicial para K ----
KChute.AR1_log <- function(y,L){
  
  
  t<-2:length(y)
  X = cbind( rep(1,length(y)), c(1,log(y[t-1]))  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(y)
  if( param[2] < -0.9 ){
    param[2] = -0.9  
  }else if( param[2] > 0.9 ){
    param[2] = 0.9
  }
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.AR1_log( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], phi = param[2]  ) )
  
}
KChute.AR2_log <- function(y,L){
  
  t<-3:length(y)
  X = cbind( rep(1,length(y)), c(1,log(y[1]),log(y[t-1])), c(1,1,log(y[t-2]))   ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(y)
  
  if( param[2] < -0.9 ){
    param[2] = -0.9  
  }else if( param[2] > 0.9 ){
    param[2] = 0.9
  }
  if( param[3] < -0.9 ){
    param[3] = -0.9  
  }else if( param[3] > 0.9 ){
    param[3] = 0.9
  }
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ), rep(  param[2],length(alpha)  ), rep(  param[3],length(alpha)  ) )
  
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.AR2_log( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], phi = param[2:3]  ) )
  
}
KChute.MA1_log <- function(y,L){
  
  param = c( mean(log(y)), 0  )
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.MA1_log( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], theta = param[2]  ) )
  
}
KChute.MA2_log <- function(y,L){
  
  param = c( mean( log(y) ), 0, 0  )
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2], length(alpha)  ),rep(  param[3],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.MA2_log( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], theta = param[2:3]  ) )
  
}
KChute.ARMA11_log <- function(y,L){
  
  t<-2:length(y)
  X = cbind( rep(1,length(y)), c(1, log(y[t-1]))  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(y)
  if( param[2] < -0.8 ){
    param[2] = -0.8  
  }else if( param[2] > 0.8 ){
    param[2] = 0.8
  }
  
  
  
  param = c( param, 0  )
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ),rep(  param[3],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.ARMA11_log( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], phi = param[2], theta = param[3]  ) )
  
}

KChute.AR1_W <- function(y,L){
  
  
  t <- 2:length(y)
  X = cbind( rep(1,length(y)), c( 1, log(lambertW0(y[t-1]))  )  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log( lambertW0(y) )
  if( param[2] < -0.8 ){
    param[2] = -0.8  
  }else if( param[2] > 0.8 ){
    param[2] = 0.8
  }
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.AR1_W( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], phi = param[2]  ) )
  
}
KChute.AR2_W <- function(y,L){
  
  t<-3:length(y)
  X = cbind( rep(1,length(y)), c(1,log(y[1]),log(lambertW0(y[t-1]))  ), c(1,1,log(lambertW0(y[t-2])) )   ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(lambertW0(y))
  if( param[2] < -0.8 ){
    param[2] = -0.8  
  }else if( param[2] > 0.8 ){
    param[2] = 0.8
  }
  if( param[2] < -0.8 ){
    param[2] = -0.8  
  }else if( param[2] > 0.8 ){
    param[2] = 0.8
  }
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ), rep(  param[2],length(alpha)  ), rep(  param[3],length(alpha)  ) )
  
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.AR2_W( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], phi = param[2:3]  ) )
  
}
KChute.MA1_W <- function(y,L){
  
  param = c( mean(log(lambertW0(y))), 0  )
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.MA1_W( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], theta = param[2]  ) )
  
}
KChute.MA2_W <- function(y,L){
  
  param = c( mean( log(lambertW0(y)) ), 0, 0  )
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2], length(alpha)  ),rep(  param[3],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.MA2_W( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], theta = param[2:3]  ) )
  
}
KChute.ARMA11_W <- function(y,L){
  
  t<-2:length(y)
  X = cbind( rep(1,length(y)), c(1, log(lambertW0(y[t-1])) )  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(lambertW0(y))
  if( param[2] < -0.8 ){
    param[2] = -0.8  
  }else if( param[2] > 0.8 ){
    param[2] = 0.8
  }
  
  param = c( param, 0  )
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ),rep(  param[3],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.ARMA11_W( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], phi = param[2], theta = param[3]  ) )
  
}

KChute.AR1_sqrt <- function(y,L){
  
  
  t<-2:length(y)
  X = cbind( rep(1,length(y)), c( 1, sqrt(y[t-1]) )  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%sqrt(y)
  if( param[2] < -0.8 ){
    param[2] = -0.8  
  }else if( param[2] > 0.8 ){
    param[2] = 0.8
  }
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ), rep(  param[2],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.AR1_sqrt( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], phi = param[2]  ) )
  
}
KChute.AR2_sqrt <- function(y,L){
  
  t<-3:length(y)
  X = cbind( rep(1,length(y)), c(1,sqrt(y[1]),sqrt(y[t-1])), c(1,1,sqrt(y[t-2]))   ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%sqrt(y)
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ), rep(  param[2],length(alpha)  ), rep(  param[3],length(alpha)  ) )
  
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.AR2_sqrt( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], phi = param[2:3]  ) )
  
}
KChute.MA1_sqrt <- function(y,L){
  
  param = c( mean(sqrt(y)), 0  )
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.MA1_sqrt( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], theta = param[2]  ) )
  
}
KChute.MA2_sqrt <- function(y,L){
  
  param = c( mean( sqrt(y) ), 0, 0  )
  
  alpha = seq(1,100,10)
  par = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2], length(alpha)  ),rep(  param[3],length(alpha)  ) )
  
  lk=c()
  for(i in 1:length(alpha)) lk[i] = lk.K.MA2_sqrt( y, par[i,], L  )
  
  a = cbind( alpha, lk)
  chute.a = subset( a, a[,2] == max(a[,2])  )[1]
  return( chute = c(alpha = chute.a, delta=param[1], theta = param[2:3]  ) )
  
}
KChute.ARMA11_sqrt <- function(y,L){
  
  t<-2:length(y)
  X = cbind( rep(1,length(y)), c(1, sqrt(y[t-1]))  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%sqrt(y)
  if( param[2] < -0.8 ){
    param[2] = -0.8  
  }else if( param[2] > 0.8 ){
    param[2] = 0.8
  }
  
  param = c( param, 0  )
  
  alpha = seq(1,100,10)
  theta = seq(-0.7,0.7,0.2)
  par1 = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ), rep(theta[1], length(alpha) ) )
  par2 = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ), rep(theta[2], length(alpha) ) )
  par3 = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ), rep(theta[3], length(alpha) ) )
  par4 = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ), rep(theta[4], length(alpha) ) )
  par5 = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ), rep(theta[5], length(alpha) ) )
  par6 = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ), rep(theta[6], length(alpha) ) )
  par7 = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ), rep(theta[7], length(alpha) ) )
  par8 = cbind( alpha, rep(  param[1],length(alpha)  ),rep(  param[2],length(alpha)  ), rep(theta[8], length(alpha) ) )
  
  par = rbind(par1, par2, par3, par4, par5, par6, par7, par8) 
  
  lk=c()
  for(i in 1:dim(par)[1] ) lk[i] = lk.K.ARMA11_sqrt( y, par[i,], L  )
  
  a = cbind( par, lk)
  chute = subset( a, a[,5] == max(a[,5])  )[1:4]
  return( chute = chute )
  
}


#### CML - K-ARMA ----


cml.K.AR1_log <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  A1 <- matrix(c(-1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(0.01,50000,0.99,0.99)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.AR1_log(z,L) 
  lk.Knew = function(param)  lk.K.AR1_log(z, param,L)
  res1 <- maxLik(lk.Knew, grad = NULL, hess = NULL, param, method = "BFGS" ) #, constraints=list(ineqA = A1, ineqB = B1)
  code = res1$code
  est = res1$estimate
  alpha = est[1]; delta=est[2];phi=est[3]
  
  mu.hat[1] = exp(  delta + phi*log( z0 ) )
  t <- 2:n 
  mu.hat[t] = exp( delta + phi*log( z[t-1] ) )
  
  AIC = -2*lk.K.AR1_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.AR1_log(z, est,L) + 2*log(n)
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )
  
}
cml.K.AR2_log <- function(z,L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(-0.1,50000,0.95,0.95,0.95,0.95)
  
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.AR2_log(z,L) 
  lk.Knew = function(param)  lk.K.AR2_log(z, param,L)
  #grad.Knew = function(param)  grad.K.AR2_log(z, param)
  res1 <- maxLik(lk.Knew, grad = NULL, hess = NULL, param, method = "BFGS" ) #, constraints=list(ineqA = A1, ineqB = B1)
  est = res1$estimate
  code = res1$code
  alpha = est[1];delta=est[2];phi=est[3:4]
  
  mu.hat[1] = exp( delta + phi[1]*( log(z0) ) + phi[2]*( log(z0) ) )
  mu.hat[2] = exp( delta + phi[1]*( log(z[1]) ) + phi[2]*( log(z0) ) )
  
  t <- 3:n 
  mu.hat[t] = exp( delta + phi[1]*log( z[t-1] ) + phi[2]*log( z[t-2] )  )
  
  AIC = -2*lk.K.AR2_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.AR2_log(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )
}
cml.K.MA1_log <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(-0.1,50000,0.95,0.95)
  
  #z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.MA1_log(z,L) 
  lk.new = function(param)  lk.K.MA1_log(z, param,L)
  #grad.new = function(param)  grad.K.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS" ) # , constraints=list(ineqA = A1, ineqB = B1)
  est = res1$estimate
  alpha = est[1];delta=est[2];theta=est[3]
  code = res1$code
  z0 = mean(z); mu0=z0
  
  mu.hat[1] = exp( delta + theta*( log(z0) - log(mu0) ) )
  
  for(t in 2:n) mu.hat[t] = exp( delta + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) )
  
  AIC = -2*lk.K.MA1_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.MA1_log(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}
cml.K.MA2_log <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(-0.1,50000,0.95,0.95,0.95,0.95)
  
  require(purrr)
  
  #A1 <- matrix(c(1, -1, 0,0,0,0,0,0 ), 2,4)
  #B1 <- c(-1,200)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.MA2_log(z,L) 
  lk.new = function(param)  lk.K.MA2_log(z, param,L)
  #grad.new = function(param)  grad.K.MA(z, param,2)
  maxlik.new = function(logLik,grad,hess,start,method)  maxLik(logLik = lk.new,grad=NULL,hess = NULL, start=param, method = "BFGS") #, constraints=list(ineqA = A1, ineqB = B1)
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS" ) # , constraints=list(ineqA = A1, ineqB = B1)
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "CG", constraints=list(ineqA = A1, ineqB = B1) )
  }else{
    res1 <- res1$result
  }
  
  if(res1$code == 100){
    est = rep(NA,length(param)); mu.hat = rep(NA,length(z))
  }else{
    
    est = res1$estimate
    alpha = est[1];delta=est[2];theta=est[3:4]
    code = res1$code
    
    z0 = mean(z); mu0=z0
    mu.hat[1:2] = exp( delta + theta[1]*( log(z0) - log(mu0) ) + theta[2]*( log(z0) - log(mu0) ) )
    
    for(t in 3:n) mu.hat[t] = exp( delta + theta[1]*( log(z[t-1]) - log(mu.hat[t-1]) ) + theta[2]*( log(z[t-2]) - log(mu.hat[t-2]) ) )
    
  }
  
  AIC = -2*lk.K.MA2_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.MA2_log(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}
cml.K.ARMA11_log <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, 0, 0,0 , -1,0,0,0, 0,0,1,0, 0,0,-1,0, 0,0,0,1,0,0,0,-1 ), 6,4, byrow = T)
  B1 <- c(-0.5,50000, 0.95, 0.95, 0.95, 0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.ARMA11_log(z,L)
  #param[1] = min(param[1], 49.999)
  lk.new = function(param)  lk.K.ARMA11_log(z, param,L)
  #grad.new = function(param)  grad.K.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS" ) # , constraints=list(ineqA = A1, ineqB = B1)
  est = res1$estimate
  alpha = est[1];delta=est[2];phi=est[3];theta=est[4]
  code = res1$code
  z0 = mean(z); mu0=z0
  
  mu.hat[1] = exp( delta + phi*log(z0) + theta*( log(z0) - log(mu0) ) )
  
  for(t in 2:n) mu.hat[t] = exp( delta + phi*log(z[t-1]) + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) )
  
  AIC = -2*lk.K.ARMA11_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.ARMA11_log(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}

auto.K_log <- function(y, criterio,L){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  
  ar1 = cml.K.AR1_log(y,L);ar2 = cml.K.AR2_log(y,L)
  ma1 = cml.K.MA1_log(y,L); ma2 = cml.K.MA2_log(y,L)
  arma11 = cml.K.ARMA11_log(y,L)
  
  if(criterio == "AIC"){
    
    criter.ar1 =  data.frame(AIC=ar1$AIC,model="AR1")
    criter.ar2 =  data.frame(AIC=ar2$AIC,model="AR2")
    criter.ma1 =  data.frame(AIC=ma1$AIC,model="MA1")
    criter.ma2 =  data.frame(AIC=ma2$AIC,model="MA2")
    criter.arma11 =  data.frame(AIC=arma11$AIC,model="ARMA11")
    
  }else if(criterio == "BIC"){
    
    criter.ar1 =  data.frame(BIC=ar1$BIC,model="AR1")
    criter.ar2 =  data.frame(BIC=ar2$BIC,model="AR2")
    criter.ma1 =  data.frame(BIC=ma1$BIC,model="MA1")
    criter.ma2 =  data.frame(BIC=ma2$BIC,model="MA2")
    criter.arma11 =  data.frame(BIC=arma11$BIC,model="ARMA11")     
    
  }
  
  CRITERIOS = rbind(criter.ar1, criter.ar2, criter.ma1, criter.ma2, criter.arma11)
  
  best.model = subset(CRITERIOS, CRITERIOS[,1] == min(CRITERIOS[,1], na.rm = T))
  
  if(best.model$model == "AR1"){
    est = ar1$est; mu.hat = ar1$mu.hat
  }else if( best.model$model == "AR2" ){
    est = ar2$est; mu.hat = ar2$mu.hat
  }else if( best.model$model == "MA1" ){
    est = ma1$est; mu.hat = ma1$mu.hat
  }else if( best.model$model == "MA2" ){
    est = ma2$est; mu.hat = ma2$mu.hat
  }else if( best.model$model == "ARMA11" ){
    est = arma11$est; mu.hat = arma11$mu.hat
  }
  
  return(list(model = best.model$model,criterio=criterio, est=est ,mu.hat=mu.hat) )  
}


cml.K.AR1_W <- function(z,L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(-0.1,50000,0.95,0.95)
  
  
  require(lamW)
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.AR1_W(z,L) 
  lk.Knew = function(param)  lk.K.AR1_W(z, param,L)
  res1 <- maxLik(lk.Knew, grad = NULL, hess = NULL, param, method = "BFGS") # , constraints=list(ineqA = A1, ineqB = B1) 
  code = res1$code
  est = res1$estimate
  delta=est[2];phi=est[3]
  eta=c()
  
  eta[1] = delta + phi*( log(lambertW0(mean(z))) ) 
  mu.hat[1] = exp( eta[1] + exp(eta[1])  )  
  
  t <- 2:n 
  eta[t] = delta + phi*( log(lambertW0(z[t-1])) ) 
  mu.hat[t] = exp( eta[t] + exp(eta[t])  )
  
  AIC = -2*lk.K.AR1_W(z, est, L) + 2*length(est) 
  BIC = -2*lk.K.AR1_W(z, est, L) + 2*log(n)
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )
  
}
cml.K.AR2_W <- function(z,L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  require(lamW)
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(-0.1,50000,0.95,0.95,0.95,0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.AR2_W(z,L) 
  lk.Knew = function(param)  lk.K.AR2_W(z, param,L)
  #grad.Knew = function(param)  grad.K.AR2_W(z, param)
  res1 <- maxLik(lk.Knew, grad = NULL, hess = NULL, param, method = "BFGS" ) # , constraints=list(ineqA = A1, ineqB = B1)
  est = res1$estimate
  delta=est[2];phi=est[3:4]
  code = res1$code
  
  eta=c()
  
  eta[1] = delta + phi[1]*( log(lambertW0(mean(z))) ) + phi[2]*( log(lambertW0(mean(z))) ) 
  eta[2] = delta + phi[1]*( log(lambertW0(z[1])) ) + phi[2]*( log(lambertW0(mean(z))) )
  mu.hat[1:2] = exp( eta[1:2] + exp(eta[1:2])  )  
  
  t <- 3:n 
  eta[t] = delta + phi[1]*( log(lambertW0(z[t-1])) ) + phi[2]*( log(lambertW0(z[t-2])) )
  mu.hat[t] = exp( eta[t] + exp(eta[t])  )  
  
  AIC = -2*lk.K.AR2_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.AR2_W(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )
}
cml.K.MA1_W <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(-0.1,50000,0.95,0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  require(lamW)
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.MA1_W(z,L) 
  lk.new = function(param)  lk.K.MA1_W(z, param,L)
  #grad.new = function(param)  grad.K.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS" ) # , constraints=list(ineqA = A1, ineqB = B1)
  est = res1$estimate
  delta=est[2];theta=est[3]
  code = res1$code
  z0 = mean(z); mu0=z0
  eta=c()
  
  eta[1] = delta + theta*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  ) 
  mu.hat[1] = exp( eta[1] + exp(eta[1])  )  
  
  for(t in 2:n){
    eta[t] = delta + theta*( log(lambertW0(mean(z[t-1]))) - log(lambertW0(mu.hat[t-1]))  ) 
    mu.hat[t] = exp( eta[t] + exp(eta[t])  ) 
  }
  
  AIC = -2*lk.K.MA1_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.MA1_W(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}
cml.K.MA2_W <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(purrr)
  require(lamW)
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0), 2,4)
  B1 <- c(50000,-1.3)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.MA2_W(z,L) 
  lk.new = function(param)  lk.K.MA2_W(z, param,L)
  #grad.new = function(param)  grad.K.MA(z, param,2)
  maxlik.new = function(logLik,grad,hess,start,method)  maxLik(logLik = lk.new,grad=NULL,hess = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = NULL, start=param, method = "CG")
  }else{
    res1 <- res1$result
  }
  
  if(res1$code == 100){
    est = rep(NA,length(param)); mu.hat = rep(NA,length(z))
  }else{
    
    
    est = res1$estimate
    code = res1$code
    
    delta=est[2];theta=est[3:4]
    z0 = mean(z); mu0=z0
    eta=c()
    
    eta[1] = delta + theta[1]*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  ) + theta[2]*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  )
    mu.hat[1] = exp( eta[1] + exp(eta[1])  ) 
    eta[2] = delta + theta[1]*( log(lambertW0(z[1])) - log(lambertW0(mu.hat[1]))  ) + theta[2]*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  )
    mu.hat[2] = exp( eta[2] + exp(eta[2])  ) 
    
    
    for(t in 3:n){
      eta[t] = delta + theta[1]*( log(lambertW0(z[t-1])) - log(lambertW0(mu.hat[t-1]))  ) + theta[2]*( log(lambertW0(z[t-2])) - log(lambertW0(mu.hat[t-2]))  )
      mu.hat[t] = exp( eta[t] + exp(eta[t])  ) 
    } 
  }
  
  AIC = -2*lk.K.MA2_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.MA2_W(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC,code = code) )  
}
cml.K.ARMA11_W <- function(z,L ){
  require(purrr)
  require(lamW)
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, 0, 0,0 , -1,0,0,0, 0,0,1,0, 0,0,-1,0, 0,0,0,1,0,0,0,-1 ), 6,4, byrow = T)
  B1 <- c(-0.5,50000, 0.95, 0.95, 0.95, 0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();
  n = length(z)
  
  param = KChute.ARMA11_W(z,L)
  lk.new = function(param)  lk.K.ARMA11_W(z, param,L)
  #grad.new = function(param)  grad.K.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS" ) # , constraints=list(ineqA = A1, ineqB = B1)
  est = res1$estimate
  alpha = est[1];delta=est[2];phi=est[3];theta=est[4]
  code = res1$code
  z0 = mean(z); mu0=z0
  eta=c()
  
  eta[1] = delta + phi*( log(lambertW0(z0))  ) + theta*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  )
  mu.hat[1] = exp( eta[1] + exp(eta[1])  ) 
  
  for(t in 2:n){
    eta[t] = delta + phi*( log(lambertW0(z[t-1]))  ) + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu.hat[t-1]))  )
    mu.hat[t] = exp( eta[t] + exp(eta[t])  ) 
  }
  
  AIC = -2*lk.K.ARMA11_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.ARMA11_W(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}

auto.K_W <- function(y, criterio,L){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  require(lamW)
  ar1 = cml.K.AR1_W(y,L);ar2 = cml.K.AR2_W(y,L)
  ma1 = cml.K.MA1_W(y,L); ma2 = cml.K.MA2_W(y,L)
  arma11 = cml.K.ARMA11_W(y,L)
  
  if(criterio == "AIC"){
    
    criter.ar1 =  data.frame(AIC=ar1$AIC,model="AR1")
    criter.ar2 =  data.frame(AIC=ar2$AIC,model="AR2")
    criter.ma1 =  data.frame(AIC=ma1$AIC,model="MA1")
    criter.ma2 =  data.frame(AIC=ma2$AIC,model="MA2")
    criter.arma11 =  data.frame(AIC=arma11$AIC,model="ARMA11")
    
  }else if(criterio == "BIC"){
    
    criter.ar1 =  data.frame(BIC=ar1$BIC,model="AR1")
    criter.ar2 =  data.frame(BIC=ar2$BIC,model="AR2")
    criter.ma1 =  data.frame(BIC=ma1$BIC,model="MA1")
    criter.ma2 =  data.frame(BIC=ma2$BIC,model="MA2")
    criter.arma11 =  data.frame(BIC=arma11$BIC,model="ARMA11")     
    
  }
  
  CRITERIOS = rbind(criter.ar1, criter.ar2, criter.ma1, criter.ma2, criter.arma11)
  
  best.model = subset(CRITERIOS, CRITERIOS[,1] == min(CRITERIOS[,1], na.rm = T))
  
  if(best.model$model == "AR1"){
    est = ar1$est; mu.hat = ar1$mu.hat
  }else if( best.model$model == "AR2" ){
    est = ar2$est; mu.hat = ar2$mu.hat
  }else if( best.model$model == "MA1" ){
    est = ma1$est; mu.hat = ma1$mu.hat
  }else if( best.model$model == "MA2" ){
    est = ma2$est; mu.hat = ma2$mu.hat
  }else if( best.model$model == "ARMA11" ){
    est = arma11$est; mu.hat = arma11$mu.hat
  }
  
  return(list(model = best.model$model,criterio=criterio, est=est ,mu.hat=mu.hat) )  
}


cml.K.AR1_sqrt <- function(z, L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(-0.1,50000,0.95,0.95)
  
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.AR1_sqrt(z,L) 
  lk.Knew = function(param)  lk.K.AR1_sqrt(z, param,L)
  res1 <- maxLik(lk.Knew, grad = NULL, hess = NULL, param, method = "BFGS" ) # , constraints=list(ineqA = A1, ineqB = B1)
  code = res1$code
  est = res1$estimate
  alpha = est[1]; delta=est[2];phi=est[3]
  
  eta=c()
  eta[1] = delta + phi*sqrt( z0 )
  mu.hat[1] = eta[1]^2
  t <- 2:n 
  eta[t] = delta + phi*sqrt( z[t-1] ) 
  mu.hat[t] = eta[t]^2
  
  AIC = -2*lk.K.AR1_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.AR1_sqrt(z, est,L) + 2*log(n)
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )
  
}
cml.K.AR2_sqrt <- function(z, L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(-0.1,50000,0.95,0.95,0.95,0.95)
  
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.AR2_sqrt(z,L) 
  lk.Knew = function(param)  lk.K.AR2_sqrt(z, param,L)
  #grad.Knew = function(param)  grad.K.AR2_sqrt(z, param)
  res1 <- maxLik(lk.Knew, grad = NULL, hess = NULL, param, method = "BFGS" ) # , constraints=list(ineqA = A1, ineqB = B1)
  est = res1$estimate
  code = res1$code
  
  alpha = est[1];delta=est[2];phi=est[3:4]
  eta=c()
  eta[1] = delta + phi[1]*( sqrt(z0) ) + phi[2]*( sqrt(z0) )
  eta[2] = delta + phi[1]*( sqrt(z[1]) ) + phi[2]*( sqrt(z0) )
  mu.hat[1] = eta[1]^2
  mu.hat[2] = eta[2]^2
  
  t <- 3:n 
  eta[t] =  delta + phi[1]*sqrt( z[t-1] ) + phi[2]*sqrt( z[t-2] )
  mu.hat[t] = eta[t]^2
  
  AIC = -2*lk.K.AR2_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.AR2_sqrt(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC,code=code) )
}
cml.K.MA1_sqrt <- function(z, L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(-0.1,50000,0.95,0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.MA1_sqrt(z,L) 
  lk.new = function(param)  lk.K.MA1_sqrt(z, param,L)
  #grad.new = function(param)  grad.K.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS" ) # , constraints=list(ineqA = A1, ineqB = B1)
  est = res1$estimate
  code = res1$code
  alpha = est[1];delta=est[2];theta=est[3]
  
  z0 = mean(z); mu0=z0
  eta=c()
  eta[1] = delta + theta*( sqrt(z0) - sqrt(mu0) )
  mu.hat[1] =eta[1]^2
  
  for(t in 2:n){
    eta[t] = delta + theta*( sqrt(z[t-1]) - sqrt(mu.hat[t-1]) )
    mu.hat[t] = eta[t]^2
  }
  
  AIC = -2*lk.K.MA1_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.MA1_sqrt(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC,code = code) )  
}
cml.K.MA2_sqrt <- function(z, L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(purrr)
  
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(-0.1,50000,0.95,0.95,0.95,0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = KChute.MA2_sqrt(z,L) 
  lk.new = function(param)  lk.K.MA2_sqrt(z, param,L)
  #grad.new = function(param)  grad.K.MA(z, param,2)
  maxlik.new = function(logLik,grad,hess,start,method)  maxLik(logLik = lk.new,grad=NULL,hess = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "CG")
  }else{
    res1 <- res1$result
  }
  
  if(res1$code == 100){
    est = rep(NA,length(param)); mu.hat = rep(NA,length(z))
  }else{
    
    est = res1$estimate
    alpha = est[1];delta=est[2];theta=est[3:4]
    code = res1$code
    
    z0 = mean(z); mu0=z0
    eta=c()
    eta[1:2] = delta + theta[1]*( sqrt(z0) - sqrt(mu0) ) + theta[2]*( sqrt(z0) - sqrt(mu0) )
    mu.hat[1:2] = eta[1:2]^2
    
    for(t in 3:n){
      eta[t] = delta + theta[1]*( sqrt(z[t-1]) - sqrt(mu.hat[t-1]) ) + theta[2]*( sqrt(z[t-2]) - sqrt(mu.hat[t-2]) ) 
      mu.hat[t] = eta[t]^2
      
    } 
  }
  
  AIC = -2*lk.K.MA2_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.MA2_sqrt(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC,code = code) )  
}
cml.K.ARMA11_sqrt <- function(z, L ){
  require(purrr)
  
  A1 <- matrix(c(1, 0, 0,0 , -1,0,0,0, 0,0,1,0, 0,0,-1,0, 0,0,0,1,0,0,0,-1 ), 6,4, byrow = T)
  B1 <- c(-0.5,50000, 0.95, 0.95, 0.95, 0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();
  n = length(z)
  
  param = KChute.ARMA11_sqrt(z,L)
  lk.new = function(param)  lk.K.ARMA11_sqrt(z, param,L)
  #grad.new = function(param)  grad.K.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS") # , constraints=list(ineqA = A1, ineqB = B1) 
  est = res1$estimate
  alpha = est[1];delta=est[2];phi=est[3];theta=est[4]
  code = res1$code
  z0 = mean(z); mu0=z0
  
  eta=c()
  eta[1] = delta + phi*sqrt( z0 ) + theta*( sqrt(z0) - sqrt(mu0) )
  mu.hat[1] = eta[1]^2
  
  for(t in 2:n){
    eta[t] = delta + phi*sqrt( z[t-1] ) + theta*( sqrt(z[t-1]) - sqrt(mu.hat[t-1]) )
    mu.hat[t] =eta[t]^2
  }
  
  AIC = -2*lk.K.ARMA11_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.K.ARMA11_sqrt(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}


auto.K_sqrt <- function(y, criterio, L){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  
  ar1 = cml.K.AR1_sqrt(y,L);ar2 = cml.K.AR2_sqrt(y,L)
  ma1 = cml.K.MA1_sqrt(y,L); ma2 = cml.K.MA2_sqrt(y,L)
  arma11 = cml.K.ARMA11_sqrt(y,L)
  
  if(criterio == "AIC"){
    
    criter.ar1 =  data.frame(AIC=ar1$AIC,model="AR1")
    criter.ar2 =  data.frame(AIC=ar2$AIC,model="AR2")
    criter.ma1 =  data.frame(AIC=ma1$AIC,model="MA1")
    criter.ma2 =  data.frame(AIC=ma2$AIC,model="MA2")
    criter.arma11 =  data.frame(AIC=arma11$AIC,model="ARMA11")
    
  }else if(criterio == "BIC"){
    
    criter.ar1 =  data.frame(BIC=ar1$BIC,model="AR1")
    criter.ar2 =  data.frame(BIC=ar2$BIC,model="AR2")
    criter.ma1 =  data.frame(BIC=ma1$BIC,model="MA1")
    criter.ma2 =  data.frame(BIC=ma2$BIC,model="MA2")
    criter.arma11 =  data.frame(BIC=arma11$BIC,model="ARMA11")     
    
  }
  
  CRITERIOS = rbind(criter.ar1, criter.ar2, criter.ma1, criter.ma2, criter.arma11)
  
  best.model = subset(CRITERIOS, CRITERIOS[,1] == min(CRITERIOS[,1], na.rm = T))
  
  if(best.model$model == "AR1"){
    est = ar1$est; mu.hat = ar1$mu.hat
  }else if( best.model$model == "AR2" ){
    est = ar2$est; mu.hat = ar2$mu.hat
  }else if( best.model$model == "MA1" ){
    est = ma1$est; mu.hat = ma1$mu.hat
  }else if( best.model$model == "MA2" ){
    est = ma2$est; mu.hat = ma2$mu.hat
  }else if( best.model$model == "ARMA11" ){
    est = arma11$est; mu.hat = arma11$mu.hat
  }
  
  return(list(model = best.model$model,criterio=criterio, est=est ,mu.hat=mu.hat) )  
}


#### Erro padrao boot ----

EP.K.AR1_log <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_AR1_log(N, par, L)
      ajuste = cml.K.AR1_log( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.AR2_log <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_AR2_log(N, par, L)
      ajuste = cml.K.AR2_log( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.MA1_log <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_MA1_log(N, par, L)
      ajuste = cml.K.MA1_log( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.MA2_log <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_MA2_log(N, par, L)
      ajuste = cml.K.MA2_log( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.ARMA11_log <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_ARMA11_log(N, par, L)
      ajuste = cml.K.ARMA11_log( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}

EP.K.AR1_W <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
 #   error = 1;
#    while( error == 1 ){  
      y <- rK_AR1_W(N, par, L)
      ajuste = cml.K.AR1_W( y, L  )
 #     error = ifelse( ajuste$code == 0, 0, 1  )
#    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.AR2_W <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_AR2_W(N, par, L)
      ajuste = cml.K.AR2_W( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.MA1_W <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_MA1_W(N, par, L)
      ajuste = cml.K.MA1_W( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.MA2_W <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_MA2_W(N, par, L)
      ajuste = cml.K.MA2_W( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.ARMA11_W <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_ARMA11_W(N, par, L)
      ajuste = cml.K.ARMA11_W( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}

EP.K.AR1_sqrt <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    #error = 1;
    #while( error == 1 ){  
      y <- rK_AR1_sqrt(N, par, L)
      ajuste = cml.K.AR1_sqrt( y, L  )
     # error = ifelse( ajuste$code == 0, 0, 1  )
    #}
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.AR2_sqrt <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_AR2_sqrt(N, par, L)
      ajuste = cml.K.AR2_sqrt( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.MA1_sqrt <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_MA1_sqrt(N, par, L)
      ajuste = cml.K.MA1_sqrt( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.MA2_sqrt <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
    while( error == 1 ){  
      y <- rK_MA2_sqrt(N, par, L)
      ajuste = cml.K.MA2_sqrt( y, L  )
      error = ifelse( ajuste$code == 0, 0, 1  )
    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.K.ARMA11_sqrt <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
   # error = 1;
  #  while( error == 1 ){  
      y <- rK_ARMA11_sqrt(N, par, L)
      ajuste = cml.K.ARMA11_sqrt( y, L  )
   #   error = ifelse( ajuste$code == 0, 0, 1  )
  #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}


EP.gamma.AR1_log <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
   # error = 1;
  #  while( error == 1 ){  
      y <- rgama_AR1_log(N, par, L)
      ajuste = cml.gamma.AR1_log( y, L  )
   #   error = ifelse( ajuste$code == 0, 0, 1  )
  #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.AR2_log <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
  #  error = 1;
  #  while( error == 1 ){  
      y <- rgama_AR2_log(N, par, L)
      ajuste = cml.gamma.AR2_log( y, L  )
   #   error = ifelse( ajuste$code == 0, 0, 1  )
  #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.MA1_log <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
  #  error = 1;
  #  while( error == 1 ){  
      y <- rgama_MA1_log(N, par, L)
      ajuste = cml.gamma.MA1_log( y, L  )
   #   error = ifelse( ajuste$code == 0, 0, 1  )
  #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.MA2_log <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
   # error = 1;
  #  while( error == 1 ){  
      y <- rgama_MA2_log(N, par, L)
      ajuste = cml.gamma.MA2_log( y, L  )
   #   error = ifelse( ajuste$code == 0, 0, 1  )
  #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.ARMA11_log <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
  #  error = 1;
  #  while( error == 1 ){  
      y <- rgama_ARMA11_log(N, par, L)
      ajuste = cml.gamma.ARMA11_log( y, L  )
   #   error = ifelse( ajuste$code == 0, 0, 1  )
  #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}

EP.gamma.AR1_W <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
 #   error = 1;
  #  while( error == 1 ){  
      y <- rgama_AR1_W(N, par, L)
      ajuste = cml.gamma.AR1_W( y, L  )
   #   error = ifelse( ajuste$code == 0, 0, 1  )
  #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.AR2_W <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
  #  error = 1;
   # while( error == 1 ){  
      y <- rgama_AR2_W(N, par, L)
      ajuste = cml.gamma.AR2_W( y, L  )
    #  error = ifelse( ajuste$code == 0, 0, 1  )
    #}
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.MA1_W <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
  #  error = 1;
  #  while( error == 1 ){  
      y <- rgama_MA1_W(N, par, L)
      ajuste = cml.gamma.MA1_W( y, L  )
   #   error = ifelse( ajuste$code == 0, 0, 1  )
  #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.MA2_W <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
  #  error = 1;
  #  while( error == 1 ){  
      y <- rgama_MA2_W(N, par, L)
      ajuste = cml.gamma.MA2_W( y, L  )
   #   error = ifelse( ajuste$code == 0, 0, 1  )
  #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.ARMA11_W <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    error = 1;
   # while( error == 1 ){  
      y <- rgama_ARMA11_W(N, par, L)
      ajuste = cml.gamma.ARMA11_W( y, L  )
    #  error = ifelse( ajuste$code == 0, 0, 1  )
    #}
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}

EP.gamma.AR1_sqrt <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
  #  error = 1;
  #  while( error == 1 ){  
      y <- rgama_AR1_sqrt(N, par, L)
      ajuste = cml.gamma.AR1_sqrt( y, L  )
   #   error = ifelse( ajuste$code == 0, 0, 1  )
  #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.AR2_sqrt <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
 #   error = 1;
  #  while( error == 1 ){  
      y <- rgama_AR2_sqrt(N, par, L)
      ajuste = cml.gamma.AR2_sqrt( y, L  )
   #   error = ifelse( ajuste$code == 0, 0, 1  )
  #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.MA1_sqrt <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
 #   error = 1;
#    while( error == 1 ){  
      y <- rgama_MA1_sqrt(N, par, L)
      ajuste = cml.gamma.MA1_sqrt( y, L  )
 #     error = ifelse( ajuste$code == 0, 0, 1  )
#    }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.MA2_sqrt <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
   # error = 1;
    #while( error == 1 ){  
      y <- rgama_MA2_sqrt(N, par, L)
      ajuste = cml.gamma.MA2_sqrt( y, L  )
     # error = ifelse( ajuste$code == 0, 0, 1  )
    #}
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}
EP.gamma.ARMA11_sqrt <- function(N, B, par, L){
  est = matrix(NA, B, length(par))
  code = c()
  
  for( i in 1:B ){
    # error = 1;
    #  while( error == 1 ){  
    y <- rgama_ARMA11_sqrt(N, par, L)
    ajuste = cml.gamma.ARMA11_sqrt( y, L  )
    #   error = ifelse( ajuste$code == 0, 0, 1  )
    #  }
    est[i,] = ajuste$est
  }
  
  Sx = apply(est , 2, sd)/sqrt(N)
  MEAN = apply(est , 2, mean)
  return( list(erro.padrao = round(Sx, 4), mean = MEAN) )
  
}



filterK_log <- function(L, Matrix.ST1,Matrix.ST2,Matrix.ST3,Matrix.ST4, cores, ...) {
  
  
  M = dim(Matrix.ST1)[1]
  iter <- function(i){
    
    x1 = Matrix.ST1[i,]
    x2 = Matrix.ST2[i,]
    x3 = Matrix.ST3[i,]
    x4 = Matrix.ST4[i,]
    
    
    K.filter_left = cml.K.ARMA11_log(x1, L) 
    K.filter_right = cml.K.ARMA11_log(x2, L) 
    K.filter_under = cml.K.ARMA11_log(x3, L) 
    K.filter_up = cml.K.ARMA11_log(x4, L) 
    
    pixel.est_left = K.filter_left$mu.hat[ length(x1)/2 + 0.5]
    pixel.est_right = K.filter_right$mu.hat[ length(x2)/2 + 0.5]
    pixel.est_under = K.filter_under$mu.hat[ length(x3)/2 + 0.5]
    pixel.est_up = K.filter_up$mu.hat[ length(x4)/2 + 0.5]
    pixel.est.mean = mean ( c(  pixel.est_left, pixel.est_right, pixel.est_under, pixel.est_up  )  )
    
    ALPHA_left = K.filter_left$est[1] 
    ALPHA_right = K.filter_right$est[1] 
    ALPHA_under = K.filter_under$est[1] 
    ALPHA_up = K.filter_up$est[1] 
    ALPHA_mean = mean ( c(  ALPHA_left, ALPHA_right, ALPHA_under, ALPHA_up  )  )
    
    DELTA_left = K.filter_left$est[2] 
    DELTA_right = K.filter_right$est[2] 
    DELTA_under = K.filter_under$est[2] 
    DELTA_up = K.filter_up$est[2] 
    DELTA_mean = mean ( c(  DELTA_left, DELTA_right, DELTA_under, DELTA_up  )  )
    
    PHI_left = K.filter_left$est[3] 
    PHI_right = K.filter_right$est[3] 
    PHI_under = K.filter_under$est[3] 
    PHI_up = K.filter_up$est[3] 
    PHI_mean = mean ( c(  PHI_left, PHI_right, PHI_under, PHI_up  )  )
    
    THETA_left = K.filter_left$est[4] 
    THETA_right = K.filter_right$est[4] 
    THETA_under = K.filter_under$est[4] 
    THETA_up = K.filter_up$est[4] 
    THETA_mean = mean ( c(  THETA_left, THETA_right, THETA_under, THETA_up  )  )
    
    res = c(pixel.est_left, pixel.est_right, pixel.est_under, pixel.est_up, pixel.est.mean,
            ALPHA_left, ALPHA_right, ALPHA_under, ALPHA_up, ALPHA_mean,
            DELTA_left, DELTA_right, DELTA_under, DELTA_up, DELTA_mean,
            PHI_left, PHI_right, PHI_under, PHI_up, PHI_mean,
            THETA_left, THETA_right, THETA_under, THETA_up, THETA_mean
    )
    return(res)  
  }
  
  
  result_mc <-
    if (.Platform$OS.type == "windows") {
      lapply( X = 1L:M, FUN = iter )
    } else {
      pbmcapply::pbmclapply( X = 1L:M, FUN = iter, mc.cores = cores )
    } %>%
    unlist %>% 
    matrix( ncol = 25, nrow = M, byrow = TRUE )  
  
  
  
  colnames(result_mc) <- c(
    "Left", "Right", "Under", "Up", "Mean",
    "ALPHA_Left", "ALPHA_Right", "ALPHA_Under", "ALPHA_Up", "ALPHA_Mean",
    "DELTA_Left", "DELTA_Right", "DELTA_Under", "DELTA_Up", "DELTA_Mean",
    "PHI_Left", "PHI_Right", "PHI_Under", "PHI_Up", "PHI_Mean",
    "THETA_Left", "THETA_Right", "THETA_Under", "THETA_Up", "THETA_Mean"
  )
  
  return(list(result_mc = result_mc))
  
}

filterK_sqrt <- function(L, Matrix.ST1,Matrix.ST2,Matrix.ST3,Matrix.ST4, cores, ...) {
  
  
  M = dim(Matrix.ST1)[1]
  res = matrix( ncol = 25, nrow = M, byrow = TRUE )  
  for(i in 1:M){
    
    x1 = Matrix.ST1[i,]
    x2 = Matrix.ST2[i,]
    x3 = Matrix.ST3[i,]
    x4 = Matrix.ST4[i,]
    
    
    K.filter_left = cml.K.ARMA11_sqrt(x1, L) 
    K.filter_right = cml.K.ARMA11_sqrt(x2, L) 
    K.filter_under = cml.K.ARMA11_sqrt(x3, L) 
    K.filter_up = cml.K.ARMA11_sqrt(x4, L) 
    
    pixel.est_left = K.filter_left$mu.hat[ length(x1)/2 + 0.5]
    pixel.est_right = K.filter_right$mu.hat[ length(x2)/2 + 0.5]
    pixel.est_under = K.filter_under$mu.hat[ length(x3)/2 + 0.5]
    pixel.est_up = K.filter_up$mu.hat[ length(x4)/2 + 0.5]
    pixel.est.mean = mean ( c(  pixel.est_left, pixel.est_right, pixel.est_under, pixel.est_up  )  )
    
    ALPHA_left = K.filter_left$est[1] 
    ALPHA_right = K.filter_right$est[1] 
    ALPHA_under = K.filter_under$est[1] 
    ALPHA_up = K.filter_up$est[1] 
    ALPHA_mean = mean ( c(  ALPHA_left, ALPHA_right, ALPHA_under, ALPHA_up  )  )
    
    DELTA_left = K.filter_left$est[2] 
    DELTA_right = K.filter_right$est[2] 
    DELTA_under = K.filter_under$est[2] 
    DELTA_up = K.filter_up$est[2] 
    DELTA_mean = mean ( c(  DELTA_left, DELTA_right, DELTA_under, DELTA_up  )  )
    
    PHI_left = K.filter_left$est[3] 
    PHI_right = K.filter_right$est[3] 
    PHI_under = K.filter_under$est[3] 
    PHI_up = K.filter_up$est[3] 
    PHI_mean = mean ( c(  PHI_left, PHI_right, PHI_under, PHI_up  )  )
    
    THETA_left = K.filter_left$est[4] 
    THETA_right = K.filter_right$est[4] 
    THETA_under = K.filter_under$est[4] 
    THETA_up = K.filter_up$est[4] 
    THETA_mean = mean ( c(  THETA_left, THETA_right, THETA_under, THETA_up  )  )
    
    res[i,] = c(pixel.est_left, pixel.est_right, pixel.est_under, pixel.est_up, pixel.est.mean,
                ALPHA_left, ALPHA_right, ALPHA_under, ALPHA_up, ALPHA_mean,
                DELTA_left, DELTA_right, DELTA_under, DELTA_up, DELTA_mean,
                PHI_left, PHI_right, PHI_under, PHI_up, PHI_mean,
                THETA_left, THETA_right, THETA_under, THETA_up, THETA_mean
    )
    print(i)
    #return(res)  
  }
  
  
  #result_mc <-
  #  pbmclapply( X = 1L:M, FUN = iter, mc.cores = cores ) %>%
  #  unlist %>% 
  #  matrix( ncol = 25, nrow = M, byrow = TRUE )  
  
  
  
  colnames(res) <- c(
    "Left", "Right", "Under", "Up", "Mean",
    "ALPHA_Left", "ALPHA_Right", "ALPHA_Under", "ALPHA_Up", "ALPHA_Mean",
    "DELTA_Left", "DELTA_Right", "DELTA_Under", "DELTA_Up", "DELTA_Mean",
    "PHI_Left", "PHI_Right", "PHI_Under", "PHI_Up", "PHI_Mean",
    "THETA_Left", "THETA_Right", "THETA_Under", "THETA_Up", "THETA_Mean"
  )
  
  return(list(result_mc = result_mc))
  
}

filter_W <- function(L, Matrix.ST1,Matrix.ST2,Matrix.ST3,Matrix.ST4, cores, ...) {
  
  
  M = dim(Matrix.ST1)[1]
  iter <- function(i){
    
    x1 = Matrix.ST1[i,]
    x2 = Matrix.ST2[i,]
    x3 = Matrix.ST3[i,]
    x4 = Matrix.ST4[i,]
    
    
    K.filter_left = cml.K.ARMA11_W(x1, L) 
    K.filter_right = cml.K.ARMA11_W(x2, L) 
    K.filter_under = cml.K.ARMA11_W(x3, L) 
    K.filter_up = cml.K.ARMA11_W(x4, L) 
    
    pixel.est_left = K.filter_left$mu.hat[ length(x1)/2 + 0.5]
    pixel.est_right = K.filter_right$mu.hat[ length(x2)/2 + 0.5]
    pixel.est_under = K.filter_under$mu.hat[ length(x3)/2 + 0.5]
    pixel.est_up = K.filter_up$mu.hat[ length(x4)/2 + 0.5]
    pixel.est.mean = mean ( c(  pixel.est_left, pixel.est_right, pixel.est_under, pixel.est_up  )  )
    
    ALPHA_left = K.filter_left$est[1] 
    ALPHA_right = K.filter_right$est[1] 
    ALPHA_under = K.filter_under$est[1] 
    ALPHA_up = K.filter_up$est[1] 
    ALPHA_mean = mean ( c(  ALPHA_left, ALPHA_right, ALPHA_under, ALPHA_up  )  )
    
    DELTA_left = K.filter_left$est[2] 
    DELTA_right = K.filter_right$est[2] 
    DELTA_under = K.filter_under$est[2] 
    DELTA_up = K.filter_up$est[2] 
    DELTA_mean = mean ( c(  DELTA_left, DELTA_right, DELTA_under, DELTA_up  )  )
    
    PHI_left = K.filter_left$est[3] 
    PHI_right = K.filter_right$est[3] 
    PHI_under = K.filter_under$est[3] 
    PHI_up = K.filter_up$est[3] 
    PHI_mean = mean ( c(  PHI_left, PHI_right, PHI_under, PHI_up  )  )
    
    THETA_left = K.filter_left$est[4] 
    THETA_right = K.filter_right$est[4] 
    THETA_under = K.filter_under$est[4] 
    THETA_up = K.filter_up$est[4] 
    THETA_mean = mean ( c(  THETA_left, THETA_right, THETA_under, THETA_up  )  )
    
    res = c(pixel.est_left, pixel.est_right, pixel.est_under, pixel.est_up, pixel.est.mean,
            ALPHA_left, ALPHA_right, ALPHA_under, ALPHA_up, ALPHA_mean,
            DELTA_left, DELTA_right, DELTA_under, DELTA_up, DELTA_mean,
            PHI_left, PHI_right, PHI_under, PHI_up, PHI_mean,
            THETA_left, THETA_right, THETA_under, THETA_up, THETA_mean
    )
    return(res)  
  }
  
  
  result_mc <-
    if (.Platform$OS.type == "windows") {
      lapply( X = 1L:M, FUN = iter )
    } else {
      pbmcapply::pbmclapply( X = 1L:M, FUN = iter, mc.cores = cores )
    } %>%
    unlist %>% 
    matrix( ncol = 25, nrow = M, byrow = TRUE )  
  
  
  
  colnames(result_mc) <- c(
    "Left", "Right", "Under", "Up", "Mean",
    "ALPHA_Left", "ALPHA_Right", "ALPHA_Under", "ALPHA_Up", "ALPHA_Mean",
    "DELTA_Left", "DELTA_Right", "DELTA_Under", "DELTA_Up", "DELTA_Mean",
    "PHI_Left", "PHI_Right", "PHI_Under", "PHI_Up", "PHI_Mean",
    "THETA_Left", "THETA_Right", "THETA_Under", "THETA_Up", "THETA_Mean"
  )
  
  return(list(result_mc = result_mc))
  
}

