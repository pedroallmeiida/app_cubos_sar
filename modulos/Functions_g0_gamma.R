#rm(list=ls())

#### Fun??es para o processo Gama
library(maxLik)
library(invgamma)
library(nortest)
library(tseries)
library(forecast)
library(ggplot2)
library(reshape2)
library(purrr)


#### Geradores ####

# Gera??o distribui??o G0
vP = par
rG0_AR1_log = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3]
  mu=c();z=c()
  z0 = 1
  
  mu[1] = exp(
    delta + phin*( log(z0) )
  )
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
  
  for(t in 2:(N+1)){
    mu[t] = exp(
      delta + phin*( log(z[t-1]) ) 
    )
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  

  return(z[2:(N+1)])
}
rG0_AR2_log = function(N,vP,L){
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
  
  for(t in 3:(N+2)){
    mu[t] = exp(
      delta + phin[1]*( log(z[t-1]) ) + phin[2]*( log(z[t-2]) ) 
    )
    z[t] = rinvgamma( 1, shape=-alpha, rate=mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look); 
  }
  
  return(z[3:(N+2)])
}
rG0_MA1_log = function(N,par,L){
  # N: number of outcomes
  # par: vector of parameters  
  
  q=1; look = L
  alpha = par[1]; delta = par[2]; theta = par[3]
  mu=c();z=c();
  
  error=rep(0, N)
  
  error0 = 0
  mu[1] = exp(
    delta + theta*error0
  )
  
  z[1] = rinvgamma( 1, shape =-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape = look, scale = 1/look); 
  
  
  for(t in 2:(N+1)){
    
    error[t-1] = log( z[t-1] ) - log( mu[t-1] )
    mu[t] = exp(
      delta + theta*error[t-1]
    )
    z[t] = rinvgamma( 1, shape = -alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape = look, scale = 1/look); 
  }
  
  return(z[2:(N+1)])
}
rG0_ARMA11_log = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3];theta = vP[4]
  mu=c();z=c()
  z0 = 1;mu0=1
  
  mu[1] = exp(
    delta + phin*( log(z0) ) + theta*( log(z0) - log(mu0) )
  )
  z[1] = rinvgamma( 1, shape=-alpha, rate=mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
  
  
  for(t in 2:(N+1)){
    mu[t] = exp(
      delta + phin*( log(z[t-1]) ) + theta*( log(z[t-1]) - log(mu[t-1]) ) 
    )
    z[t] = rinvgamma( 1, shape=-alpha, rate=mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  return(z[2:(N+1)])
}

rG0_AR1_W = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3]
  mu=c();z=c()
  z0 = 1
  
  require(lamW)
  # to find optimal lambda
  eta=c()
  
  eta[1] = delta + phin*( log(lambertW0(1)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )  
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
  
  for(t in 2:(N+1)){
    eta[t] = delta + phin*( log(lambertW0( z[t-1] )) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  return(z[2:(N+1)])
}
rG0_AR2_W = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p = 2
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  z0 = 1
  
  require(lamW)
  eta=c()
  eta[1] = delta + phin[1]*( log(lambertW0(1)) ) + phin[2]*( log(lambertW0( 1 )) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
  eta[2] = delta + phin[1]*( log(lambertW0(z[1])) ) + phin[2]*( log(lambertW0( 1 )) ) 
  mu[2] = exp( eta[2] + exp(eta[2])  )  
  z[2] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
  
  for(t in 3:(N+2)){
    eta[t] = delta + phin[1]*( log(lambertW0( z[t-1] )) ) + phin[2]*( log(lambertW0( z[t-2] )) )
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }

  return(z[3:(N+2)])
}
rG0_MA1_W = function(N,par,L){
  # N: number of outcomes
  # par: vector of parameters  
  
  q=1; look = L
  alpha = par[1]; delta = par[2]; theta = par[3]
  mu=c();z=c();
  
  
  require(lamW)
  eta=c()
  eta[1] = delta 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
 
  for(t in 2:(N+1)){
    eta[t] = delta + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  return(z[2:(N+1)])
}
rG0_ARMA11_W = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3];theta = vP[4]
  mu=c();z=c()
  z0 = 1;mu0=1
  eta=c()
  
  require(lamW)
  
  eta[1] = delta + phin*( log(lambertW0(1)) ) + theta*( log(lambertW0(1)) - log(lambertW0(1)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
 
  for(t in 2:(N+1)){
    eta[t] = delta + phin*( log(lambertW0( z[t-1] )) ) + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  return(z[2:(N+1)])
}

rG0_AR1_boxcox = function(N,vP,lambda,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  look=L
  alpha = vP[1]; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  eta=c()
  
  eta[1] = delta + phin*( ( 1^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  
  for(t in 2:N){
    eta[t] = delta + phin*( ( z[t-1]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda) 
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  return(z)
  
}
rG0_AR2_boxcox = function(N,vP,lambda,L){
  # N: number of outcomes
  # vP: vector of parameters  
  look=L
  alpha = vP[1]; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  eta=c()
  
  z0 = rinvgamma( 1, shape=-alpha, rate = 1*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  
  eta[1] = delta + phin[1]*( ( z0^(1/lambda) - 1 )/lambda ) + phin[2]*( ( z0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  
  eta[2] = delta + phin[1]*( ( z[1]^(1/lambda) - 1 )/lambda ) + phin[2]*( ( z0^(1/lambda) - 1 )/lambda ) 
  mu[2] = (lambda*eta[2] + 1)^(1/lambda)
  z[2] = rinvgamma( 1, shape=-alpha, rate = mu[2]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  
  for(t in 3:N){
    eta[t] = delta + phin[1]*( ( z[t-1]^(1/lambda) - 1 )/lambda ) + phin[2]*( ( z[t-2]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda)
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
    
  }
  
  
  return(z)
  
}
rG0_MA1_boxcox = function(N,vP,lambda,L){
  # N: number of outcomes
  # vP: vector of parameters  
  look=L
  alpha = vP[1]; delta = vP[2]; theta = vP[3:length(vP)]
  mu=c();z=c();  eta=c()
  
  mu0 = 1
  z0 = rinvgamma( 1, shape=-alpha, rate = 1*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  
  
  eta[1] = delta + theta*( ( z0^(1/lambda) - 1 )/lambda - ( mu0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)  
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  
  for(t in 2:N){
    eta[t] = delta + theta*( ( z[t-1]^(1/lambda) - 1 )/lambda - ( mu[t-1]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda) 
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  
  return(z)
  
}
rG0_ARMA11_boxcox = function(N,vP,lambda,L){
  # N: number of outcomes
  # vP: vector of parameters  
  look=L
  alpha = vP[1]; delta = vP[2]; phin = vP[3:length(vP)]; theta = vP[3:length(vP)]
  mu=c();z=c()
  eta=c()
  
  mu0 = 1
  z0 = rinvgamma( 1, shape=-alpha, rate = 1*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  
  eta[1] = delta + phin*( ( z0^(1/lambda) - 1 )/lambda )  + theta*( ( z0^(1/lambda) - 1 )/lambda - ( mu0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)  
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  
  for(t in 2:N){
    eta[t] = delta + phin*( ( z[t-1]^(1/lambda) - 1 )/lambda )  + theta*( ( z[t-1]^(1/lambda) - 1 )/lambda - ( mu[t-1]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda) 
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  return(z)
  
}


rG0_AR1_inv = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3]
  mu=c();z=c()
  z0 = 1
  eta=c()
  eta[1] = delta + phin*( 1/(z0) )
  mu[1] = 1/eta[1]
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
  
  
  for(t in 2:N){
    eta[t] = delta + phin*( 1/(z[t-1]) ) 
    mu[t] = 1/eta[t]
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
    }
  
  return(z)
}
rG0_AR2_inv = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p = 2
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  z0 = 1
  eta=c()
  
  eta[1] = delta + phin[1]*( 1/(z0) ) + phin[2]*( 1/(z0) )
  mu[1] = 1/eta[1]
  z[1] = rinvgamma( 1, shape=-alpha, rate=mu[1]*(-alpha-1) ) 
  eta[2] = delta + phin[1]*( 1/(z[1]) ) + phin[2]*( 1/(z0) ) 
  mu[2] = 1/eta[2]
  z[2] = rinvgamma( 1, shape=-alpha, rate=mu[2]*(-alpha-1) ) 
  
  for(t in 3:N){
    eta[t] = delta + phin[1]*( 1/(z[t-1]) ) + phin[2]*( 1/(z[t-2]) ) 
    mu[t] = 1/eta[t]
    z[t] = rinvgamma( 1, shape=-alpha, rate=mu[t]*(-alpha-1) )
  }
  
  au = rgamma(N, shape=look, scale=1/look); 
  x = au*z
  return(x)
}
rG0_MA1_inv = function(N,par,L){
  # N: number of outcomes
  # par: vector of parameters  
  
  q=1; look = L
  alpha = par[1]; delta = par[2]; theta = par[3]
  mu=c();z=c();
  eta=c()
  
  eta[1] = delta + theta*(  1/( 1 ) - 1/( 1 )  )
  mu[1] = 1
  z[1] = rinvgamma( 1, shape =-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape = look, scale = 1/look); 
  
  
  for(t in 2:N){
    eta[t] = delta + theta*( 1/( z[t-1] ) - 1/( mu[t-1] )  )
    mu[t] = 1/eta[t]
    z[t] = rinvgamma( 1, shape = -alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape = look, scale = 1/look); 
  }
  
  return(z)
}
rG0_ARMA11_inv = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3];theta = vP[4]
  mu=c();z=c()
  z0 = 1;mu0=1
  eta=c()
  
  eta[1] = delta + phin*( 1/(z0) ) + theta*( 1/(z0) - 1/(mu0) )
  mu[1] = 1/eta[1]
  z[1] = rinvgamma( 1, shape=-alpha, rate=mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
  
  
  for(t in 2:N){
    eta[t] = delta + phin*( 1/(z[t-1]) ) + theta*( 1/(z[t-1]) - 1/(mu[t-1]) ) 
    mu[t] = 1/eta[t]
    z[t] = rinvgamma( 1, shape=-alpha, rate=mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  return(z)
}


rG0_AR1_sqrt = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3]
  mu=c();z=c()
  z0 = 1
  eta=c()
  eta[1] = delta + phin*( sqrt(z0) )
  mu[1] = eta[1]^2
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
  
  
  for(t in 2:(N+1)){
    eta[t] = delta + phin*( sqrt(z[t-1]) ) 
    mu[t] = eta[t]^2
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  return(z[2:(N+1)])
}
rG0_AR2_sqrt = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p = 2
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  z0 = 1
  eta=c()
  
  eta[1] = delta + phin[1]*( sqrt(z0) ) + phin[2]*( sqrt(z0) )
  mu[1] = eta[1]^2
  z[1] = rinvgamma( 1, shape=-alpha, rate=mu[1]*(-alpha-1) ) 
  eta[2] = delta + phin[1]*( sqrt(z[1]) ) + phin[2]*( sqrt(z0) ) 
  mu[2] = eta[2]^2
  z[2] = rinvgamma( 1, shape=-alpha, rate=mu[2]*(-alpha-1) ) 
  
  for(t in 3:(N+2)){
    eta[t] = delta + phin[1]*( sqrt(z[t-1]) ) + phin[2]*( sqrt(z[t-2]) ) 
    mu[t] = eta[t]^2
    z[t] = rinvgamma( 1, shape=-alpha, rate=mu[t]*(-alpha-1) )* rgamma(1, shape=look, scale=1/look); 
  }
  
  return(z[3:(N+2)])
}
rG0_MA1_sqrt = function(N,par,L){
  # N: number of outcomes
  # par: vector of parameters  
  
  q=1; look = L
  alpha = par[1]; delta = par[2]; theta = par[3]
  mu=c();z=c();
  eta=c()
  
  eta[1] = delta + theta*(  sqrt( 1 ) - sqrt( 1 )  )
  mu[1] = 1
  z[1] = rinvgamma( 1, shape =-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape = look, scale = 1/look); 
  
  
  for(t in 2:(N+1)){
    eta[t] = delta + theta*( sqrt( z[t-1] ) - sqrt( mu[t-1] )  )
    mu[t] = eta[t]^2
    z[t] = rinvgamma( 1, shape = -alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape = look, scale = 1/look); 
  }
  
  return(z[2:(N+1)])
}
rG0_ARMA11_sqrt = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3];theta = vP[4]
  mu=c();z=c()
  z0 = 1;mu0=1
  eta=c()
  
  eta[1] = delta + phin*( sqrt(z0) ) + theta*( sqrt(z0) - sqrt(mu0) )
  mu[1] =eta[1]^2
  z[1] = rinvgamma( 1, shape=-alpha, rate=mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
  
  
  for(t in 2:(N+1)){
    eta[t] = delta + phin*( sqrt(z[t-1]) ) + theta*( sqrt(z[t-1]) - sqrt(mu[t-1]) ) 
    mu[t] = eta[t]^2
    z[t] = rinvgamma( 1, shape=-alpha, rate=mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  return(z[2:(N+1)])
}

## Gera??o da distribui??o gamma


rgama_AR1_log = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  
  delta = vP[1]; phin = vP[2:length(vP)]
  mu=c();z=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = L, scale = 1/(L/1) )
  
  mu[1] = exp(
    delta +  phin*( log(z0) )    
  )
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )
  
  
  for(t in 2:N){
    mu[t] = exp(
      delta + phin*( log(z[t-1]) )   
    )
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
    
  }
  
  return(z)
  
}
rgama_AR2_log = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  delta = vP[1]; phin = vP[2:length(vP)]
  mu=c();z=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = L, scale = 1/(L/1) )
  
  mu[1] = exp(
    delta +  phin[1]*( log(z0) ) +  phin[2]*( log(z0) )  
  )
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )
  
  mu[2] = exp(
    delta +  phin[1]*( log(z[1]) ) +  phin[2]*( log(z0) )   
  )
  z[2] = rgamma( 1, shape = L, scale = 1/(L/mu[2]) )
  
  
  for(t in 3:N){
    mu[t] = exp(
      delta +  phin[1]*( log(z[t-1]) ) +  phin[2]*( log(z[t-2]) )   
    )
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
    
  }
  
  return(z)
  
}
rgama_MA1_log = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  delta = vP[1]; theta = vP[2]
  mu=c();z=c();
  
  mu0 = 1
  x0 = rgamma( 1, shape = L, scale = 1/(L/1) )
  z0 = rgamma( 1, shape = L, scale = 1/(L/1) )
  
  mu[1] = exp(
    delta +  theta*( log(z0) - log(mu0)  )  
  )
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )
  
  for(t in 2:N){
    mu[t] = exp(
      delta +  theta*( log(z[t-1]) - log(mu[t-1]) )  
    )
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
    
  }
  
  return(z)
  
}
rgama_ARMA11_log = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  delta = vP[1]; phin = vP[2]; theta = vP[3]
  mu=c();z=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = L, scale = 1/(L/1) )
  
  mu[1] = exp(
    delta + phin*( log(z0) ) +  theta*( log(z0) - log(mu0)  )  
  )
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )
  
  for(t in 2:N){
    mu[t] = exp(
      delta + phin*( log(z[t-1]) ) +  theta*( log(z[t-1]) - log(mu[t-1]) )  
    )
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
    
  }
  
  return(z)
  
}

rgama_AR1_W = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  
  delta = vP[1]; phin = vP[2:length(vP)]
  mu=c();z=c()
  eta=c()
  
  require(lamW)
  # to find optimal lambda
  
  eta[1] = delta + phin*( log(lambertW0(1)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )  
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )
  
  for(t in 2:N){
    eta[t] = delta + phin*( log(lambertW0( z[t-1] )) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
  }
  
  return(z)
  
}
rgama_AR2_W = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  delta = vP[1]; phin = vP[2:length(vP)]
  mu=c();z=c()
  eta=c()
  
  
  require(lamW)
  z0 = rgamma( 1, shape = L, scale = 1/(L/1) )
  
  eta[1] = delta + phin[1]*( log(lambertW0(z0)) ) + phin[2]*( log(lambertW0( z0 )) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )
  eta[2] = delta + phin[1]*( log(lambertW0(z[1])) ) + phin[2]*( log(lambertW0( z0 )) ) 
  mu[2] = exp( eta[2] + exp(eta[2])  )  
  z[2] = rgamma( 1, shape = L, scale = 1/(L/mu[2]) )
  
  for(t in 3:N){
    eta[t] = delta + phin[1]*( log(lambertW0( z[t-1] )) ) + phin[2]*( log(lambertW0( z[t-2] )) )
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
  }
  
  
  return(z)
  
}
rgama_MA1_W = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  delta = vP[1]; theta = vP[2]
  mu=c();z=c();  eta=c()

  require(lamW)
  
  mu0 = 1
  z0 = rgamma( 1, shape = L, scale = 1/(L/1) )
  
  
  eta[1] = delta + theta*( log(lambertW0(z0)) - log(lambertW0(mu0)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) ) 
  
  for(t in 2:N){
    eta[t] = delta + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
  }
  
  
  return(z)
  
}
rgama_ARMA11_W = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  delta = vP[1]; phin = vP[2]; theta = vP[3]
  mu=c();z=c()
  eta=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = L, scale = 1/(L/1) )

  require(lamW)
  
  eta[1] = delta + phin*( log(lambertW0(z0)) ) + theta*( log(lambertW0(z0)) - log(lambertW0(mu0)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )
  
  for(t in 2:N){
    eta[t] = delta + phin*( log(lambertW0( z[t-1] )) ) + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    z[t] = rgamma( 1, shape = L, scale = 1/(L/mu[t]) )
  }
  
  return(z)
  
}

rgama_AR1_boxcox = function(N,vP,lambda,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  
  delta = vP[1]; phin = vP[2:length(vP)]
  mu=c();z=c()
  eta=c()
  
  eta[1] = delta + phin*( ( 1^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )
  
  for(t in 2:N){
    eta[t] = delta + phin*( ( z[t-1]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda) 
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
  }
  
  return(z)
  
}
rgama_AR2_boxcox = function(N,vP,lambda,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  delta = vP[1]; phin = vP[2:length(vP)]
  mu=c();z=c()
  eta=c()
  
  z0 = rgamma( 1, shape = L, scale = 1/(L/1) )
  
  eta[1] = delta + phin[1]*( ( z0^(1/lambda) - 1 )/lambda ) + phin[2]*( ( z0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )
  
  eta[2] = delta + phin[1]*( ( z[1]^(1/lambda) - 1 )/lambda ) + phin[2]*( ( z0^(1/lambda) - 1 )/lambda ) 
  mu[2] = (lambda*eta[2] + 1)^(1/lambda)
  z[2] = rgamma( 1, shape = L, scale = 1/(L/mu[2]) )
  
  for(t in 3:N){
    eta[t] = delta + phin[1]*( ( z[t-1]^(1/lambda) - 1 )/lambda ) + phin[2]*( ( z[t-2]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda)
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
  }
  
  
  return(z)
  
}
rgama_MA1_boxcox = function(N,vP,lambda,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  delta = vP[1]; theta = vP[2]
  mu=c();z=c();  eta=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = L, scale = 1/(L/1) )
  
  
  eta[1] = delta + theta*( ( z0^(1/lambda) - 1 )/lambda - ( mu0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)  
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) ) 
  
  for(t in 2:N){
    eta[t] = delta + theta*( ( z[t-1]^(1/lambda) - 1 )/lambda - ( mu[t-1]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda) 
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
  }
  
  
  return(z)
  
}
rgama_ARMA11_boxcox = function(N,vP,lambda,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  delta = vP[1]; phin = vP[2]; theta = vP[3]
  mu=c();z=c()
  eta=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = L, scale = 1/(L/1) )
  
  eta[1] = delta + phin*( ( z0^(1/lambda) - 1 )/lambda )  + theta*( ( z0^(1/lambda) - 1 )/lambda - ( mu0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)  
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )
  
  for(t in 2:N){
    eta[t] = delta + phin*( ( z[t-1]^(1/lambda) - 1 )/lambda )  + theta*( ( z[t-1]^(1/lambda) - 1 )/lambda - ( mu[t-1]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda) 
    z[t] = rgamma( 1, shape = L, scale = 1/(L/mu[t]) )
  }
  
  return(z)
  
}

rgama_AR1_inv = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  look = L; delta = vP[1]; phin = vP[2]
  mu=c();z=c()
  z0 = 1
  eta=c()
  eta[1] = delta + phin*( 1/(z0) )
  mu[1] = 1/eta[1]
  z[1] = rgamma(1,shape = L, scale = 1/(L/mu[1]))
  
  for(t in 2:N){
    eta[t] = delta + phin*( 1/(z[t-1]) ) 
    mu[t] = 1/eta[t]
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
  }
  
  
  return(z)
}
rgama_AR2_inv = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p = 2
  look = L; delta = vP[1]; phin = vP[2:length(vP)]
  mu=c();z=c()
  z0 = 1
  eta=c()
  
  eta[1] = delta + phin[1]*( 1/(z0) ) + phin[2]*( 1/(z0) )
  mu[1] = 1/eta[1]
  z[1] = rgamma(1,shape = L, scale = 1/(L/mu[1]))
  eta[2] = delta + phin[1]*( 1/(z[1]) ) + phin[2]*( 1/(z0) ) 
  mu[2] = 1/eta[2]
  z[2] = rgamma(1,shape = L, scale = 1/(L/mu[2]))
  
  for(t in 3:N){
    eta[t] = delta + phin[1]*( 1/(z[t-1]) ) + phin[2]*( 1/(z[t-2]) ) 
    mu[t] = 1/eta[t]
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
  }
  
  return(z)
}
rgama_MA1_inv = function(N,vP,L){
  # N: number of outcomes
  # par: vector of parameters  
  
  q=1; look = L
  delta = vP[1]; theta = vP[2]
  mu=c();z=c();
  eta=c()
  
  eta[1] = 1
  mu[1] = 1/eta[1]
  z[1] = rgamma( 1,shape = L, scale = 1/(L/mu[1]) )  
  
  for(t in 2:N){
    eta[t] = delta + theta*(  1/( z[t-1] ) - 1/( mu[t-1] )  )
    mu[t] = 1/eta[t]
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
    }
  
  return(z)
}
rgama_ARMA11_inv = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  look = L; delta = vP[1]; phin = vP[2];theta = vP[3]
  mu=c();z=c()
  z0 = 1;mu0=1
  eta=c()
  
  eta[1] = delta + phin*( 1/(z0) ) + theta*( 1/(z0) - 1/(mu0) )
  mu[1] = 1/eta[1]
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )  
  
  for(t in 2:N){
    eta[t] = delta + phin*( 1/(z[t-1]) ) + theta*( 1/(z[t-1]) - 1/(mu[t-1]) ) 
    mu[t] = 1/eta[t]
    z[t] = rgamma( 1, shape = L, scale = 1/(L/mu[t]) )
    }
  
  return(z)
}


rgama_AR1_sqrt = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  look = L; delta = vP[1]; phin = vP[2]
  mu=c();z=c()
  z0 = 1
  eta=c()
  eta[1] = delta + phin*( sqrt(z0) )
  mu[1] = eta[1]^2
  z[1] = rgamma(1,shape = L, scale = 1/(L/mu[1]))
  
  for(t in 2:N){
    eta[t] = delta + phin*( sqrt(z[t-1]) ) 
    mu[t] = eta[t]^2
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
  }
  
  
  return(z)
}
rgama_AR2_sqrt = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p = 2
  look = L; delta = vP[1]; phin = vP[2:length(vP)]
  mu=c();z=c()
  z0 = 1
  eta=c()
  
  eta[1] = delta + phin[1]*( sqrt(z0) ) + phin[2]*( sqrt(z0) )
  mu[1] = eta[1]^2
  z[1] = rgamma(1,shape = L, scale = 1/(L/mu[1]))
  eta[2] = delta + phin[1]*( sqrt(z[1]) ) + phin[2]*( sqrt(z0) ) 
  mu[2] = eta[2]^2
  z[2] = rgamma(1,shape = L, scale = 1/(L/mu[2]))
  
  for(t in 3:N){
    eta[t] = delta + phin[1]*( sqrt(z[t-1]) ) + phin[2]*( sqrt(z[t-2]) ) 
    mu[t] = eta[t]^2
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
  }
  
  return(z)
}
rgama_MA1_sqrt = function(N,vP,L){
  # N: number of outcomes
  # par: vector of parameters  
  
  q=1; look = L
  delta = vP[1]; theta = vP[2]
  mu=c();z=c();
  eta=c()
  
  eta[1] = 1
  mu[1] = eta[1]^2
  z[1] = rgamma( 1,shape = L, scale = 1/(L/mu[1]) )  
  
  for(t in 2:N){
    eta[t] = delta + theta*(  sqrt( z[t-1] ) - sqrt( mu[t-1] )  )
    mu[t] = eta[t]^2
    z[t] = rgamma(1,shape = L, scale = 1/(L/mu[t]))
  }
  
  return(z)
}
rgama_ARMA11_sqrt = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  look = L; delta = vP[1]; phin = vP[2];theta = vP[3]
  mu=c();z=c()
  z0 = 1;mu0=1
  eta=c()
  
  eta[1] = delta + phin*( sqrt(z0) ) + theta*( sqrt(z0) - sqrt(mu0) )
  mu[1] = eta[1]^2
  z[1] = rgamma( 1, shape = L, scale = 1/(L/mu[1]) )  
  
  for(t in 2:N){
    eta[t] = delta + phin*( sqrt(z[t-1]) ) + theta*( sqrt(z[t-1]) - sqrt(mu[t-1]) ) 
    mu[t] = eta[t]^2
    z[t] = rgamma( 1, shape = L, scale = 1/(L/mu[t]) )
  }
  
  return(z)
}

##### end geradores ####



####    Log verossimilhan?as     #####

# Distribui??o G0

#AR(1)
lk.g0.AR1_log = function(z,par,L){
  
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
  
  mu[1] = exp(
    delta + phin*( log(z0) )
  )
  
  v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
  lk[1] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1]) + (alpha-L)*log(v[1])
  
  
  for(t in (p+1):n){
    mu[t] = exp(
      delta + phin*( log(z[t-1]) ) 
    )
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
    
  }
  
  return(sum(lk[2:n]))
  
}
lk.g0.AR2_log = function(z,par,L){
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

  mu[1] = exp(
    delta + phin[1]*( log(z0) ) + phin[2]*( log(z0) ) 
  )
  mu[2] = exp(
    delta + phin[1]*( log(z[1]) ) + phin[2]*( log(z0) ) 
  )
  
  v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
  lk[1] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1]) - alpha*log(-alpha-1) -log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1]) + (alpha-L)*log(v[1])
  v[2] = ( mu[2]*(-alpha-1) + L*z[2] )
  lk[2] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[2]) - alpha*log(-alpha-1) -log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[2]) + (alpha-L)*log(v[2])
  
  for(t in (p+1):n){
    mu[t] = exp(
      delta + phin[1]*( log(z[t-1]) ) + phin[2]*( log(z[t-2]) ) 
    )
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1)-log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
    
  }
  
  return(sum(lk[3:n]))
  
}
lk.g0.MA1_log <- function(z, par,L){
  
  n=length(z)
  alpha <- par[1]
  delta<- par[2]
  theta <- par[3]
  
  v=c()
  m=1;q=1
  
  eta <- mu <- v <- rep(NA,n)
  error<-rep(0,n)
  
  
  for(t in (m+1):n){
    eta[t] <- delta + theta*error[t-q] 
    mu[t]   <- exp( eta[t] )
    error[t] <- log(z[t]) - log(mu[t]) 
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 
lk.g0.MA2_log <- function(z, par,L){
  
  n=length(z)
  alpha = par[1]
  delta<- par[2]
  theta <- par[3:4]
  z0=mean(z);mu0=1
  v=c()
  m=2;q=2
  eta <- mu <- v <- c()
  eta[1] <- delta + theta[1]*( log(z0) - log(mu0)  ) + theta[2]*( log(z0) - log(mu0)  ) 
  mu[1]  <- exp( eta[1] )
  eta[2] <- delta + theta[1]*( log(z[1]) - log(mu[1])  ) + theta[2]*( log(z0) - log(mu0)  ) 
  mu[2]   <- exp( eta[2] )
  
  for(t in (m+1):n){
    eta[t] <- delta + theta[1]*( log(z[t-1]) - log(mu[t-1])  ) + theta[2]*( log(z[t-2]) - log(mu[t-2])  ) 
    mu[t]  <- exp( eta[t] )
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 
lk.g0.ARMA11_log <- function(z, par,L){
  
  n=length(z)
  alpha <- par[1]
  delta <- par[2]
  phi <- par[3]
  theta <- par[4]
  
  m=1;q=1;p=1
  
  eta <- mu <- v <- rep(NA,n)
  error<-rep(0,n)
  
  
  for(t in (m+1):n){
    eta[t] <- delta + phi*log(z[t-p]) + theta*error[t-q] 
    mu[t]   <- exp( eta[t] )
    error[t] <- log(z[t]) - log(mu[t]) 
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 


lk.g0.AR1_W = function(z,par,L){
  
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
  
  v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
  lk[1] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1]) + (alpha-L)*log(v[1])
  
  
  for(t in (p+1):n){
    eta[t] = delta + phin*( log(lambertW0(z[t-1])) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
    
  }
  
  return(sum(lk[2:n]))
  
}
lk.g0.AR2_W = function(z,par,L){
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


  v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
  lk[1] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1]) - alpha*log(-alpha-1) -log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1]) + (alpha-L)*log(v[1])
  v[2] = ( mu[2]*(-alpha-1) + L*z[2] )
  lk[2] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[2]) - alpha*log(-alpha-1) -log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[2]) + (alpha-L)*log(v[2])
  
  for(t in (p+1):n){
    eta[t] = delta + phin[1]*( log(lambertW0(z[t-1])) ) + phin[2]*( log(lambertW0( z[t-2] )) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1)-log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
    
  }
  
  return(sum(lk[3:n]))
  
}
lk.g0.MA1_W <- function(z, par,L){
  
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
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 
lk.g0.MA2_W <- function(z, par,L){
  
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
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 
lk.g0.ARMA11_W <- function(z, par,L){
  
  n=length(z)
  alpha <- par[1];delta <- par[2];phi <- par[3];theta <- par[4]
  
  m=1;q=1;p=1
  
  z0 = mean(z)
  
  eta <- mu <- v <- c()
  eta[1] = delta + phi*( log(lambertW0(z0)) ) + theta*( log(lambertW0(z0)) - log(lambertW0(z0)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  v[1] <- ( mu[1]*(-alpha-1) + L*z[1] )
  
  for(t in (m+1):n){
    eta[t] <- delta + phi*( log(lambertW0(z[t-1])) ) + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) 
    mu[t]   <- exp( eta[t] + exp(eta[t])  )   
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 

lk.g0.AR1_boxcox = function(z,par,lambda,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(z)
  delta = par[1]; phin=par[2]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  z0 = rep(mean(z),p)
  
  #lambda <- BoxCox.lambda(z)
  
  # to find optimal lambda
  eta=c()
  
  eta[1] = delta + phin*( ( z0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda) 
  v[1] = ( mu[t]*(-alpha-1) + L*z[1] )
  lk[1] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1]) + (alpha-L)*log(v[1])
  
  t<-2:n
  eta[t] = delta + phin*( ( z[t-1]^(1/lambda) - 1 )/lambda ) 
  mu[t] = (lambda*eta[t] + 1)^(1/lambda)
  v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
  
  
  
  return(sum(lk[2:n]))
  
}
lk.g0.AR2_boxcox = function(z,par,lambda,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(z)
  delta = par[1]; phin=par[2:length(par)]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  z0 = rep(mean(z),p)
  
  #lambda <- BoxCox.lambda(z)
  
  # to find optimal lambda
  eta=c()
  
  eta[1] = delta + phin[1]*( ( z0^(1/lambda) - 1 )/lambda ) + phin[2]*( ( z0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda) 
  eta[2] = delta + phin[1]*( ( z[1]^(1/lambda) - 1 )/lambda ) + phin[2]*( ( z0^(1/lambda) - 1 )/lambda ) 
  mu[2] = (lambda*eta[2] + 1)^(1/lambda) 
  
  v[1:2] = ( mu[t]*(-alpha-1) + L*z[1:2] )
  lk[1:2] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1:2]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1:2]) + (alpha-L)*log(v[1:2])
  
  t<-2:n
  eta[t] = delta + phin[1]*( ( z[t-1]^(1/lambda) - 1 )/lambda ) + phin[2]*( ( z[t-2]^(1/lambda) - 1 )/lambda ) 
  mu[t] = (lambda*eta[t] + 1)^(1/lambda)
  v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
  
  
  
  return(sum(lk[2:n]))
  
}
lk.g0.MA1_boxcox = function(z,par,lambda,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 1
  n = length(z)
  delta = par[1]; thetan=par[2]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,q)
  z0 = rep(mean(z),q)
  
  #lambda <- BoxCox.lambda(z)
  
  eta <- mu <- v <- c()
  
  eta[1] = delta + theta*( ( z0^(1/lambda) - 1 )/lambda - ( mu0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)
  v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
  lk[1] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1]) + (alpha-L)*log(v[1])
  
  for(t in 2:n){
    eta[t] = delta + theta*( ( z[t-1]^(1/lambda) - 1 )/lambda - ( mu[t-1]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda) 
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  }
  t <- 2:n
  lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
  
  
  return( sum(lk[2:n]) )
  
}
lk.g0.ARMA11_boxcox = function(z,par,lambda,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(z)
  delta = par[1]; phin = par[2]; thetan=par[3]
  mu=c();v=c();lk=c()
  
  #lambda <- BoxCox.lambda(z)
  
  mu0 = rep(1,p)
  z0 = rep( mean(z), p )
  
  eta <- mu <- v <- c()
  eta[1] = delta + phin*( ( z0^(1/lambda) - 1 )/lambda  ) + thetan*( ( z0^(1/lambda) - 1 )/lambda - ( mu0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)
  v[1] = ( mu[t]*(-alpha-1) + L*z[1] )
  lk[1] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1]) + (alpha-L)*log(v[1])
  
  for(t in (p+1):n){
    eta[t] = delta + phin*( ( z[t-1]^(1/lambda) - 1 )/lambda  ) + thetan*( ( z[t-1]^(1/lambda) - 1 )/lambda - ( mu[t-1]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda)
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  
  t <- (p+1):n
  lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
  
  
  return(sum(lk[2:n]))
  
}


lk.g0.AR1_inv = function(z,par,L){
  
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
  eta[1] = delta + phin*( 1/(z0) )
  mu[1] = 1/eta[1]
  
  v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
  lk[1] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1]) + (alpha-L)*log(v[1])
  
  
  for(t in (p+1):n){
    eta[t] = delta + phin*( 1/(z[t-1]) ) 
    mu[t] = 1/eta[t]
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
    
  }
  
  return(sum(lk[2:n]))
  
}
lk.g0.AR2_inv = function(z,par,L){
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
  
  eta[1] = delta + phin[1]*( 1/(z0) ) + phin[2]*( 1/(z0) ) 
  mu[1] = 1/eta[1]
  eta[2] = delta + phin[1]*( 1/(z[1]) ) + phin[2]*( 1/(z0) ) 
  mu[2] = 1/eta[2]
  
  v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
  lk[1] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1]) - alpha*log(-alpha-1) -log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1]) + (alpha-L)*log(v[1])
  v[2] = ( mu[2]*(-alpha-1) + L*z[2] )
  lk[2] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[2]) - alpha*log(-alpha-1) -log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[2]) + (alpha-L)*log(v[2])
  
  for(t in (p+1):n){
    eta[t] = delta + phin[1]*( 1/(z[t-1]) ) + phin[2]*( 1/(z[t-2]) ) 
    mu[t] = 1/eta[t]
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1)-log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
    
  }
  
  return(sum(lk[3:n]))
  
}
lk.g0.MA1_inv <- function(z, par,L){
  
  n=length(z)
  alpha <- par[1]
  delta<- par[2]
  theta <- par[3]

  m=1;q=1
  
  eta <- mu <- v <- c()
  eta[1] = delta + theta*(  1/(mean(z)) - 1/(1)    )
  mu[1] = 1/eta[1]
  v[1] <- ( mu[1]*(-alpha-1) + L*z[1] )

  for(t in (m+1):n){
    eta[t] <- delta + theta*( 1/(z[t-1]) - 1/(mu[t-1])    )
    mu[t] <- 1/eta[t]
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 
lk.g0.MA2_inv <- function(z, par,L){
  
  n=length(z)
  alpha = par[1]
  delta<- par[2]
  theta <- par[3:4]
  z0=mean(z);mu0=1
  v=c()
  m=2;q=2
  eta <- mu <- v <- c()
  eta[1] <- delta + theta[1]*( 1/(z0) - 1/(mu0)  ) + theta[2]*( 1/(z0) - 1/(mu0)  ) 
  mu[1]  <- 1/eta[1]
  eta[2] <- delta + theta[1]*( 1/(z[1]) - 1/(mu[1])  ) + theta[2]*( 1/(z0) - 1/(mu0)  ) 
  mu[2]   <- 1/eta[2]
  
  for(t in (m+1):n){
    eta[t] <- delta + theta[1]*( 1/(z[t-1]) - 1/(mu[t-1])  ) + theta[2]*( 1/(z[t-2]) - 1/(mu[t-2])  ) 
    mu[t]  <- 1/eta[t]
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 
lk.g0.ARMA11_inv <- function(z, par,L){
  
  n=length(z)
  alpha <- par[1]
  delta <- par[2]
  phi <- par[3]
  theta <- par[4]
  
  m=1;q=1;p=1
  
  eta <- mu <- v <- c()
  eta[1] = delta + phi*1/( mean(z)  ) + theta*(  1/(mean(z)) - 1/(1)   ) 
  mu[1] = 1/eta[1]
  
  for(t in (m+1):n){
    eta[t] <- delta + phi*1/(z[t-1]) + theta*(  1/(z[t-1]) - 1/(mu[t-1])   )
    mu[t]   <- 1/eta[t]
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 


lk.g0.AR1_sqrt = function(z,par,L){
  
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
  
  v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
  lk[1] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1]) + (alpha-L)*log(v[1])
  
  
  for(t in (p+1):n){
    eta[t] = delta + phin*( sqrt(z[t-1]) ) 
    mu[t] = eta[t]^2
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1) - log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
    
  }
  
  return(sum(lk[2:n]))
  
}
lk.g0.AR2_sqrt = function(z,par,L){
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
  
  v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
  lk[1] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[1]) - alpha*log(-alpha-1) -log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[1]) + (alpha-L)*log(v[1])
  v[2] = ( mu[2]*(-alpha-1) + L*z[2] )
  lk[2] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[2]) - alpha*log(-alpha-1) -log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[2]) + (alpha-L)*log(v[2])
  
  for(t in (p+1):n){
    eta[t] = delta + phin[1]*( sqrt(z[t-1]) ) + phin[2]*( sqrt(z[t-2]) ) 
    mu[t] = eta[t]^2
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1)-log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
    
  }
  
  return(sum(lk[3:n]))
  
}


lk.g0.MA1_sqrt <- function(z, par,L){
  
  n=length(z)
  alpha <- par[1]
  delta<- par[2]
  theta <- par[3]
  
  m=1;q=1
  
  eta <- mu <- v <- c()
  eta[1] = delta + theta*(  sqrt(mean(z)) - sqrt(1)    )
  mu[1] = eta[1]^2
  v[1] <- ( mu[1]*(-alpha-1) + L*z[1] )
  
  for(t in (m+1):n){
    eta[t] <- delta + theta*( sqrt(z[t-1]) - sqrt(mu[t-1])    )
    mu[t] <- eta[t]^2
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 
lk.g0.MA2_sqrt <- function(z, par,L){
  
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
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 
lk.g0.ARMA11_sqrt <- function(z, par,L){
  
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
    v[t] <- ( mu[t]*(-alpha-1) + L*z[t] )
  }
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 

#### Log verossimilhan?a do processo Gama ----

lk.gamma.AR1_log = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  delta = par[1]; phin=par[2]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  y0 = rep(mean(y),p)
  
  mu[1] = exp(
    delta + phin*( log(y0) )
  )
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  t<-2:n
    mu[t] = exp(
      delta + phin*( log(y[t-1]) ) 
    )
    lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
    
  
  
  return(sum(lk[2:n]))
  
}
lk.gamma.AR2_log = function(y,par,L){
  
  ## ---> Input
  ## y: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=2
  n=length(y)
  delta = par[1]; phi=par[2:3]
  mu=c();v=c();lk=c()
  y0 = mean(y)
  
  mu[1] = exp(
    delta + phi[1]*( log(y0) ) + phi[2]*( log(y0) )
  )
  mu[2] = exp(
    delta + phi[1]*( log(y[1]) ) + phi[2]*( log(y0) )
  )
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1]) 
  lk[2] = (L-1)*log(y[2]) + L*log(L) - log( gamma(L)  ) - L*(  y[2]/mu[2] )  - L*log(mu[2]) 
  
  t <- (p+1):n
    mu[t] = exp(
      delta + phi[1]*( log(y[t-1]) ) + phi[2]*( log(y[t-2]) )
    )
    lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] )  - L*log(mu[t]) 
    
  
  
  return(sum(lk[3:n]))
  
}
lk.gamma.MA1_log = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 1
  n = length(y)
  delta = par[1]; thetan=par[2]
  mu=c();lk=c()
  
  mu0 = rep(1,q)
  y0 = rep(mean(y),q)
  
  mu[1] = exp(
    delta + thetan*( log(y0) - log(mu0) )
  )
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in 2:n){
    mu[t] = exp( delta + thetan*( log(y[t-1]) - log(mu[t-1]) ) )
  }
    t <- 2:n
    lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L) ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return( sum(lk[2:n]) )
  
}
lk.gamma.MA2_log = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 2
  n = length(y)
  delta = par[1]; thetan=par[2:3]
  mu=c();v=c();lk=c()
  
  mu0 = 1
  y0 = mean(y)
  
  mu[1:2] = exp(
    delta + thetan[1]*( log(y0) - log(mu0) ) + thetan[2]*( log(y0) - log(mu0) )
  )
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  lk[2] = (L-1)*log(y[2]) + L*log(L) - log( gamma(L)  ) - L*(  y[2]/mu[2] )  - L*log(mu[2])
  
  for(t in (q+1):n){
    mu[t] = exp(
      delta + thetan[1]*( log(y[t-1]) - log(mu[t-1]) )  + thetan[2]*( log(y[t-2]) - log(mu[t-2]) )
    )
  }
t<-(q+1):n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L) ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return(sum(lk[(q+1):n]))
  
}
lk.gamma.ARMA11_log = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  delta = par[1]; phin = par[2]; thetan=par[3]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  y0 = rep( mean(y), p )
  
  mu[1] = exp(
    delta + phin*log(y0) + thetan*( log(y0) - log(mu0) )
  )
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in (p+1):n){
    mu[t] = exp(
      delta + phin*log(y[t-1]) + thetan*( log(y[t-1]) - log(mu[t-1]) ) 
    )
  }
  
  
  t <- (p+1):n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return(sum(lk[2:n]))
  
}

lk.gamma.AR1_W = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  delta = par[1]; phin=par[2]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  y0 = rep(mean(y),p)
  
  require(lamW)
  # to find optimal lambda
  eta=c()
  
  eta[1] = delta + phin*( log(lambertW0(y0)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  t<-2:n
  eta[t] = delta + phin*( log(lambertW0(y[t-1])) ) 
  mu[t] = exp( eta[t] + exp(eta[t])  ) 
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  
  return(sum(lk[2:n]))
  
}
lk.gamma.AR2_W = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  delta = par[1]; phin=par[2:3]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  y0 = rep(mean(y),p)
  
  require(lamW)
  # to find optimal lambda
  eta=c()
  
  eta[1] = delta + phin[1]*( log(lambertW0(y0)) ) + phin[2]*( log(lambertW0(y0)) ) 
  eta[2] = delta + phin[1]*( log(lambertW0(y[1])) ) + phin[2]*( log(lambertW0(y0)) ) 
  
  mu[1:2] = exp( eta[1:2] + exp(eta[1:2])  )  
  lk[1:2] = (L-1)*log(y[1:2]) + L*log(L) - log( gamma(L)  ) - L*(  y[1:2]/mu[1:2] )  - L*log(mu[1:2])
  
  t<-3:n
  eta[t] = delta + phin[1]*( log(lambertW0(y[t-1])) ) + phin[2]*( log(lambertW0(y[t-2])) ) 
  mu[t] = exp( eta[t] + exp(eta[t])  ) 
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  
  return(sum(lk[3:n]))
  
}
lk.gamma.MA1_W = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 1
  n = length(y)
  delta = par[1]; theta=par[2]
  mu=c();lk=c()
  
  mu0 = rep(1,q)
  y0 = rep(mean(y),q)
  
  require(lamW)
  eta <- mu <- v <- c()
  
  eta[1] = delta + theta*( log(lambertW0(1)) - log(lambertW0(1)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in 2:n){
    eta[t] = delta + theta*( log(lambertW0(y[t-1])) - log(lambertW0(mu[t-1])) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
  }
  t <- 2:n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L) ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return( sum(lk[2:n]) )
  
}
lk.gamma.MA2_W = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 1
  n = length(y)
  delta = par[1]; theta=par[2:3]
  mu=c();lk=c()
  
  mu0 = rep(1,q)
  y0 = rep(mean(y),q)
  
  require(lamW)
  eta <- mu <- v <- c()
  
  eta[1] = delta + theta[1]*( log(lambertW0(y0)) - log(lambertW0(mu0)) ) + theta[1]*( log(lambertW0(y0)) - log(lambertW0(mu0)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )  
  eta[2] = delta + theta[1]*( log(lambertW0(y[1])) - log(lambertW0(mu[1])) ) + theta[2]*( log(lambertW0(y0)) - log(lambertW0(mu0)) ) 
  mu[2] = exp( eta[2] + exp(eta[2])  )  
  lk[1:2] = (L-1)*log(y[1:2]) + L*log(L) - log( gamma(L)  ) - L*(  y[1:2]/mu[1:2] )  - L*log(mu[1:2])
  
  for(t in 3:n){
    eta[t] = delta + theta[1]*( log(lambertW0(y[t-1])) - log(lambertW0(mu[t-1])) ) + theta[1]*( log(lambertW0(y[t-2])) - log(lambertW0(mu[t-2])) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )  
  }
  t <- 2:n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L) ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return( sum(lk[2:n]) )
  
}
lk.gamma.ARMA11_W = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  delta = par[1]; phin = par[2]; theta=par[3]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  y0 = rep( mean(y), p )
  
  eta <- mu <- v <- c()
  eta[1] = delta + phin*( log(lambertW0(1)) ) + theta*( log(lambertW0(1)) - log(lambertW0(1)) ) 
  mu[1] = exp( eta[1] + exp(eta[1])  )   
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in (p+1):n){
    eta[t] = delta + phin*( log(lambertW0(y[t-1])) ) + theta*( log(lambertW0(y[t-1])) - log(lambertW0(mu[t-1])) ) 
    mu[t] = exp( eta[t] + exp(eta[t])  )
  }
  
  
  t <- (p+1):n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return(sum(lk[2:n]))
  
}

lk.gamma.AR1_boxcox = function(y,par,lambda,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  delta = par[1]; phin=par[2]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  y0 = rep(mean(y),p)
  
  #lambda <- BoxCox.lambda(y)
  
  # to find optimal lambda
  eta=c()
  
  eta[1] = delta + phin*( ( y0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda) 
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  t<-2:n
  eta[t] = delta + phin*( ( y[t-1]^(1/lambda) - 1 )/lambda ) 
  mu[t] = (lambda*eta[t] + 1)^(1/lambda)
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  
  return(sum(lk[2:n]))
  
}
lk.gamma.MA1_boxcox = function(y,par,lambda,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 1
  n = length(y)
  delta = par[1]; theta=par[2]
  mu=c();lk=c()
  
  mu0 = rep(1,q)
  y0 = rep(mean(y),q)
  
  #lambda <- BoxCox.lambda(y)
  
  eta <- mu <- v <- c()
  
  eta[1] = delta + theta*( ( y0^(1/lambda) - 1 )/lambda - ( mu0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in 2:n){
    eta[t] = delta + theta*( ( y[t-1]^(1/lambda) - 1 )/lambda - ( mu[t-1]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda) 
  }
  t <- 2:n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L) ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return( sum(lk[2:n]) )
  
}
lk.gamma.ARMA11_boxcox = function(y,par,lambda,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  delta = par[1]; phin = par[2]; thetan=par[3]
  mu=c();v=c();lk=c()
  
  #lambda <- BoxCox.lambda(y)
  
  mu0 = rep(1,p)
  y0 = rep( mean(y), p )
  
  eta <- mu <- v <- c()
  eta[1] = delta + phin*( ( y0^(1/lambda) - 1 )/lambda  ) + thetan*( ( y0^(1/lambda) - 1 )/lambda - ( mu0^(1/lambda) - 1 )/lambda ) 
  mu[1] = (lambda*eta[1] + 1)^(1/lambda)
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in (p+1):n){
    eta[t] = delta + phin*( ( y[t-1]^(1/lambda) - 1 )/lambda  ) + thetan*( ( y[t-1]^(1/lambda) - 1 )/lambda - ( mu[t-1]^(1/lambda) - 1 )/lambda ) 
    mu[t] = (lambda*eta[t] + 1)^(1/lambda)
  }
  
  
  t <- (p+1):n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return(sum(lk[2:n]))
  
}



lk.gamma.AR1_inv = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  p=1
  n=length(y)
  delta = par[1]; phin=par[2]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  y0 = rep(mean(y),p)
  eta=c()
  eta[1] = delta + phin*( 1/(y0) )
  mu[1] = 1/eta[1]
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  t<-2:n
  eta[t] = delta + phin*( 1/(y[t-1]) ) 
  mu[t] = 1/eta[t]
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  
  return(sum(lk[2:n]))
  
}
lk.gamma.AR2_inv = function(y,par,L){
  
  ## ---> Input
  ## y: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=2
  n=length(y)
  delta = par[1]; phi=par[2:3]
  mu=c();v=c();lk=c()
  y0 = mean(y)
  
  eta=c()
  eta[1] = delta + phi[1]*( 1/(y0) ) + phi[2]*( 1/(y0) )
  eta[2] = delta + phi[1]*( 1/(y[1]) ) + phi[2]*( 1/(y0) )
  mu[1] = 1/eta[1]
  mu[2] = 1/eta[2]
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1]) 
  lk[2] = (L-1)*log(y[2]) + L*log(L) - log( gamma(L)  ) - L*(  y[2]/mu[2] )  - L*log(mu[2]) 
  
  t <- (p+1):n
  eta[t] = delta + phi[1]*( 1/(y[t-1]) ) + phi[2]*( 1/(y[t-2]) )
  mu[t] = 1/eta[t]
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] )  - L*log(mu[t]) 
  
  
  return(sum(lk[3:n]))
  
}
lk.gamma.MA1_inv = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 1
  n = length(y)
  delta = par[1]; thetan=par[2]
  mu=c();lk=c()
  
  mu0 = rep(1,q)
  y0 = rep(mean(y),q)
  eta=c()
  eta[1] = delta + thetan*( 1/(y0) - 1/(mu0) )
  mu[1] = 1/eta[1]
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in 2:n){
    eta[t] = delta + thetan*( 1/(y[t-1]) - 1/(mu[t-1]) )
    mu[t] = 1/eta[t]
  }
  t <- 2:n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L) ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return( sum(lk[2:n]) )
  
}
lk.gamma.MA2_inv = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 2
  n = length(y)
  delta = par[1]; thetan=par[2:3]
  mu=c();v=c();lk=c()
  
  mu0 = 1
  y0 = mean(y)
  eta=c()
  eta[1] = delta + thetan[1]*( 1/(y0) - 1/(mu0) ) + thetan[2]*( 1/(y0) - 1/(mu0) )
  mu[1] = 1/eta[1]
  eta[2] = delta + thetan[1]*( 1/(y[1]) - 1/(mu[1]) ) + thetan[2]*( 1/(y0) - 1/(mu0) )
  mu[2] = 1/eta[2]
  
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  lk[2] = (L-1)*log(y[2]) + L*log(L) - log( gamma(L)  ) - L*(  y[2]/mu[2] )  - L*log(mu[2])
  
  for(t in (q+1):n){
    eta[t] = delta + thetan[1]*( 1/(y[t-1]) - 1/(mu[t-1]) )  + thetan[2]*( 1/(y[t-2]) - 1/(mu[t-2]) )
    mu[t] = 1/eta[t]
  }
  t<-(q+1):n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L) ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return(sum(lk[(q+1):n]))
  
}
lk.gamma.ARMA11_inv = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  delta = par[1]; phin = par[2]; thetan=par[3]
  mu=c();v=c();lk=c();eta=c()
  
  mu0 = rep(1,p)
  y0 = rep( mean(y), p )
  eta[1] = delta + phin*1/(y0) + thetan*( 1/(y0) - 1/(mu0) )
  mu[1] = 1/eta[1]
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in (p+1):n){
    eta[t] = delta + phin*(1/(y[t-1])) + thetan*( 1/(y[t-1]) - 1/(mu[t-1]) ) 
    mu[t] = 1/eta[1]
  }
  
  
  t <- (p+1):n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return(sum(lk[2:n]))
  
}

lk.gamma.AR1_sqrt = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  p=1
  n=length(y)
  delta = par[1]; phin=par[2]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  y0 = rep(mean(y),p)
  eta=c()
  eta[1] = delta + phin*( sqrt(y0) )
  mu[1] = eta[1]^2
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  t<-2:n
  eta[t] = delta + phin*( sqrt(y[t-1]) ) 
  mu[t] = eta[t]^2
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  
  return(sum(lk[2:n]))
  
}
lk.gamma.AR2_sqrt = function(y,par,L){
  
  ## ---> Input
  ## y: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=2
  n=length(y)
  delta = par[1]; phi=par[2:3]
  mu=c();v=c();lk=c()
  y0 = mean(y)
  
  eta=c()
  eta[1] = delta + phi[1]*( sqrt(y0) ) + phi[2]*( sqrt(y0) )
  eta[2] = delta + phi[1]*( sqrt(y[1]) ) + phi[2]*( sqrt(y0) )
  mu[1] = eta[1]^2
  mu[2] = eta[2]^2
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1]) 
  lk[2] = (L-1)*log(y[2]) + L*log(L) - log( gamma(L)  ) - L*(  y[2]/mu[2] )  - L*log(mu[2]) 
  
  t <- (p+1):n
  eta[t] = delta + phi[1]*( sqrt(y[t-1]) ) + phi[2]*( sqrt(y[t-2]) )
  mu[t] = eta[t]^2
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] )  - L*log(mu[t]) 
  
  
  return(sum(lk[3:n]))
  
}
lk.gamma.MA1_sqrt = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 1
  n = length(y)
  delta = par[1]; thetan=par[2]
  mu=c();lk=c()
  
  mu0 = rep(1,q)
  y0 = rep(mean(y),q)
  eta=c()
  eta[1] = delta + thetan*( sqrt(y0) - sqrt(mu0) )
  mu[1] = eta[1]^2
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in 2:n){
    eta[t] = delta + thetan*( sqrt(y[t-1]) - sqrt(mu[t-1]) )
    mu[t] = eta[t]^2
  }
  
  t <- 2:n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L) ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return( sum(lk[2:n]) )
  
}
lk.gamma.MA2_sqrt = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 2
  n = length(y)
  delta = par[1]; thetan=par[2:3]
  mu=c();v=c();lk=c()
  
  mu0 = 1
  y0 = mean(y)
  eta=c()
  eta[1] = delta + thetan[1]*( sqrt(y0) - sqrt(mu0) ) + thetan[2]*( sqrt(y0) - sqrt(mu0) )
  mu[1] = eta[1]^2
  eta[2] = delta + thetan[1]*( sqrt(y[1]) - sqrt(mu[1]) ) + thetan[2]*( sqrt(y0) - sqrt(mu0) )
  mu[2] = eta[2]^2
  
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  lk[2] = (L-1)*log(y[2]) + L*log(L) - log( gamma(L)  ) - L*(  y[2]/mu[2] )  - L*log(mu[2])
  
  for(t in (q+1):n){
    eta[t] = delta + thetan[1]*( sqrt(y[t-1]) - sqrt(mu[t-1]) )  + thetan[2]*( sqrt(y[t-2]) - sqrt(mu[t-2]) )
    mu[t] = eta[t]^2
  }
  t<-(q+1):n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L) ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return(sum(lk[(q+1):n]))
  
}
lk.gamma.ARMA11_sqrt = function(y,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  p=1
  n=length(y)
  delta = par[1]; phin = par[2]; thetan=par[3]
  mu=c();v=c();lk=c();eta=c()
  
  mu0 = rep(1,p)
  y0 = rep( mean(y), p )
  eta[1] = delta + phin*sqrt(y0) + thetan*( sqrt(y0) - sqrt(mu0) )
  mu[1] = eta[1]^2
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in (p+1):n){
    eta[t] = delta + phin*(sqrt(y[t-1])) + thetan*( sqrt(y[t-1]) - sqrt(mu[t-1]) ) 
    mu[t] = eta[t]^2
  }
  
  
  t <- (p+1):n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return(sum(lk[2:n]))
  
}

#### chute inicial para GAMA ----

gChute.AR1_log <- function(z,L){
    t<-2:length(z)
    X = cbind( rep(1,length(z)), c(1,log(z[t-1]))  ) 
    param = ( solve(t(X)%*%X) )%*%t(X)%*%log(z)

  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  chute = c(delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
gChute.AR2_log <- function(z,L){
  
  t<-3:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),log(z[1]),log(z[t-1])), c(mean(z),mean(z),log(z[t-2]))   ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(z)
  
  chute = c(delta.hat = param[1], phi.hat = param[2:3] )
  return(chute)
  
}
gChute.MA1_log <- function(z,L){

  param = c( mean(log(z)), 0  )
  
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  
  chute = c( delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
gChute.MA2_log <- function(z,L){

  param = c( mean(log(z)), 0, 0  )
  
  
  chute = c(delta.hat = param[1], theta.hat = param[2:3] )
  return(chute)
  
}
gChute.ARMA11_log <- function(z,L){
  
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(1,log(z[t-1]))  ) 
  betan = ( solve(t(X)%*%X) )%*%t(X)%*%log(z)
  
  param = c(betan,0)
  
  if( param[2] < -1 ){
    param[2] = -0.8  
  }else if( param[2] > 1 ){
    param[2] = 0.8
  }
  if( param[3] < -1 ){
    param[3] = -0.8  
  }else if( param[3] > 1 ){
    param[3] = 0.8
  }
  
  chute = c( delta.hat = param[1], phi.hat = param[2], theta.hat=param[3] )
  return(chute)
  
} 

gChute.AR1_W <- function(z,L){
  require(lamW)
  # to find optimal lambda
  eta=c();z0=mean(z)
  
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(1,log(lambertW0(z[t-1])))  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(lambertW0(z))
  
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  chute = c(delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
gChute.AR2_W <- function(z,L){
  
  t<-3:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),log(lambertW0(z[1])),log(lambertW0(z[t-1]))), c(mean(z),mean(z),log(lambertW0(z[t-2])))   ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(lambertW0(z))
  
  chute = c(delta.hat = param[1], phi.hat = param[2:3] )
  return(chute)
  
}
gChute.MA1_W <- function(z,L){
  
  param = c( mean(log(z)), 0  )
  
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  
  chute = c( delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
gChute.MA2_W <- function(z,L){
  
  param = c( mean(log(z)), 0, 0  )
  
  
  chute = c(delta.hat = param[1], theta.hat = param[2:3] )
  return(chute)
  
}
gChute.ARMA11_W <- function(z,L){
  
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),log(lambertW0(z[t-1])))  ) 
  betan =( solve(t(X)%*%X) )%*%t(X)%*%log(lambertW0(z))
  
  param = c(betan,0)
  
  if( param[2] < -1 ){
    param[2] = -0.8  
  }else if( param[2] > 1 ){
    param[2] = 0.8
  }
  if( param[3] < -1 ){
    param[3] = -0.8  
  }else if( param[3] > 1 ){
    param[3] = 0.8
  }
  
  chute = c( delta.hat = param[1], phi.hat = param[2], theta.hat=param[3] )
  return(chute)
  
} 

gChute.AR1_inv <- function(z,L){
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(1,1/(z[t-1]))  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%(1/(z))
  
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  chute = c(delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
gChute.AR2_inv <- function(z,L){
  
  t<-3:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),1/(z[1]),1/(z[t-1])), c(mean(z),mean(z),1/(z[t-2]))   ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%(1/(z))
  
  chute = c(delta.hat = param[1], phi.hat = param[2:3] )
  return(chute)
  
}
gChute.MA1_inv <- function(z,L){
  
  param = c( mean(1/(z)), 0  )
  
  chute = c( delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
gChute.MA2_inv <- function(z,L){
  
  param = c( mean(1/(z)), 0, 0  )
  
  
  chute = c(delta.hat = param[1], theta.hat = param[2:3] )
  return(chute)
  
}
gChute.ARMA11_inv <- function(z,L){
  
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(1,1/(z[t-1]))  ) 
  betan = ( solve(t(X)%*%X) )%*%t(X)%*%(1/(z))
  
  param = c(betan,0)
  
  if( param[2] < -0.8 ){
    param[2] = -0.8  
  }else if( param[2] > 0.8 ){
    param[2] = 0.8
  }
  if( param[3] < -0.8 ){
    param[3] = -0.8  
  }else if( param[3] > 0.8 ){
    param[3] = 0.8
  }
  
  chute = c( delta.hat = param[1], phi.hat = param[2], theta.hat=param[3] )
  return(chute)
  
} 


gChute.AR1_sqrt <- function(z,L){
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(1,sqrt(z[t-1]))  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%(sqrt(z))
  
  
  if( param[2] < -0.9 ){
    param[2] = -0.9  
  }else if( param[2] > 0.9 ){
    param[2] = 0.9
  }
  
  chute = c( delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
gChute.AR2_sqrt <- function(z,L){
  
  t<-3:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),sqrt(z[1]),sqrt(z[t-1])), c(mean(z),mean(z),sqrt(z[t-2]))   ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%(sqrt(z))
  
  chute = c(delta.hat = param[1], phi.hat = param[2:3] )
  return(chute)
  
}
gChute.MA1_sqrt <- function(z,L){
  
  param = c( mean(sqrt(z)), 0  )
  
  chute = c( delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
gChute.MA2_sqrt <- function(z,L){
  
  param = c( mean(sqrt(z)), 0, 0  )
  
  
  chute = c(delta.hat = param[1], theta.hat = param[2:3] )
  return(chute)
  
}
gChute.ARMA11_sqrt <- function(z,L){
  
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(1,sqrt(z[t-1]))  ) 
  betan = ( solve(t(X)%*%X) )%*%t(X)%*%(sqrt(z))
  
  param = c(betan,0)
  
  if( param[2] < -0.8 ){
    param[2] = -0.8  
  }else if( param[2] > 0.8 ){
    param[2] = 0.8
  }
  if( param[3] < -0.8 ){
    param[3] = -0.8  
  }else if( param[3] > 0.8 ){
    param[3] = 0.8
  }
  
  chute = c( delta.hat = param[1], phi.hat = param[2], theta.hat=param[3] )
  return(chute)
  
} 


### chute inicial para G0 ----


g0Chute.AR1_log <- function(z,L){
  

  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),log(z[t-1]))  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(z)
  
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) 
    
    mu[1] = exp(  delta + phi*log( z0 ) )
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      mu[t] = exp( delta + phi*log( z[t-1] ) )
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  alpha.hat = alpha.hat$root
  
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
g0Chute.AR2_log <- function(z,L){
  
  t<-3:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),log(z[1]),log(z[t-1])), c(mean(z),mean(z),log(z[t-2]))   ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(z)
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2:length(param)]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) 
    
    mu[1] = exp(  delta + phi[1]*log( z0 ) + phi[2]*log(z0) )
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    mu[2] = exp(  delta + phi[1]*log( z[1] ) + phi[2]*log(z0) )
    v[2] = ( mu[2]*(-alpha-1) + L*z[2] )
    
    
    for(t in 3:n){
      mu[t] = exp( delta + phi[1]*log( z[t-1] ) + phi[2]*log( z[t-2] ) )
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, c(-50,-1.1),tol = 0.00001, z=z, param = param,L=L,extendInt = "yes")
  alpha.hat = alpha.hat$root
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  if( param[3] < -1 ){
    param[3] = -0.9  
  }else if( param[3] > 1 ){
    param[3] = 0.9
  }
  
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], phi.hat = param[2:length(param)] )
  return(chute)
  
}
g0Chute.MA1_log <- function(z,L){
  
  
  param = c( mean(log(z)), 0  )
  
  f <- function(alpha,z,param,L){
    delta = param[1]; theta = param[2]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    
    mu[1] = exp(  delta + theta*( log( z0 ) - log(mu0)  ) )
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      mu[t] = exp(  delta + theta*( log( z[t-1] ) - log(mu[t-1])  ) )
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  alpha.hat = alpha.hat$root
  
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat =alpha.hat, delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
g0Chute.MA2_log <- function(z,L){
  
  
  param = c( mean(log(z)), 0, 0  )

  f <- function(alpha,z,param,L){
    delta = param[1]; theta = param[2:3]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    
    mu[1] = exp(  delta + theta[1]*( log( z0 ) - log(1)  ) + theta[2]*( log( z0 ) - log(1)  ) )
    mu[2] = exp(  delta + theta[1]*( log( z[1] ) - log(mu[1])  ) + theta[2]*( log( z0 ) - log(1)  ) )
    
    v[1] = mu[1]*(-alpha-1) + L*z[1] 
    v[2] = mu[2]*(-alpha-1) + L*z[2] 
    
    for(t in 3:n){
      mu[t] = exp(  delta + theta[1]*( log( z[t-1] ) - log(mu[t-1])  ) + theta[2]*( log( z[t-2] ) - log(mu[t-2])  ) )
      v[t] = mu[t]*(-alpha-1) + L*z[t]
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  alpha.hat = alpha.hat$root
  
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], theta.hat = param[2:3] )
  return(chute)
  
}
g0Chute.ARMA11_log <- function(z,L){

  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(median(z),log(z[t-1]))  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(z)
  param=c(param,0)
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2];theta = param[3]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    
    mu[1] = exp(  delta + phi*log(z0) + theta*( log( z0 ) - log(mu0)  ) )
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      mu[t] = exp(  delta + + phi*log(z[t-1])  + theta*( log( z[t-1] ) - log(mu[t-1])  ) )
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  alpha.hat = alpha.hat$root
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], phi.hat = param[2],theta.hat = param[3] )
  return(chute)
  
}

g0Chute.AR1_W <- function(z,L){
  
  require(lamW)
  # to find optimal lambda
  z0=mean(z)
  
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(z0,log(lambertW0(z[t-1])))  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(lambertW0(z))
  
  
  if( param[2] <= -1 ){
    param[2] = -0.9  
  }else if( param[2] >= 1 ){
    param[2] = 0.9
  }
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) 
    
    eta=c()
    
    eta[1] = delta + phi*( log(lambertW0(z0)) ) 
    mu[1] = exp( eta[1] + exp(eta[1])  )  
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      eta[t] = delta + phi*( log(lambertW0(z[t-1])) ) 
      mu[t] = exp( eta[t] + exp(eta[t])  )  
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  alpha.hat = alpha.hat$root
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat =  alpha.hat   , delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
g0Chute.AR2_W <- function(z,L){
  require(lamW)
  
  t<-3:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),log(lambertW0(z[1])),log(lambertW0(z[t-1]))), c(mean(z),mean(z),log(lambertW0(z[t-2])))   ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%log(lambertW0(z))
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2:length(param)]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) 
    eta=c()
    
    eta[1] = delta + phi[1]*( log(lambertW0(z0)) ) + phi[2]*( log(lambertW0(z0)) ) 
    mu[1] = exp( eta[1] + exp(eta[1])  )  
    eta[2] = delta + phi[1]*( log(lambertW0(z[1])) ) + phi[2]*( log(lambertW0(z0)) )
    mu[2] = exp( eta[2] + exp(eta[2])  )  
    
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    v[2] = ( mu[2]*(-alpha-1) + L*z[2] )
    
    
    for(t in 3:n){
      eta[t] = delta + phi[1]*( log(lambertW0(z[t-1])) ) + phi[2]*( log(lambertW0(z[t-2])) ) 
      mu[t] = exp( eta[t] + exp(eta[t])  ) 
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, c(-50,-1.1),tol = 0.00001, z=z, param = param,L=L,extendInt = "yes")
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  if( param[3] < -1 ){
    param[3] = -0.9  
  }else if( param[3] > 1 ){
    param[3] = 0.9
  }
  alpha.hat = alpha.hat$root
  
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat =  alpha.hat  , delta.hat = param[1], phi.hat = param[2:length(param)] )
  return(chute)
  
}
g0Chute.MA1_W <- function(z,L){
  require(lamW)
  
  
  param = c( mean(log(lambertW0(z))), 0  )
  
  f <- function(alpha,z,param,L){
    delta = param[1]; theta = param[2]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    eta=c()
    
    eta[1] = delta + theta*( log(lambertW0(z0)) - log(lambertW0(mu0)) ) 
    mu[1] = exp( eta[1] + exp(eta[1])  ) 
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      eta[t] = delta + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) 
      mu[t] = exp( eta[t] + exp(eta[t])  ) 
      v[t] = mu[t]*(-alpha-1) + L*z[t] 
    }
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  alpha.hat = alpha.hat$root
  
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
g0Chute.MA2_W <- function(z,L){
  require(lamW)
  
  
  param = c( mean(log(lambertW0(z))), 0, 0  )
  
  f <- function(alpha,z,param,L){
    delta = param[1]; theta = param[2:3]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    eta=c()
    
    eta[1] = delta + theta[1]*( log(lambertW0(z0)) - log(lambertW0(mu0)) ) + theta[2]*( log(lambertW0(z0)) - log(lambertW0(mu0)) ) 
    mu[1] = exp( eta[1] + exp(eta[1])  ) 
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    eta[2] = delta + theta[1]*( log(lambertW0(z[1])) - log(lambertW0(mu[1])) ) + theta[2]*( log(lambertW0(z0)) - log(lambertW0(mu0)) ) 
    mu[2] = exp( eta[2] + exp(eta[2])  ) 
    v[2] = ( mu[2]*(-alpha-1) + L*z[2] )
    
    for(t in 3:n){
      eta[t] = delta + theta[1]*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) + theta[2]*( log(lambertW0(z[t-2])) - log(lambertW0(mu[t-2])) ) 
      mu[t] = exp( eta[t] + exp(eta[t])  ) 
      v[t] = mu[t]*(-alpha-1) + L*z[t] 
    }
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  alpha.hat = alpha.hat$root
  
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], theta.hat = param[2:3] )
  return(chute)
  
}
g0Chute.ARMA11_W <- function(z,L){
  require(lamW)
  
  t<-2:length(z)
  z0=median(z)
  
  X = cbind( rep(1,length(z)), c(z0,log(lambertW0(z[t-1])))  ) 
  ar = ( solve(t(X)%*%X) )%*%t(X)%*%log(lambertW0(z))
  param = c(ar, 0)
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2];theta = param[3]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    eta=c()
    eta[1] = delta + phi*( log(lambertW0(z0))  ) + theta*( log(lambertW0(z0)) - log(lambertW0(mu0)) ) 
    mu[1] = exp( eta[1] + exp(eta[1])  ) 
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      eta[t] = delta + phi*( log(lambertW0(z[t-1]))  ) + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu[t-1])) ) 
      mu[t] = exp( eta[t] + exp(eta[t])  ) 
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  alpha.hat = alpha.hat$root
  
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], phi.hat = param[2],theta.hat = param[3] )
  return(chute)
  
}


g0Chute.AR1_inv <- function(z,L){
  
  
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),1/(z[t-1]))  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%( 1/(z) )
  
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) 
    eta=c()
    eta[1] = delta + phi*1/( z0 )
    mu[1] = 1/eta[1]
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      eta[t] = delta + phi*1/( z[t-1] )
      mu[t] = 1/eta[t]
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  
  
  chute = c( alpha.hat = ifelse(alpha.hat$root < -50, -50, alpha.hat$root   ), delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
g0Chute.AR2_inv <- function(z,L){
  
  t<-3:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),1/(z[1]),1/(z[t-1])), c(mean(z),mean(z),1/(z[t-2]))   ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%(1/(z))
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2:length(param)]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) 
    
    
    eta=c()
    eta[1] = delta + phi[1]*1/( z0 ) + phi[2]*1/(z0)
    eta[2] = delta + phi[1]*1/( z[1] ) + phi[2]*1/(z0)
    
    mu[1] = 1/eta[1]
    mu[2] = 1/eta[2]
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    v[2] = ( mu[2]*(-alpha-1) + L*z[2] )
    
    
    for(t in 3:n){
      eta[t] = delta + phi[1]*1/( z[t-1] ) + phi[2]*1/( z[t-2] )
      mu[t] = 1/eta[t]
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, c(-50,-1.1),tol = 0.00001, z=z, param = param,L=L,extendInt = "yes")
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  if( param[3] < -1 ){
    param[3] = -0.9  
  }else if( param[3] > 1 ){
    param[3] = 0.9
  }
  
  chute = c( alpha.hat = ifelse(alpha.hat$root < -50, -40, alpha.hat$root   ), delta.hat = param[1], phi.hat = param[2:length(param)] )
  return(chute)
  
}
g0Chute.MA1_inv <- function(z,L){
  
  
  param = c( mean(1/(z)), 0  )
  
  f <- function(alpha,z,param,L){
    delta = param[1]; theta = param[2]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    eta=c()
    
    eta[1] = delta + theta*( 1/( z0 ) - 1/(mu0)  )
    mu[1] = 1/eta[1]
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      eta[t] = delta + theta*( 1/( z[t-1] ) - 1/(mu[t-1])  )
      mu[t] =1/eta[t]
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  
  
  chute = c( alpha.hat = ifelse(alpha.hat$root < -50, -50,alpha.hat$root), delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
g0Chute.MA2_inv <- function(z,L){
  
  
  param = c( mean(1/(z)), 0, 0  )
  
  f <- function(alpha,z,param,L){
    delta = param[1]; theta = param[2:3]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    eta=c()
    eta[1] = delta + theta[1]*( 1/( z0 ) - 1/(1)  ) + theta[2]*( 1/( z0 ) - 1/(1)  )
    mu[1] = 1/eta[1]
    eta[2] = delta + theta[1]*( 1/( z[1] ) - 1/(mu[1])  ) + theta[2]*( 1/( z0 ) - 1/(1)  )
    mu[2] = 1/eta[2]
    
    v[1] = mu[1]*(-alpha-1) + L*z[1] 
    v[2] = mu[2]*(-alpha-1) + L*z[2] 
    
    for(t in 3:n){
      eta[t] = delta + theta[1]*( 1/( z[t-1] ) - 1/(mu[t-1])  ) + theta[2]*( 1/( z[t-2] ) - 1/(mu[t-2])  )
      mu[t] = 1/eta[t]
      v[t] = mu[t]*(-alpha-1) + L*z[t]
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  
  
  chute = c( alpha.hat = ifelse(alpha.hat$root < -50, -50,alpha.hat$root), delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
g0Chute.ARMA11_inv <- function(z,L){
  
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(median(z),1/(z[t-1]))  ) 
  betan = ( solve(t(X)%*%X) )%*%t(X)%*%(1/(z))
  param=c(param,0)
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2];theta = param[3]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    eta=c()
    eta[1] = delta + phi*1/(z0) + theta*( 1/( z0 ) - 1/(mu0)  )
    mu[1] = 1/eta[1]
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      eta[t] = delta + phi*1/(z[t-1])  + theta*( 1/( z[t-1] ) - 1/(mu[t-1])  )
      mu[t] = 1/eta[t]
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  chute = c( alpha.hat = ifelse(alpha.hat$root < -50, -50,alpha.hat$root), delta.hat = param[1], phi.hat = param[2],theta.hat = param[3] )
  return(chute)
  
}


g0Chute.AR1_sqrt <- function(z,L){
  
  
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),sqrt(z[t-1]))  ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%( sqrt(z) )
  
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) 
    eta=c()
    eta[1] = delta + phi*sqrt( z0 )
    mu[1] = eta[1]^2
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      eta[t] = delta + phi*sqrt( z[t-1] )
      mu[t] = eta[t]^2
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  alpha.hat = alpha.hat$root
  
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat =  alpha.hat, delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
g0Chute.AR2_sqrt <- function(z,L){
  
  t<-3:length(z)
  X = cbind( rep(1,length(z)), c(mean(z),sqrt(z[1]),sqrt(z[t-1])), c(mean(z),mean(z),sqrt(z[t-2]))   ) 
  param = ( solve(t(X)%*%X) )%*%t(X)%*%(sqrt(z))
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2:length(param)]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) 
    
    
    eta=c()
    eta[1] = delta + phi[1]*sqrt( z0 ) + phi[2]*sqrt(z0)
    eta[2] = delta + phi[1]*sqrt( z[1] ) + phi[2]*sqrt(z0)
    
    mu[1] =eta[1]^2
    mu[2] = eta[2]^2
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    v[2] = ( mu[2]*(-alpha-1) + L*z[2] )
    
    
    for(t in 3:n){
      eta[t] = delta + phi[1]*sqrt( z[t-1] ) + phi[2]*sqrt( z[t-2] )
      mu[t] = eta[t]^2
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, c(-50,-1.1),tol = 0.00001, z=z, param = param,L=L,extendInt = "yes")
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  if( param[3] < -1 ){
    param[3] = -0.9  
  }else if( param[3] > 1 ){
    param[3] = 0.9
  }
  
  alpha.hat = alpha.hat$root
  
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat = alpha.hat , delta.hat = param[1], phi.hat = param[2:length(param)] )
  return(chute)
  
}
g0Chute.MA1_sqrt <- function(z,L){
  
  
  param = c( mean(sqrt(z)), 0  )
  
  f <- function(alpha,z,param,L){
    delta = param[1]; theta = param[2]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    eta=c()
    
    eta[1] = delta + theta*( sqrt( z0 ) - sqrt(mu0)  )
    mu[1] = eta[1]^2
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      eta[t] = delta + theta*( sqrt( z[t-1] ) - sqrt(mu[t-1])  )
      mu[t] =eta[t]^2
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  alpha.hat = alpha.hat$root
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
g0Chute.MA2_sqrt <- function(z,L){
  
  
  param = c( mean(sqrt(z)), 0, 0  )
  
  f <- function(alpha,z,param,L){
    delta = param[1]; theta = param[2:3]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    eta=c()
    eta[1] = delta + theta[1]*( sqrt( z0 ) - sqrt(1)  ) + theta[2]*( sqrt( z0 ) - sqrt(1)  )
    mu[1] = eta[1]^2
    eta[2] = delta + theta[1]*( sqrt( z[1] ) - sqrt(mu[1])  ) + theta[2]*( sqrt( z0 ) - sqrt(1)  )
    mu[2] = eta[2]^2
    
    v[1] = mu[1]*(-alpha-1) + L*z[1] 
    v[2] = mu[2]*(-alpha-1) + L*z[2] 
    
    for(t in 3:n){
      eta[t] = delta + theta[1]*( sqrt( z[t-1] ) - sqrt(mu[t-1])  ) + theta[2]*( sqrt( z[t-2] ) - sqrt(mu[t-2])  )
      mu[t] = eta[t]^2
      v[t] = mu[t]*(-alpha-1) + L*z[t]
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  alpha.hat = alpha.hat$root
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], theta.hat = param[2:3] )
  return(chute)
  
}
g0Chute.ARMA11_sqrt <- function(z,L){
  
  t<-2:length(z)
  X = cbind( rep(1,length(z)), c(median(z),sqrt(z[t-1]))  ) 
  betan = ( solve(t(X)%*%X) )%*%t(X)%*%(sqrt(z))
  param=c(betan,0)
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2];theta = param[3]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    eta=c()
    eta[1] = delta + phi*sqrt(z0) + theta*( sqrt( z0 ) - sqrt(mu0)  )
    mu[1] = eta[1]^2
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      eta[t] = delta + phi*sqrt(z[t-1])  + theta*( sqrt( z[t-1] ) - sqrt(mu[t-1])  )
      mu[t] = eta[t]^2
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  alpha.hat = alpha.hat$root
  if( alpha.hat <= -50 ){
    alpha.hat = -40
  }else if( alpha.hat > -1.4 ){
    alpha.hat = -1.4
  }
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], phi.hat = param[2],theta.hat = param[3] )
  return(chute)
  
}


#### CML - gamma-ARMA ----

cml.gamma.AR1_log <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR1_log(z,L) 
  lk.gammanew = function(param)  lk.gamma.AR1_log(z, param,L)
  res1 <- maxLik(lk.gammanew, grad = NULL, hess = NULL, param, method = "BFGS")

  est = res1$estimate
  delta=est[1];phi=est[2]
  
  mu.hat[1] = exp(  delta + phi*log( z0 ) )
  t <- 2:n 
  mu.hat[t] = exp( delta + phi*log( z[t-1] ) )
  
  AIC = -2*lk.gamma.AR1_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.AR1_log(z, est,L) + 2*log(n)
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )

}
cml.gamma.AR2_log <- function(z,L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  

  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR2_log(z,L) 
  lk.gammanew = function(param)  lk.gamma.AR2_log(z, param,L)
  #grad.gammanew = function(param)  grad.gamma.AR2_log(z, param)
  res1 <- maxLik(lk.gammanew, grad = NULL, hess = NULL, param, method = "BFGS")
  est = res1$estimate
  delta=est[1];phi=est[2:3]
  
  mu.hat[1] = exp( delta + phi[1]*( log(z0) ) + phi[2]*( log(z0) ) )
  mu.hat[2] = exp( delta + phi[1]*( log(z[1]) ) + phi[2]*( log(z0) ) )
  
  t <- 3:n 
  mu.hat[t] = exp( delta + phi[1]*log( z[t-1] ) + phi[2]*log( z[t-2] )  )
  
  AIC = -2*lk.gamma.AR2_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.AR2_log(z, est,L) + 2*log(n)
  
 
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )
}
cml.gamma.MA1_log <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  #A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  #B1 <- c(-1,100,1,1)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA1_log(z,L) 
  lk.new = function(param)  lk.gamma.MA1_log(z, param,L)
  #grad.new = function(param)  grad.gamma.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS")
  est = res1$estimate
  delta=est[1];theta=est[2]
  
  z0 = mean(z); mu0=z0
  
  mu.hat[1] = exp( delta + theta*( log(z0) - log(mu0) ) )
  
  for(t in 2:n) mu.hat[t] = exp( delta + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) )
  
  AIC = -2*lk.gamma.MA1_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.MA1_log(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
cml.gamma.MA2_log <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(purrr)
  
  #A1 <- matrix(c(1, -1, 0,0,0,0,0,0 ), 2,4)
  #B1 <- c(-1,100)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA2_log(z,L) 
  lk.new = function(param)  lk.gamma.MA2_log(z, param,L)
  #grad.new = function(param)  grad.gamma.MA(z, param,2)
  maxlik.new = function(logLik,grad,hess,start,method)  maxLik(logLik = lk.new,grad=NULL,hess = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
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
    delta=est[1];theta=est[2:3]
    
    z0 = mean(z); mu0=z0
    mu.hat[1:2] = exp( delta + theta[1]*( log(z0) - log(mu0) ) + theta[2]*( log(z0) - log(mu0) ) )
    
    for(t in 3:n) mu.hat[t] = exp( delta + theta[1]*( log(z[t-1]) - log(mu.hat[t-1]) ) + theta[2]*( log(z[t-2]) - log(mu.hat[t-2]) ) )
    
  }
  
  AIC = -2*lk.gamma.MA2_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.MA2_log(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
cml.gamma.ARMA11_log <- function(z,L ){
  require(purrr)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  #A1 <- matrix(c(1, -1, 0,0,0,0, 0,0,0,0,0,0, 0,0,1,-1,0,0, 0,0,0,0,1,-1 ), 6,4)
  #B1 <- c(-1,100,1,1,1,1)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.ARMA11_log(z,L) 
  lk.new = function(param)  lk.gamma.ARMA11_log(z, param,L)
  #grad.new = function(param)  grad.gamma.ARMA11_log(z, param)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "CG")
  }else{
    res1 <- res1$result
  }
  
  
  est = res1$estimate
  delta=est[1];phi=est[2];theta=est[3] 
  
  z0 = mean(z); mu0=z0
  
  mu.hat[1] = exp( delta + phi*log( z0 ) + theta*( log(z0) - log(mu0) ) )
  
  for(t in 2:n) mu.hat[t] = exp( delta + phi*log( z[t-1] ) + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) )
  
  AIC = -2*lk.gamma.ARMA11_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.ARMA11_log(z, est,L) + 2*log(n)
  

  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
auto.gamma_log <- function(y, criterio,L){
  # y: serie
  # criterio: criterio de informação AIC e BIC

  ar1 = cml.gamma.AR1_log(y,L);ar2 = cml.gamma.AR2_log(y,L)
  ma1 = cml.gamma.MA1_log(y,L); ma2 = cml.gamma.MA2_log(y,L)
  arma11 = cml.gamma.ARMA11_log(y,L)
  
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

cml.gamma.AR1_W <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  require(lamW)
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR1_W(z,L) 
  lk.gammanew = function(param)  lk.gamma.AR1_W(z, param,L)
  res1 <- maxLik(lk.gammanew, grad = NULL, hess = NULL, param, method = "BFGS")
  
  est = res1$estimate
  delta=est[1];phi=est[2]
  eta=c()
  
  eta[1] = delta + phi*( log(lambertW0(mean(z))) ) 
  mu.hat[1] = exp( eta[1] + exp(eta[1])  )  
  
  t <- 2:n 
  eta[t] = delta + phi*( log(lambertW0(z[t-1])) ) 
  mu.hat[t] = exp( eta[t] + exp(eta[t])  )
  
  AIC = -2*lk.gamma.AR1_W(z, est, L) + 2*length(est) 
  BIC = -2*lk.gamma.AR1_W(z, est, L) + 2*log(n)
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )
  
}
cml.gamma.AR2_W <- function(z,L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  require(lamW)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR2_W(z,L) 
  lk.gammanew = function(param)  lk.gamma.AR2_W(z, param,L)
  #grad.gammanew = function(param)  grad.gamma.AR2_W(z, param)
  res1 <- maxLik(lk.gammanew, grad = NULL, hess = NULL, param, method = "BFGS")
  est = res1$estimate
  delta=est[1];phi=est[2:3]
  
  eta=c()
  
  eta[1] = delta + phi[1]*( log(lambertW0(mean(z))) ) + phi[2]*( log(lambertW0(mean(z))) ) 
  eta[2] = delta + phi[1]*( log(lambertW0(z[1])) ) + phi[2]*( log(lambertW0(mean(z))) )
  mu.hat[1:2] = exp( eta[1:2] + exp(eta[1:2])  )  
  
  t <- 3:n 
  eta[t] = delta + phi[1]*( log(lambertW0(z[t-1])) ) + phi[2]*( log(lambertW0(z[t-2])) )
  mu.hat[t] = exp( eta[t] + exp(eta[t])  )  
  
  AIC = -2*lk.gamma.AR2_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.AR2_W(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )
}
cml.gamma.MA1_W <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  #A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  #B1 <- c(-1,100,1,1)
  z = ifelse(z == 0, summary(z)[2], z)
  require(lamW)
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA1_W(z,L) 
  lk.new = function(param)  lk.gamma.MA1_W(z, param,L)
  #grad.new = function(param)  grad.gamma.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS")
  est = res1$estimate
  delta=est[1];theta=est[2]
  
  z0 = mean(z); mu0=z0
  eta=c()
  
  eta[1] = delta + theta*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  ) 
  mu.hat[1] = exp( eta[1] + exp(eta[1])  )  
  
  for(t in 2:n){
    eta[t] = delta + theta*( log(lambertW0(mean(z[t-1]))) - log(lambertW0(mu.hat[t-1]))  ) 
    mu.hat[t] = exp( eta[t] + exp(eta[t])  ) 
  }
  
  AIC = -2*lk.gamma.MA1_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.MA1_W(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
cml.gamma.MA2_W <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(purrr)
  require(lamW)
  #A1 <- matrix(c(1, -1, 0,0,0,0,0,0 ), 2,4)
  #B1 <- c(-1,100)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA2_W(z,L) 
  lk.new = function(param)  lk.gamma.MA2_W(z, param,L)
  #grad.new = function(param)  grad.gamma.MA(z, param,2)
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
    delta=est[1];theta=est[2:3]
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
  
  AIC = -2*lk.gamma.MA2_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.MA2_W(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
cml.gamma.ARMA11_W <- function(z,L ){
  require(purrr)
  require(lamW)
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  #A1 <- matrix(c(1, -1, 0,0,0,0, 0,0,0,0,0,0, 0,0,1,-1,0,0, 0,0,0,0,1,-1 ), 6,4)
  #B1 <- c(-1,100,1,1,1,1)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.ARMA11_W(z,L) 
  lk.new = function(param)  lk.gamma.ARMA11_W(z, param,L)
  #grad.new = function(param)  grad.gamma.ARMA11_W(z, param)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "CG")
  }else{
    res1 <- res1$result
  }
  
  
  est = res1$estimate
  delta=est[1];phi=est[2];theta=est[3] 
  z0 = mean(z); mu0=z0
  eta=c()
  
  eta[1] = delta + phi*( log(lambertW0(z0))  ) + theta*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  )
  mu.hat[1] = exp( eta[1] + exp(eta[1])  ) 
  
  for(t in 2:n){
    eta[t] = delta + phi*( log(lambertW0(z[t-1]))  ) + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu.hat[t-1]))  )
    mu.hat[t] = exp( eta[t] + exp(eta[t])  ) 
  }
  
  AIC = -2*lk.gamma.ARMA11_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.ARMA11_W(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
auto.gamma_W <- function(y, criterio,L){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  require(lamW)
  ar1 = cml.gamma.AR1_W(y,L);ar2 = cml.gamma.AR2_W(y,L)
  ma1 = cml.gamma.MA1_W(y,L); ma2 = cml.gamma.MA2_W(y,L)
  arma11 = cml.gamma.ARMA11_W(y,L)
  
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


cml.gamma.AR1_inv <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR1_inv(z,L) 
  lk.gammanew = function(param)  lk.gamma.AR1_inv(z, param,L)
  res1 <- maxLik(lk.gammanew, grad = NULL, hess = NULL, param, method = "BFGS")
  
  est = res1$estimate
  delta=est[1];phi=est[2]
  eta=c()
  eta[1] = delta + phi*1/( z0 ) 
  mu.hat[1] = 1/eta[1]
  t <- 2:n 
  eta[t] = delta + phi*1/( z[t-1] )
  mu.hat[t] = 1/eta[t]
  
  AIC = -2*lk.gamma.AR1_inv(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.AR1_inv(z, est,L) + 2*log(n)
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )
  
}
cml.gamma.AR2_inv <- function(z,L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR2_inv(z,L) 
  lk.gammanew = function(param)  lk.gamma.AR2_inv(z, param,L)
  #grad.gammanew = function(param)  grad.gamma.AR2_inv(z, param)
  res1 <- maxLik(lk.gammanew, grad = NULL, hess = NULL, param, method = "BFGS")
  est = res1$estimate
  delta=est[1];phi=est[2:3]
  eta=c()
  eta[1] = delta + phi[1]*( 1/(z0) ) + phi[2]*( 1/(z0) )
  eta[2] = delta + phi[1]*( 1/(z[1]) ) + phi[2]*( 1/(z0) )
  mu.hat[1] = 1/eta[1]
  mu.hat[2] = 1/eta[2]
  
  t <- 3:n 
  eta[t] = delta + phi[1]*1/( z[t-1] ) + phi[2]*1/( z[t-2] )
  mu.hat[t] = 1/eta[t]
  
  AIC = -2*lk.gamma.AR2_inv(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.AR2_inv(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )
}
cml.gamma.MA1_inv <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  #A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  #B1 <- c(-1,100,1,1)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA1_inv(z,L) 
  lk.new = function(param)  lk.gamma.MA1_inv(z, param,L)
  #grad.new = function(param)  grad.gamma.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS")
  est = res1$estimate
  delta=est[1];theta=est[2]
  
  z0 = mean(z); mu0=z0
  eta=c()
  eta[1] = delta + theta*( 1/(z0) - 1/(mu0) )
  mu.hat[1] = 1/eta[1]
  
  for(t in 2:n){ 
  eta[t] = delta + theta*( 1/(z[t-1]) - 1/(mu.hat[t-1]) );
  mu.hat[t] = 1/eta[t]
  }
  
  AIC = -2*lk.gamma.MA1_inv(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.MA1_inv(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
cml.gamma.MA2_inv <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(purrr)
  
  #A1 <- matrix(c(1, -1, 0,0,0,0,0,0 ), 2,4)
  #B1 <- c(-1,100)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA2_inv(z,L) 
  lk.new = function(param)  lk.gamma.MA2_inv(z, param,L)
  #grad.new = function(param)  grad.gamma.MA(z, param,2)
  maxlik.new = function(logLik,grad,hess,start,method)  maxLik(logLik = lk.new,grad=NULL,hess = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
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
    delta=est[1];theta=est[2:3]
    z0 = mean(z); mu0=z0
    eta=c()
    eta[1] = delta + theta[1]*( 1/(z0) - 1/(mu0) ) + theta[2]*( 1/(z0) - 1/(mu0) )
    mu.hat[1] = 1/eta[1]
    eta[2] = delta + theta[1]*( 1/(z[1]) - 1/(mu[1]) ) + theta[2]*( 1/(z0) - 1/(mu0) )
    mu.hat[2] = 1/eta[2]
    
    for(t in 3:n){ 
      eta[t] = delta + theta[1]*( 1/(z[t-1]) - 1/(mu.hat[t-1]) ) + theta[2]*( 1/(z[t-2]) - 1/(mu.hat[t-2]) )
      mu.hat[t] = 1/eta[t]
      } 
    
  }
  
  AIC = -2*lk.gamma.MA2_inv(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.MA2_inv(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
cml.gamma.ARMA11_inv <- function(z,L ){
  require(purrr)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  #A1 <- matrix(c(1, -1, 0,0,0,0, 0,0,0,0,0,0, 0,0,1,-1,0,0, 0,0,0,0,1,-1 ), 6,4)
  #B1 <- c(-1,100,1,1,1,1)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.ARMA11_inv(z,L) 
  lk.new = function(param)  lk.gamma.ARMA11_inv(z, param,L)
  #grad.new = function(param)  grad.gamma.ARMA11_inv(z, param)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "CG")
  }else{
    res1 <- res1$result
  }
  
  
  est = res1$estimate
  delta=est[1];phi=est[2];theta=est[3] 
  
  z0 = mean(z); mu0=z0
  eta=c()
  eta[1] = delta + phi*1/( z0 ) + theta*( 1/(z0) - 1/(mu0) )
  mu.hat[1] = 1/eta[1]
  
  for(t in 2:n){
    eta[t] = delta + phi*1/( z[t-1] ) + theta*( 1/(z[t-1]) - 1/(mu.hat[t-1]) )
    mu.hat[t] = 1/eta[t]
  } 
    
  
  AIC = -2*lk.gamma.ARMA11_inv(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.ARMA11_inv(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
auto.gamma_inv <- function(y, criterio,L){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  
  ar1 = cml.gamma.AR1_inv(y,L);ar2 = cml.gamma.AR2_inv(y,L)
  ma1 = cml.gamma.MA1_inv(y,L); ma2 = cml.gamma.MA2_inv(y,L)
  arma11 = cml.gamma.ARMA11_inv(y,L)
  
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



cml.gamma.AR1_sqrt <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR1_sqrt(z,L) 
  lk.gammanew = function(param)  lk.gamma.AR1_sqrt(z, param,L)
  res1 <- maxLik(lk.gammanew, grad = NULL, hess = NULL, param, method = "BFGS")
  
  est = res1$estimate
  delta=est[1];phi=est[2]
  eta=c()
  eta[1] = delta + phi*sqrt( z0 ) 
  mu.hat[1] = eta[1]^2
  t <- 2:n 
  eta[t] = delta + phi*sqrt( z[t-1] )
  mu.hat[t] = eta[t]^2
  
  AIC = -2*lk.gamma.AR1_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.AR1_sqrt(z, est,L) + 2*log(n)
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )
  
}
cml.gamma.AR2_sqrt <- function(z,L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR2_sqrt(z,L) 
  lk.gammanew = function(param)  lk.gamma.AR2_sqrt(z, param,L)
  #grad.gammanew = function(param)  grad.gamma.AR2_sqrt(z, param)
  res1 <- maxLik(lk.gammanew, grad = NULL, hess = NULL, param, method = "BFGS")
  est = res1$estimate
  delta=est[1];phi=est[2:3]
  eta=c()
  eta[1] = delta + phi[1]*( sqrt(z0) ) + phi[2]*( sqrt(z0) )
  eta[2] = delta + phi[1]*( sqrt(z[1]) ) + phi[2]*( sqrt(z0) )
  mu.hat[1] = eta[1]^2
  mu.hat[2] = eta[2]^2
  
  t <- 3:n 
  eta[t] = delta + phi[1]*sqrt( z[t-1] ) + phi[2]*sqrt( z[t-2] )
  mu.hat[t] = eta[t]^2
  
  AIC = -2*lk.gamma.AR2_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.AR2_sqrt(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )
}
cml.gamma.MA1_sqrt <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  #A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  #B1 <- c(-1,100,1,1)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA1_sqrt(z,L) 
  lk.new = function(param)  lk.gamma.MA1_sqrt(z, param,L)
  #grad.new = function(param)  grad.gamma.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS")
  est = res1$estimate
  delta=est[1];theta=est[2]
  
  z0 = mean(z); mu0=z0
  eta=c()
  eta[1] = delta + theta*( sqrt(z0) - sqrt(mu0) )
  mu.hat[1] = eta[1]^2
  
  for(t in 2:n){ 
    eta[t] = delta + theta*( sqrt(z[t-1]) - sqrt(mu.hat[t-1]) );
    mu.hat[t] = eta[t]^2
  }
  
  AIC = -2*lk.gamma.MA1_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.MA1_sqrt(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
cml.gamma.MA2_sqrt <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(purrr)
  
  #A1 <- matrix(c(1, -1, 0,0,0,0,0,0 ), 2,4)
  #B1 <- c(-1,100)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA2_sqrt(z,L) 
  lk.new = function(param)  lk.gamma.MA2_sqrt(z, param,L)
  #grad.new = function(param)  grad.gamma.MA(z, param,2)
  maxlik.new = function(logLik,grad,hess,start,method)  maxLik(logLik = lk.new,grad=NULL,hess = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
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
    delta=est[1];theta=est[2:3]
    z0 = mean(z); mu0=z0
    eta=c()
    eta[1] = delta + theta[1]*( sqrt(z0) - sqrt(mu0) ) + theta[2]*( sqrt(z0) - sqrt(mu0) )
    mu.hat[1] = eta[1]^2
    eta[2] = delta + theta[1]*( sqrt(z[1]) - sqrt(mu.hat[1]) ) + theta[2]*( sqrt(z0) - sqrt(mu0) )
    mu.hat[2] = eta[2]^2
    
    for(t in 3:n){ 
      eta[t] = delta + theta[1]*( sqrt(z[t-1]) - sqrt(mu.hat[t-1]) ) + theta[2]*( sqrt(z[t-2]) - sqrt(mu.hat[t-2]) )
      mu.hat[t] = eta[t]^2
    } 
    
  }
  
  AIC = -2*lk.gamma.MA2_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.MA2_sqrt(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
cml.gamma.ARMA11_sqrt <- function(z,L ){
  require(purrr)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  #A1 <- matrix(c(1, -1, 0,0,0,0, 0,0,0,0,0,0, 0,0,1,-1,0,0, 0,0,0,0,1,-1 ), 6,4)
  #B1 <- c(-1,100,1,1,1,1)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.ARMA11_sqrt(z,L) 
  lk.new = function(param)  lk.gamma.ARMA11_sqrt(z, param,L)
  #grad.new = function(param)  grad.gamma.ARMA11_sqrt(z, param)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "CG")
  }else{
    res1 <- res1$result
  }
  
  
  est = res1$estimate
  delta=est[1];phi=est[2];theta=est[3] 
  
  z0 = mean(z); mu0=z0
  eta=c()
  eta[1] = delta + phi*sqrt( z0 ) + theta*( sqrt(z0) - sqrt(mu0) )
  mu.hat[1] = eta[1]^2
  
  for(t in 2:n){
    eta[t] = delta + phi*sqrt( z[t-1] ) + theta*( sqrt(z[t-1]) - sqrt(mu.hat[t-1]) )
    mu.hat[t] = eta[t]^2
  } 
  
  
  AIC = -2*lk.gamma.ARMA11_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.gamma.ARMA11_sqrt(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
auto.gamma_sqrt <- function(y, criterio,L){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  
  ar1 = cml.gamma.AR1_sqrt(y,L);ar2 = cml.gamma.AR2_sqrt(y,L)
  ma1 = cml.gamma.MA1_sqrt(y,L); ma2 = cml.gamma.MA2_sqrt(y,L)
  arma11 = cml.gamma.ARMA11_sqrt(y,L)
  
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


#### CML - g0-ARMA ----

cml.g0.AR1_log <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(100,-1.3,0.95,0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR1_log(z,L) 
  lk.g0new = function(param)  lk.g0.AR1_log(z, param,L)
  res1 <- maxLik(lk.g0new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  code = res1$code
  est = res1$estimate
  alpha = est[1]; delta=est[2];phi=est[3]
  
  
  
  mu.hat[1] = exp(  delta + phi*log( z0 ) )
  t <- 2:n 
  mu.hat[t] = exp( delta + phi*log( z[t-1] ) )
  
  AIC = -2*lk.g0.AR1_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.AR1_log(z, est,L) + 2*log(n)
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )
  
}
cml.g0.AR2_log <- function(z,L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(100,-1.3,0.95,0.95,0.95,0.95)
  
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR2_log(z,L) 
  lk.g0new = function(param)  lk.g0.AR2_log(z, param,L)
  #grad.g0new = function(param)  grad.g0.AR2_log(z, param)
  res1 <- maxLik(lk.g0new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  est = res1$estimate
  code = res1$code
  alpha = est[1];delta=est[2];phi=est[3:4]
  
  mu.hat[1] = exp( delta + phi[1]*( log(z0) ) + phi[2]*( log(z0) ) )
  mu.hat[2] = exp( delta + phi[1]*( log(z[1]) ) + phi[2]*( log(z0) ) )
  
  t <- 3:n 
  mu.hat[t] = exp( delta + phi[1]*log( z[t-1] ) + phi[2]*log( z[t-2] )  )
  
  AIC = -2*lk.g0.AR2_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.AR2_log(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )
}
cml.g0.MA1_log <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(100,-1.3,0.95,0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA1_log(z,L) 
  lk.new = function(param)  lk.g0.MA1_log(z, param,L)
  #grad.new = function(param)  grad.g0.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  est = res1$estimate
  alpha = est[1];delta=est[2];theta=est[3]
  code = res1$code
  z0 = mean(z); mu0=z0
  
  mu0 = exp( delta + theta*( log(z0) - log(mu0) ) )
  
  mu.hat[1] = exp( delta + theta*( log(z0) - log(mu0) ) )
  
  for(t in 2:n) mu.hat[t] = exp( delta + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) )
  
  AIC = -2*lk.g0.MA1_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA1_log(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}
cml.g0.MA2_log <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0), 2,4)
  B1 <- c(100,-1.3)
  
  require(purrr)
  
  #A1 <- matrix(c(1, -1, 0,0,0,0,0,0 ), 2,4)
  #B1 <- c(-1,100)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA2_log(z,L) 
  lk.new = function(param)  lk.g0.MA2_log(z, param,L)
  #grad.new = function(param)  grad.g0.MA(z, param,2)
  maxlik.new = function(logLik,grad,hess,start,method)  maxLik(logLik = lk.new,grad=NULL,hess = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
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
    mu0 = exp( delta + theta[1]*( log(z0) - log(mu0) ) + theta[2]*( log(z0) - log(mu0) ) )
    
    mu.hat[1] = exp( delta + theta[1]*( log(z0) - log(mu0) ) + theta[2]*( log(z0) - log(mu0) ) )
    mu.hat[2] = exp( delta + theta[1]*( log(z[1]) - log(mu.hat[1]) ) + theta[2]*( log(z0) - log(mu0) ) )
    
    for(t in 3:n) mu.hat[t] = exp( delta + theta[1]*( log(z[t-1]) - log(mu.hat[t-1]) ) + theta[2]*( log(z[t-2]) - log(mu.hat[t-2]) ) )
    
  }
  
  AIC = -2*lk.g0.MA2_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA2_log(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}
cml.g0.ARMA11_log <- function(z,L ){
  require(purrr)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(100,-1.3,0.95,0.95,0.95,0.95)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.ARMA11_log(z,L) 
  lk.new = function(param)  lk.g0.ARMA11_log(z, param,L)
  #grad.new = function(param)  grad.g0.ARMA11_log(z, param)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "CG", constraints=list(ineqA = A1, ineqB = B1) )
  }else{
    res1 <- res1$result
  }
  
  code = res1$code
  est = res1$estimate
  delta=est[2];phi=est[3];theta=est[4] 
  
  z0 = mean(z); mu0=z0
  mu0 = exp( delta + phi*log( z0 ) + theta*( log(z0) - log(mu0) ) )
  
  mu.hat[1] = exp( delta + phi*log( z0 ) + theta*( log(z0) - log(mu0) ) )
  
  for(t in 2:n) mu.hat[t] = exp( delta + phi*log( z[t-1] ) + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) )
  
  AIC = -2*lk.g0.ARMA11_log(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.ARMA11_log(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}
auto.g0_log <- function(y, criterio,L){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  
  ar1 = cml.g0.AR1_log(y,L);ar2 = cml.g0.AR2_log(y,L)
  ma1 = cml.g0.MA1_log(y,L); ma2 = cml.g0.MA2_log(y,L)
  arma11 = cml.g0.ARMA11_log(y,L)
  
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

cml.g0.AR1_W <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(100,-1.3,0.95,0.95)
  
  
  require(lamW)
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR1_W(z,L) 
  lk.g0new = function(param)  lk.g0.AR1_W(z, param,L)
  res1 <- maxLik(lk.g0new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  code = res1$code
  est = res1$estimate
  delta=est[2];phi=est[3]
  eta=c()
  
  eta[1] = delta + phi*( log(lambertW0(mean(z))) ) 
  mu.hat[1] = exp( eta[1] + exp(eta[1])  )  
  
  t <- 2:n 
  eta[t] = delta + phi*( log(lambertW0(z[t-1])) ) 
  mu.hat[t] = exp( eta[t] + exp(eta[t])  )
  
  AIC = -2*lk.g0.AR1_W(z, est, L) + 2*length(est) 
  BIC = -2*lk.g0.AR1_W(z, est, L) + 2*log(n)
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )
  
}
cml.g0.AR2_W <- function(z,L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  require(lamW)
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(100,-1.3,0.95,0.95,0.95,0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR2_W(z,L) 
  lk.g0new = function(param)  lk.g0.AR2_W(z, param,L)
  #grad.g0new = function(param)  grad.g0.AR2_W(z, param)
  res1 <- maxLik(lk.g0new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
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
  
  AIC = -2*lk.g0.AR2_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.AR2_W(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )
}
cml.g0.MA1_W <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(100,-1.3,0.95,0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  require(lamW)
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA1_W(z,L) 
  lk.new = function(param)  lk.g0.MA1_W(z, param,L)
  #grad.new = function(param)  grad.g0.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  est = res1$estimate
  delta=est[2];theta=est[3]
  code = res1$code
  z0 = mean(z); mu0=z0
  eta=c()
  
  eta0 = delta + theta*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  ) 
  mu0 = exp( eta0 + exp(eta0)  )  
  
  eta[1] = delta + theta*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  ) 
  mu.hat[1] = exp( eta[1] + exp(eta[1])  )  
  
  for(t in 2:n){
    eta[t] = delta + theta*( log(lambertW0(mean(z[t-1]))) - log(lambertW0(mu.hat[t-1]))  ) 
    mu.hat[t] = exp( eta[t] + exp(eta[t])  ) 
  }
  
  AIC = -2*lk.g0.MA1_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA1_W(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}
cml.g0.MA2_W <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(purrr)
  require(lamW)
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0), 2,4)
  B1 <- c(100,-1.3)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA2_W(z,L) 
  lk.new = function(param)  lk.g0.MA2_W(z, param,L)
  #grad.new = function(param)  grad.g0.MA(z, param,2)
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
    
    eta0 = delta + theta[1]*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  ) + theta[2]*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  )
    mu0 = exp( eta0 + exp(eta0)  ) 
    
    eta[1] = delta + theta[1]*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  ) + theta[2]*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  )
    mu.hat[1] = exp( eta[1] + exp(eta[1])  ) 
    eta[2] = delta + theta[1]*( log(lambertW0(z[1])) - log(lambertW0(mu.hat[1]))  ) + theta[2]*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  )
    mu.hat[2] = exp( eta[2] + exp(eta[2])  ) 
    
    
    for(t in 3:n){
      eta[t] = delta + theta[1]*( log(lambertW0(z[t-1])) - log(lambertW0(mu.hat[t-1]))  ) + theta[2]*( log(lambertW0(z[t-2])) - log(lambertW0(mu.hat[t-2]))  )
      mu.hat[t] = exp( eta[t] + exp(eta[t])  ) 
    } 
  }
  
  AIC = -2*lk.g0.MA2_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA2_W(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC,code = code) )  
}
cml.g0.ARMA11_W <- function(z,L ){
  require(purrr)
  require(lamW)
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(100,-1.3,0.95,0.95,0.95,0.95)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = median(z)
  n = length(z)
  
  param = g0Chute.ARMA11_W(z,L) 
  lk.new = function(param)  lk.g0.ARMA11_W(z, param,L)
  #grad.new = function(param)  grad.g0.ARMA11_W(z, param)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "CG", constraints=list(ineqA = A1, ineqB = B1) )
  }else{
    res1 <- res1$result
  }
  
  code = res1$code
  est = res1$estimate
  delta=est[2];phi=est[3];theta=est[4] 
  z0 = mean(z); mu0=z0
  eta=c()
  
  eta0 = delta + phi*( log(lambertW0(z0))  ) + theta*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  )
  mu0 = exp( eta0 + exp(eta0)  ) 
  
  eta[1] = delta + phi*( log(lambertW0(z0))  ) + theta*( log(lambertW0(mean(z))) - log(lambertW0(mu0))  )
  mu.hat[1] = exp( eta[1] + exp(eta[1])  ) 
  
  for(t in 2:n){
    eta[t] = delta + phi*( log(lambertW0(z[t-1]))  ) + theta*( log(lambertW0(z[t-1])) - log(lambertW0(mu.hat[t-1]))  )
    mu.hat[t] = exp( eta[t] + exp(eta[t])  ) 
  }
  
  AIC = -2*lk.g0.ARMA11_W(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.ARMA11_W(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}
auto.g0_W <- function(y, criterio,L){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  require(lamW)
  ar1 = cml.g0.AR1_W(y,L);ar2 = cml.g0.AR2_W(y,L)
  ma1 = cml.g0.MA1_W(y,L); ma2 = cml.g0.MA2_W(y,L)
  arma11 = cml.g0.ARMA11_W(y,L)
  
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

      
cml.g0.AR1_inv <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(100,-1.3,1,1)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR1_inv(z,L) 
  lk.g0new = function(param)  lk.g0.AR1_inv(z, param,L)
  res1 <- maxLik(lk.g0new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  
  est = res1$estimate
  alpha = est[1]; delta=est[2];phi=est[3]
  
  eta=c()
  eta[1] = delta + phi*1/( z0 )
  mu.hat[1] = 1/eta[1]
  t <- 2:n 
  eta[t] = delta + phi*1/( z[t-1] ) 
  mu.hat[t] = 1/eta[t]
  
  AIC = -2*lk.g0.AR1_inv(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.AR1_inv(z, est,L) + 2*log(n)
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )
  
}
cml.g0.AR2_inv <- function(z,L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(100,-1.3,1,1,1,1)
  
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR2_inv(z,L) 
  lk.g0new = function(param)  lk.g0.AR2_inv(z, param,L)
  #grad.g0new = function(param)  grad.g0.AR2_inv(z, param)
  res1 <- maxLik(lk.g0new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  est = res1$estimate
  alpha = est[1];delta=est[2];phi=est[3:4]
  eta=c()
  eta[1] = delta + phi[1]*( 1/(z0) ) + phi[2]*( 1/(z0) )
  eta[2] = delta + phi[1]*( 1/(z[1]) ) + phi[2]*( 1/(z0) )
  mu.hat[1] = 1/eta[1]
  mu.hat[2] = 1/eta[2]
  
  t <- 3:n 
  eta[t] =  delta + phi[1]*1/( z[t-1] ) + phi[2]*1/( z[t-2] )
  mu.hat[t] = 1/eta[t]
  
  AIC = -2*lk.g0.AR2_inv(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.AR2_inv(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )
}
cml.g0.MA1_inv <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(100,-1.3,1,1)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA1_inv(z,L) 
  lk.new = function(param)  lk.g0.MA1_inv(z, param,L)
  #grad.new = function(param)  grad.g0.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  est = res1$estimate
  alpha = est[1];delta=est[2];theta=est[3]
  
  z0 = mean(z); mu0=z0
  eta=c()
  eta[1] = delta + theta*( 1/(z0) - 1/(mu0) )
  mu.hat[1] =1/eta[1]
  
  for(t in 2:n){
    eta[t] = delta + theta*( 1/(z[t-1]) - 1/(mu.hat[t-1]) )
    mu.hat[t] = 1/eta[t]
  }
  
  AIC = -2*lk.g0.MA1_inv(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA1_inv(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
cml.g0.MA2_inv <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0), 2,4)
  B1 <- c(100,-1.3)
  
  require(purrr)
  
  #A1 <- matrix(c(1, -1, 0,0,0,0,0,0 ), 2,4)
  #B1 <- c(-1,100)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA2_inv(z,L) 
  lk.new = function(param)  lk.g0.MA2_inv(z, param,L)
  #grad.new = function(param)  grad.g0.MA(z, param,2)
  maxlik.new = function(logLik,grad,hess,start,method)  maxLik(logLik = lk.new,grad=NULL,hess = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
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
    
    z0 = mean(z); mu0=z0
    eta=c()
    eta[1:2] = delta + theta[1]*( 1/(z0) - 1/(mu0) ) + theta[2]*( 1/(z0) - 1/(mu0) )
    mu.hat[1:2] = 1/eta[1:2]
    
    for(t in 3:n){
      eta[t] = delta + theta[1]*( 1/(z[t-1]) - 1/(mu.hat[t-1]) ) + theta[2]*( 1/(z[t-2]) - 1/(mu.hat[t-2]) ) 
      mu.hat[t] = 1/eta[t]
      
    } 
  }
  
  AIC = -2*lk.g0.MA2_inv(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA2_inv(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
cml.g0.ARMA11_inv <- function(z,L ){
  require(purrr)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(300,-1.3,1,1,1,1)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = median(z)
  n = length(z)
  
  param = g0Chute.ARMA11_inv(z,L) 
  lk.new = function(param)  lk.g0.ARMA11_inv(z, param,L)
  #grad.new = function(param)  grad.g0.ARMA11_inv(z, param)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "CG", constraints=list(ineqA = A1, ineqB = B1) )
  }else{
    res1 <- res1$result
  }
  
  
  est = res1$estimate
  delta=est[2];phi=est[3];theta=est[4] 
  
  z0 = median(z); mu0=1
  eta=c()
  eta[1] = delta + phi*1/( z0 ) + theta*( 1/(z0) - 1/(mu0) )
  mu.hat[1] = 1/eta[1]
  
  for(t in 2:n){
    eta[t] = delta + phi*1/( z[t-1] ) + theta*( 1/(z[t-1]) - 1/(mu.hat[t-1]) )
    mu.hat[t] =1/eta[t]
  }
  
  AIC = -2*lk.g0.ARMA11_inv(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.ARMA11_inv(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC) )  
}
auto.g0_inv <- function(y, criterio,L){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  
  ar1 = cml.g0.AR1_inv(y,L);ar2 = cml.g0.AR2_inv(y,L)
  ma1 = cml.g0.MA1_inv(y,L); ma2 = cml.g0.MA2_inv(y,L)
  arma11 = cml.g0.ARMA11_inv(y,L)
  
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

cml.g0.AR1_sqrt <- function(z, L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(100,-1.3,0.95,0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR1_sqrt(z,L) 
  lk.g0new = function(param)  lk.g0.AR1_sqrt(z, param,L)
  res1 <- maxLik(lk.g0new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  code = res1$code
  est = res1$estimate
  alpha = est[1]; delta=est[2];phi=est[3]

  
  eta=c()
  eta[1] = delta + phi*sqrt( z0 )
  mu.hat[1] = eta[1]^2
  t <- 2:n 
  eta[t] = delta + phi*sqrt( z[t-1] ) 
  mu.hat[t] = eta[t]^2
  
  AIC = -2*lk.g0.AR1_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.AR1_sqrt(z, est,L) + 2*log(n)
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )
  
}
cml.g0.AR2_sqrt <- function(z, L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(100,-1.3,0.95,0.95,0.95,0.95)
  
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR2_sqrt(z,L) 
  lk.g0new = function(param)  lk.g0.AR2_sqrt(z, param,L)
  #grad.g0new = function(param)  grad.g0.AR2_sqrt(z, param)
  res1 <- maxLik(lk.g0new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
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
  
  AIC = -2*lk.g0.AR2_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.AR2_sqrt(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC,code=code) )
}
cml.g0.MA1_sqrt <- function(z, L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(100,-1.3,0.95,0.95)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA1_sqrt(z,L) 
  lk.new = function(param)  lk.g0.MA1_sqrt(z, param,L)
  #grad.new = function(param)  grad.g0.MA1(z, param)
  res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  est = res1$estimate
  code = res1$code
  alpha = est[1];delta=est[2];theta=est[3]
  
  z0 = mean(z); mu0=z0
  eta=c()
  eta0 = delta + theta*( sqrt(z0) - sqrt(mu0) )
  mu0 =eta0^2
  
  eta[1] = delta + theta*( sqrt(z0) - sqrt(mu0) )
  mu.hat[1] =eta[1]^2
  
  for(t in 2:n){
    eta[t] = delta + theta*( sqrt(z[t-1]) - sqrt(mu.hat[t-1]) )
    mu.hat[t] = eta[t]^2
  }
  
  AIC = -2*lk.g0.MA1_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA1_sqrt(z, est,L) + 2*log(n)
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC,code = code) )  
}
cml.g0.MA2_sqrt <- function(z, L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(purrr)
  
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0), 2,4)
  B1 <- c(100,-1.3)
  
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA2_sqrt(z,L) 
  lk.new = function(param)  lk.g0.MA2_sqrt(z, param,L)
  #grad.new = function(param)  grad.g0.MA(z, param,2)
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
    #eta0 = delta + theta[1]*( sqrt(z0) - sqrt(mu0) ) + theta[2]*( sqrt(z0) - sqrt(mu0) )
    #mu0 = eta0^2
    
    eta[1] = delta + theta[1]*( sqrt(z0) - sqrt(mu0) ) + theta[2]*( sqrt(z0) - sqrt(mu0) )
    mu.hat[1] = eta[1]^2
    
    eta[2] = delta + theta[1]*( sqrt(z[1]) - sqrt(mu.hat[1]) ) + theta[2]*( sqrt(z0) - sqrt(mu0) )
    mu.hat[2] = eta[2]^2
    
    for(t in 3:n){
      eta[t] = delta + theta[1]*( sqrt(z[t-1]) - sqrt(mu.hat[t-1]) ) + theta[2]*( sqrt(z[t-2]) - sqrt(mu.hat[t-2]) ) 
      mu.hat[t] = eta[t]^2
      
    } 
  }
  
  AIC = -2*lk.g0.MA2_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA2_sqrt(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC,code = code) )  
}
cml.g0.ARMA11_sqrt <- function(z, L ){
  require(purrr)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0, 0,0, 0,0, 0,0, 0,0,1,-1,0,0,0,0,0,0,1,-1), 6,4)
  B1 <- c(100,-1.3,0.95,0.95,0.95,0.95)
  z = ifelse(z == 0, summary(z)[2], z)
  
  mu.hat = c();z0 = median(z)
  n = length(z)
  
  param = g0Chute.ARMA11_sqrt(z,L) 
  lk.new = function(param)  lk.g0.ARMA11_sqrt(z, param,L)
  #grad.new = function(param)  grad.g0.ARMA11_sqrt(z, param)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS")
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = NULL, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1) )
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "CG")
  }else{
    res1 <- res1$result
  }
  
  code = res1$code
  est = res1$estimate
  delta=est[2];phi=est[3];theta=est[4] 
  
  z0 = mean(z); mu0=z0
  eta=c()
  eta0 = delta + phi*sqrt( z0 ) + theta*( sqrt(z0) - sqrt(mu0) )
  mu0 = eta0^2
  
  eta[1] = delta + phi*sqrt( z0 ) + theta*( sqrt(z0) - sqrt(mu0) )
  mu.hat[1] = eta[1]^2
  
  for(t in 2:n){
    eta[t] = delta + phi*sqrt( z[t-1] ) + theta*( sqrt(z[t-1]) - sqrt(mu.hat[t-1]) )
    mu.hat[t] =eta[t]^2
  }
  
  AIC = -2*lk.g0.ARMA11_sqrt(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.ARMA11_sqrt(z, est,L) + 2*log(n)
  
  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, code = code) )  
}
auto.g0_sqrt <- function(y, criterio, L){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  
  ar1 = cml.g0.AR1_sqrt(y,L);ar2 = cml.g0.AR2_sqrt(y,L)
  ma1 = cml.g0.MA1_sqrt(y,L); ma2 = cml.g0.MA2_sqrt(y,L)
  arma11 = cml.g0.ARMA11_sqrt(y,L)
  
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
  

forecast.AR1 <- function(y,par,h,funcLink){
  #y: vetor
  #par: parametros estimados
  #h: passos para previsao
  #funcLink: funcao de ligacao
  
  k=1:h
  delta=par[1];phi = par[2]
  n=length(y)
  
  if( funcLink == "log" ){
  eta = delta + phi*log(y[n+k-1])
  muh = exp(eta)
  }else if(funcLink == "Wlamb"){
    require(lamW)
    eta= delta + phi*( log( lambertW0( y[n+k-1] )  ) ) 
    muh = exp( eta + exp(eta)  )  
  }else if(funcLink == "sqrt"){
    eta = delta + phi*(sqrt(y[n+k-1]) )
    muh = eta^2
  }
  
  return(muh = muh)
  
}
forecast.AR2 <- function(y,par,h,funcLink){
  #y: vetor
  #par: parametros estimados
  #h: passos para previsao
  #funcLink: funcao de ligacao
  
  k=1:h
  delta=par[1];phi = par[2:3]
  n=length(y)
  if( funcLink == "log" ){
    eta = delta + phi[1]*log(y[n+k-1]) + phi[2]*log(y[n+k-2])
    muh = exp(eta)
  }else if(funcLink == "Wlamb"){
    require(lamW)
    eta= delta + phi[1]*( log( lambertW0( y[n+k-1] )  ) ) + phi[2]*( log( lambertW0( y[n+k-2] )  ) )  
    muh = exp( eta + exp(eta)  )  
  }else if(funcLink == "sqrt"){
    eta = delta + phi[1]*(sqrt(y[n+k-1]) ) + phi[2]*(sqrt(y[n+k-2]) )
    muh = eta^2
  }
  
  return(muh = muh)
  
}
forecast.MA1 <- function(y,par,h,funcLink){
  #y: vetor
  #par: parametros estimados
  #h: passos para previsao
  #funcLink: funcao de ligacao
  
  k=1:h
  delta=par[1];theta = par[2]
  n=length(y)
  mu = c(1); eta=c(0)
  
  if( funcLink == "log" ){
    for(t in 2:n){ eta[t] = delta + theta*( log(y[t-1])  - log( mu[t-1] ) ) ; mu[t] = exp( eta[t] )}
    etah = delta + theta*(log(y[n+k-1])  - log( mu[n+k-1] ) )
    muh = exp(etah)
  }else if(funcLink == "Wlamb"){
    require(lamW)
    for(t in 2:n) {eta[t] = delta + theta*( log( lambertW0( y[t-1] )  ) - log( lambertW0( mu[t-1] )  ) )   ; mu[t] = exp( eta[t] + exp(eta[t]) )}
    etah = delta + theta*( log( lambertW0( y[n+k-1] )  ) - log( lambertW0( mu[n+k-1] )  ) ) 
    muh = exp( etah + exp(etah)  )  
  }else if(funcLink == "sqrt"){
    for(t in 2:n){ eta[t] =  delta + theta*( sqrt(y[t-1])  - sqrt( mu[t-1]  )   ) ;mu[t] = eta[t]^2}
    etah = delta + theta*( sqrt(y[n+k-1]) - sqrt(mu[n+k-1]) )
    muh = etah^2
  }
  
  return(muh = muh)
  
}
forecast.MA2 <- function(y,par,h,funcLink){
  #y: vetor
  #par: parametros estimados
  #h: passos para previsao
  #funcLink: funcao de ligacao
  
  k=1:h
  delta=par[1];theta = par[2:3]
  n=length(y)
  mu=c(1) ; eta=c()
  if( funcLink == "log" ){
    for(t in 3:n){ eta[t] = delta + theta[1]*( sqrt(y[t-1]) - sqrt(mu[t-1])    ) + theta[2]*( sqrt(y[t-2]) - sqrt(mu[t-2])    )  ; mu[t] = exp(eta[t])}
    etah = delta + theta[1]*( sqrt(y[n+k-1]) - sqrt(mu[n+k-1])    ) + theta[2]*( sqrt(y[n+k-2]) - sqrt(mu[n+k-2])    )
    muh = exp(etah)
  }else if(funcLink == "Wlamb"){
    require(lamW)
    for(t in 3:n) {eta[t] = delta + theta[1]*( log( lambertW0( y[t-1] )  ) - log( lambertW0( mu[t-1] )  ) ) + theta[2]*( log( lambertW0( y[t-2] )  ) - log( lambertW0( mu[t-2] )  ) )   ; 
    mu[t] = exp( eta[t] + exp(eta[t]) )}
    etah= theta[1]*( log( lambertW0( y[n+k-1] )  ) - log( lambertW0( mu[n+k-1] )  ) ) + theta[2]*( log( lambertW0( y[n+k-2] )  ) - log( lambertW0( mu[n+k-2] )  ) )
    muh = exp( etah + exp(etah)  )  
  }else if(funcLink == "sqrt"){
    for(t in 3:n){ eta[t] = delta + theta[1]*( sqrt(y[t-1]) - sqrt(mu[t-1])    ) + theta[2]*( sqrt(y[t-2]) - sqrt(mu[t-2])    )  ; 
    mu[t] = eta[t]^2}
    etah = delta + theta[1]*(  sqrt(y[n+k-1]) - sqrt(mu[n+k-1])    ) + theta[2]*(  sqrt(y[n+k-2]) - sqrt(mu[n+k-2])    )
    muh = etah^2
  }
  
  return(muh = muh)
  
}
forecast.ARMA11 <- function(y,par,h,funcLink){
  #y: vetor
  #par: parametros estimados
  #h: passos para previsao
  #funcLink: funcao de ligacao
  
  k=1:h
  delta=par[1];phi = par[2]; theta=par[3]
  n=length(y)
  mu=c(1);eta=c()
  
  if( funcLink == "log" ){
    for(t in 2:n) {eta[t] = delta + phi*( log(y[t-1]) ) + theta*( log(y[t-1])  - log( mu[t-1] ) ) ; 
    mu[t] = exp(eta[t])}
    etah = delta + phi*log(y[n+k-1]) + theta*( log(y[n+k-1])  - log( mu[n+k-1] ) ) ; 
    muh = exp(etah)
  }else if(funcLink == "Wlamb"){
    require(lamW)
    for(t in 2:n) {eta[t] = delta + phi*( log(y[t-1]) ) + theta*( log(y[t-1])  - log( mu[t-1] ) ) ; 
    mu[t] = exp( eta[t] + exp(eta[t])  )} 
    etah= delta + phi*( log( lambertW0( y[n+k-1] )  ) )  +theta*( log( lambertW0( y[n+k-1] )  )  - log( lambertW0( mu[n+k-1] )  )  )
    muh = exp( etah + exp(etah)  )  
  }else if(funcLink == "sqrt"){
    for(t in 2:n) {eta[t] = delta + phi*( log(y[t-1]) ) + theta*( log(y[t-1])  - log( mu[t-1] ) ) ; 
    mu[t] = eta[t]^2} 
    etah = delta + phi*(sqrt(y[n+k-1]) ) + theta*( sqrt(y[n+k-1])  - sqrt( mu[n+k-1] ) )
    muh = etah^2
  }
  
  return(muh = muh)
  
}


PRESS.AR1 <- function(x,k,L, funcLink,model){
  #x: serie 
  #model: gamma ARMA or G0-ARMA
  #k: quantidade da serie minima para o kfolds
  n=length(x)
  press=c()
  h=1;
  
  if( model == "gamma" ){
  for(i in 0:(n-(k+1))) 
  {
    y=x[1:(k+i)];yi = x[(k+i+1)]
    if(funcLink == "log"){
    est = cml.gamma.AR1_log(y,L)$est
    prev = forecast.AR1(y,est,h,"log")  
    press[i] = (yi - prev)^2
    }else if(funcLink == "Wlamb"){
      est = cml.gamma.AR1_W(y,L)$est
      prev = forecast.AR1(y,est,h,"Wlamb")  
      press[i] = (yi - prev)^2
    }else if(funcLink == "sqrt"){
      est = cml.gamma.AR1_sqrt(y,L)$est
      prev = forecast.AR1(y,est,h,"sqrt")  
      press[i] = (yi - prev)^2
    }
    
  }
  }else if( model == "g0"){
    for(i in 0:(n-(k+1))) 
    {
      y=x[1:(k+i)];yi = x[(k+i+1)]
      if(funcLink == "log"){
        est = cml.g0.AR1_log(y,L)$est[-1]
        prev = forecast.AR1(y,est,h,"log")  
        press[i] = (yi - prev)^2
      }else if(funcLink == "Wlamb"){
        est = cml.g0.AR1_W(y,L)$est[-1]
        prev = forecast.AR1(y,est,h,"Wlamb")  
        press[i] = (yi - prev)^2
      }else if(funcLink == "sqrt"){
        est = cml.g0.AR1_sqrt(y,L)$est[-1]
        prev = forecast.AR1(y,est,h,"sqrt")  
        press[i] = (yi - prev)^2
      }
      
    }
    
  }
  
  
  SST = sum((x - mean(x))^2)
  R2 = 1-sum(press)/SST
  MSE = mean(press)
  return(list(PRESS = sum(press), MSE=MSE))

}
PRESS.AR2 <- function(x,k,L, funcLink,model){
  #y: serie 
  #y.hat: serie estimada
  #k: quantidade da serie minima para o kfolds
  n=length(x)
  press=c()
  h=1;
  
  if( model == "gamma" ){
  for(i in 0:(n-(k+1))) 
  {
    y=x[1:(k+i)];yi = x[(k+i+1)]
    if(funcLink == "log"){
      est = cml.gamma.AR2_log(y,L)$est
      prev = forecast.AR2(y,est,h,"log")  
      press[i] = (yi - prev)^2
    }else if(funcLink == "Wlamb"){
      est = cml.gamma.AR2_W(y,L)$est
      prev = forecast.AR2(y,est,h,"Wlamb")  
      press[i] = (yi - prev)^2
    }else if(funcLink == "sqrt"){
      est = cml.gamma.AR2_sqrt(y,L)$est
      prev = forecast.AR2(y,est,h,"sqrt")  
      press[i] = (yi - prev)^2
    }
  }
  }else if(model == "g0"){
    for(i in 0:(n-(k+1))) 
    {
      y=x[1:(k+i)];yi = x[(k+i+1)]
      if(funcLink == "log"){
        est = cml.g0.AR2_log(y,L)$est
        prev = forecast.AR2(y,est,h,"log")  
        press[i] = (yi - prev)^2
      }else if(funcLink == "Wlamb"){
        est = cml.g0.AR2_W(y,L)$est
        prev = forecast.AR2(y,est,h,"Wlamb")  
        press[i] = (yi - prev)^2
      }else if(funcLink == "sqrt"){
        est = cml.g0.AR2_sqrt(y,L)$est
        prev = forecast.AR2(y,est,h,"sqrt")  
        press[i] = (yi - prev)^2
      }
    }
    
  }
  SST = sum((x - mean(x))^2)
  R2 = 1-sum(press)/SST
  MSE = mean(press)
  return(list(PRESS = sum(press), MSE=MSE))
  
}
PRESS.MA1 <- function(x,k,L, funcLink,model){
  #y: serie 
  #y.hat: serie estimada
  #k: quantidade da serie minima para o kfolds
  n=length(x)
  press=c()
  h=1;
  
  if(model == "gamma"){
  for(i in 0:(n-(k+1))) 
  {
    y=x[1:(k+i)];yi = x[(k+i+1)]
    if(funcLink == "log"){
      est = cml.gamma.MA1_log(y,L)$est
      prev = forecast.MA1(y,est,h,"log")  
      press[i] = (yi - prev)^2
    }else if(funcLink == "Wlamb"){
      est = cml.gamma.MA1_W(y,L)$est
      prev = forecast.MA1(y,est,h,"Wlamb")  
      press[i] = (yi - prev)^2
    }else if(funcLink == "sqrt"){
      est = cml.gamma.MA1_sqrt(y,L)$est
      prev = forecast.MA1(y,est,h,"sqrt")  
      press[i] = (yi - prev)^2
    }
    
  }
}else if(model == "g0"){
  for(i in 0:(n-(k+1))) 
  {
    y=x[1:(k+i)];yi = x[(k+i+1)]
    if(funcLink == "log"){
      est = cml.gamma.MA1_log(y,L)$est
      prev = forecast.MA1(y,est,h,"log")  
      press[i] = (yi - prev)^2
    }else if(funcLink == "Wlamb"){
      est = cml.g0.MA1_W(y,L)$est
      prev = forecast.MA1(y,est,h,"Wlamb")  
      press[i] = (yi - prev)^2
    }else if(funcLink == "sqrt"){
      est = cml.g0.MA1_sqrt(y,L)$est
      prev = forecast.MA1(y,est,h,"sqrt")  
      press[i] = (yi - prev)^2
    }
  }
}
  SST = sum((x - mean(x))^2)
  R2 = 1-sum(press)/SST
  MSE = mean(press)
  return(list(PRESS = sum(press), MSE=MSE))
  
}
PRESS.MA2 <- function(x,k,L, funcLink,model){
  #y: serie 
  #y.hat: serie estimada
  #k: quantidade da serie minima para o kfolds
  n=length(x)
  press=c()
  h=1;
  
  if(model == "gamma"){
    
  for(i in 0:(n-(k+1))) 
  {
    y=x[1:(k+i)];yi = x[(k+i+1)]
    if(funcLink == "log"){
      est = cml.gamma.MA2_log(y,L)$est
      prev = forecast.MA2(y,est,h,"log")  
      press[i] = (yi - prev)^2
    }else if(funcLink == "Wlamb"){
      est = cml.gamma.MA2_W(y,L)$est
      prev = forecast.MA2(y,est,h,"Wlamb")  
      press[i] = (yi - prev)^2
    }else if(funcLink == "sqrt"){
      est = cml.gamma.MA2_sqrt(y,L)$est
      prev = forecast.MA2(y,est,h,"sqrt")  
      press[i] = (yi - prev)^2
    }
    
  }
  }else if(model == "g0"){
    
    for(i in 0:(n-(k+1))) 
    {
      y=x[1:(k+i)];yi = x[(k+i+1)]
      if(funcLink == "log"){
        est = cml.g0.MA2_log(y,L)$est
        prev = forecast.MA2(y,est,h,"log")  
        press[i] = (yi - prev)^2
      }else if(funcLink == "Wlamb"){
        est = cml.g0.MA2_W(y,L)$est
        prev = forecast.MA2(y,est,h,"Wlamb")  
        press[i] = (yi - prev)^2
      }else if(funcLink == "sqrt"){
        est = cml.g0.MA2_sqrt(y,L)$est
        prev = forecast.MA2(y,est,h,"sqrt")  
        press[i] = (yi - prev)^2
      }
      
    }
  }    
  SST = sum((x - mean(x))^2)
  R2 = 1-sum(press)/SST
  MSE = mean(press)
  return(list(PRESS = sum(press), MSE=MSE))
  
}
PRESS.ARMA11 <- function(x,k,L, funcLink,model){
  #y: serie 
  #y.hat: serie estimada
  #k: quantidade da serie minima para o kfolds
  n=length(x)
  press=c()
  h=1;
  
  if(model == "gamma"){
  for(i in 0:(n-(k+1))) 
  {
    y=x[1:(k+i)];yi = x[(k+i+1)]
    if(funcLink == "log"){
      est = cml.gamma.ARMA11_log(y,L)$est
      prev = forecast.ARMA11(y,est,h,"log")  
      press[i] = (yi - prev)^2
    }else if(funcLink == "Wlamb"){
      est = cml.gamma.ARMA11_W(y,L)$est
      prev = forecast.ARMA11(y,est,h,"Wlamb")  
      press[i] = (yi - prev)^2
    }else if(funcLink == "sqrt"){
      est = cml.gamma.ARMA11_sqrt(y,L)$est
      prev = forecast.ARMA11(y,est,h,"sqrt")  
      press[i] = (yi - prev)^2
    }
    
  }
  }else if( model == "g0"  ){
    for(i in 0:(n-(k+1))) 
    {
      y=x[1:(k+i)];yi = x[(k+i+1)]
      if(funcLink == "log"){
        est = cml.g0.ARMA11_log(y,L)$est
        prev = forecast.ARMA11(y,est,h,"log")  
        press[i] = (yi - prev)^2
      }else if(funcLink == "Wlamb"){
        est = cml.g0.ARMA11_W(y,L)$est
        prev = forecast.ARMA11(y,est,h,"Wlamb")  
        press[i] = (yi - prev)^2
      }else if(funcLink == "sqrt"){
        est = cml.g0.ARMA11_sqrt(y,L)$est
        prev = forecast.ARMA11(y,est,h,"sqrt")  
        press[i] = (yi - prev)^2
      }
      
    }
  }
  SST = sum((x - mean(x))^2)
  R2 = 1-sum(press)/SST
  
  MSE = mean(press)
  return(list(PRESS = sum(press), MSE=MSE))  
}









