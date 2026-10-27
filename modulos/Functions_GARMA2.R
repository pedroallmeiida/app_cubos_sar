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

# Geracao distribuicao G0



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
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) ) 
  
  for(t in 2:N){
    mu[t] = exp(
      delta + phin*( log(z[t-1]) ) 
    )
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )
  }
  
  au = rgamma(N, shape=look, scale=1/look); 
  x = au*z
  return(x)
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
  
  for(t in 3:N){
    mu[t] = exp(
      delta + phin[1]*( log(z[t-1]) ) + phin[2]*( log(z[t-2]) ) 
    )
    z[t] = rinvgamma( 1, shape=-alpha, rate=mu[t]*(-alpha-1) )
  }
  
  au = rgamma(N, shape=look, scale=1/look); 
  x = au*z
  return(x)
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
  
  
  for(t in 2:N){
    
    error[t-1] = log( z[t-1] ) - log( mu[t-1] )
    mu[t] = exp(
      delta + theta*error[t-1]
    )
    z[t] = rinvgamma( 1, shape = -alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape = look, scale = 1/look); 
  }
  
  return(z)
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
  
  
  for(t in 2:N){
    mu[t] = exp(
      delta + phin*( log(z[t-1]) ) + theta*( log(z[t-1]) - log(mu[t-1]) ) 
    )
    z[t] = rinvgamma( 1, shape=-alpha, rate=mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  return(z)
}


## Gera??o da distribui??o gamma



rgama_AR1_log = function(N,vP){
  # N: number of outcomes
  # vP: vector of parameters  
  
  
  alpha = vP[1]; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = alpha, scale = 1/(alpha/1) )
  
  mu[1] = exp(
    delta +  phin*( log(z0) )    
  )
  z[1] = rgamma( 1, shape = alpha, scale = 1/(alpha/mu[1]) )
  
  
  for(t in 2:N){
    mu[t] = exp(
      delta + phin*( log(z[t-1]) )   
    )
    z[t] = rgamma(1,shape = alpha, scale = 1/(alpha/mu[t]))
    
  }
  
  return(z)
  
}

rgama_AR2_log = function(N,vP){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = alpha, scale = 1/(alpha/1) )
  
  mu[1] = exp(
    delta +  phin[1]*( log(z0) ) +  phin[2]*( log(z0) )  
  )
  z[1] = rgamma( 1, shape = alpha, scale = 1/(alpha/mu[1]) )
  
  mu[2] = exp(
    delta +  phin[1]*( log(z[1]) ) +  phin[2]*( log(z0) )   
  )
  z[2] = rgamma( 1, shape = alpha, scale = 1/(alpha/mu[2]) )
  
  
  for(t in 3:N){
    mu[t] = exp(
      delta +  phin[1]*( log(z[t-1]) ) +  phin[2]*( log(z[t-2]) )   
    )
    z[t] = rgamma(1,shape = alpha, scale = 1/(alpha/mu[t]))
    
  }
  
  return(z)
  
}

rgama_MA1_log = function(N,vP){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; delta = vP[2]; theta = vP[3]
  mu=c();z=c();
  
  mu0 = 1
  x0 = rgamma( 1, shape = alpha, scale = 1/(alpha/1) )
  z0 = rgamma( 1, shape = alpha, scale = 1/(alpha/1) )
  
  mu[1] = exp(
    delta +  theta*( log(z0) - log(mu0)  )  
  )
  z[1] = rgamma( 1, shape = alpha, scale = 1/(alpha/mu[1]) )
  
  for(t in 2:N){
    mu[t] = exp(
      delta +  theta*( log(z[t-1]) - log(mu[t-1]) )  
    )
    z[t] = rgamma(1,shape = alpha, scale = 1/(alpha/mu[t]))
    
  }
  
  return(z)
  
}

rgama_ARMA11_log = function(N,vP){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; delta = vP[2]; phin = vP[3]; theta = vP[4]
  mu=c();z=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = alpha, scale = 1/(alpha/1) )
  
  mu[1] = exp(
    delta + phin*( log(z0) ) +  theta*( log(z0) - log(mu0)  )  
  )
  z[1] = rgamma( 1, shape = alpha, scale = 1/(alpha/mu[1]) )
  
  for(t in 2:N){
    mu[t] = exp(
      delta + phin*( log(z[t-1]) ) +  theta*( log(z[t-1]) - log(mu[t-1]) )  
    )
    z[t] = rgamma(1,shape = alpha, scale = 1/(alpha/mu[t]))
    
  }
  
  return(z)
  
}


rG0_AR1_boxcox = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p=1
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3]
  mu=c();z=c()
  z0 = 1
  
  mu[1] = ( lambda*(delta + phin*( log(z0) ) ) + 1 )^(1/lambda)
  z[1] = rinvgamma( 1, shape=-alpha, rate = mu[1]*(-alpha-1) ) 
  
  for(t in 2:N){
    mu[t] = ( lambda*(delta + phin*( log(z[t-1]) ) ) + 1 )^(1/lambda)
    z[t] = rinvgamma( 1, shape=-alpha, rate = mu[t]*(-alpha-1) )
  }
  
  au = rgamma(N, shape=look, scale=1/look); 
  x = au*z
  return(x)
}

rG0_AR2_boxcox = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  p = 2
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  z0 = 1
  
  mu[1] = ( lambda*(  delta + phin[1]*( log(z0) ) + phin[2]*( log(z0) ) ) + 1)^(1/lambda)
  z[1] = rinvgamma( 1, shape=-alpha, rate=mu[1]*(-alpha-1) ) 
  
  mu[2] = ( lambda*(  delta + phin[1]*( log(z[1]) ) + phin[2]*( log(z0) ) ) + 1)^(1/lambda)
  z[2] = rinvgamma( 1, shape=-alpha, rate=mu[2]*(-alpha-1) ) 
  
  for(t in 3:N){
    mu[t] = ( lambda*(  delta + phin[1]*( log(z[t-1]) ) + phin[2]*( log(z[t-2]) ) ) + 1)^(1/lambda)
    z[t] = rinvgamma( 1, shape=-alpha, rate=mu[t]*(-alpha-1) )
  }
  
  au = rgamma(N, shape=look, scale=1/look); 
  x = au*z
  return(x)
}

rG0_MA1_boxcox = function(N,par,L){
  # N: number of outcomes
  # par: vector of parameters  
  
  q=1; look = L
  alpha = par[1]; delta = par[2]; theta = par[3]
  mu=c();z=c();
  
  error=rep(0, N)
  
  error0 = 0
  mu[1] = (  lambda*(delta + theta*error0 ) + 1)^(1/lambda) 
  
  z[1] = rinvgamma( 1, shape =-alpha, rate = mu[1]*(-alpha-1) )*rgamma(1, shape = look, scale = 1/look); 
  
  
  for(t in 2:N){
    
    error[t-1] = log( z[t-1] ) - log( mu[t-1] )
    mu[t] = ( lambda*(delta + theta*error[t-1] ) + 1)^(1/lambda) 
    z[t] = rinvgamma( 1, shape = -alpha, rate = mu[t]*(-alpha-1) )*rgamma(1, shape = look, scale = 1/look); 
  }
  
  return(z)
}

rG0_ARMA11_boxcox = function(N,vP,L){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; look = L; delta = vP[2]; phin = vP[3];theta = vP[4]
  mu=c();z=c()
  z0 = 1;mu0=1
  
  mu[1] =   ( lambda*( delta + phin*( log(z0) ) + theta*( log(z0) - log(mu0) ) ) + 1)^(1/lambda)
  z[1] = rinvgamma( 1, shape=-alpha, rate=mu[1]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look) 
  
  
  for(t in 2:N){
    mu[t] = ( lambda*( delta + phin*( log(z[t-1]) ) + theta*( log(z[t-1]) - log(mu[t-1]) ) ) + 1)^(1/lambda)
    z[t] = rinvgamma( 1, shape=-alpha, rate=mu[t]*(-alpha-1) )*rgamma(1, shape=look, scale=1/look)
  }
  
  return(z)
}


## Gera??o da distribui??o gamma



rgama_AR1_boxcox = function(N,vP){
  # N: number of outcomes
  # vP: vector of parameters  
  
  
  alpha = vP[1]; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = alpha, scale = 1/(alpha/1) )
  
  mu[1] = (
    lambda*(delta +  phin*( log(z0) ) ) + 1   
  )^(1/lambda)
  z[1] = rgamma( 1, shape = alpha, scale = 1/(alpha/mu[1]) )
  
  
  for(t in 2:N){
    mu[t] = (
      lambda*(delta + phin*( log(z[t-1]) ) ) + 1 
    )^(1/lambda)
    z[t] = rgamma(1,shape = alpha, scale = 1/(alpha/mu[t]))
    
  }
  
  return(z)
  
}

rgama_AR2_boxcox = function(N,vP){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; delta = vP[2]; phin = vP[3:length(vP)]
  mu=c();z=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = alpha, scale = 1/(alpha/1) )
  
  mu[1] = (
    lambda*(delta +  phin[1]*( log(z0) ) +  phin[2]*( log(z0) )   )  +1
  )^(1/lambda)
  z[1] = rgamma( 1, shape = alpha, scale = 1/(alpha/mu[1]) )
  
  mu[2] = (
    lambda*(delta +  phin[1]*( log(z[1]) ) +  phin[2]*( log(z0) ) )  +1   
  )^(1/lambda)
  z[2] = rgamma( 1, shape = alpha, scale = 1/(alpha/mu[2]) )
  
  
  for(t in 3:N){
    mu[t] = (
      lambda*(delta +  phin[1]*( log(z[t-1]) ) +  phin[2]*( log(z[t-2]) ) )  +1  
    )^(1/lambda)
    z[t] = rgamma(1,shape = alpha, scale = 1/(alpha/mu[t]))
    
  }
  
  return(z)
  
}

rgama_MA1_boxcox = function(N,vP){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; delta = vP[2]; theta = vP[3]
  mu=c();z=c();
  
  mu0 = 1
  x0 = rgamma( 1, shape = alpha, scale = 1/(alpha/1) )
  z0 = rgamma( 1, shape = alpha, scale = 1/(alpha/1) )
  
  mu[1] = (
    lambda*(delta +  theta*( log(z0) - log(mu0)  )  )  +1
  )^(1/lambda)
  z[1] = rgamma( 1, shape = alpha, scale = 1/(alpha/mu[1]) )
  
  for(t in 2:N){
    mu[t] = (
      lambda*(delta +  theta*( log(z[t-1]) - log(mu[t-1]) )  )  +1
    )^(1/lambda)
    z[t] = rgamma(1,shape = alpha, scale = 1/(alpha/mu[t]))
    
  }
  
  return(z)
  
}

rgama_ARMA11_boxcox = function(N,vP){
  # N: number of outcomes
  # vP: vector of parameters  
  
  alpha = vP[1]; delta = vP[2]; phin = vP[3]; theta = vP[4]
  mu=c();z=c()
  
  mu0 = 1
  z0 = rgamma( 1, shape = alpha, scale = 1/(alpha/1) )
  
  mu[1] = (
    lambda*(delta + phin*( log(z0) ) +  theta*( log(z0) - log(mu0)  )  )  +1
  )^(1/lambda)
  z[1] = rgamma( 1, shape = alpha, scale = 1/(alpha/mu[1]) )
  
  for(t in 2:N){
    mu[t] = (
      lambda*(delta + phin*( log(z[t-1]) ) +  theta*( log(z[t-1]) - log(mu[t-1]) ) )  +1
    )^(1/lambda)
    z[t] = rgamma(1,shape = alpha, scale = 1/(alpha/mu[t]))
    
  }
  
  return(z)
  
}





##### end geradores ####



####    Log verossimilhan?as     #####

#### Distribuicao G0 ####

## funcao  de ligacao log() ----

lk.g0.AR1.log = function(z,par,L){
  
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
  
  mu=c(); v = c(); 
  
  mu[1] <- exp( delta + phin*( log(mean(z)) )  )
  v[1] <- ( mu[1]*(-alpha-1) + L*z[1] )
  
 t <- (p+1):n
 mu[t] = exp( delta + phin*( log(z[t-1]) )  )
 v[t] = ( mu[t]*(-alpha-1) + L*z[t] )


  mun <- mu[(p+1):n]
  zn <- z[(p+1):n]
  vn <- v[(p+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
  
}

lk.g0.AR2.log = function(z,par,L){
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## <-- Output
  ## log.verossimilhan?a 
  
  p = 2
  n = length(z)
  alpha = par[1]; delta = par[2]; phin=par[3:length(par)]
  
  z0 <- rep(mean(z),p)
  mu <- v <- rep(NA,n)
  
  mu[1] = exp( delta + phin[1]*( log(z0[1]) ) + phin[2]*( log(z0[2]) ) )
  mu[2] = exp( delta + phin[1]*( log(z[1]) ) + phin[2]*( log(z0[2]) ) )
  
  t <- (p+1):n
    mu[t] = exp( delta + phin[1]*( log(z[t-1]) ) + phin[2]*( log(z[t-2]) ) )
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  
  
  mun <- mu[(p+1):n]
  zn <- z[(p+1):n]
  vn <- v[(p+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
  
}

lk.g0.ARP.log = function(z,par,p,L){
  
  n=length(z)
  alpha = par[1]; delta = par[2]; phin=par[3:length(par)]
  mu=c();v=c();lk=c()
  
  z0 = rep(mean(z),p)
  zz = matrix(NA,n,p)
  zz[1:p,] =  log(z0)  
  
  
  for(j in 1:p){
    if(j < p){
      mu[j] = exp( delta + phin%*%zz[j,])
    }else if(j == p & j == 1){
      mu[j] = exp(delta + phin%*%( c(log(z0[1]), log(z0[1:(p-1)]) )   ))
    }else if(j == p & j != 1){
      mu[j] = exp( delta + phin%*%( c(log(z[j-1]), log(z0[1:(p-1)]) ) ) )
    }
    v[j] = ( mu[j]*(-alpha-1) + L*z[j] )
    lk[j] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[j]) - alpha*log(-alpha-1) -log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[j]) + (alpha-L)*log(v[j])
  }
  
  t<-(p+1):n

    for(i in t:(t+p-(2*p-1))) zz[t,i-(t-1) + 2*(t-i)] = log( z[i- 1] )   
  
  mu[t] = exp(
    delta + phin%*%zz[t,] 
  )
  v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  lk[t] = L*log(L) + log( gamma(L - alpha) ) - alpha*log(mu[t]) - alpha*log(-alpha-1)-log( gamma(-alpha) ) -log( gamma(L) ) + (L-1)*log(z[t]) + (alpha-L)*log(v[t])
  
  
  
  return(sum(lk[(p+1):n]))
  
}

lk.g0.MA1.log <- function(z, par,L){
  q=1
  n=length(z)
  alpha = par[1]; delta<- par[2]; theta <- par[3]
  
  v <- rep(NA,n)
  mu=c()
  mu[1] <- exp( delta + theta*( log( mean(z) ) - log( 1 )  ) ) 
  
  for(t in (q+1):n){
    mu[t] <- exp( delta + theta*( log(z[t-1]) - log(mu[t-1])  ) )
  }
  
  t<-(q+1):n
  v[t] <- mu[t]*(-alpha-1) + L*z[t] 
  
  mun <- mu[(q+1):n]
  zn <- z[(q+1):n]
  vn <- v[(q+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 

lk.g0.MA2.log = function(z,par,L){
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 2
  n = length(z)
  alpha = par[1]; delta = par[2]; thetan=par[3:length(par)]
  
  mu = c();
  z0 = rep(mean(z),q)
  mu0 = 1
  v <- rep(NA,n)
  
  mu[1] = exp( delta + thetan[1]*( log(z0[1]) - log(mu0) ) + thetan[2]*( log(z0[2]) - log(mu0) ) )
  mu[2] = exp( delta + thetan[1]*( log(z[1]) - log(mu[1]) ) + thetan[2]*( log(z0[2]) - log(mu0) ) )
  
  
  for(t in (q+1):n){
    mu[t] = exp( delta + thetan[1]*( log(z[t-1]) - log(mu[t-1]) ) + thetan[2]*( log(z[t-2]) - log(mu[t-2]) ) )
  }
  
  t<-(q+1):n
  v[t] <- mu[t]*(-alpha-1) + L*z[t] 
  
  mun <- mu[(q+1):n]
  zn <- z[(q+1):n]
  vn <- v[(q+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
  
}

lk.g0.ARMA11.log <- function(z, par,L){
  
  n=length(z)
  alpha = par[1]; delta<- par[2]
  phi <- par[3]; theta <- par[4]
  
  q=1;p=1; m = max(p,q)
  
  v <- rep(NA,n)
  mu=c()
  mu[1] <- exp( delta + phi*log( mean(z) ) + theta*( log( mean(z) ) - log( 1 )  ) ) 
  
  
  
  for(t in (m+1):n){
    mu[t]   <- exp( delta + phi*log(z[t-p]) + theta*( log(z[t-q]) - log(mu[t-q]) ) )
  }
  
  t<-(m+1):n
  v[t] <- mu[t]*(-alpha-1) + L*z[t] 
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 


## funcao  de ligacao log() ----

lk.g0.AR1.boxcox = function(z,par,L){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  p=1
  n=length(z)
  alpha = par[1]; delta = par[2]; phin=par[3]
  
  
  
  mu=c(); v = c(); 
  
  mu[1] <- ( lambda*( delta + phin*( log(mean(z)) ) ) + 1  )^(1/lambda)
  v[1] <- (  mu[1]*(-alpha-1) + L*z[1] )
  
  
  
  t <- (p+1):n
  mu[t] = ( lambda*( delta + phin*( log(z[t-1]) ) ) + 1  )^(1/lambda)
  v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  
  
  mun <- mu[(p+1):n]
  zn <- z[(p+1):n]
  vn <- v[(p+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
  
}

lk.g0.AR2.boxcox = function(z,par,L){
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## <-- Output
  ## log.verossimilhan?a 
  
  
  
  p = 2
  n = length(z)
  alpha = par[1]; delta = par[2]; phin=par[3:length(par)]
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  
  z0 <- rep(mean(z),p)
  mu <- v <- rep(NA,n)
  
  mu[1] = (lambda*( delta + phin[1]*( log(mean(z)) ) + phin[2]*( log(mean(z)) ) ) + 1  )^(1/lambda)
  mu[2] = (lambda*( delta + phin[1]*( log(z[1]) ) + phin[2]*( log(mean(z)) ) ) + 1 )^(1/lambda)
  
  
  t <- (p+1):n
  mu[t] = (lambda*( delta + phin[1]*( log(z[t-1]) ) + phin[2]*( log(z[t-2]) ) ) + 1 )^(1/lambda)
  v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  
  
  mun <- mu[(p+1):n]
  zn <- z[(p+1):n]
  vn <- v[(p+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
  
}

lk.g0.MA1.boxcox <- function(z, par,L){
  q=1
  n=length(z)
  alpha = par[1]; delta<- par[2]; theta <- par[3]
  
  v <- rep(NA,n)
  mu=c()
  mu[1] <- ( lambda*( delta + theta*( log( mean(z) ) - log( 1 )  ) ) + 1)^(1/lambda) 
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  
  for(t in (q+1):n){
    mu[t] <- ( lambda*( delta + theta*( log(z[t-1]) - log(mu[t-1])  ) ) + 1)^(1/lambda) 
  }
  
  t<-(q+1):n
  v[t] <- mu[t]*(-alpha-1) + L*z[t] 
  
  mun <- mu[(q+1):n]
  zn <- z[(q+1):n]
  vn <- v[(q+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 

lk.g0.MA2.boxcox = function(z,par,L){
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(alpha,delta,Phi)
  ## <-- Output
  ## log.verossimilhan?a 
  
  q = 2
  n = length(z)
  alpha = par[1]; delta = par[2]; thetan=par[3:length(par)]
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  
  mu = c();
  z0 = rep(mean(z),q)
  mu0 = 1
  v <- rep(NA,n)
  
  mu[1] = ( lambda*(delta + thetan[1]*( log(z0[1]) - log(mu0) ) + thetan[2]*( log(z0[2]) - log(mu0) ) ) + 1 )^(1/lambda)
  mu[2] = ( lambda*(delta + thetan[1]*( log(z[1]) - log(mu[1]) ) + thetan[2]*( log(z0[2]) - log(mu0) ) ) + 1 )^(1/lambda)
  
  
  for(t in (q+1):n){
    mu[t] = ( lambda*( delta + thetan[1]*( log(z[t-1]) - log(mu[t-1]) ) + thetan[2]*( log(z[t-2]) - log(mu[t-2]) ) ) +1 )^(1/lambda)
  }
  
  t<-(q+1):n
  v[t] <- mu[t]*(-alpha-1) + L*z[t] 
  
  mun <- mu[(q+1):n]
  zn <- z[(q+1):n]
  vn <- v[(q+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
  
}

lk.g0.ARMA11.boxcox <- function(z, par,L){
  
  n=length(z)
  alpha = par[1]; delta<- par[2]
  phi <- par[3]; theta <- par[4]
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  
  q=1;p=1; m = max(p,q)
  
  v <- rep(NA,n)
  mu=c()
  mu[1] <- ( lambda*(delta + phi*log( mean(z) ) + theta*( log( mean(z) ) - log( 1 )  ) ) + 1 )^(1/lambda) 
  
  
  
  for(t in (m+1):n){
    mu[t]   <- ( lambda*(delta + phi*log(z[t-p]) + theta*( log(z[t-q]) - log(mu[t-q]) ) ) + 1 )^(1/lambda)
  }
  
  t<-(m+1):n
  v[t] <- mu[t]*(-alpha-1) + L*z[t] 
  
  mun <- mu[(m+1):n]
  zn<-z[(m+1):n]
  vn <-v[(m+1):n]
  
  
  lv = L*log(L) + log(gamma(L-alpha)) - alpha*log(mun) - alpha*log(-alpha-1) - log(gamma(-alpha)) - log(gamma(L)) + (L-1)*log(zn) + (alpha-L)*log(vn) 
  return(sum(lv))
} 




### Log verossimilhan?a do processo Gama

lk.gamma.AR1.boxcox = function(y,par){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  require(forecast)
  lambda = BoxCox.lambda(y)
  
  p=1
  n=length(y)
  L = par[1]; delta = par[2]; phin=par[3]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)
  y0 = rep(mean(y),p)
  
  mu[1] <- ( lambda*( delta + phin*( log(mean(y)) ) ) + 1  )^(1/lambda)
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  t<-2:n
  mu[t] <- ( lambda*( delta + phin*( log(y[t-1]) ) ) + 1  )^(1/lambda)
  
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  
  return(sum(lk[2:n]))
  
}

lk.gamma.AR2.boxcox = function(y,par){
  
  ## ---> Input
  ## y: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  require(forecast)
  lambda = BoxCox.lambda(y)
  
  p=2
  n=length(y)
  L = par[1]; delta = par[2]; phin=par[3:4]
  mu=c();v=c();lk=c()
  y0 = mean(y)
  
  mu[1] = (lambda*( delta + phin[1]*( log(mean(y)) ) + phin[2]*( log(mean(y)) ) ) + 1  )^(1/lambda)
  mu[2] = (lambda*( delta + phin[1]*( log(y[1]) ) + phin[2]*( log(mean(y)) ) ) + 1 )^(1/lambda)
  
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1]) 
  lk[2] = (L-1)*log(y[2]) + L*log(L) - log( gamma(L)  ) - L*(  y[2]/mu[2] )  - L*log(mu[2]) 
  
  t <- (p+1):n
  mu[t] = (lambda*( delta + phin[1]*( log( y[t-1] ) ) + phin[2]*( log( y[t-2] ) ) ) + 1 )^(1/lambda)
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] )  - L*log(mu[t]) 
  
  
  
  return(sum(lk[3:n]))
  
}

lk.gamma.MA1.boxcox = function(y,par){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  require(forecast)
  lambda = BoxCox.lambda(y)
  
  q = 1
  n = length(y)
  L = par[1]; delta = par[2]; thetan=par[3]
  mu=c();lk=c()
  
  mu0 = rep(1,q)
  y0 = rep(mean(y),q)
  
  mu[1] <- ( lambda*( delta + thetan*( log( mean(y) ) - log( 1 )  ) ) + 1)^(1/lambda) 
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in 2:n){
    mu[t] = ( lambda*( delta + thetan*( log(y[t-1]) - log(mu[t-1]) ) ) + 1)^(1/lambda) 
  }
  
  t <- 2:n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L) ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return( sum(lk[2:n]) )
  
}

lk.gamma.MA2.boxcox= function(y,par){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  require(forecast)
  lambda = BoxCox.lambda(y)
  
  q = 2
  n = length(y)
  L = par[1]; delta = par[2]; thetan=par[3:4]
  mu=c();v=c();lk=c()
  
  mu0 = 1
  y0 = mean(y)
  
  mu[1] = ( lambda*(delta + thetan[1]*( log(y0) - log(mu0) ) + thetan[2]*( log(y0) - log(mu0) ) ) + 1 )^(1/lambda)
  mu[2] = ( lambda*(delta + thetan[1]*( log(y[1]) - log(mu[1]) ) + thetan[2]*( log(y0) - log(mu0) ) ) + 1 )^(1/lambda)
  
  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  lk[2] = (L-1)*log(y[2]) + L*log(L) - log( gamma(L)  ) - L*(  y[2]/mu[2] )  - L*log(mu[2])
  

  
  for(t in (q+1):n){
    mu[t] = ( lambda*( delta + thetan[1]*( log(y[t-1]) - log(mu[t-1]) ) + thetan[2]*( log(y[t-2]) - log(mu[t-2]) ) ) +1 )^(1/lambda)
  }
  
  
  t<-(q+1):n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L) ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return(sum(lk[(q+1):n]))
  
}

lk.gamma.ARMA11.boxcox= function(y,par){
  
  ## ---> Input
  ## z: s?rie resposta
  ## x: covari?vel 
  ## par: vetor de par?metros c(L,delta,Phi)
  ## L: numero de looks
  ## <-- Output
  ## log.verossimilhan?a 
  
  require(forecast)
  lambda = BoxCox.lambda(y)
  
  p=1
  n=length(y)
  L = par[1]; delta = par[2]; phin = par[3]; thetan=par[4]
  mu=c();v=c();lk=c()
  
  mu0 = rep(1,p)

  mu[1] = ( lambda*(delta + phin*log(mean(y))  + thetan*( log(mean(y)) - log(mu0) ) ) + 1 )^(1/lambda)

  
  lk[1] = (L-1)*log(y[1]) + L*log(L) - log( gamma(L)  ) - L*(  y[1]/mu[1] )  - L*log(mu[1])
  
  for(t in (p+1):n){
    mu[t] = ( lambda*(delta + phin*log(mean(y[t-1]))  + thetan*( log( mean(y) ) - log(mu[t-1]) ) ) + 1 )^(1/lambda)
  }
  
  
  t <- (p+1):n
  lk[t] = (L-1)*log(y[t]) + L*log(L) - log( gamma(L)  ) - L*(  y[t]/mu[t] ) - L*log(mu[t]) 
  
  
  return(sum(lk[2:n]))
  
}




#### end log verossimilhan?as ####









#### Entropias ####


entropy.gamma <- function(y, par){
  
  n=length(y)
  L = par[1]; delta = par[2]; phi = par[3];
  mu=c();
  y0 = mean(y)
  
  H=c()
  
  mu[1] = exp( delta + phi*log(y0) )
  
  t <- 2:n
  mu[t] = exp( delta + phi*log(y[t-1]) )
  
  H[t] = L + log( gamma(L) )  + (1 - L)*digamma(L) - log(L) + log(mu[t])
  sum(H)
}



#### Vetor Gradiente #### 


## Gamma

grad.gamma.AR1.log <- function(y,param) {
  L <- param[1]
  delta <- param[2]
  phi <- param[3]
  n = length(y)
  mu=c()
  
  mu[1] = exp( delta + phi*log( mean(y) ) )
  t<-2:n
  mu[t] = exp( delta + phi*log( y[t-1] ) )
  
  
  
  g <- matrix(0, length(y)-1, 3)
  g[,1] <- log(L) - digamma(L) + 1 + log(y[2:n]) - y[2:n]/mu[2:n] - log(mu[2:n])
  g[,2] <- L*( y[2:n]/mu[2:n] -1 )
  g[,3] <- L*( ( y[2:n]/mu[2:n] -1 )*log(y[1:(n-1)]) )
  colSums(g)
  g
}
grad.gamma.AR2.log <- function(y,param) {
  L <- param[1]
  delta <- param[2]
  phi <- param[3:4]
  n = length(y)
  mu=c();p=2
  
  mu0 = 1
  y0 = mean(y)
  
  mu[1] = exp(
    delta + phi[1]*( log(y0) ) + phi[2]*( log(y0) )
  )
  mu[2] = exp(
    delta + phi[1]*( log(y[1]) ) + phi[2]*( log(y0) )
  )
  
  t<-(p+1):n
  mu[t] = exp( delta + phi[1]*log(y[t-1]) + phi[2]*( log(y[t-2]) ) )
  
  g <- matrix(0, n-p, 4)
  g[,1] <- log(L) - digamma(L) + 1 + log(y[t]) - y[t]/mu[t] - log(mu[t])
  g[,2] <- L*( y[t]/mu[t] -1 )
  g[,3] <- L*( ( y[t]/mu[t] -1 )*log(y[t-1]) )
  g[,4] <- L*( ( y[t]/mu[t] -1 )*log(y[t-2]) )
  colSums(g)
  g
}
grad.gamma.MA1.log <- function(y,param) {
  m=q = 1
  L <- param[1]
  delta <- param[2]
  theta <- param[3]
  
  n = length(y)
  mu = c();mu0 = 1; y0 = 1
  error = c()
  
  mu[1] = exp( delta + theta*( log( y0 ) - log(mu0) ) )
  #error[1] = log( y0 ) - log( mu0 )
  
  
  for(t in (m+1):n){ 
    
    mu[t] = exp( delta + theta*(log( y[t-q] ) - log( mu[t-q] )) )
    #error[t] = log( y[t-q] ) - log( mu[t-q] )
    
  }
  
  
  deta.ddelta <- rep(0, n)
  deta.dtheta <- rep(0, n)
  
  for(t in (m+1):n)
  {
    deta.ddelta[t] <- 1 - theta*deta.ddelta[t-q]
    deta.dtheta[t] <- log( y[t-q] ) - log( mu[t-q] ) - theta*deta.dtheta[t-q]
  }
  
  t <- (m+1):n
  g <- matrix(0, length(y)-1, 3)
  g[,1] <- log(L) - digamma(L) + 1 + log(y[t]) - y[t]/mu[t] - log(mu[t])
  g[,2] <- L*( y[t]/mu[t] -1 )*(deta.ddelta[t])
  g[,3] <- L*( y[t]/mu[t] -1 )*(deta.dtheta[t])
  colSums(g)
  g
}
grad.gamma.MA.log <- function(y,par,q){
  
  m=q
  L = par[1]; delta = par[2]; theta = par[3:length(par)]
  n = length(y)
  mu <- rep(1,m);
  g <- matrix(0, n, q+2)
  
  deta.ddelta <- rep(0,n)
  deta.dtheta <- matrix(0, ncol=q,nrow=n)
  
  
  for(t in (m+1):n){
    mu[t] <- exp( 
      delta + theta%*%( log(y[t-1:q]) - log(mu[t-1:q]) ) 
    )
    deta.ddelta[t] <- 1 - theta%*%deta.ddelta[t-1:q]
    deta.dtheta[t,] <- log( y[t-1:q] ) - log( mu[t-1:q] ) - theta%*%deta.dtheta[t-1:q,]
  }
  
  t <- (m+1):n
  
  g[t,1] <- log(L) - digamma(L) + 1 + log(y[t]) - y[t]/mu[t] - log(mu[t])
  g[t,2] <- L*( y[t]/mu[t] -1 )%*%(deta.ddelta[t])
  g[t,3:length(par)] <- L*( y[t]/mu[t] -1 )%*%(deta.dtheta[t,])
  
  
  U=g[(m+1):n,]
  return(U)
}
grad.gamma.ARMA11.log <- function(y,param) {
  m = q = p = 1
  L <- param[1]
  delta <- param[2]
  phi <- param[3]
  theta <- param[4]
  
  n = length(y)
  mu = c();mu0 = 1; y0 = 1
  error = c()
  
  mu[1] = exp( delta + phi*log(y0) + theta*( log( y0 ) - log(mu0) ) )
  #error[1] = log( y0 ) - log( mu0 )
  
  deta.ddelta <- rep(0, n)
  deta.dphi <- rep(0, n)
  deta.dtheta <- rep(0, n)
  
  
  for(t in (m+1):n){ 
    
    mu[t] = exp( delta + phi*log(y[t-p]) + theta*(log( y[t-q] ) - log( mu[t-q] )) )
    deta.ddelta[t] <- 1 - theta*deta.ddelta[t-q]
    deta.dphi[t] <- log( y[t-p] ) - theta*deta.dphi[t-p]
    deta.dtheta[t] <- log( y[t-q] ) - log( mu[t-q] ) - theta*deta.dtheta[t-q]
    
  }
  
  
  t <- (m+1):n
  g <- matrix(0, length(y)-1, 4)
  g[,1] <- log(L) - digamma(L) + 1 + log(y[t]) - y[t]/mu[t] - log(mu[t])
  g[,2] <- L*( y[t]/mu[t] -1 )*(deta.ddelta[t])
  g[,3] <- L*( y[t]/mu[t] -1 )*(deta.dphi[t])
  g[,4] <- L*( y[t]/mu[t] -1 )*(deta.dtheta[t])
  colSums(g)
  g
}

### G0I 

grad.g0.AR1.log <- function(z,par,L){
  
  alpha <- par[1]; delta <- par[2]; phi <- par[3]
  
  n = length(z)
  mu = c()
  v = rep(NA,n)
  
  
  mu[1] = exp( delta + phi*log( mean(z) ) )
  v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
  
  t<-2:n
  mu[t] = exp( delta + phi*log( z[t-1] ) )
  v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  
  
  c1 = (alpha-L)*(-alpha-1)
  g <- matrix(0, length(z)-1, 3)
  g[,1] <- digamma(-alpha) - digamma(L - alpha) - log( -alpha-1 ) + alpha/(-alpha-1) + log(v[t]/mu[t]) - (alpha-L)*(mu[t]/v[t])
  g[,2] <- c1*(mu[t]/v[t]) - alpha
  g[,3] <- (c1*(mu[t]/v[t]) - alpha)*log(z[t-1])
  g
}
grad.g0.AR.log <- function(z,par,L) {
  
  alpha <- par[1]; delta <- par[2]; phi <- par[3:length(par)]
  
  n = length(z)
  mu <- v <- rep(NA,n)
  p = length(phi)
  
  c1 = (alpha-L)*(-alpha-1)
  g <- matrix(0, n, length(par) )
  
  for(t in (p+1):n){
    mu[t] = exp( delta + phi%*%log( z[t -1:p] ) )
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    g[t,1] <- digamma(-alpha) - digamma(L - alpha) - log( -alpha-1 ) + alpha/(-alpha-1) + log(v[t]/mu[t]) - (alpha-L)*(mu[t]/v[t])
    g[t,2] <- c1*(mu[t]/v[t]) - alpha
    g[t,3:length(par)] <- (c1*(mu[t]/v[t]) - alpha)%*%log(z[t-1:p])
  }
  U = g[(p+1):n,]
  U
}
grad.g0.MA1.log <- function(z,par,L){
  
  m=q=1
  alpha <- par[1]
  delta <- par[2]
  theta <- par[3]
  
  n = length(z)
  mu = c();mu0 = 1;
  v = rep(NA,n)
  
  mu[1] = exp( delta + theta*( log( mean(z) ) - log(mu0) ) )
  
  deta.ddelta <- rep(0, n)
  deta.dtheta <- rep(0, n)
  
  for(t in (m+1):n){ 
    mu[t] = exp( delta + theta*(log( z[t-q] ) - log( mu[t-q] )) )
    deta.ddelta[t] <- 1 - theta*deta.ddelta[t-q]
    deta.dtheta[t] <- log( z[t-q] ) - log( mu[t-q] ) - theta*deta.dtheta[t-q]
  }
  
  
  t <- (m+1):n
  v[t] =  mu[t]*(-alpha-1) + L*z[t] 
  c1 = (alpha-L)*(-alpha-1)
  g <- matrix(0, length(z)-1, length(par))
  g[,1] <- digamma(-alpha) - digamma(L - alpha) - log( -alpha-1 ) + alpha/(-alpha-1) + log(v[t]/mu[t]) - (alpha-L)*( mu[t]/v[t] )
  g[,2] <- (c1*(mu[t]/v[t]) - alpha)%*%deta.ddelta[t]
  g[,3] <- (c1*(mu[t]/v[t]) - alpha)%*%deta.dtheta[t]
  g
}
grad.g0.MA.log <- function(z,par,L){
  
  alpha <- par[1]; delta <- par[2]; theta <- par[3:length(par)]
  m = q = length(theta)
  
  n = length(z)
  mu = c();mu0 = 1;
  v = rep(NA,n)
  
  mu[1:q] = exp( delta + theta[1]*( log( mean(z) ) - log(mu0) )  + theta[2]*( log( mean(z) ) - log(mu0) ) )

  deta.ddelta <- rep(0, n)
  deta.dtheta <- matrix(0, ncol=q,nrow=n)
  
  c1 = (alpha-L)*(-alpha-1)
  g <- matrix(0, n, length(par))
  
  for(t in (m+1):n){ 
    mu[t] = exp( delta + theta*(log( z[t-q] ) - log( mu[t-q] )) )
    deta.ddelta[t] <- 1 - theta*deta.ddelta[t-q]
    deta.dtheta[t,] <- log( z[t-1:q] ) - log( mu[t-1:q] ) - theta%*%deta.dtheta[t-1:q,] 
  }
  
  t <- (m+1):n
  v[t] =  mu[t]*(-alpha-1) + L*z[t] 
  g[t,1] <- digamma(-alpha) - digamma(L - alpha) - log( -alpha-1 ) + alpha/(-alpha-1) + log(v[t]/mu[t]) - (alpha-L)*( mu[t]/v[t] )
  g[t,2] <- (c1*(mu[t]/v[t]) - alpha)*deta.ddelta[t]
  g[t,3:length(par)] <- (c1*(mu[t]/v[t]) - alpha)%*%deta.dtheta[t,]
  
  U = g[(m+1):n,]
  U
}
grad.g0.ARMA11.log <- function(z,par,L) {
  
  alpha <- par[1]; delta <- par[2]; phi <- par[3]; theta <- par[4]
  
  n = length(z)
  mu=c();mu0 = 1; z0=mean(z)
  v=rep(NA,n)
  p=q=1; m = max(p,q)
  
  mu[1] = exp( delta + phi*log(z0) + theta*(log( z0 ) - log(mu0) ) )
  
  deta.ddelta <- rep(0, n)
  deta.dphi <- rep(0, n)
  deta.dtheta <- rep(0, n)
  
  for(t in (m+1):n)
  {
    mu[t] = exp( delta + phi*log(z[t-1]) + theta*(log( z[t-1] ) - log(mu[t-1]) ) )
    deta.ddelta[t] <- 1 - theta*deta.ddelta[t-q]
    deta.dphi[t] <- log( z[t-p] ) - theta*deta.dphi[t-p]
    deta.dtheta[t] <- log( z[t-q] ) - log( mu[t-q] ) - theta*deta.dtheta[t-q]
  }
  
  t<-(m+1)
  v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  c1 = (alpha-L)*(-alpha-1)
  g <- matrix(0, length(z)-1, 4)
  g[,1] <- digamma(-alpha) - digamma(L - alpha) - log( -alpha-1 ) + alpha/(-alpha-1) + log(v[t]/mu[t]) - (alpha-L)*( mu[t]/v[t] )
  g[,2] <- (c1*(mu[t]/v[t]) - alpha)*deta.ddelta[t]
  g[,3] <- (c1*(mu[t]/v[t]) - alpha)*deta.dphi[t]
  g[,4] <- (c1*(mu[t]/v[t]) - alpha)*deta.dtheta[t]
  g
}


## Gamma

grad.gamma.AR1.boxcox <- function(y,param) {
  L <- param[1]
  delta <- param[2]
  phi <- param[3]
  n = length(y)
  mu=c()
  
  require(forecast)
  lambda = BoxCox.lambda(y)
  
  mu[1] = ( lambda*( delta + phi*log( mean(y) ) ) + 1)^(1/lambda)
  t<-2:n
  mu[t] = ( lambda*( delta + phi*log( y[t-1] ) ) + 1)^(1/lambda)
  
  
  
  g <- matrix(0, length(y)-1, 3)
  g[,1] <- log(L) - digamma(L) + 1 + log(y[2:n]) - y[2:n]/mu[2:n] - log(mu[2:n])
  g[,2] <- L*( y[2:n]/mu[2:n] -1 )
  g[,3] <- L*( ( y[2:n]/mu[2:n] -1 )*log(y[1:(n-1)]) )
  colSums(g)
  g
}
grad.gamma.AR2.boxcox <- function(y,param) {
  L <- param[1]
  delta <- param[2]
  phi <- param[3:4]
  n = length(y)
  mu=c();p=2
  
  require(forecast)
  lambda = BoxCox.lambda(y)
  
  mu0 = 1
  y0 = mean(y)
  
  mu[1] = ( lambda*( delta + phi[1]*( log(y0) ) + phi[2]*( log(y0) ) ) + 1)^(1/lambda)
  mu[2] = ( lambda*( delta + phi[1]*( log(y[1]) ) + phi[2]*( log(y0) ) ) + 1)^(1/lambda)
  
  t<-(p+1):n
  mu[t] = ( lambda*( delta + phi[1]*log(y[t-1]) + phi[2]*( log(y[t-2]) ) ) + 1)^(1/lambda)
  
  g <- matrix(0, n-p, 4)
  g[,1] <- log(L) - digamma(L) + 1 + log(y[t]) - y[t]/mu[t] - log(mu[t])
  g[,2] <- L*( y[t]/mu[t] -1 )
  g[,3] <- L*( ( y[t]/mu[t] -1 )*log(y[t-1]) )
  g[,4] <- L*( ( y[t]/mu[t] -1 )*log(y[t-2]) )
  colSums(g)
  g
}
grad.gamma.MA1.boxcox <- function(y,param) {
  m=q = 1
  L <- param[1]
  delta <- param[2]
  theta <- param[3]
  
  require(forecast)
  lambda = BoxCox.lambda(y)
  
  n = length(y)
  mu = c();mu0 = 1; y0 = 1
  error = c()
  
  mu[1] = ( lambda*( delta + theta*( log( y0 ) - log(mu0) ) ) + 1)^(1/lambda)
  #error[1] = log( y0 ) - log( mu0 )
  
  
  for(t in (m+1):n){ 
    mu[t] = ( lambda*( delta + theta*(log( y[t-q] ) - log( mu[t-q] )) ) + 1)^(1/lambda)
  }
  
  
  deta.ddelta <- rep(0, n)
  deta.dtheta <- rep(0, n)
  
  for(t in (m+1):n)
  {
    deta.ddelta[t] <- 1 - theta*deta.ddelta[t-q]
    deta.dtheta[t] <- log( y[t-q] ) - log( mu[t-q] ) - theta*deta.dtheta[t-q]
  }
  
  t <- (m+1):n
  g <- matrix(0, length(y)-1, 3)
  g[,1] <- log(L) - digamma(L) + 1 + log(y[t]) - y[t]/mu[t] - log(mu[t])
  g[,2] <- L*( y[t]/mu[t] -1 )*(deta.ddelta[t])
  g[,3] <- L*( y[t]/mu[t] -1 )*(deta.dtheta[t])
  colSums(g)
  g
}
grad.gamma.MA.boxcox <- function(y,par,q){
  
  require(forecast)
  lambda = BoxCox.lambda(y)
  
  m=q
  L = par[1]; delta = par[2]; theta = par[3:length(par)]
  n = length(y)
  mu <- rep(1,m);
  g <- matrix(0, n, q+2)
  
  deta.ddelta <- rep(0,n)
  deta.dtheta <- matrix(0, ncol=q,nrow=n)
  
  
  for(t in (m+1):n){
    mu[t] <- ( lambda*( delta + theta%*%( log(y[t-1:q]) - log(mu[t-1:q]) ) ) + 1)^(1/lambda)
    deta.ddelta[t] <- 1 - theta%*%deta.ddelta[t-1:q]
    deta.dtheta[t,] <- log( y[t-1:q] ) - log( mu[t-1:q] ) - theta%*%deta.dtheta[t-1:q,]
  }
  
  t <- (m+1):n
  
  g[t,1] <- log(L) - digamma(L) + 1 + log(y[t]) - y[t]/mu[t] - log(mu[t])
  g[t,2] <- L*( y[t]/mu[t] -1 )%*%(deta.ddelta[t])
  g[t,3:length(par)] <- L*( y[t]/mu[t] -1 )%*%(deta.dtheta[t,])
  
  
  U=g[(m+1):n,]
  return(U)
}
grad.gamma.ARMA11.boxcox <- function(y,param) {
  m = q = p = 1
  L <- param[1]
  delta <- param[2]
  phi <- param[3]
  theta <- param[4]
  
  require(forecast)
  lambda = BoxCox.lambda(y)
  
  n = length(y)
  mu = c();mu0 = 1; y0 = 1
  error = c()
  
  mu[1] = ( lambda*( delta + phi*log(y0) + theta*( log( y0 ) - log(mu0) ) ) + 1)^(1/lambda)
  #error[1] = log( y0 ) - log( mu0 )
  
  deta.ddelta <- rep(0, n)
  deta.dphi <- rep(0, n)
  deta.dtheta <- rep(0, n)
  
  
  for(t in (m+1):n){ 
    
    mu[t] = ( lambda*( delta + phi*log(y[t-p]) + theta*(log( y[t-q] ) - log( mu[t-q] )) ) + 1)^(1/lambda)
    deta.ddelta[t] <- 1 - theta*deta.ddelta[t-q]
    deta.dphi[t] <- log( y[t-p] ) - theta*deta.dphi[t-p]
    deta.dtheta[t] <- log( y[t-q] ) - log( mu[t-q] ) - theta*deta.dtheta[t-q]
    
  }
  
  
  t <- (m+1):n
  g <- matrix(0, length(y)-1, 4)
  g[,1] <- log(L) - digamma(L) + 1 + log(y[t]) - y[t]/mu[t] - log(mu[t])
  g[,2] <- L*( y[t]/mu[t] -1 )*(deta.ddelta[t])
  g[,3] <- L*( y[t]/mu[t] -1 )*(deta.dphi[t])
  g[,4] <- L*( y[t]/mu[t] -1 )*(deta.dtheta[t])
  colSums(g)
  g
}

grad.g0.AR1.boxcox <- function(z,par,L){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  alpha <- par[1]; delta <- par[2]; phi <- par[3]
  
  n = length(z)
  mu = c()
  v = rep(NA,n)
  

  
  mu[1] = ( lambda*( delta + phi*log( mean(z) ) ) + 1)^(1/lambda)
  v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
  
  t<-2:n
  mu[t] = ( lambda*( delta + phi*log( z[t-1] ) ) + 1)^(1/lambda)
  v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  
  
  c1 = (alpha-L)*(-alpha-1)
  g <- matrix(0, length(z)-1, 3)
  g[,1] <- digamma(-alpha) - digamma(L - alpha) - log( -alpha-1 ) + alpha/(-alpha-1) + log(v[t]/mu[t]) - (alpha-L)*(mu[t]/v[t])
  g[,2] <- c1*(mu[t]/v[t]) - alpha
  g[,3] <- (c1*(mu[t]/v[t]) - alpha)*log(z[t-1])
  g
}
grad.g0.AR.boxcox <- function(z,par,L) {
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  alpha <- par[1]; delta <- par[2]; phi <- par[3:length(par)]
  
  n = length(z)
  mu <- v <- rep(NA,n)
  p = length(phi)
  
  c1 = (alpha-L)*(-alpha-1)
  g <- matrix(0, n, length(par) )
  
  for(t in (p+1):n){
    mu[t] = ( lambda*( delta + phi%*%log( z[t -1:p] ) ) + 1)^(1/lambda)
    v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    g[t,1] <- digamma(-alpha) - digamma(L - alpha) - log( -alpha-1 ) + alpha/(-alpha-1) + log(v[t]/mu[t]) - (alpha-L)*(mu[t]/v[t])
    g[t,2] <- c1*(mu[t]/v[t]) - alpha
    g[t,3:length(par)] <- (c1*(mu[t]/v[t]) - alpha)%*%log(z[t-1:p])
  }
  U = g[(p+1):n,]
  U
}
grad.g0.MA1.boxcox <- function(z,par,L){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  m=q=1
  alpha <- par[1]
  delta <- par[2]
  theta <- par[3]
  
  n = length(z)
  mu = c();mu0 = 1;
  v = rep(NA,n)
  
  mu[1] = ( lambda*( delta + theta*( log( mean(z) ) - log(mu0) ) ) + 1)^(1/lambda)
  
  deta.ddelta <- rep(0, n)
  deta.dtheta <- rep(0, n)
  
  for(t in (m+1):n){ 
    mu[t] = ( lambda*( delta + theta*(log( z[t-q] ) - log( mu[t-q] )) ) + 1)^(1/lambda)
    deta.ddelta[t] <- 1 - theta*deta.ddelta[t-q]
    deta.dtheta[t] <- log( z[t-q] ) - log( mu[t-q] ) - theta*deta.dtheta[t-q]
  }
  
  
  t <- (m+1):n
  v[t] =  mu[t]*(-alpha-1) + L*z[t] 
  c1 = (alpha-L)*(-alpha-1)
  g <- matrix(0, length(z)-1, length(par))
  g[,1] <- digamma(-alpha) - digamma(L - alpha) - log( -alpha-1 ) + alpha/(-alpha-1) + log(v[t]/mu[t]) - (alpha-L)*( mu[t]/v[t] )
  g[,2] <- (c1*(mu[t]/v[t]) - alpha)%*%deta.ddelta[t]
  g[,3] <- (c1*(mu[t]/v[t]) - alpha)%*%deta.dtheta[t]
  g
}
grad.g0.MA.boxcox <- function(z,par,L){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  alpha <- par[1]; delta <- par[2]; theta <- par[3:length(par)]
  m = q = length(theta)
  
  n = length(z)
  mu = c();mu0 = 1;
  v = rep(NA,n)
  
  mu[1:q] = ( lambda*( delta + theta[1]*( log( mean(z) ) - log(mu0) )  + theta[2]*( log( mean(z) ) - log(mu0) ) ) + 1)^(1/lambda)
  
  deta.ddelta <- rep(0, n)
  deta.dtheta <- matrix(0, ncol=q,nrow=n)
  
  c1 = (alpha-L)*(-alpha-1)
  g <- matrix(0, n, length(par))
  
  for(t in (m+1):n){ 
    mu[t] = ( lambda*( delta + theta*(log( z[t-q] ) - log( mu[t-q] )) ) + 1)^(1/lambda)
    deta.ddelta[t] <- 1 - theta*deta.ddelta[t-q]
    deta.dtheta[t,] <- log( z[t-1:q] ) - log( mu[t-1:q] ) - theta%*%deta.dtheta[t-1:q,] 
  }
  
  t <- (m+1):n
  v[t] =  mu[t]*(-alpha-1) + L*z[t] 
  g[t,1] <- digamma(-alpha) - digamma(L - alpha) - log( -alpha-1 ) + alpha/(-alpha-1) + log(v[t]/mu[t]) - (alpha-L)*( mu[t]/v[t] )
  g[t,2] <- (c1*(mu[t]/v[t]) - alpha)*deta.ddelta[t]
  g[t,3:length(par)] <- (c1*(mu[t]/v[t]) - alpha)%*%deta.dtheta[t,]
  
  U = g[(m+1):n,]
  U
}
grad.g0.ARMA11.boxcox <- function(z,par,L) {
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  alpha <- par[1]; delta <- par[2]; phi <- par[3]; theta <- par[4]
  
  n = length(z)
  mu=c();mu0 = 1; z0=mean(z)
  v=rep(NA,n)
  p=q=1; m = max(p,q)
  
  mu[1] = ( lambda*( delta + phi*log(z0) + theta*(log( z0 ) - log(mu0) ) ) + 1)^(1/lambda)
  
  deta.ddelta <- rep(0, n)
  deta.dphi <- rep(0, n)
  deta.dtheta <- rep(0, n)
  
  for(t in (m+1):n)
  {
    mu[t] = ( lambda*( delta + phi*log(z[t-1]) + theta*(log( z[t-1] ) - log(mu[t-1]) ) ) + 1)^(1/lambda)
    deta.ddelta[t] <- 1 - theta*deta.ddelta[t-q]
    deta.dphi[t] <- log( z[t-p] ) - theta*deta.dphi[t-p]
    deta.dtheta[t] <- log( z[t-q] ) - log( mu[t-q] ) - theta*deta.dtheta[t-q]
  }
  
  t<-(m+1)
  v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
  c1 = (alpha-L)*(-alpha-1)
  g <- matrix(0, length(z)-1, 4)
  g[,1] <- digamma(-alpha) - digamma(L - alpha) - log( -alpha-1 ) + alpha/(-alpha-1) + log(v[t]/mu[t]) - (alpha-L)*( mu[t]/v[t] )
  g[,2] <- (c1*(mu[t]/v[t]) - alpha)*deta.ddelta[t]
  g[,3] <- (c1*(mu[t]/v[t]) - alpha)*deta.dphi[t]
  g[,4] <- (c1*(mu[t]/v[t]) - alpha)*deta.dtheta[t]
  g
}

#### end vetor gradiente ####


#### Estima??o preliminar #### 


## chute inicial para GAMA

gChute.AR1.log <- function(z){
  
  ols.AR1 <- function(z){
    
    z0 = mean(z);t<-2:length(z)
    Z = data.frame(z1 = c(z0,z[t-1]))
    
    model <- lm(log(z) ~  log(z1), data = Z)
    chute = model$coefficients
    return(chute)
  }
  
  param = ols.AR1(z)
  
  f <- function(L,z, param){
    delta = param[1]; phi = param[2]
    n = length(z)
    mu = c()
    z0 = mean(z) 
    mu[1] = exp(  delta + phi*log( z0 ) )
    
    t <- 2:n 
    mu[t] = exp( delta + phi*log( z[t-1] ) )
    
    U.L = (n-1)*( - digamma(L) + 1 + log(L) ) + sum( log(z/mu) ) - sum(z/mu)   
    
    return(U.L)
    
  }
  
  L.hat = uniroot(f, c(1,100), tol = 0.0001,z=z, param = param,extendInt = "yes")
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  if( L.hat$root < 2 ){
    L.hat$root= 2  
  }else if( L.hat$root > 100 ){
    L.hat$root = 100
  }
  
  chute = c( L.hat = L.hat$root, delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
gChute.AR2.log <- function(z){
  
  ols.AR2 <- function(z){
    
    z0 = mean(z); t<-2:length(z); k<-3:length(z)
    Z = data.frame(z1 = c(z0,z[t-1]),z2 = c(z0,z0,z[k-2]) )
    
    model <- lm(log(z) ~  log(z1) + log(z2), data = Z)
    chute = model$coefficients
    return(chute)
  }
  
  
  
  param = ols.AR2(z)
  
  f <- function(L,z, param){
    delta = param[1]; phin = param[2:length(param)]
    n = length(z)
    mu = c()
    z0 = mean(z) 
    mu[1] = exp(  delta + phin[1]*log( z0 ) + phin[2]*log(z0))
    mu[2] = exp(  delta + phin[1]*log( z[1] ) + phin[2]*log(z0) )
    
    t <- 3:n 
    mu[t] = exp( delta + phin[1]*log( z[t-1] ) + phin[2]*log(z[t-2]) )
    
    U.L = n*( - digamma(L) + 1 + log(L) ) + sum( log(z/mu) - z/mu ) 
    
    return(U.L)
    
  }
  
  L.hat = uniroot(f, c(1,100), tol = 0.0001,z=z, param = param,extendInt = "yes")
  
  if( L.hat$root < 2 ){
    L.hat$root= 2  
  }else if( L.hat$root > 100 ){
    L.hat$root = 100
  }
  
  chute = c( L.hat =L.hat$root, delta.hat = param[1], phi.hat = param[2:3] )
  return(chute)
  
}
gChute.ARP.log <- function(z,p){
  
  ols.ARP <- function(z,p){
    
    n = length(z)
    z = matriz(NA, n, p)
    
    for(i in 1:p){
      t<-(i+1):n
      Z[,i] = c( rep(z0,i), z[t-i]  )
    }
    
    
    model <- lm(log(z) ~  log(Z) )
    chute = model$coefficients
    return(chute)
  }
  
  param = ols.ARP(z,p)
  
  f <- function(L,z, param){
    delta = param[1]; phin = param[2:length(param)]
    n = length(z)
    mu = c()
    z0 = mean(z) 
    
    U.L = n*( - digamma(L) + 1 + log(L) ) + sum( log(z/mu) - z/mu ) 
    
    return(U.L)
    
  }
  
  L.hat = uniroot(f, c(0,100), tol = 0.0001,z=z, param = param,extendInt = "yes")
  
  if( L.hat$root < 2 ){
    L.hat$root= 2  
  }else if( L.hat$root > 100 ){
    L.hat$root = 100
  }
  
  chute = c( L.hat = L.hat$root, delta.hat = param[1], phi.hat = param[2:length(param)] )
  return(chute)
  
} 
gChute.MA1.log <- function(z){
  
  param = c( mean(log(z)), 0  )
  
  f <- function(L,z, param){
    
    delta = param[1]; theta = param[2]
    n = length(z)
    mu = c()
    z0 = mean(z) ; mu0=1
    mu[1] = exp(  delta + theta*( log(z0)  - log(mu0) )  )
    
    
    for(t in 2:n) mu[t] = exp( delta + theta*( log(z[t-1])  - log(mu[t-1]) )  )
    
    U.L = n*( - digamma(L) + 1 + log(L) ) + sum( log(z/mu) ) - sum(z/mu)   
    
    return(U.L)
    
  }
  
  L.hat = uniroot(f, c(1,100), tol = 0.0001,z=z, param = param,extendInt = "yes")
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  if( L.hat$root < 2 ){
    L.hat$root= 2  
  }else if( L.hat$root > 100 ){
    L.hat$root = 100
  }
  
  chute = c( L.hat = L.hat$root, delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
gChute.MA2.log <- function(z){
  
  param = c( mean(log(z)), 0, 0  )
  
  f <- function(L,z, param){
    
    delta = param[1]; theta = param[2:3]
    n = length(z)
    mu = c()
    z0 = mean(z) ; mu0=1
    mu[1:2] = exp(  delta + theta[1]*( log(z0)  - log(mu0) ) + theta[2]*( log(z0)  - log(mu0) )  )
    
    
    for(t in 3:n)
    {
      mu[t] = exp( delta + theta[1]*( log(z[t-1])  - log(mu[t-1]) ) + theta[2]*( log(z[t-2])  - log(mu[t-2]) ) )
    }     
    
    U.L = n*( - digamma(L) + 1 + log(L) ) + sum( log(z/mu) ) - sum(z/mu)   
    
    return(U.L)
    
  }
  
  L.hat = uniroot(f, c(1.5,100), tol = 0.0001,z=z, param = param,extendInt = "yes")
  
  
  if( L.hat$root < 2 ){
    L.hat$root= 2 
  }else if( L.hat$root > 100 ){
    L.hat$root = 100
  }
  
  chute = c( L.hat = L.hat$root, delta.hat = param[1], theta.hat = param[2:3] )
  return(chute)
  
}
gChute.ARMA11.log <- function(z,p){
  
  ols.AR1 <- function(z){
    
    z0 = mean(z);t<-2:length(z)
    Z = data.frame(z1 = c(z0,z[t-1]))
    
    model <- lm(log(z) ~  log(z1), data = Z)
    chute = model$coefficients
    return(chute)
  }
  
  param = c(ols.AR1(z),0)
  
  f <- function(L,z, param){
    delta = param[1]; phi = param[2]; theta = param[3]
    n = length(z);m=1
    mu = c()
    z0 = mean(z) ; mu0 = 1
    mu[1] = exp(  delta + phi*log( z0 ) + theta*( log(z0) - log(mu0)  ) )
    
    for(t in 2:n) mu[t] = exp(  delta + phi*log( z[t-1] ) + theta*( log(z[t-1]) - log(mu[t-1])  ) )
    
    U.L = (n-m)*( - digamma(L) + 1 + log(L) ) + sum( log(z[-1]/mu[-1]) ) - sum( z[-1]/mu[-1] )    
    
    return(U.L)
    
  }
  
  L.hat = uniroot(f, c(1,100), tol = 0.0001,z=z, param = param,extendInt = "yes")
  
  
  if( param[2] < -1 ){
    param[2] = -0.8  
  }else if( param[2] > 1 ){
    param[2] = 0.8
  }
  
  if( L.hat$root < 2 ){
    L.hat$root= 2  
  }else if( L.hat$root > 100 ){
    L.hat$root = 100
  }
  
  chute = c( L.hat = L.hat$root, delta.hat = param[1], phi.hat = param[2],theta.hat=param[3] )
  return(chute)
  
} 

## chute inicial para GAMA

gChute.AR1.boxcox <- function(z){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  ols.AR1 <- function(z){
    
    z0 = mean(z);t<-2:length(z)
    Z = data.frame(z1 = c(z0,z[t-1]))
    
    model <- lm( (z^lambda - 1)/lambda ~  (z1^lambda - 1)/lambda, data = Z)
    chute = model$coefficients
    return(chute)
  }
  
  param = ols.AR1(z)
  
  f <- function(L,z, param){
    delta = param[1]; phi = param[2]
    n = length(z)
    mu = c()
    z0 = mean(z) 
    mu[1] = ( lambda*(  delta + phi*log( z0 ) ) + 1)^(1/lambda)
    
    t <- 2:n 
    mu[t] = ( lambda*( delta + phi*log( z[t-1] ) ) + 1)^(1/lambda)
    
    U.L = (n-1)*( - digamma(L) + 1 + log(L) ) + sum( log(z/mu) ) - sum(z/mu)   
    
    return(U.L)
    
  }
  
  L.hat = uniroot(f, c(1,100), tol = 0.0001,z=z, param = param,extendInt = "yes")
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  if( L.hat$root < 2 ){
    L.hat$root= 2  
  }else if( L.hat$root > 100 ){
    L.hat$root = 100
  }
  
  chute = c( L.hat = L.hat$root, delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
gChute.AR2.boxcox <- function(z){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  ols.AR2 <- function(z){
    
    z0 = mean(z); t<-2:length(z); k<-3:length(z)
    Z = data.frame(z1 = c(z0,z[t-1]),z2 = c(z0,z0,z[k-2]) )
    
    model <- lm((z^lambda - 1)/lambda~  (z1^lambda - 1)/lambda + (z2^lambda - 1)/lambda, data = Z)
    chute = model$coefficients
    return(chute)
  }
  
  
  
  param = ols.AR2(z)
  
  f <- function(L,z, param){
    delta = param[1]; phin = param[2:length(param)]
    n = length(z)
    mu = c()
    z0 = mean(z) 
    mu[1] = ( lambda*(  delta + phin[1]*log( z0 ) + phin[2]*log(z0)) + 1)^(1/lambda)
    mu[2] = ( lambda*(  delta + phin[1]*log( z[1] ) + phin[2]*log(z0) ) + 1)^(1/lambda)
    
    t <- 3:n 
    mu[t] = ( lambda*( delta + phin[1]*log( z[t-1] ) + phin[2]*log(z[t-2]) ) + 1)^(1/lambda)
    
    U.L = n*( - digamma(L) + 1 + log(L) ) + sum( log(z/mu) - z/mu ) 
    
    return(U.L)
    
  }
  
  L.hat = uniroot(f, c(1,100), tol = 0.0001,z=z, param = param,extendInt = "yes")
  
  if( L.hat$root < 2 ){
    L.hat$root= 2  
  }else if( L.hat$root > 100 ){
    L.hat$root = 100
  }
  
  chute = c( L.hat =L.hat$root, delta.hat = param[1], phi.hat = param[2:3] )
  return(chute)
  
}
gChute.MA1.boxcox <- function(z){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  param = c( mean( (z^lambda - 1)/lambda ), 0  )
  
  f <- function(L,z, param){
    
    delta = param[1]; theta = param[2]
    n = length(z)
    mu = c()
    z0 = mean(z) ; mu0=1
    mu[1] = ( lambda*(  delta + theta*( log(z0)  - log(mu0) )  ) + 1)^(1/lambda)
    
    
    for(t in 2:n) mu[t] = ( lambda*( delta + theta*( log(z[t-1])  - log(mu[t-1]) )  ) + 1)^(1/lambda)
    
    U.L = n*( - digamma(L) + 1 + log(L) ) + sum( log(z/mu) ) - sum(z/mu)   
    
    return(U.L)
    
  }
  
  L.hat = uniroot(f, c(1,100), tol = 0.0001,z=z, param = param,extendInt = "yes")
  
  if( param[2] < -1 ){
    param[2] = -0.9  
  }else if( param[2] > 1 ){
    param[2] = 0.9
  }
  
  if( L.hat$root < 2 ){
    L.hat$root= 2  
  }else if( L.hat$root > 100 ){
    L.hat$root = 100
  }
  
  chute = c( L.hat = L.hat$root, delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
gChute.MA2.boxcox<- function(z){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  param = c( mean((z^lambda - 1)/lambda), 0, 0  )
  
  f <- function(L,z, param){
    
    delta = param[1]; theta = param[2:3]
    n = length(z)
    mu = c()
    z0 = mean(z) ; mu0=1
    mu[1:2] = ( lambda*(  delta + theta[1]*( log(z0)  - log(mu0) ) + theta[2]*( log(z0)  - log(mu0) )  ) + 1)^(1/lambda)
    
    
    for(t in 3:n)
    {
      mu[t] = ( lambda*( delta + theta[1]*( log(z[t-1])  - log(mu[t-1]) ) + theta[2]*( log(z[t-2])  - log(mu[t-2]) ) ) + 1)^(1/lambda)
    }     
    
    U.L = n*( - digamma(L) + 1 + log(L) ) + sum( log(z/mu) ) - sum(z/mu)   
    
    return(U.L)
    
  }
  
  L.hat = uniroot(f, c(1.5,100), tol = 0.0001,z=z, param = param,extendInt = "yes")
  
  
  if( L.hat$root < 2 ){
    L.hat$root= 2 
  }else if( L.hat$root > 100 ){
    L.hat$root = 100
  }
  
  chute = c( L.hat = L.hat$root, delta.hat = param[1], theta.hat = param[2:3] )
  return(chute)
  
}
gChute.ARMA11.boxcox <- function(z,p){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  ols.AR1 <- function(z){
    
    z0 = mean(z);t<-2:length(z)
    Z = data.frame(z1 = c(z0,z[t-1]))
    
    model <- lm((z^lambda - 1)/lambda ~  (z1^lambda - 1)/lambda, data = Z)
    chute = model$coefficients
    return(chute)
  }
  
  param = c(ols.AR1(z),0)
  
  f <- function(L,z, param){
    delta = param[1]; phi = param[2]; theta = param[3]
    n = length(z);m=1
    mu = c()
    z0 = mean(z) ; mu0 = 1
    mu[1] = ( lambda*(  delta + phi*log( z0 ) + theta*( log(z0) - log(mu0)  ) ) + 1)^(1/lambda)
    
    for(t in 2:n) mu[t] = ( lambda*(  delta + phi*log( z[t-1] ) + theta*( log(z[t-1]) - log(mu[t-1])  ) ) + 1)^(1/lambda)
    
    U.L = (n-m)*( - digamma(L) + 1 + log(L) ) + sum( log(z[-1]/mu[-1]) ) - sum( z[-1]/mu[-1] )    
    
    return(U.L)
    
  }
  
  L.hat = uniroot(f, c(1,100), tol = 0.0001,z=z, param = param,extendInt = "yes")
  
  
  if( param[2] < -1 ){
    param[2] = -0.8  
  }else if( param[2] > 1 ){
    param[2] = 0.8
  }
  
  if( L.hat$root < 2 ){
    L.hat$root= 2  
  }else if( L.hat$root > 100 ){
    L.hat$root = 100
  }
  
  chute = c( L.hat = L.hat$root, delta.hat = param[1], phi.hat = param[2],theta.hat=param[3] )
  return(chute)
  
} 



ols.AR <- function(z,p){
  n = length(z)
  z0 = mean(z);
  t <- (p+1):length(z)
  Z = matrix(z0,nrow = length(z), ncol = p )
  
  for(t in (p+1):n){
    Z[t,1:p] = z[t-1:p]
  }
  
  X = cbind(rep(1,n), Z)
  
  chute = solve( t(X)%*%X )%*%t(X)%*%log(z)   
  
  return(t(chute))
}
g0Chute.AR1.log <- function(z,L){
  
  p=1
  n = length(z)
  
  z0 = mean(z);t<-2:length(z)
  Z = data.frame(z1 = c(z0,z[t-1]))
  
  model <- lm(log(z) ~  log(z1), data = Z)
  param = model$coefficients
  
  if( param[2] < -1 ){
    param[2] = -0.95  
  }else if( param[2] > 1 ){
    param[2] = 0.95
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
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00002, z=z, param = param, L=L,extendInt = "yes")
  
  alpha.hat = alpha.hat$root
  
  if( alpha.hat < -50 ){
    alpha.hat = -50
  }else if( alpha.hat > -1.1 ){
    alpha.hat = -1.2
  }
  
  
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
g0Chute.AR2.log <- function(z,L){
  
  p=2
  n = length(z)
  
  param = ols.AR(z,p)
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2:length(param)]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) 
    
    mu[1] = exp(  delta + phi[1]*log( z0 ) + phi[2]*log(z0) )
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    mu[2] = exp(  delta + phi[1]*log( z[1] ) + phi[2]*log(z0) )
    v[2] = ( mu[2]*(-alpha-1) + L*z[2] )
    
    
    for(t in (p+1):n){
      mu[t] = exp( delta + phi[1]*log( z[t-1] ) + phi[2]*log( z[t-2] ) )
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    #vn = v[(p+1):n]
    #mun = mu[(p+1):n]
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, c(-50,-1.1),tol = 0.00001, z=z, param = param,L=L,extendInt = "yes")
  
  
  chute = c( alpha.hat = alpha.hat$root, delta.hat = param[1], phi.hat = param[2:length(param)] )
  return(chute)
  
}
g0Chute.MA1.log <- function(z,L){
  
  n = length(z)
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
  
  
  chute = c( alpha.hat = ifelse(alpha.hat$root < -50, -50,alpha.hat$root), delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
g0Chute.MA2.log <- function(z,L){
  
  n = length(z)
  param = c( mean(log(z)), 0 , 0 )
  
  f <- function(alpha,z,param,L){
    delta = param[1]; theta = param[2:length(param)]
    n = length(z); q = length(theta)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    
    mu[1:q] = exp(  delta + theta[1]*( log( z0 ) - log(mu0)  ) + theta[2]*( log( z0 ) - log(mu0)  ) )
    v[1:q] = ( mu[1:q]*(-alpha-1) + L*z[1:q] )
    
    for(t in 2:n){
      mu[t] = exp(  delta + theta[1]*( log( z[t-1] ) - log(mu[t-1])  ) + theta[2]*( log( z[t-2] ) - log(mu[t-2])  ) )
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  
  
  chute = c( alpha.hat = ifelse(alpha.hat$root < -50, -50,alpha.hat$root), delta.hat = param[1], theta.hat = param[2:3] )
  return(chute)
  
}
g0Chute.ARMA11.log <- function(z,L){
  
  p=1
  n = length(z)
  
  param = c(ols.AR(z,p),0)
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2];theta = param[3]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    
    mu[1] = exp(  delta + phi*log(z0) + theta*( log( z0 ) - log(mu0)  ) )
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      mu[t] = exp(  delta + phi*log(z[t-1])  + theta*( log( z[t-1] ) - log(mu[t-1])  ) )
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  
  if(alpha.hat$root < -50){
    alpha.hat$root = -50
  }else if(alpha.hat$root > -1.1){
    alpha.hat$root = -1.1
  }
  
  if(param[2] < -0.98){
    param[2] = -0.9
  }else if(param[2] > 0.98){
    param[2] = 0.9
  }
  
  chute = c( alpha.hat = alpha.hat$root, delta.hat = param[1], phi.hat = param[2],theta.hat = param[3] )
  return(chute)
  
}

ols.AR.boxcox <- function(z,p){
  n = length(z)
  z0 = mean(z);
  t <- (p+1):length(z)
  Z = matrix(z0,nrow = length(z), ncol = p )
  
  for(t in (p+1):n){
    Z[t,1:p] = z[t-1:p]
  }
  
  X = cbind(rep(1,n), Z)
  
  chute = solve( t(X)%*%X )%*%t(X)%*%((z^lambda - 1)/lambda)
  
  return(t(chute))
}
g0Chute.AR1.boxcox <- function(z,L){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  p=1
  n = length(z)
  
  z0 = mean(z);t<-2:length(z)
  Z = data.frame(z1 = c(z0,z[t-1]))
  
  model <- lm((z^lambda - 1)/lambda ~  (z1^lambda - 1)/lambda, data = Z)
  param = model$coefficients
  
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
     
    mu[1] = ( lambda*(  delta + phi*log( z0 ) ) + 1)^(1/lambda)
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      mu[t] = ( lambda*( delta + phi*log( z[t-1] ) ) + 1)^(1/lambda)
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-100,-1.1), tol = 0.00002, z=z, param = param, L=L,extendInt = "yes")
  
  alpha.hat = alpha.hat$root
  
  if( alpha.hat < -100 ){
    alpha.hat = -100
  }else if( alpha.hat > -1.1 ){
    alpha.hat = -1.2
  }
  
  
  
  chute = c( alpha.hat = alpha.hat, delta.hat = param[1], phi.hat = param[2] )
  return(chute)
  
}
g0Chute.AR2.boxcox <- function(z,L){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  p=2
  n = length(z)
  
  param = ols.AR(z,p)
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2:length(param)]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) 
    
    mu[1] = ( lambda*(  delta + phi[1]*log( z0 ) + phi[2]*log(z0) ) + 1)^(1/lambda)
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    mu[2] = ( lambda*(  delta + phi[1]*log( z[1] ) + phi[2]*log(z0) ) + 1)^(1/lambda)
    v[2] = ( mu[2]*(-alpha-1) + L*z[2] )
    
    
    for(t in (p+1):n){
      mu[t] = ( lambda*( delta + phi[1]*log( z[t-1] ) + phi[2]*log( z[t-2] ) ) + 1)^(1/lambda)
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    #vn = v[(p+1):n]
    #mun = mu[(p+1):n]
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, c(-50,-1.1),tol = 0.00001, z=z, param = param,L=L,extendInt = "yes")
  
  
  chute = c( alpha.hat = alpha.hat$root, delta.hat = param[1], phi.hat = param[2:length(param)] )
  return(chute)
  
}
g0Chute.MA1.boxcox <- function(z,L){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  n = length(z)
  param = c( mean((z^lambda - 1)/lambda), 0  )
  
  f <- function(alpha,z,param,L){
    delta = param[1]; theta = param[2]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    
    mu[1] = ( lambda*(  delta + theta*( log( z0 ) - log(mu0)  ) ) + 1)^(1/lambda)
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      mu[t] = ( lambda*(  delta + theta*( log( z[t-1] ) - log(mu[t-1])  ) ) + 1)^(1/lambda)
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  
  
  chute = c( alpha.hat = ifelse(alpha.hat$root < -50, -50,alpha.hat$root), delta.hat = param[1], theta.hat = param[2] )
  return(chute)
  
}
g0Chute.MA2.boxcox <- function(z,L){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  n = length(z)
  param = c( mean((z^lambda - 1)/lambda), 0 , 0 )
  
  f <- function(alpha,z,param,L){
    delta = param[1]; theta = param[2:length(param)]
    n = length(z); q = length(theta)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    
    mu[1:q] = ( lambda*(  delta + theta[1]*( log( z0 ) - log(mu0)  ) + theta[2]*( log( z0 ) - log(mu0)  ) ) + 1)^(1/lambda)
    v[1:q] = ( mu[1:q]*(-alpha-1) + L*z[1:q] )
    
    for(t in 2:n){
      mu[t] = ( lambda*(  delta + theta[1]*( log( z[t-1] ) - log(mu[t-1])  ) + theta[2]*( log( z[t-2] ) - log(mu[t-2])  ) ) + 1)^(1/lambda)
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  
  
  chute = c( alpha.hat = ifelse(alpha.hat$root < -50, -50,alpha.hat$root), delta.hat = param[1], theta.hat = param[2:3] )
  return(chute)
  
}
g0Chute.ARMA11.boxcox <- function(z,L){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  p=1
  n = length(z)
  
  param = c(ols.AR(z,p),0)
  
  f <- function(alpha,z,param,L){
    delta = param[1]; phi = param[2];theta = param[3]
    n = length(z)
    mu = c();v=c()
    z0 = mean(z) ; mu0 = 1
    
    mu[1] = ( lambda*(  delta + phi*log(z0) + theta*( log( z0 ) - log(mu0)  ) ) + 1)^(1/lambda)
    v[1] = ( mu[1]*(-alpha-1) + L*z[1] )
    
    for(t in 2:n){
      mu[t] = ( lambda*(  delta + phi*log(z[t-1])  + theta*( log( z[t-1] ) - log(mu[t-1])  ) ) + 1)^(1/lambda)
      v[t] = ( mu[t]*(-alpha-1) + L*z[t] )
    }
    
    
    U.alpha = n*( digamma(-alpha) - digamma(L-alpha) + alpha/(-alpha-1) - log( -alpha-1 ) ) + sum( log(v/mu) ) - (alpha - L)*sum( mu/v )
    
    return(U.alpha)
    
  }
  
  alpha.hat = uniroot(f, interval = c(-50,-1.1), tol = 0.00001, z=z, param = param, L=L,extendInt = "yes")
  
  if(alpha.hat$root < -50){
    alpha.hat$root = -50
  }else if(alpha.hat$root > -1.1){
    alpha.hat$root = -1.1
  }
  
  if(param[2] < -0.98){
    param[2] = -0.9
  }else if(param[2] > 0.98){
    param[2] = 0.9
  }
  
  chute = c( alpha.hat = alpha.hat$root, delta.hat = param[1], phi.hat = param[2],theta.hat = param[3] )
  return(chute)
  
}

## chute inicial dos param?tros para a distribui??o g0





#### end estima??o preliminar ####

### Criterios de Seleção 

SSIM <- function(x,y,k1,k2,L){
  C1 = (k1*L)^2
  C2 = (k2*L)^2  
  N = length( c(na.remove(y)) )
  
  mu.x = mean(x, na.rm =  TRUE)
  mu.y = mean(y, na.rm =  TRUE)
  sigma.x = sqrt( (1/(N-1))*sum( (x - mu.x)^2, na.rm = TRUE) ) 
  sigma.y = sqrt( (1/(N-1))*sum( (y - mu.y)^2, na.rm = TRUE ) )
  
  med.x = median(x, na.rm = TRUE) # mediana
  med.y = median(y, na.rm = TRUE) # mediana 
  mad.x = median( abs( x - med.x ), na.rm = TRUE ) # MAD - equivalente ao sd
  mad.y =  median( abs( y - med.y ), na.rm = TRUE )  # MAD - equivalente ao sd
  
  sigma.xy =  (1/(N-1))*sum( (x - mu.x)*(y - mu.y), na.rm = TRUE ) 
  mad.xy =  median( abs( (x - med.x)*(y - med.y) ), na.rm = TRUE ) 
  
  
  SSIM = ( (2*mu.x*mu.y + C1)*(2*sigma.xy + C2)  )/( (mu.x^2 + mu.y^2 + C1)*(sigma.x^2 + sigma.y^2 + C2) )
  
  SSIM.robust = ( (2*med.x*med.y + C1)*(2*mad.xy + C2)  )/( (med.x^2 + med.y^2 + C1)*(mad.x^2 + mad.y^2 + C2) )
  
  list(SSIM = SSIM, SSIM.robust=SSIM.robust)
  
}

acuracia <- function(y, y.hat){
  n=length(c(na.remove(y.hat)) )
  
  MSE = mean((y - y.hat)^2, na.rm = TRUE)
  MAE = mean(abs(y - y.hat), na.rm = TRUE)
  RMSE = sqrt(mean((y - y.hat)^2, na.rm = TRUE))
  MSE.rob = median((y - y.hat)^2, na.rm = TRUE)
  MAE.rob = median(abs(y - y.hat), na.rm = TRUE)
  RMSE.rob = sqrt(median((y - y.hat)^2, na.rm = TRUE))
  SSIM = SSIM(y,y.hat,k1=0,k2=0,L=1)$SSIM
  SSIM.robust = SSIM(y,y.hat,k1=0,k2=0,L=1)$SSIM.robust
  beta = sum( (y - mean(y, na.rm = TRUE))*(y.hat - mean(y.hat, na.rm = TRUE)), na.rm = TRUE )/sqrt( sum( (y - mean(y, na.rm = TRUE))^2, na.rm = TRUE )*sum( (y.hat - mean(y.hat, na.rm = TRUE))^2, na.rm = TRUE )  )
  
  
  acuracia = c(SSIM = SSIM, beta=beta, MSE= MSE , MAE = MAE, RMSE = RMSE, SSIM.robust = SSIM.robust, MSE.rob= MSE.rob , MAE.rob = MAE.rob, RMSE.rob = RMSE.rob  )
  round(acuracia, 4)
}


### Estimacao 


cml.AR1.log <- function(z ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(-1,200,1,1)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR1(z) 
  lk.gammanew = function(param)  lk.gamma.AR1(z, param)
  grad.gammanew = function(param)  grad.gamma.AR1(z, param)
  try(
    res1 <- maxLik(lk.gammanew, grad = grad.gammanew, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  )
  est = res1$estimate
  L = est[1];delta=est[2];phi=est[3]
  
  mu.hat[1] = exp(  delta + phi*log( z0 ) )
  
  t <- 2:n 
  mu.hat[t] = exp( delta + phi*log( z[t-1] ) )
  
  AIC = -2*lk.gamma.AR1(z, est) + 2*length(est) 
  BIC = -2*lk.gamma.AR1(z, est) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )
  
}
cml.AR2.log <- function(z){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  
  A1 <- matrix(c(1, -1, 0,0,0,0,0,0 ), 2,4)
  B1 <- c(-1,200)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR2(z) 
  lk.gammanew = function(param)  lk.gamma.AR2(z, param)
  grad.gammanew = function(param)  grad.gamma.AR2(z, param)
  res1 <- maxLik(lk.gammanew, grad = grad.gammanew, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  est = res1$estimate
  L = est[1];delta=est[2];phi=est[3:4]
  
  mu.hat[1] = exp( delta + phi[1]*( log(z0) ) + phi[2]*( log(z0) ) )
  mu.hat[2] = exp( delta + phi[1]*( log(z[1]) ) + phi[2]*( log(z0) ) )
  
  t <- 3:n 
  mu.hat[t] = exp( delta + phi[1]*log( z[t-1] ) + phi[2]*log( z[t-2] )  )
  
  AIC = -2*lk.gamma.AR2(z, est) + 2*length(est) 
  BIC = -2*lk.gamma.AR2(z, est) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )
}
cml.MA1.log <- function(z ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(-1,200,1,1)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA1(z) 
  lk.new = function(param)  lk.gamma.MA1(z, param)
  grad.new = function(param)  grad.gamma.MA1(z, param)
  try(
    res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  )
  est = res1$estimate
  L = est[1];delta=est[2];theta=est[3]
  
  z0 = mean(z); mu0=1
  
  mu.hat[1] = exp( delta + theta*( log(z0) - log(mu0) ) )
  
  for(t in 2:n) mu.hat[t] = exp( delta + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) )
  
  AIC = -2*lk.gamma.MA1(z, est) + 2*length(est) 
  BIC = -2*lk.gamma.MA1(z, est) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
cml.MA2.log <- function(z ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(purrr)
  
  A1 <- matrix(c(1, -1, 0,0,0,0,0,0 ), 2,4)
  B1 <- c(-1,200)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA2(z) 
  lk.new = function(param)  lk.gamma.MA2(z, param)
  grad.new = function(param)  grad.gamma.MA(z, param,2)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = grad.new, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = grad.new, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "CG", constraints=list(ineqA = A1, ineqB = B1))
  }else{
    res1 <- res1$result
  }
  
  if(res1$code == 100){
    est = rep(NA,length(param)); mu.hat = rep(NA,length(z))
  }else{
    
    est = res1$estimate
    L = est[1];delta=est[2];theta=est[3:4]
    
    z0 = mean(z); mu0=1
    mu.hat[1:2] = exp( delta + theta[1]*( log(z0) - log(mu0) ) + theta[2]*( log(z0) - log(mu0) ) )
    
    for(t in 3:n) mu.hat[t] = exp( delta + theta[1]*( log(z[t-1]) - log(mu.hat[t-1]) ) + theta[2]*( log(z[t-2]) - log(mu.hat[t-2]) ) )
    
  }
  
  AIC = -2*lk.gamma.MA2(z, est) + 2*length(est) 
  BIC = -2*lk.gamma.MA2(z, est) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
cml.ARMA11.log <- function(z ){
  require(purrr)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1, 0,0,0,0, 0,0,0,0,0,0, 0,0,1,-1,0,0, 0,0,0,0,1,-1 ), 6,4)
  B1 <- c(-1,200,1,1,1,1)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.ARMA11(z) 
  lk.new = function(param)  lk.gamma.ARMA11(z, param)
  grad.new = function(param)  grad.gamma.ARMA11(z, param)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = grad.new, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = grad.new, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "CG", constraints=list(ineqA = A1, ineqB = B1))
  }else{
    res1 <- res1$result
  }
  
  
  est = res1$estimate
  L = est[1];delta=est[2];phi=est[3];theta=est[4] 
  
  z0 = mean(z); mu0=1
  
  mu.hat[1] = exp( delta + phi*log( z0 ) + theta*( log(z0) - log(mu0) ) )
  
  for(t in 2:n) mu.hat[t] = exp( delta + phi*log( z[t-1] ) + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) )
  
  AIC = -2*lk.gamma.ARMA11(z, est) + 2*length(est) 
  BIC = -2*lk.gamma.ARMA11(z, est) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
auto.gamma.log <- function(y, criterio){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  
  ar1 = cml.AR1(y);ar2 = cml.AR2(y)
  ma1 = cml.MA1(y); ma2 = cml.MA2(y)
  arma11 = cml.ARMA11(y)
  
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
    
  }else if(criterio == "UQI"){
    
    criter.ar1 =  data.frame(UQI=ar1$UQI,model="AR1")
    criter.ar2 =  data.frame(UQI=ar2$UQI,model="AR2")
    criter.ma1 =  data.frame(UQI=ma1$UQI,model="MA1")
    criter.ma2 =  data.frame(UQI=ma2$UQI,model="MA2")
    criter.arma11 =  data.frame(UQI=arma11$UQI,model="ARMA11")     
    
  }else if(criterio == "beta"){
    
    criter.ar1 =  data.frame(beta=ar1$beta,model="AR1")
    criter.ar2 =  data.frame(beta=ar2$beta,model="AR2")
    criter.ma1 =  data.frame(beta=ma1$beta,model="MA1")
    criter.ma2 =  data.frame(beta=ma2$beta,model="MA2")
    criter.arma11 =  data.frame(beta=arma11$beta,model="ARMA11")     
    
  }else if(criterio == "MSE"){
    
    criter.ar1 =  data.frame(MSE=ar1$MSE,model="AR1")
    criter.ar2 =  data.frame(MSE=ar2$MSE,model="AR2")
    criter.ma1 =  data.frame(MSE=ma1$MSE,model="MA1")
    criter.ma2 =  data.frame(MSE=ma2$MSE,model="MA2")
    criter.arma11 =  data.frame(MSE=arma11$MSE,model="ARMA11")     
    
  }else if(criterio == "MAE"){
    
    criter.ar1 =  data.frame(MAE=ar1$MAE,model="AR1")
    criter.ar2 =  data.frame(MAE=ar2$MAE,model="AR2")
    criter.ma1 =  data.frame(MAE=ma1$MAE,model="MA1")
    criter.ma2 =  data.frame(MAE=ma2$MAE,model="MA2")
    criter.arma11 =  data.frame(MAE=arma11$MAE,model="ARMA11")     
    
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

cml.AR1.boxcox <- function(z ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(-1,200,1,1)
  
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR1.boxcox(z) 
  lk.gammanew = function(param)  lk.gamma.AR1.boxcox(z, param)
  grad.gammanew = function(param)  grad.gamma.AR1.boxcox(z, param)
  try(
    res1 <- maxLik(lk.gammanew, grad = grad.gammanew, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  )
  est = res1$estimate
  L = est[1];delta=est[2];phi=est[3]
  
  mu.hat[1] = ( lambda*(  delta + phi*log( z0 ) ) + 1)^(1/lambda)
  
  t <- 2:n 
  mu.hat[t] = ( lambda*( delta + phi*log( z[t-1] ) ) + 1)^(1/lambda)
  
  AIC = -2*lk.gamma.AR1.boxcox(z, est) + 2*length(est) 
  BIC = -2*lk.gamma.AR1.boxcox(z, est) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )
  
}
cml.AR2.boxcox <- function(z){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  A1 <- matrix(c(1, -1, 0,0,0,0,0,0 ), 2,4)
  B1 <- c(-1,200)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.AR2(z) 
  lk.gammanew = function(param)  lk.gamma.AR2.boxcox(z, param)
  grad.gammanew = function(param)  grad.gamma.AR2.boxcox(z, param)
  res1 <- maxLik(lk.gammanew, grad = grad.gammanew, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  est = res1$estimate
  L = est[1];delta=est[2];phi=est[3:4]
  
  mu.hat[1] = ( lambda*( delta + phi[1]*( log(z0) ) + phi[2]*( log(z0) ) ) + 1)^(1/lambda)
  mu.hat[2] = ( lambda*( delta + phi[1]*( log(z[1]) ) + phi[2]*( log(z0) ) ) + 1)^(1/lambda)
  
  t <- 3:n 
  mu.hat[t] = ( lambda*( delta + phi[1]*log( z[t-1] ) + phi[2]*log( z[t-2] )  ) + 1)^(1/lambda)
  
  AIC = -2*lk.gamma.AR2.boxcox(z, est) + 2*length(est) 
  BIC = -2*lk.gamma.AR2.boxcox(z, est) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )
}
cml.MA1.boxcox <- function(z ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(-1,200,1,1)
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA1.boxcox(z) 
  lk.new = function(param)  lk.gamma.MA1.boxcox(z, param)
  grad.new = function(param)  grad.gamma.MA1.boxcox(z, param)
  try(
    res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  )
  est = res1$estimate
  L = est[1];delta=est[2];theta=est[3]
  
  z0 = mean(z); mu0=1
  
  mu.hat[1] = ( lambda*( delta + theta*( log(z0) - log(mu0) ) ) + 1)^(1/lambda)
  
  for(t in 2:n) mu.hat[t] = ( lambda*( delta + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) ) + 1)^(1/lambda)
  
  AIC = -2*lk.gamma.MA1.boxcox(z, est) + 2*length(est) 
  BIC = -2*lk.gamma.MA1.boxcox(z, est) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
cml.MA2.boxcox <- function(z ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  require(purrr)
  
  A1 <- matrix(c(1, -1, 0,0,0,0,0,0 ), 2,4)
  B1 <- c(-1,200)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.MA2.boxcox(z) 
  lk.new = function(param)  lk.gamma.MA2.boxcox(z, param)
  grad.new = function(param)  grad.gamma.MA.boxcox(z, param,2)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = grad.new, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = grad.new, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "CG", constraints=list(ineqA = A1, ineqB = B1))
  }else{
    res1 <- res1$result
  }
  
  if(res1$code == 100){
    est = rep(NA,length(param)); mu.hat = rep(NA,length(z))
  }else{
    
    est = res1$estimate
    L = est[1];delta=est[2];theta=est[3:4]
    
    z0 = mean(z); mu0=1
    mu.hat[1:2] = ( lambda*( delta + theta[1]*( log(z0) - log(mu0) ) + theta[2]*( log(z0) - log(mu0) ) ) + 1)^(1/lambda)
    
    for(t in 3:n) mu.hat[t] = ( lambda*( delta + theta[1]*( log(z[t-1]) - log(mu.hat[t-1]) ) + theta[2]*( log(z[t-2]) - log(mu.hat[t-2]) ) ) + 1)^(1/lambda)
    
  }
  
  AIC = -2*lk.gamma.MA2.boxcox(z, est) + 2*length(est) 
  BIC = -2*lk.gamma.MA2.boxcox(z, est) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
cml.ARMA11.boxcox <- function(z ){
  require(purrr)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  A1 <- matrix(c(1, -1, 0,0,0,0, 0,0,0,0,0,0, 0,0,1,-1,0,0, 0,0,0,0,1,-1 ), 6,4)
  B1 <- c(-1,200,1,1,1,1)
  
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = gChute.ARMA11.boxcox(z) 
  lk.new = function(param)  lk.gamma.ARMA11.boxcox(z, param)
  grad.new = function(param)  grad.gamma.ARMA11.boxcox(z, param)
  maxlik.new = function(logLik,grad,hess,start,method,constraints)  maxLik(logLik = lk.new, grad = grad.new, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  
  safe <- safely(maxlik.new, otherwise = "ERROR")
  
  res1 <- safe(
    maxLik(logLik = lk.new, grad = grad.new, start=param, method = "BFGS", constraints=list(ineqA = A1, ineqB = B1))
  )
  
  
  if(res1$result[1] == "ERROR"){
    res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "CG", constraints=list(ineqA = A1, ineqB = B1))
  }else{
    res1 <- res1$result
  }
  
  
  est = res1$estimate
  L = est[1];delta=est[2];phi=est[3];theta=est[4] 
  
  z0 = mean(z); mu0=1
  
  mu.hat[1] = ( lambda*( delta + phi*log( z0 ) + theta*( log(z0) - log(mu0) ) ) + 1)^(1/lambda)
  
  for(t in 2:n) mu.hat[t] = ( lambda*( delta + phi*log( z[t-1] ) + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) ) + 1)^(1/lambda)
  
  AIC = -2*lk.gamma.ARMA11.boxcox(z, est) + 2*length(est) 
  BIC = -2*lk.gamma.ARMA11.boxcox(z, est) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
auto.gamma.boxcox <- function(y, criterio){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  
  ar1 = cml.AR1.boxcox(y);ar2 = cml.AR2.boxcox(y)
  ma1 = cml.MA1.boxcox(y); ma2 = cml.MA2.boxcox(y)
  arma11 = cml.ARMA11.boxcox(y)
  
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
    
  }else if(criterio == "UQI"){
    
    criter.ar1 =  data.frame(UQI=ar1$UQI,model="AR1")
    criter.ar2 =  data.frame(UQI=ar2$UQI,model="AR2")
    criter.ma1 =  data.frame(UQI=ma1$UQI,model="MA1")
    criter.ma2 =  data.frame(UQI=ma2$UQI,model="MA2")
    criter.arma11 =  data.frame(UQI=arma11$UQI,model="ARMA11")     
    
  }else if(criterio == "beta"){
    
    criter.ar1 =  data.frame(beta=ar1$beta,model="AR1")
    criter.ar2 =  data.frame(beta=ar2$beta,model="AR2")
    criter.ma1 =  data.frame(beta=ma1$beta,model="MA1")
    criter.ma2 =  data.frame(beta=ma2$beta,model="MA2")
    criter.arma11 =  data.frame(beta=arma11$beta,model="ARMA11")     
    
  }else if(criterio == "MSE"){
    
    criter.ar1 =  data.frame(MSE=ar1$MSE,model="AR1")
    criter.ar2 =  data.frame(MSE=ar2$MSE,model="AR2")
    criter.ma1 =  data.frame(MSE=ma1$MSE,model="MA1")
    criter.ma2 =  data.frame(MSE=ma2$MSE,model="MA2")
    criter.arma11 =  data.frame(MSE=arma11$MSE,model="ARMA11")     
    
  }else if(criterio == "MAE"){
    
    criter.ar1 =  data.frame(MAE=ar1$MAE,model="AR1")
    criter.ar2 =  data.frame(MAE=ar2$MAE,model="AR2")
    criter.ma1 =  data.frame(MAE=ma1$MAE,model="MA1")
    criter.ma2 =  data.frame(MAE=ma2$MAE,model="MA2")
    criter.arma11 =  data.frame(MAE=arma11$MAE,model="ARMA11")     
    
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


cml.g0.AR1.log <- function(z, L ){
  
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(2000,-1,1,1)
  # constraints=list(ineqA = A1, ineqB = B1)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR1.log(z,L) 
  lk.new = function(param)  lk.g0.AR1.log(z, param,L)
  grad.new = function(param)  grad.g0.AR1.log(z, param,L)
  (res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "BFGS" , constraints=list(ineqA = A1, ineqB = B1)))
  est = res1$estimate
  code = res1$code
  alpha = est[1]; delta=est[2]; phi=est[3]
  
  mu.hat[1] = exp(  delta + phi*log( z0 ) )
  
  for(t in 2:n) mu.hat[t] = exp( delta + phi*log( z[t-1] ) )
  
  AIC = -2*lk.g0.AR1.log(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.AR1.log(z, est,L) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE, convergence = code) )
  
}
 cml.g0.AR2.log <- function(z,L){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR2.log(z,L) 
  lk.new = function(param)  lk.g0.AR2.log(z, param,L)
  grad.new = function(param)  grad.g0.AR.log(z, param,L)
  (res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "BFGS"))
  est = res1$estimate
  alpha = est[1];delta=est[2];phi=est[3:4]
  
  mu.hat[1] = exp( delta + phi[1]*( log(z0) ) + phi[2]*( log(z0) ) )
  mu.hat[2] = exp( delta + phi[1]*( log(z[1]) ) + phi[2]*( log(z0) ) )
  
  for(t in 3:n) mu.hat[t] = exp( delta + phi[1]*log( z[t-1] ) + phi[2]*log( z[t-2] )  )
  
  AIC = -2*lk.g0.AR2.log(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.AR2.log(z, est,L) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )
}
cml.g0.MA1.log <- function(z,L ){
  
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(2000,-1.01,1,1)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA1.log(z, L) 
  lk.new = function(param)  lk.g0.MA1.log(z, param, L)
  grad.new = function(param)  grad.g0.MA1.log(z, param, L)
  (res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "BFGS",constraints=list(ineqA = A1, ineqB = B1)))
  est = res1$estimate
  alpha = est[1];delta=est[2];theta=est[3]
  
  z0 = mean(z); mu0=1
  
  mu.hat[1] = exp( delta + theta*( log(z0) - log(mu0) ) )
  
  for(t in 2:n) mu.hat[t] = exp( delta + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) )
  
  AIC = -2*lk.g0.MA1.log(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA1.log(z, est,L) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4]  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
cml.g0.MA2.log <- function(z,L ){
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA2.log(z,L) 
  lk.new = function(param)  lk.g0.MA2.log(z, param,L)
  grad.new = function(param)  grad.g0.MA.log(z, param,L)
  try(
    (res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS"))
  )
  
  if(res1$code == 100){
    est = rep(NA,length(param)); mu.hat = rep(NA,length(z))
  }else{
    
    est = res1$estimate
    alpha = est[1];delta=est[2];theta=est[3:4]
    
    z0 = mean(z); mu0=1
    mu.hat[1:2] = exp( delta + theta[1]*( log(z0) - log(mu0) ) + theta[2]*( log(z0) - log(mu0) ) )
    
    for(t in 3:n) mu.hat[t] = exp( delta + theta[1]*( log(z[t-1]) - log(mu.hat[t-1]) ) + theta[2]*( log(z[t-2]) - log(mu.hat[t-2]) ) )
    
  }
  
  AIC = -2*lk.g0.MA2.log(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA2.log(z, est,L) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4]  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
cml.g0.ARMA11.log <- function(z,L ){
  
  
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.ARMA11.log(z,L) 
  lk.new = function(param)  lk.g0.ARMA11.log(z, param,L)
  grad.new = function(param)  grad.g0.ARMA11.log(z, param,L)
  (res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "BFGS"))
  est = res1$estimate
  alpha = est[1];delta=est[2];phi=est[3];theta=est[4] 
  
  z0 = mean(z); mu0=1
  
  mu.hat[1] = exp( delta + phi*log( z0 ) + theta*( log(z0) - log(mu0) ) )
  
  for(t in 2:n) mu.hat[t] = exp( delta + phi*log( z[t-1] ) + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) )
  
  AIC = -2*lk.g0.ARMA11.log(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.ARMA11.log(z, est,L) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
auto.g0.log <- function(y,L, criterio){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  
  ar1 = cml.g0.AR1.log(y,L);ar2 = cml.g0.AR2.log(y,L)
  ma1 = cml.g0.MA1.log(y,L); ma2 = cml.g0.MA2.log(y,L)
  arma11 = cml.g0.ARMA11.log(y,L)
  
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
    
  }else if(criterio == "UQI"){
    
    criter.ar1 =  data.frame(UQI=ar1$UQI,model="AR1")
    criter.ar2 =  data.frame(UQI=ar2$UQI,model="AR2")
    criter.ma1 =  data.frame(UQI=ma1$UQI,model="MA1")
    criter.ma2 =  data.frame(UQI=ma2$UQI,model="MA2")
    criter.arma11 =  data.frame(UQI=arma11$UQI,model="ARMA11")     
    
  }else if(criterio == "beta"){
    
    criter.ar1 =  data.frame(beta=ar1$beta,model="AR1")
    criter.ar2 =  data.frame(beta=ar2$beta,model="AR2")
    criter.ma1 =  data.frame(beta=ma1$beta,model="MA1")
    criter.ma2 =  data.frame(beta=ma2$beta,model="MA2")
    criter.arma11 =  data.frame(beta=arma11$beta,model="ARMA11")     
    
  }else if(criterio == "MSE"){
    
    criter.ar1 =  data.frame(MSE=ar1$MSE,model="AR1")
    criter.ar2 =  data.frame(MSE=ar2$MSE,model="AR2")
    criter.ma1 =  data.frame(MSE=ma1$MSE,model="MA1")
    criter.ma2 =  data.frame(MSE=ma2$MSE,model="MA2")
    criter.arma11 =  data.frame(MSE=arma11$MSE,model="ARMA11")     
    
  }else if(criterio == "MAE"){
    
    criter.ar1 =  data.frame(MAE=ar1$MAE,model="AR1")
    criter.ar2 =  data.frame(MAE=ar2$MAE,model="AR2")
    criter.ma1 =  data.frame(MAE=ma1$MAE,model="MA1")
    criter.ma2 =  data.frame(MAE=ma2$MAE,model="MA2")
    criter.arma11 =  data.frame(MAE=arma11$MAE,model="ARMA11")     
    
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

cml.g0.AR1.boxcox <- function(z, L ){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(200,-1.01,1,1)
  # constraints=list(ineqA = A1, ineqB = B1)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR1.boxcox(z,L) 
  lk.new = function(param)  lk.g0.AR1.boxcox(z, param,L)
  grad.new = function(param)  grad.g0.AR1.boxcox(z, param,L)
  (res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "BFGS" , constraints=list(ineqA = A1, ineqB = B1)))
  est = res1$estimate
  code = res1$code
  alpha = est[1]; delta=est[2]; phi=est[3]
  
  mu.hat[1] = ( lambda*(  delta + phi*log( z0 ) ) + 1)^(1/lambda)
  
  for(t in 2:n) mu.hat[t] = ( lambda*( delta + phi*log( z[t-1] ) ) + 1)^(1/lambda)
  
  AIC = -2*lk.g0.AR1.boxcox(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.AR1.boxcox(z, est,L) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE, convergence = code) )
  
}
cml.g0.AR2.boxcox <- function(z,L){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.AR2.boxcox(z,L) 
  lk.new = function(param)  lk.g0.AR2.boxcox(z, param,L)
  grad.new = function(param)  grad.g0.AR.boxcox(z, param,L)
  (res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "BFGS"))
  est = res1$estimate
  alpha = est[1];delta=est[2];phi=est[3:4]
  
  mu.hat[1] = ( lambda*( delta + phi[1]*( log(z0) ) + phi[2]*( log(z0) ) ) + 1)^(1/lambda)
  mu.hat[2] = ( lambda*( delta + phi[1]*( log(z[1]) ) + phi[2]*( log(z0) ) ) + 1)^(1/lambda)
  
  for(t in 3:n) mu.hat[t] = ( lambda*( delta + phi[1]*log( z[t-1] ) + phi[2]*log( z[t-2] )  ) + 1)^(1/lambda)
  
  AIC = -2*lk.g0.AR2.boxcox(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.AR2.boxcox(z, est,L) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )
}
cml.g0.MA1.boxcox <- function(z,L ){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  A1 <- matrix(c(1, -1,0,0, 0,0,0,0,0,0,1,-1), 4,3)
  B1 <- c(200,-1.01,1,1)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA1.boxcox(z, L) 
  lk.new = function(param)  lk.g0.MA1.boxcox(z, param, L)
  grad.new = function(param)  grad.g0.MA1.boxcox(z, param, L)
  (res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "BFGS",constraints=list(ineqA = A1, ineqB = B1)))
  est = res1$estimate
  alpha = est[1];delta=est[2];theta=est[3]
  
  z0 = mean(z); mu0=1
  
  mu.hat[1] = ( lambda*( delta + theta*( log(z0) - log(mu0) ) ) + 1)^(1/lambda)
  
  for(t in 2:n) mu.hat[t] = ( lambda*( delta + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) ) + 1)^(1/lambda)
  
  AIC = -2*lk.g0.MA1.boxcox(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA1.boxcox(z, est,L) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4]  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
cml.g0.MA2.boxcox <- function(z,L ){
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.MA2.boxcox(z,L) 
  lk.new = function(param)  lk.g0.MA2.boxcox(z, param,L)
  grad.new = function(param)  grad.g0.MA.boxcox(z, param,L)
  try(
    (res1 <- maxLik(lk.new, grad = NULL, hess = NULL, param, method = "BFGS"))
  )
  
  if(res1$code == 100){
    est = rep(NA,length(param)); mu.hat = rep(NA,length(z))
  }else{
    
    est = res1$estimate
    alpha = est[1];delta=est[2];theta=est[3:4]
    
    z0 = mean(z); mu0=1
    mu.hat[1:2] = ( lambda*( delta + theta[1]*( log(z0) - log(mu0) ) + theta[2]*( log(z0) - log(mu0) ) ) + 1)^(1/lambda)
    
    for(t in 3:n) mu.hat[t] = ( lambda*( delta + theta[1]*( log(z[t-1]) - log(mu.hat[t-1]) ) + theta[2]*( log(z[t-2]) - log(mu.hat[t-2]) ) ) + 1)^(1/lambda)
    
  }
  
  AIC = -2*lk.g0.MA2.boxcox(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.MA2.boxcox(z, est,L) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4]  
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
cml.g0.ARMA11.boxcox <- function(z,L ){
  
  
  require(forecast)
  lambda = BoxCox.lambda(z)
  
  #y : série
  # metodo: m?todo de otimiza??o: BFGS, CG, SANN, NM e NR
  mu.hat = c();z0 = mean(z)
  n = length(z)
  
  param = g0Chute.ARMA11.boxcox(z,L) 
  lk.new = function(param)  lk.g0.ARMA11.boxcox(z, param,L)
  grad.new = function(param)  grad.g0.ARMA11.boxcox(z, param,L)
  (res1 <- maxLik(lk.new, grad = grad.new, hess = NULL, param, method = "BFGS"))
  est = res1$estimate
  alpha = est[1];delta=est[2];phi=est[3];theta=est[4] 
  
  z0 = mean(z); mu0=1
  
  mu.hat[1] = ( lambda*( delta + phi*log( z0 ) + theta*( log(z0) - log(mu0) ) ) + 1)^(1/lambda)
  
  for(t in 2:n) mu.hat[t] = ( lambda*( delta + phi*log( z[t-1] ) + theta*( log(z[t-1]) - log(mu.hat[t-1]) ) ) + 1)^(1/lambda)
  
  AIC = -2*lk.g0.ARMA11.boxcox(z, est,L) + 2*length(est) 
  BIC = -2*lk.g0.ARMA11.boxcox(z, est,L) + 2*log(n)
  
  medidas = acuracia(z, mu.hat)
  
  UQI = medidas[1]
  beta = medidas[2]
  MSE = medidas[3]
  MAE = medidas[4] 
  
  return( list(est=est, mu.hat=mu.hat, AIC = AIC, BIC = BIC, UQI = UQI, beta = beta, MSE = MSE, MAE = MAE) )  
}
auto.g0.boxcox <- function(y,L, criterio){
  # y: serie
  # criterio: criterio de informação AIC e BIC
  
  ar1 = cml.g0.AR1.boxcox(y,L);ar2 = cml.g0.AR2.boxcox(y,L)
  ma1 = cml.g0.MA1.boxcox(y,L); ma2 = cml.g0.MA2.boxcox(y,L)
  arma11 = cml.g0.ARMA11.boxcox(y,L)
  
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
    
  }else if(criterio == "UQI"){
    
    criter.ar1 =  data.frame(UQI=ar1$UQI,model="AR1")
    criter.ar2 =  data.frame(UQI=ar2$UQI,model="AR2")
    criter.ma1 =  data.frame(UQI=ma1$UQI,model="MA1")
    criter.ma2 =  data.frame(UQI=ma2$UQI,model="MA2")
    criter.arma11 =  data.frame(UQI=arma11$UQI,model="ARMA11")     
    
  }else if(criterio == "beta"){
    
    criter.ar1 =  data.frame(beta=ar1$beta,model="AR1")
    criter.ar2 =  data.frame(beta=ar2$beta,model="AR2")
    criter.ma1 =  data.frame(beta=ma1$beta,model="MA1")
    criter.ma2 =  data.frame(beta=ma2$beta,model="MA2")
    criter.arma11 =  data.frame(beta=arma11$beta,model="ARMA11")     
    
  }else if(criterio == "MSE"){
    
    criter.ar1 =  data.frame(MSE=ar1$MSE,model="AR1")
    criter.ar2 =  data.frame(MSE=ar2$MSE,model="AR2")
    criter.ma1 =  data.frame(MSE=ma1$MSE,model="MA1")
    criter.ma2 =  data.frame(MSE=ma2$MSE,model="MA2")
    criter.arma11 =  data.frame(MSE=arma11$MSE,model="ARMA11")     
    
  }else if(criterio == "MAE"){
    
    criter.ar1 =  data.frame(MAE=ar1$MAE,model="AR1")
    criter.ar2 =  data.frame(MAE=ar2$MAE,model="AR2")
    criter.ma1 =  data.frame(MAE=ma1$MAE,model="MA1")
    criter.ma2 =  data.frame(MAE=ma2$MAE,model="MA2")
    criter.arma11 =  data.frame(MAE=arma11$MAE,model="ARMA11")     
    
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




