## from http://www.medepi.org/epitools/rfunctions/cibinom.html

uci.binomial <- function(x,n,a=0,b=1,alpha=.05){
  ## Finds the exact upper confidence limit for binomial count
  ## R v0.90.1, Tomas Aragon, created 01/01/2000, edited
  ## x = observed count (non-negative integer)
  ## n = number of Bernoulli trials
  ## a = default lower bound for bisection method
  ## b = default upper bound for bisection method
  ## alpha = .05 (default), for 95% confidence interval
  ## This function is used by ci.binomial()
  
  if(x==0) answer <- 1-(alpha)^(1/n)
  else{
    answer <- (a+b)/2
    if(abs(a-b) >= 0.00000001){
      lefta <- pbinom(x,n,a)-alpha/2
      centera <- pbinom(x,n,answer)-alpha/2
      righta <- pbinom(x,n,b)-alpha/2
      if(lefta*centera < 0){
        answer <- uci.binomial(x,n,a,answer)
      }
      else{
        answer <- uci.binomial(x,n,answer,b)
      }
    }
  }
  answer
}

lci.binomial <- function(x,n,a=0,b=1,alpha=.05){
  ## Finds the exact lower confidence limit for binomial count
  ## R v0.90.1, Tomas Aragon, created 01/01/2000, edited
  ## x = observed count (non-negative integer)
  ## n = number of Bernoulli trials
  ## a = default lower bound for bisection method
  ## b = default upper bound for bisection method
  ## alpha = .05 (default), for 95% confidence interval
  ## This function is used by ci.binomial()

  if(x==0) answer <- 0
  else{
    answer <- (a+b)/2
    if(abs(a-b) >= 0.00000001){
      lefta <- 1-pbinom(x,n,a)+dbinom(x,n,a)-alpha/2
      centera <- 1-pbinom(x,n,answer)+dbinom(x,n,answer)-alpha/2
      righta <- 1-pbinom(x,n,b)+dbinom(x,n,b)-alpha/2
      if(lefta*centera < 0){
        answer <- lci.binomial(x,n,a,answer)
      }
      else{
        answer <- lci.binomial(x,n,answer,b)
      }
    }
  }
  answer
}

ci.binomial <- function(x,n,alpha=0.05){
  ## Finds the exact 95% confidence interval for binomial count
  ## R v0.90.1, Tomas Aragon, created 01/01/2000, edited
  ## x = observed count (non-negative integer)
  ## number of Bernoulli trials
  ## alpha = .05 (default), for 95% confidence interval
  ## Output includes x and CIs, and p=x/n and CIs
  ## Uses lci.binomial() and uci.binomial() for bisection method
  xn <- cbind(x,n)
  ci.pois <- function(xx,aa=alpha){
    lci <- lci.binomial(xx[1],xx[2],alpha=aa)
    uci <- uci.binomial(xx[1],xx[2],alpha=aa)
    c(lci=lci,uci=uci)
  }
  prob <- t(apply(xn,1,ci.pois))
  ans <- cbind(as.vector(x),n,prob*n,as.vector(x/n),prob)
  dimnames(ans)[[2]] <- c("x","n","x.lci","x.uci","p=x/n","p.lci","p.uci")
  ans
}
