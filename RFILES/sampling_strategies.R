# Different strategies to estimate parameters from median and higher quantile
# and sample

# Sample normal distribution
sample.normal <- function(q50, qh, ph, pp) {
  mean <- q50                  # median equals mean
  dist <- qnorm(ph)            # qh / sd
  sd   <- abs(qh-q50) / dist
  return(qnorm(pp, mean=mean, sd=sd)) # apply sampling
}

# Sample lognormal distribution
sample.lnormal <- function(ql, q50, qh, ph, pl, pp) {
  
  skew = (qh - q50) - (q50 - ql)
  
  if(skew > 0) {     # skewed to higher values
    fit     <- fit.lnormal(ql, q50, qh, pl, ph)
    samples <- fit[3] + 
        exp( qnorm(pp, mean=fit[1], sd=fit[2]) )
    } else {
    
    if(skew == 0) {  # not skewed
      samples <- sample.normal(q50, qh, ph, pp)
    } else {         # skewed to lower values
      samples <- sample.normal(q50, qh, ph, pp)
    }
  }

  return(samples)
}

# Sample beta distribution
sample.beta <- function(ql, q50, qh, ph, pl, pp, mn.l, mx) {
  fit <- fit.beta.lowerb(ql, q50, qh, pl, ph, mn.l,mx)
  
  samples <- fit[1] + (fit[2] - fit[1]) * qbeta(pp,fit[3],fit[4])
  return(samples)
}

#sample.lnormal <- function(ql, q50, qh, ph, pl, pp) {

#    if(ql < 0) {
#      q50.star <- log(q50 + abs(ql)) # apply offset and take logarithm
#      qh.star  <- log(qh  + abs(ql)) # apply offset and take logarithm
#    } else {
#      q50.star <- log(q50) # take logarithm
#      qh.star  <- log(qh)  # take logarithm
#   } 

#  samples.star <- sample.normal(q50.star, qh.star, ph, pp)
#  samples <- exp(samples.star) - ql

#  return(samples)
#}


# fit lognormal distribution with offsetd
fit.lnormal <- function(ql, q50, qh, pl, ph) {

  f <- function(sigma) {
    e.mu = (qh - q50) / (exp(qhsn * sigma) -1)
    tau  = q50 - e.mu
    return( ql - (tau + e.mu * exp(qlsn * sigma)) )
  }
  
  qlsn = qnorm(pl)
  qhsn = qnorm(ph)
  
  lower = 0.001
  upper = 1000*(qh-ql)
  
  f(lower)
  
  rc <- uniroot(f, lower=lower, upper=upper, tol=0.0001)
  
  sigma = rc$root
  e.mu  = (qh - q50) / (exp(qhsn * sigma) -1)
  tau   = q50 - e.mu
  mu    = log(e.mu)
  
  return(c(mu=mu,sigma=sigma,tau=tau))
} 

# fit beta distribution
fit.beta.lowerb <- function(ql, q50, qh, pl, ph, mn.l, mx) {
  require(rriskDistributions)
  
  f <- function(mn) {
    q    <- (ql  - mn)/(mx - mn)
    q[2] <- (q50 - mn)/(mx - mn)
    q[3] <- (qh  - mn)/(mx - mn)
    p    <- c(pl,0.5,ph)

    sh <- get.beta.par(p=p[1:2], q=q[1:2], plot=F, show.output=F)
    
    return(sum(q - qbeta(p,sh[1],sh[2])))
  }
  
  lower=mn.l
  upper=ql
  
  if(ql==q50) {
    mn = ql
    
    q    <- (q50 - mn)/(mx - mn) + 0.00001
    q[2] <- (qh  - mn)/(mx - mn)
    p    <- c(0.5,ph)

  } else {
    if(f(lower)*f(upper)<=0) {
      rc <- uniroot(f, lower=lower, upper=upper, tol=0.0001)
      mn <- rc$root
    } else {
      mn <- lower
    }

    q    <- (ql  - mn)/(mx - mn)
    q[2] <- (q50 - mn)/(mx - mn)
    q[3] <- (qh  - mn)/(mx - mn)
    p    <- c(pl,0.5,ph)
  }
  
  fit.weights <- c(rep(1,length(p)-1),3)
  
  sh <- 
    get.beta.par(p=p, q=q, plot=F, show.output=F, fit.weights = fit.weights)
  
  return(c(mn,mx,sh))
} 


