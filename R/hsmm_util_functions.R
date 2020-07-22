
.rlog <- function(n,p,M=10000) sample(1:M,n,TRUE,.dlog(1:M,p))

.dlog <- function(x,p) (-1 / log(1-p)) * p^x / x




.logdistrfit <- function(wt) {
  xbar = sum(wt*(1:length(wt)))
  fn <- function(p) xbar + p/((1-p)*log(1-p))
  uniroot(fn,c(1e-10,1-1e-10))$root
}      


#rng for a shifted poisson distrubtion
.rpois.hsmm.sojourn <- function(n,lambda,shift) rpois(n,lambda)+shift

#density function for a shifted poisson distribution
.dpois.hsmm.sojourn <- function(x=NULL,lambda,shift,log=FALSE)   {
  if(shift<0) stop(".dpois.hsmm.sojourn: shift must be > 0")
  if(log) dpois(x-shift,lambda,log=TRUE)
  else dpois(x-shift,lambda)
}
